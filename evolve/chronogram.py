#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""
Real dated branch lengths from the Fish Tree of Life chronogram.

Replaces the Grafen (topology-derived) branch lengths in the neutral-Mk test with
REAL divergence times: the time-calibrated phylogeny of ray-finned fishes
(Rabosky et al. 2018, treePL penalized-likelihood dating; branch lengths in millions
of years). For our all-fish taxa this is the right, high-coverage chronogram -- and,
crucially, the divergence times are INDEPENDENT of our trait data (estimated from
molecular clocks + fossils), so using them cannot fit the convergence result.

We download the full ~11.6k-species chronogram once, prune to our taxa (collapsing
suppressed nodes while SUMMING their branch lengths to preserve distances), and
normalise so the deepest root-to-tip path = 1 (the absolute timescale is absorbed
into the Mk rate; relative divergence-time structure is what informs the null).
"""

from __future__ import annotations

import lzma
import re
import sys
from pathlib import Path

import requests

from polarity import Node, parse_newick, _leaves

CACHE = Path(__file__).resolve().parent / "cache"
URL = "https://fishtreeoflife.org/downloads/actinopt_12k_treePL.tre.xz"
TRE = CACHE / "actinopt_12k_treePL.tre"


def load_chronogram() -> Node:
    if not TRE.exists():
        CACHE.mkdir(parents=True, exist_ok=True)
        data = lzma.decompress(requests.get(URL, timeout=180).content)
        TRE.write_bytes(data)
    sys.setrecursionlimit(200000)
    return parse_newick(TRE.read_text())


def all_tip_labels(tree: Node) -> set:
    out = set()
    stack = [tree]
    while stack:
        n = stack.pop()
        if n.is_leaf:
            if n.label:
                out.add(n.label)
        else:
            stack.extend(n.children)
    return out


def prune_to(node: Node, keep: set):
    """Return a copy of `node` pruned to leaves in `keep`; degree-1 nodes collapsed
    with branch lengths summed (distances preserved). None if nothing kept."""
    if node.is_leaf:
        return Node(label=node.label, brlen=node.brlen) if node.label in keep else None
    kept = []
    for c in node.children:
        pc = prune_to(c, keep)
        if pc is not None:
            kept.append(pc)
    if not kept:
        return None
    if len(kept) == 1:
        child = kept[0]
        child.brlen = (child.brlen or 0.0) + (node.brlen or 0.0)
        return child
    new = Node(label=node.label, brlen=node.brlen)
    new.children = kept
    return new


def _max_depth(node: Node, acc: float = 0.0) -> float:
    d = acc + (node.brlen or 0.0)
    if node.is_leaf:
        return d
    return max(_max_depth(c, d) for c in node.children)


def normalize_unit_depth(tree: Node) -> None:
    """Scale all branch lengths so the deepest root-to-tip path = 1 (in place)."""
    depth = _max_depth(tree) or 1.0

    def scale(n: Node):
        if n.brlen is not None:
            n.brlen = n.brlen / depth
        for c in n.children:
            scale(c)
    scale(tree)


def branch_length_dict(tree: Node) -> dict:
    """id(node) -> branch length, for mk.* (root -> 0.0)."""
    out = {}

    def walk(n: Node):
        out[id(n)] = n.brlen if n.brlen is not None else 0.0
        for c in n.children:
            walk(c)
    walk(tree)
    return out


def to_newick(node: Node) -> str:
    """Serialise a tree (with branch lengths) back to a Newick string."""
    if node.is_leaf:
        s = node.label or ""
    else:
        s = "(" + ",".join(to_newick(c) for c in node.children) + ")" + (node.label or "")
    if node.brlen is not None:
        s += ":" + repr(node.brlen)
    return s


# Mammal chronogram (PHYLACINE / Upham 2019, dated; Nexus with integer translate).
MAMMAL_NEX = CACHE / "mammal_small.nex"
MAMMAL_NWK = CACHE / "mammal_mcc.nwk"


def load_mammal_chronogram() -> Node:
    """First dated tree from the PHYLACINE mammal Nexus, tips relabelled via the
    translate block. Cached as a single relabelled Newick after first parse."""
    sys.setrecursionlimit(200000)
    if MAMMAL_NWK.exists():
        return parse_newick(MAMMAL_NWK.read_text())
    transl: dict = {}
    in_tr = False
    first_tree = None
    with open(MAMMAL_NEX, "r", errors="replace") as f:
        for line in f:
            s = line.strip()
            low = s.lower()
            if low.startswith("translate"):
                in_tr = True
                continue
            if in_tr:
                ended = s.startswith(";")
                body = s.rstrip(";").rstrip(",").strip()
                if body:
                    parts = body.split(None, 1)
                    if len(parts) == 2 and parts[0].isdigit():
                        transl[parts[0]] = parts[1].strip().rstrip(",")
                if ended or low.startswith("tree"):
                    in_tr = False
            if low.startswith("tree") and first_tree is None:
                rest = s[s.find("=") + 1:].strip()
                rest = re.sub(r"^\[&[^]]*\]", "", rest).strip()  # strip [&R]
                first_tree = rest
                break
    tree = parse_newick(first_tree)
    for lf in _leaves(tree):
        lf.label = transl.get(lf.label, lf.label)
    MAMMAL_NWK.write_text(to_newick(tree) + ";")
    return tree


# Bird chronogram (Jetz et al. 2012, dated; VertLife Stage2 Hackett full trees).
# The download is a 1000-tree pseudo-posterior archive; we stream just the FIRST
# tree out of the zip (decompress until the first newline) and cache it, so we
# never have to pull the whole ~190 MB file.
BIRD_ZIP_URL = ("https://data.vertlife.org/birdtree/Stage2/"
                "HackettStage2_0001_1000.zip")
BIRD_NWK = CACHE / "bird_tree.nwk"


def _stream_first_tree(url: str) -> str:
    """Stream the first member of a zip and return its first line (one Newick tree).

    Only enough compressed bytes to produce one tree are downloaded; the rest of
    the archive is never fetched. The member is a plain `one tree per line` .tre."""
    import struct
    import zlib

    r = requests.get(url, stream=True, timeout=300)
    try:
        it = r.iter_content(chunk_size=65536)
        head = next(it)
        if head[:4] != b"PK\x03\x04":
            raise ValueError("not a zip local file header")
        method = struct.unpack("<H", head[8:10])[0]
        name_len = struct.unpack("<H", head[26:28])[0]
        extra_len = struct.unpack("<H", head[28:30])[0]
        if method != 8:
            raise ValueError(f"unexpected zip compression method {method}")
        data_off = 30 + name_len + extra_len
        dec = zlib.decompressobj(-15)  # raw DEFLATE
        out = dec.decompress(head[data_off:])
        while b"\n" not in out:
            try:
                out += dec.decompress(next(it))
            except StopIteration:
                break
        return out.split(b"\n", 1)[0].decode("latin1").strip()
    finally:
        r.close()


def load_bird_chronogram() -> Node:
    """First posterior sample of the Jetz et al. 2012 dated bird tree (Hackett
    backbone, 9993 OTUs), streamed from the VertLife Stage2 archive and cached as a
    single relabelled Newick.

    Honest provenance: this is ONE pseudo-posterior tree, not the MCC consensus.
    Divergence times are still molecular-clock + fossil estimates, independent of
    the trait data, so using them cannot fit the convergence result (a different
    posterior sample would shift node ages slightly, not the trait codings)."""
    sys.setrecursionlimit(200000)
    if BIRD_NWK.exists():
        return parse_newick(BIRD_NWK.read_text())
    CACHE.mkdir(parents=True, exist_ok=True)
    newick = _stream_first_tree(BIRD_ZIP_URL)
    BIRD_NWK.write_text(newick if newick.endswith(";") else newick + ";")
    return parse_newick(newick)


def prune_dated(full: Node, taxa: list[str]):
    """Prune a chronogram to `taxa` (space-separated). Returns
    (pruned_tree, brlen_dict, matched_underscore_labels) or (None, None, set())."""
    labels = all_tip_labels(full)
    want = {t.replace(" ", "_") for t in taxa}
    keep = want & labels
    if len(keep) < 4:
        return None, None, keep
    pruned = prune_to(full, keep)
    normalize_unit_depth(pruned)
    return pruned, branch_length_dict(pruned), keep


def dated_subtree(taxa: list[str]):
    """Fish chronogram subtree for `taxa`."""
    return prune_dated(load_chronogram(), taxa)


def dated_subtree_mammal(taxa: list[str]):
    """Mammal chronogram subtree for `taxa`."""
    return prune_dated(load_mammal_chronogram(), taxa)


def dated_subtree_bird(taxa: list[str]):
    """Bird chronogram subtree for `taxa`."""
    return prune_dated(load_bird_chronogram(), taxa)


if __name__ == "__main__":
    tree = load_chronogram()
    labels = all_tip_labels(tree)
    print(f"chronogram tips: {len(labels)}")
    test = ["Danio rerio", "Gambusia marshi", "Tetraodon nigroviridis",
            "Apteronotus albifrons", "Gymnotus carapo", "Amia calva"]
    pruned, bl, keep = dated_subtree(test)
    print(f"requested {len(test)}, matched {len(keep)}: {sorted(keep)}")
    if pruned:
        from polarity import _leaves
        leaves = list(_leaves(pruned))
        print(f"pruned leaves: {len(leaves)}")
        print("normalised branch lengths (sample):",
              [round(l.brlen, 3) if l.brlen else 0 for l in leaves[:6]])
