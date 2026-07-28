#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: 2026 James Ray Hawkins

"""
Ancestral-state polarity via Fitch parsimony  (PARAMETER-FREE -- un-fittable).

The recurring trap across this whole project (cow~horse, "fin present"): similarity
from a RETAINED ANCESTRAL trait (plesiomorphy) is NOT convergence. Convergence is
the SAME DERIVED state arising INDEPENDENTLY in different lineages (homoplasy).

To tell them apart you must reconstruct ancestral states on a phylogeny and COUNT
independent origins of a state. We use Fitch parsimony because it has ZERO free
parameters -- there is nothing to tune toward a desired answer. The tree is supplied
externally (Open Tree), independent of the trait data. So the convergence count is a
deterministic function of (tree, tip states): we cannot fit it.

Outputs per character:
  * parsimony score (steps)      = number of independent state changes on the tree
  * consistency index CI = 1/steps (binary): CI=1 -> single origin, NOT convergence;
                                              CI<1 -> homoplasy = CONVERGENCE
  * independent origins of each state (down-pass), and the inferred ROOT (=ancestral)

We never assume which state is derived; the reconstruction tells us.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple


# ── Tree ────────────────────────────────────────────────────────────────────
@dataclass
class Node:
    label: Optional[str] = None
    children: List["Node"] = field(default_factory=list)
    brlen: Optional[float] = None   # branch length above this node (e.g. chronogram MY)
    # filled during reconstruction
    state_set: Set[str] = field(default_factory=set)
    state: Optional[str] = None

    @property
    def is_leaf(self) -> bool:
        return not self.children


def parse_newick(s: str) -> Node:
    """Minimal Newick parser (handles internal labels, nested single-child nodes,
    branch lengths, and the Open Tree `name_ott12345` tip style)."""
    s = s.strip().rstrip(";")
    pos = 0

    def parse_node() -> Node:
        nonlocal pos
        node = Node()
        if s[pos] == "(":
            pos += 1  # consume '('
            while True:
                node.children.append(parse_node())
                if s[pos] == ",":
                    pos += 1
                    continue
                if s[pos] == ")":
                    pos += 1
                    break
        # read token (name and/or :branchlength) until , ) ;
        m = re.match(r"[^,():;]*", s[pos:])
        token = m.group(0)
        pos += len(token)
        if ":" in token:
            name, bl = token.split(":", 1)
            node.label = name or None
            try:
                node.brlen = float(bl)
            except ValueError:
                node.brlen = None
        else:
            node.label = token or None
        # a branch length may follow a ')' as a separate ':...'
        if pos < len(s) and s[pos] == ":":
            m2 = re.match(r":[^,();]*", s[pos:])
            pos += len(m2.group(0))
            try:
                node.brlen = float(m2.group(0)[1:])
            except ValueError:
                pass
        return node

    return parse_node()


def ott_of(label: Optional[str]) -> Optional[str]:
    if not label:
        return None
    m = re.search(r"ott(\d+)", label)
    return m.group(1) if m else None


# ── Fitch parsimony ─────────────────────────────────────────────────────────
def fitch_up(node: Node, tip_state: Dict[str, str], all_states: Set[str]) -> int:
    """Up-pass: assign state_set to every node; return number of changes (steps)."""
    if node.is_leaf:
        key = ott_of(node.label) or node.label
        st = tip_state.get(key)
        node.state_set = {st} if st is not None else set(all_states)  # missing = any
        return 0
    steps = 0
    child_sets = []
    for c in node.children:
        steps += fitch_up(c, tip_state, all_states)
        child_sets.append(c.state_set)
    inter = set.intersection(*child_sets) if child_sets else set()
    if inter:
        node.state_set = inter
    else:
        node.state_set = set.union(*child_sets)
        steps += 1
    return steps


def fitch_down(node: Node, parent_state: Optional[str]) -> None:
    """Down-pass: pick a definite state per node (Fitch tie-break toward parent)."""
    if parent_state is not None and parent_state in node.state_set:
        node.state = parent_state
    else:
        # deterministic pick: sorted, so reruns are identical (no fitting)
        node.state = sorted(node.state_set)[0]
    for c in node.children:
        fitch_down(c, node.state)


def count_origins(node: Node, target: str, parent_state: Optional[str]) -> int:
    """Count edges parent!=target -> child==target (independent origins of `target`)."""
    n = 0
    if node.state == target and parent_state != target:
        n += 1
    for c in node.children:
        n += count_origins(c, target, node.state)
    return n


@dataclass
class PolarityResult:
    character: str
    steps: int
    consistency_index: float
    root_state: Optional[str]
    origins: Dict[str, int]
    n_tips_scored: int
    states: List[str]

    @property
    def is_convergent(self) -> bool:
        # the derived (non-root) state arose >= 2 times independently
        derived = [s for s in self.states if s != self.root_state]
        return any(self.origins.get(s, 0) >= 2 for s in derived)


def analyze(character: str, tree: Node, tip_state: Dict[str, str]) -> PolarityResult:
    """Full polarity analysis of one discrete character on a tree."""
    states = sorted(set(tip_state.values()))
    scored = sum(1 for n_ in _leaves(tree)
                 if (ott_of(n_.label) or n_.label) in tip_state)
    steps = fitch_up(tree, tip_state, set(states))
    # root state: if ambiguous, pick the majority tip state (reported, deterministic)
    if len(tree.state_set) == 1:
        root_state = next(iter(tree.state_set))
    else:
        from collections import Counter
        c = Counter(tip_state.values())
        root_state = c.most_common(1)[0][0]
    fitch_down(tree, root_state)
    origins = {s: count_origins(tree, s, None) for s in states}
    ci = (1.0 / steps) if steps > 0 else 1.0   # binary: m=1
    return PolarityResult(character, steps, ci, root_state, origins, scored, states)


def _leaves(node: Node):
    if node.is_leaf:
        yield node
    for c in node.children:
        yield from _leaves(c)


# ── Homology / independence via MRCA (parameter-free) ─────────────────────────
def _root_state(tree: Node, tip_state: Dict[str, str]) -> str:
    if len(tree.state_set) == 1:
        return next(iter(tree.state_set))
    from collections import Counter
    return Counter(tip_state.values()).most_common(1)[0][0]


def _index_parents(tree: Node):
    """Return (parent_by_id, leaf_node_by_key)."""
    parent: Dict[int, Node] = {}
    leaves: Dict[str, Node] = {}

    def walk(n: Node):
        if n.is_leaf:
            leaves[ott_of(n.label) or n.label] = n
        for c in n.children:
            parent[id(c)] = n
            walk(c)
    walk(tree)
    return parent, leaves


def _mrca(a: Node, b: Node, parent: Dict[int, Node]) -> Node:
    chain = []
    n = a
    while True:
        chain.append(id(n))
        if id(n) not in parent:
            break
        n = parent[id(n)]
    seen = set(chain)
    n = b
    while id(n) not in seen and id(n) in parent:
        n = parent[id(n)]
    return n


def _mrca_of(keys, leaves: Dict[str, Node], parent: Dict[int, Node]) -> Optional[Node]:
    nodes = [leaves[k] for k in keys if k in leaves]
    if not nodes:
        return None
    m = nodes[0]
    for n in nodes[1:]:
        m = _mrca(m, n, parent)
    return m


def _is_descendant(tip: Node, anc: Node, parent: Dict[int, Node]) -> bool:
    n = tip
    while True:
        if n is anc:
            return True
        if id(n) not in parent:
            return False
        n = parent[id(n)]


def dollo_recode(tree: Node, tip_state: Dict[str, str], present_state: str = "present"):
    """Recode 'absent' as INAPPLICABLE where the structure never existed.

    DOLLO assumption (uniform model choice, stated not tuned): a complex structure
    originates ONCE and is only ever lost. So its applicable region is the clade
    descended from its single origin = MRCA of all taxa that HAVE it. A taxon coding
    the non-present state OUTSIDE that clade never had the structure -> inapplicable
    (its 'absent' is primitive, not a loss). One INSIDE the clade is a true loss.

    Parameter-free and assumption-explicit: the applicable clade is pure topology
    (MRCA of present taxa); no external taxonomy is injected. Returns
    (recoded_tip_state with inapplicable taxa REMOVED, set_of_inapplicable_keys).
    """
    parent, leaves = _index_parents(tree)
    present_keys = [k for k, s in tip_state.items() if s == present_state]
    if not present_keys:
        return dict(tip_state), set()
    applicable = _mrca_of(present_keys, leaves, parent)
    inapplicable, recoded = set(), dict(tip_state)
    for k, s in tip_state.items():
        if s != present_state and k in leaves:
            if not _is_descendant(leaves[k], applicable, parent):
                inapplicable.add(k)
                del recoded[k]          # treat as missing/inapplicable for this character
    return recoded, inapplicable


def randomization_test(tree: Node, tip_state: Dict[str, str],
                       n_perm: int = 999, seed: int = 0) -> dict:
    """Phylogenetic randomization (PTP) test of homoplasy significance.

    Keeps the tree and the state frequencies FIXED, randomly permutes WHICH tips
    carry each state, and recomputes the parsimony step count (independent changes)
    via Fitch. Builds a null distribution and compares the observed count.

    Parameter-free (only n_perm, the Monte-Carlo budget). Addresses "maybe the
    state changed a lot just because change is easy / the tree is big":
      observed << null  -> the changes are CLUSTERED (phylogenetic signal); the
                           independent origins are constrained/real, not random scatter.
      observed ~= null  -> the state is as scattered as chance = labile, the apparent
                           convergence is not more than random lability.
    Returns observed_steps, null mean/sd, and p_fewer / p_more (one-sided, +1 smoothed).

    LIMIT (stated, not faked): the randomization null is the NO-STRUCTURE extreme.
    Testing excess homoplasy beyond a NEUTRAL Mk process would need branch lengths +
    a rate -- the Open Tree topology has neither, so we do not claim it.
    """
    import random
    import statistics
    keys = list(tip_state)
    vals = list(tip_state.values())
    all_states = set(vals)
    observed = fitch_up(tree, tip_state, all_states)
    rng = random.Random(seed)
    null = []
    for _ in range(n_perm):
        rng.shuffle(vals)
        null.append(fitch_up(tree, dict(zip(keys, vals)), all_states))
    le = sum(1 for s in null if s <= observed)
    ge = sum(1 for s in null if s >= observed)
    return {
        "observed_steps": observed,
        "null_mean": round(statistics.mean(null), 2),
        "null_sd": round(statistics.pstdev(null), 2),
        "p_fewer": round((le + 1) / (n_perm + 1), 4),   # clustered / signal
        "p_more": round((ge + 1) / (n_perm + 1), 4),    # over-dispersed
        "n_perm": n_perm,
    }


class TreeHomology:
    """Per-character Fitch reconstruction + MRCA-state lookup, for folding homology
    into a pairwise similarity metric.

    Build once from a tree + per-character tip states; then mrca_state(a, b, char)
    returns the reconstructed state at MRCA(a, b) for that character. A shared
    derived state between a and b is INDEPENDENT (convergence) iff mrca_state != it.
    Parameter-free (Fitch); the tree is supplied externally and never sees nothing
    it can fit.
    """

    def __init__(self, tree: Node, char_tip_states: Dict[str, Dict[str, str]]):
        self.tree = tree
        self.parent, self.leaves = _index_parents(tree)
        self.node_state: Dict[str, Dict[int, Optional[str]]] = {}
        for char, tip in char_tip_states.items():
            states = sorted(set(tip.values()))
            if len(states) < 1:
                continue
            fitch_up(tree, tip, set(states))
            fitch_down(tree, _root_state(tree, tip))
            snap: Dict[int, Optional[str]] = {}
            self._snapshot(tree, snap)
            self.node_state[char] = snap

    def _snapshot(self, node: Node, snap: Dict[int, Optional[str]]) -> None:
        snap[id(node)] = node.state
        for c in node.children:
            self._snapshot(c, snap)

    def mrca_state(self, a_key: str, b_key: str, char: str) -> Optional[str]:
        if char not in self.node_state:
            return None
        if a_key not in self.leaves or b_key not in self.leaves:
            return None
        m = _mrca(self.leaves[a_key], self.leaves[b_key], self.parent)
        return self.node_state[char].get(id(m))


def homology_filter(tree: Node, tip_state: Dict[str, str], derived_state: str):
    """Classify every pair of taxa sharing `derived_state` as either an INDEPENDENT
    acquisition (convergence) or SHARED-BY-DESCENT (homology, one event).

    A shared derived state is convergence ONLY if the pair's MRCA did NOT already
    have that state. If the MRCA is reconstructed in `derived_state`, the two taxa
    inherited it from a common ancestor -> not independent -> NOT convergence.
    Parameter-free: uses only tree topology + Fitch reconstruction.

    Returns (independent_pairs, shared_pairs, carriers).
    """
    states = sorted(set(tip_state.values()))
    fitch_up(tree, tip_state, set(states))
    fitch_down(tree, _root_state(tree, tip_state))
    parent, leaves = _index_parents(tree)
    carriers = [k for k, n in leaves.items()
                if n.state == derived_state and tip_state.get(k) == derived_state]
    independent, shared = [], []
    for i in range(len(carriers)):
        for j in range(i + 1, len(carriers)):
            a, b = carriers[i], carriers[j]
            m = _mrca(leaves[a], leaves[b], parent)
            (independent if m.state != derived_state else shared).append((a, b))
    return independent, shared, carriers


# ── Self-test on a hand-checkable toy (tests the ALGORITHM, not biology) ──────
if __name__ == "__main__":
    print("Fitch parsimony self-test (hand-checkable):\n")

    # Toy 1: ((A,B),(C,D)) with absent at A and C -> 'absent' arose TWICE
    t1 = parse_newick("((A,B),(C,D));")
    tips1 = {"A": "absent", "B": "present", "C": "absent", "D": "present"}
    r1 = analyze("toy_two_losses", t1, tips1)
    print(f"  ((A,B),(C,D)) A,C=absent B,D=present")
    print(f"    steps={r1.steps} (expect 2)  CI={r1.consistency_index:.2f} (expect 0.50)")
    print(f"    root={r1.root_state}  origins={r1.origins}  convergent={r1.is_convergent}")
    assert r1.steps == 2, "expected 2 steps"

    # Toy 2: (((A,B),C),D) with absent only at A,B (one clade) -> single origin
    t2 = parse_newick("(((A,B),C),D);")
    tips2 = {"A": "absent", "B": "absent", "C": "present", "D": "present"}
    r2 = analyze("toy_single_loss", t2, tips2)
    print(f"\n  (((A,B),C),D) A,B=absent C,D=present")
    print(f"    steps={r2.steps} (expect 1)  CI={r2.consistency_index:.2f} (expect 1.00)")
    print(f"    root={r2.root_state}  origins={r2.origins}  convergent={r2.is_convergent}")
    assert r2.steps == 1, "expected 1 step"
    assert not r2.is_convergent, "single origin is NOT convergence"

    # Open Tree style labels with ott ids
    t3 = parse_newick("((Gus_ott1,Hus_ott2),(Ius_ott3,Jus_ott4));")
    tips3 = {"1": "absent", "2": "present", "3": "absent", "4": "present"}
    r3 = analyze("toy_ott", t3, tips3)
    print(f"\n  ott-id tips, absent at 1,3 -> steps={r3.steps} (expect 2), "
          f"convergent={r3.is_convergent}")
    assert r3.steps == 2

    # Homology/MRCA: independent losses vs shared-by-descent
    # (((A,B),C),(D,E)) with absent at A and D -> two INDEPENDENT losses
    th = parse_newick("(((A,B),C),(D,E));")
    indep, shared, carr = homology_filter(
        th, {"A": "absent", "B": "present", "C": "present",
             "D": "absent", "E": "present"}, "absent")
    print(f"\n  homology: A,D absent in separate clades -> independent={indep} shared={shared}")
    assert ("A", "D") in indep and not shared, "A,D should be an independent (convergent) pair"

    # same tree, absent at A and B (one clade) -> shared by descent, NOT independent
    th2 = parse_newick("(((A,B),C),(D,E));")
    indep2, shared2, _ = homology_filter(
        th2, {"A": "absent", "B": "absent", "C": "present",
              "D": "present", "E": "present"}, "absent")
    print(f"  homology: A,B absent in one clade   -> independent={indep2} shared={shared2}")
    assert ("A", "B") in shared2 and not indep2, "A,B share one loss -> not independent"

    # Dollo inapplicability: a basal 'absent' taxon OUTSIDE the present-clade = never had it
    # (A,((B,C),(D,E))): present B,C,E ; absent A (basal, sister) and D (nested loss)
    td = parse_newick("(A,((B,C),(D,E)));")
    recoded, inappl = dollo_recode(
        td, {"A": "absent", "B": "present", "C": "present",
             "D": "absent", "E": "present"}, present_state="present")
    print(f"\n  dollo: present B,C,E; absent A(basal),D(nested) -> inapplicable={sorted(inappl)}")
    assert inappl == {"A"}, "A is basal to the present-clade -> inapplicable (never had it)"
    assert recoded.get("D") == "absent", "D is nested inside -> a true loss, kept"

    # Randomization test mechanics: a clustered trait has FEWER steps than random
    tc = parse_newick("((((A,B),C),(D,E)),((F,G),(H,I)));")
    clustered = {k: ("absent" if k in ("A", "B", "C") else "present")
                 for k in "ABCDEFGHI"}
    rt = randomization_test(tc, clustered, n_perm=499, seed=1)
    print(f"\n  randomization: clustered absent(A,B,C) obs_steps={rt['observed_steps']} "
          f"null_mean={rt['null_mean']} p_fewer={rt['p_fewer']}")
    assert rt["observed_steps"] <= rt["null_mean"], "clustered trait should have <= random steps"

    print("\nAll self-tests passed. Parameter-free: nothing here can be tuned.")
