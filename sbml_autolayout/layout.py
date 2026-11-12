"""Deterministic, rule-based layout for single metabolic pathways.

Pipeline
  1. classify participants into main / side compounds (side compounds become per-reaction clones)
  2. build the main-metabolite skeleton (directed, one edge per main substrate->product pair)
  3. detect cycles (drawn as circles) and contract them; pick the backbone as the longest path
  4. place backbone on a straight line, cycles as circles with a canonical entry point/direction
  5. apply pins / previous layout (rigid translation of the backbone block + exact pins)
  6. attach remaining metabolites as orthogonal branches (greedy, collision-checked)
  7. place side-compound clones with a fixed per-reaction template
  8. build edge geometry (lines, Bezier arcs on cycles) and place labels in discrete slots
"""
from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field, replace

import networkx as nx

from .known import find_main_path, norm_name, select_known_cycles
from .model import DRAIN, EMPTY, Metabolite, Pathway, Roles, canonical, classify, looks_pseudo, sha

ALGORITHM_VERSION = "sbml-autolayout-0.10"
Pt = tuple[float, float]


# ----------------------------------------------------------------------------- vector helpers
def add(a: Pt, b: Pt) -> Pt: return (a[0] + b[0], a[1] + b[1])
def sub(a: Pt, b: Pt) -> Pt: return (a[0] - b[0], a[1] - b[1])
def mul(a: Pt, k: float) -> Pt: return (a[0] * k, a[1] * k)
def neg(a: Pt) -> Pt: return (-a[0], -a[1])
def dist(a: Pt, b: Pt) -> float: return math.hypot(a[0] - b[0], a[1] - b[1])
def dot(a: Pt, b: Pt) -> float: return a[0] * b[0] + a[1] * b[1]
def lerp(a: Pt, b: Pt, t: float) -> Pt: return (a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t)
def mean(ps): return (sum(p[0] for p in ps) / len(ps), sum(p[1] for p in ps) / len(ps))


def unit(a: Pt) -> Pt:
    n = math.hypot(*a)
    return (a[0] / n, a[1] / n) if n > 1e-9 else (0.0, 0.0)


def seg_point_dist(p: Pt, a: Pt, b: Pt) -> float:
    ab = sub(b, a)
    L2 = dot(ab, ab)
    if L2 < 1e-12:
        return dist(p, a)
    t = max(0.0, min(1.0, dot(sub(p, a), ab) / L2))
    return dist(p, add(a, mul(ab, t)))


def segs_cross(p1: Pt, p2: Pt, q1: Pt, q2: Pt) -> bool:
    def orient(a, b, c):
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    d1, d2 = orient(q1, q2, p1), orient(q1, q2, p2)
    d3, d4 = orient(p1, p2, q1), orient(p1, p2, q2)
    return (d1 * d2 < -1e-9) and (d3 * d4 < -1e-9)


def bezier_points(seg, n=8):
    if seg[0] == "L":
        return [lerp(seg[1], seg[2], i / n) for i in range(n + 1)]
    _, p0, c1, c2, p3 = seg
    out = []
    for i in range(n + 1):
        t = i / n
        u = 1 - t
        out.append((u**3 * p0[0] + 3 * u * u * t * c1[0] + 3 * u * t * t * c2[0] + t**3 * p3[0],
                    u**3 * p0[1] + 3 * u * u * t * c1[1] + 3 * u * t * t * c2[1] + t**3 * p3[1]))
    return out


def arc_bezier(c: Pt, r: float, t0: float, t1: float):
    """Cubic Bezier approximating the circular arc from angle t0 to t1 (|t1-t0| < pi)."""
    d = t1 - t0
    k = 4 / 3 * math.tan(d / 4)
    p0 = (c[0] + r * math.cos(t0), c[1] + r * math.sin(t0))
    p3 = (c[0] + r * math.cos(t1), c[1] + r * math.sin(t1))
    c1 = (p0[0] - k * r * math.sin(t0), p0[1] + k * r * math.cos(t0))
    c2 = (p3[0] + k * r * math.sin(t1), p3[1] - k * r * math.cos(t1))
    return ("C", p0, c1, c2, p3)


def wrap(a: float) -> float:
    return (a + math.pi) % (2 * math.pi) - math.pi


# ----------------------------------------------------------------------------- parameters/results
@dataclass
class Params:
    orientation: str = "TB"          # TB: main flow top->bottom (portrait), LR: left->right (landscape)
    cycle_direction: str = "cw"      # cw / ccw (as seen on screen)
    branch_side: str = "right"       # side of the backbone used for branches (side compounds go opposite)
    spacing: float = 110.0           # metabolite-to-metabolite distance on straight segments
    cycle_spacing: float = 125.0     # arc length per metabolite on cycles
    min_cycle_radius: float = 90.0
    max_cycle_length: int = 14
    met_radius: float = 11.0
    clone_radius: float = 7.0
    rxn_half: float = 6.0
    side_offset: float = 28.0        # distance of side compounds from the reaction axis (short hooks)
    side_along: float = 16.0         # offset along the axis (upstream for substrates, downstream for products)
    side_step: float = 19.0          # extra offset for 2nd, 3rd ... side compound on the same side
    margin: float = 40.0
    reaction_labels: bool = True
    hub_degree: int = 4              # hub demotion (see model.classify)
    hub_ratio: float = 2.0
    fold: int = 0                    # max metabolites per straight backbone run (0 = auto)
    compact: bool = True             # remove empty rows/columns (down to the base spacing) at the end
    modifiers: bool = False          # draw catalysts/regulators (SBML modifiers); off by default
    pseudo: str = "drain"            # biomass/pseudo-reactions: "drain" (arrow to a small symbol next to
                                     # each precursor), "hide", or "draw" (a normal reaction node)
    pseudo_min_partners: int = 8     # a reaction with this many main partners counts as pseudo-reaction
    base_spacing: float = 0.0        # spacing before automatic enlarging (set internally)
    main_path: list = field(default_factory=list)  # ids/names of the main line; [] = canonical table; ["off"]
    auto_fold: int = 12              # auto: fold backbones longer than this into serpentine columns
    fold_gap: float = 0.0            # distance between folded columns (x spacing; 0 = as wide as labels need)
    parallel_step: float = 30.0      # spacing of parallel reactions stacked inside a circle (px)
    length_penalty: float = 3.0      # cost per extra spacing of a branch edge (same as one conflict)
    max_link: float = 0.0            # (off by default) a reaction->metabolite line longer than this (x spacing) is replaced by a
                                     # copy ONLY if the map stays connected without it (0 = never)
    split_chords: bool = False       # (off by default) draw lines through a circle as copies instead
    align_tol: float = 0.3           # metabolites closer than this (x spacing) in x or y are aligned
    grid: float = 0.5                # grid cell (x spacing) for metabolites; reactions use half of it (0 = off)
    density: float = 1.0             # >1: more compact (less white space), <1: more air
    sa_backbone: bool = False        # may the annealing step move main-line metabolites?
    sa_max_step: float = 1.0         # largest annealing move (x spacing)
    length_weight: float = 3.0       # cost per spacing of line length in the crossing search
    partner_pull: float = 1.0        # cost per spacing of distance to a node's other, already placed partners
    cycles: str = "known"            # circles for: "known" curated cycles (TCA, urea, ...), "all" loops, "none"
    label_max: int = 34              # longer names are shortened in labels
    label_mode: str = "name"         # "name" or "id"


@dataclass
class Clone:
    glyph: str
    species: str
    rxn: str
    pos: Pt
    role: str     # sidesubstrate / sideproduct
    far: bool = False


@dataclass
class Edge:
    id: str
    rxn: str
    species: str
    glyph: str
    role: str     # substrate / product / sidesubstrate / sideproduct
    segs: list
    weight: float
    reversible: bool
    drawn: str = ""     # role as drawn (product = arrow end); differs from role for flipped reactions
    far: bool = False   # main partner drawn as clone because it is far away


@dataclass
class Label:
    glyph: str        # text glyph id
    target: str       # graphical object id
    species: str | None
    rxn: str | None
    text: str
    x: float          # top-left
    y: float
    w: float
    h: float
    size: float
    kind: str = "main"      # main (bold), copy (copy of a main metabolite, bold), side (cofactor, light), reaction
    anchor: str = "start"   # text alignment within the box: start / middle / end (towards the node)
    leader: bool = False    # placed away from its node: draw a leader line
    node: Pt = (0.0, 0.0)   # centre of the labelled node


@dataclass
class Cycle:
    center: Pt
    radius: float
    direction: int
    slots: list


@dataclass
class Layout:
    pathway: Pathway
    params: Params
    roles: dict[str, Roles]
    met_pos: dict[str, Pt] = field(default_factory=dict)
    rxn_pos: dict[str, Pt] = field(default_factory=dict)
    rxn_axis: dict[str, Pt] = field(default_factory=dict)
    rxn_side: dict[str, Pt] = field(default_factory=dict)
    clones: list[Clone] = field(default_factory=list)
    edges: list[Edge] = field(default_factory=list)
    labels: dict[str, Label] = field(default_factory=dict)
    cycles: list[Cycle] = field(default_factory=list)
    ring_met: dict[str, tuple[int, float]] = field(default_factory=dict)
    ring_rxn: dict[str, tuple[int, float]] = field(default_factory=dict)
    far_links: dict[str, set] = field(default_factory=dict)   # rid -> main species drawn as clone (too far away)
    flipped: set = field(default_factory=set)                  # reversible reactions drawn against written direction
    backbone: list[str] = field(default_factory=list)
    width: float = 0.0
    height: float = 0.0
    meta: dict = field(default_factory=dict)

    # glyph id helpers -------------------------------------------------------
    @staticmethod
    def sg(mid: str, rid: str | None = None) -> str:
        return f"SG_M_{mid}" if rid is None else f"SG_M_{mid}__R_{rid}"

    @staticmethod
    def rg(rid: str) -> str:
        return f"RG_R_{rid}"


@dataclass
class Previous:
    """Positions taken from an earlier layout (or user edits); all are honoured exactly."""
    met: dict[str, Pt] = field(default_factory=dict)
    rxn: dict[str, Pt] = field(default_factory=dict)
    labels: dict[str, Pt] = field(default_factory=dict)   # text glyph id -> top-left
    clones: dict[str, Pt] = field(default_factory=dict)   # clone glyph id -> center (keeps their side)
    layout_hash: str | None = None


# ----------------------------------------------------------------------------- helpers
def canonical_orientation(pw: Pathway, roles: dict[str, Roles]) -> dict[str, int]:
    """Find the canonical cycles (known.CANONICAL_CYCLES) regardless of the written direction of
    reversible reactions, and return the drawing direction of the reactions on them:
    rid -> +1 (as written) or -1 (reversed), so the cycle runs in canonical order
    (e.g. oxaloacetate -> citrate -> ... -> malate)."""
    from .known import CANONICAL_CYCLES, match_members
    edges: dict[tuple, list] = {}
    adj: dict[str, set] = {}
    for rid in sorted(roles):
        ro, rev = roles[rid], pw.reactions[rid].reversible
        for s_ in ro.main_sub:
            for p_ in ro.main_prod:
                if s_ == p_:
                    continue
                edges.setdefault((s_, p_), []).append((rid, 1))
                if rev:
                    edges.setdefault((p_, s_), []).append((rid, -1))
                adj.setdefault(s_, set()).add(p_)
                adj.setdefault(p_, set()).add(s_)
    forced: dict[str, int] = {}
    used: set[str] = set()
    for name, members in CANONICAL_CYCLES.items():
        member_of = match_members(pw.metabolites, members)
        need = {b for b, _, _ in members}
        if set(member_of.values()) != need:
            continue
        order = {b: i for i, (b, _, _) in enumerate(members)}
        mem = set(member_of)
        # members plus intermediates between two members (e.g. cis-aconitate)
        nodes = (mem | {n for n in adj if len(adj[n] & mem) >= 2}) - used
        Gb = nx.DiGraph()
        Gb.add_edges_from(e for e in sorted(edges) if e[0] in nodes and e[1] in nodes)
        best = None
        for c in nx.simple_cycles(Gb, length_bound=len(members) + 4):
            idx = [order[member_of[m]] for m in c if m in member_of]
            if {member_of[m] for m in c if m in member_of} != need:
                continue
            n = len(members)
            forward = sum(1 for i in range(len(idx)) if (idx[(i + 1) % len(idx)] - idx[i]) % n == 1)
            key = (-forward, len(c) - len(idx), tuple(sorted(c)))
            if best is None or key < best[0]:
                best = (key, c)
        if best is None:
            continue
        c = best[1]
        for i, u in enumerate(c):
            v = c[(i + 1) % len(c)]
            for rid, d in edges.get((u, v), []):
                forced.setdefault(rid, d)
        used |= set(c)
    return forced


def path_orientation(pw: Pathway, roles: dict[str, Roles], user_spec=None) -> dict[str, int]:
    """Reversible reactions along a known main pathway (known.CANONICAL_PATHS, e.g. glycolysis) or a
    user-given --main-path are drawn in pathway order (D-glucose -> pyruvate), whatever direction
    the file writes them in (BRENDA writes e.g. G6P isomerase as F6P -> G6P, pyruvate kinase as
    pyruvate -> PEP). rid -> +1 (as written) / -1 (reversed)."""
    from .known import CANONICAL_PATHS, MIN_PATH_MEMBERS, match_members, user_path_members
    U = nx.Graph()
    for rid in sorted(roles):
        ro = roles[rid]
        for s_ in ro.main_sub:
            for p_ in ro.main_prod:
                if s_ != p_:
                    U.add_edge(s_, p_)
                    U[s_][p_].setdefault("rxns", []).append(rid)
    tables = {"main path": user_path_members(user_spec)} if user_spec else CANONICAL_PATHS
    need = 2 if user_spec else MIN_PATH_MEMBERS
    forced: dict[str, int] = {}
    for name, members in tables.items():
        member_of = match_members(pw.metabolites, members)
        cands = [[m for m, b in sorted(member_of.items()) if b == bb and m in U] for bb, _, _ in members]
        cands = [c for c in cands if c]
        if len(cands) < need:
            continue
        rank: dict[str, int] = {}
        if user_spec:
            # only some stations given: the steps in between follow the shortest route
            k = 0
            for a_, b_ in zip(cands, cands[1:]):
                best = None
                for u in a_:
                    for v in b_:
                        if nx.has_path(U, u, v):
                            pth = nx.shortest_path(U, u, v)
                            if best is None or len(pth) < len(best):
                                best = pth
                for x in best or []:
                    if x not in rank:
                        rank[x] = k
                        k += 1
        else:
            # known pathway: only a run of >= MIN_PATH_MEMBERS consecutive members that are directly
            # linked counts (a pentose phosphate map with glucose, G6P, F6P and GAP is no glycolysis);
            # the same compound in several compartments shares one rank (transport is not oriented)
            linked = [any(U.has_edge(u, v) for u in a_ for v in b_) for a_, b_ in zip(cands, cands[1:])]
            runs, start = [], 0
            for i_, ok in enumerate(linked + [False]):
                if not ok:
                    runs.append((start, i_))
                    start = i_ + 1
            lo, hi = max(runs, key=lambda r: (r[1] - r[0], -r[0]))
            if hi - lo + 1 < need:
                continue
            for i_ in range(lo, hi + 1):
                for m in cands[i_]:
                    rank[m] = i_
        if len(set(rank.values())) < need:
            continue
        for rid in sorted(roles):
            if not pw.reactions[rid].reversible:
                continue
            ro = roles[rid]
            si = [rank[m] for m in ro.main_sub if m in rank]
            pi = [rank[m] for m in ro.main_prod if m in rank]
            if si and pi and min(si) != min(pi):
                forced.setdefault(rid, 1 if min(si) < min(pi) else -1)
    return forced


def orient_reversible(pw: Pathway, roles: dict[str, Roles], forced: dict[str, int] | None = None) -> set[str]:
    """Reversible reactions are written in arbitrary direction (e.g. PGK as 3PG -> 1,3BPG).
    Choose their drawing direction by flow conservation:
      1. chains of reversible reactions through unbranched metabolites are oriented consistently
         (each inner metabolite is produced by one step and consumed by the next)
      2. whole chains are flipped while that reduces the number of dead ends, i.e. metabolites
         with >= 2 reactions that are only produced or only consumed
    Ties keep the written direction."""
    inc: dict[str, list[str]] = {}
    for rid in sorted(roles):
        ro = roles[rid]
        for m in dict.fromkeys(ro.main_sub + ro.main_prod):
            inc.setdefault(m, []).append(rid)
    forced = forced or {}
    rev = {rid for rid in roles if pw.reactions[rid].reversible and roles[rid].main_sub and roles[rid].main_prod
           and rid not in forced}
    o = {rid: forced.get(rid, 1) for rid in roles}     # canonical cycle reactions are fixed anchors

    def produces(rid, m):
        ro = roles[rid]
        return (m in ro.main_prod) if o[rid] == 1 else (m in ro.main_sub)

    parent = {r: r for r in rev}

    def find(r):
        while parent[r] != r:
            parent[r] = parent[parent[r]]
            r = parent[r]
        return r
    links = []
    for m in sorted(inc):
        rids = inc[m]
        if len(rids) == 2 and all(r in rev for r in rids):
            parent[find(rids[0])] = find(rids[1])
            links.append((m, rids[0], rids[1]))
    units: dict[str, list[str]] = {}
    for r in sorted(rev):
        units.setdefault(find(r), []).append(r)
    # 1. consistent orientation inside each chain (BFS from its first reaction, as written)
    adj: dict[str, list[tuple[str, str]]] = {}
    for m, r1, r2 in links:
        adj.setdefault(r1, []).append((m, r2))
        adj.setdefault(r2, []).append((m, r1))
    for members in units.values():
        seen = {members[0]}
        queue = [members[0]]
        while queue:
            r = queue.pop(0)
            for m, r2 in sorted(adj.get(r, [])):
                if r2 in seen:
                    continue
                if produces(r2, m) == produces(r, m):   # both produce or both consume m -> flip r2
                    o[r2] = -o[r2]
                seen.add(r2)
                queue.append(r2)

    def cost():
        dead = 0
        for m, rids in inc.items():
            if len(rids) >= 2:
                n_in = sum(produces(r, m) for r in rids)
                if n_in == 0 or n_in == len(rids):
                    dead += 1
        return dead * 1000 + sum(1 for r in rev if o[r] == -1)
    # 2. flip whole chains while it helps
    best = cost()
    improved = True
    while improved:
        improved = False
        for key in sorted(units):
            for r in units[key]:
                o[r] = -o[r]
            c = cost()
            if c < best:
                best, improved = c, True
            else:
                for r in units[key]:
                    o[r] = -o[r]
    return {r for r in roles if o[r] == -1}


def short_text(text: str, n: int) -> str:
    return text if len(text) <= n else text[: n - 1].rstrip(" ,-") + "…"


def display_name(pw: Pathway, mid: str, P: "Params", clone: bool = False) -> str:
    m = pw.metabolites[mid]
    if clone and pw.is_currency(mid):
        text = m.short()
    elif P.label_mode == "id":
        text = m.base or mid
    else:
        text = m.name
    if len(text) > P.label_max:
        text = text[: P.label_max - 1].rstrip(" ,-") + "…"
    return text + (f" [{m.comp_label}]" if m.comp_label else "")


# ----------------------------------------------------------------------------- main entry
KIND_RANK = {"stimulation": 0, "modulation": 1, "catalysis": 2, "inhibition": 3}


def is_regulator_link(pw: Pathway, rid: str, m: str) -> bool:
    """main partner that is really a modifier standing in for an empty side (see add_regulators)"""
    r = pw.reactions.get(rid)
    return bool(r and m in r.modifiers)


def add_regulators(pw: Pathway, roles: dict[str, Roles]):
    """Synthesis / degradation with an empty side ("-> VSP" switched on by MYC2, "JAZ ->" by SCF).

    The empty side gets the SBGN source/sink symbol (placed like a cofactor). If the protein made
    or degraded takes part in nothing else, its regulator stands in for the empty side, so all genes
    a transcription factor switches on fan out around it instead of floating as islands; the link
    is drawn as a modulation arc. A protein that is also used elsewhere stays there, and its
    regulator gets a line or a clone next to the reaction."""
    one_sided = [rid for rid in roles if not pw.reactions[rid].substrates or not pw.reactions[rid].products]
    if not one_sided:
        return
    if EMPTY not in pw.metabolites:
        pw.metabolites[EMPTY] = Metabolite(EMPTY, "\u2205", "c", EMPTY)
    uses: dict[str, int] = {}
    for rid, ro in roles.items():
        for m in set(ro.main_sub + ro.main_prod + ro.side_sub + ro.side_prod):
            uses[m] = uses.get(m, 0) + 1
    main_somewhere = {m for ro in roles.values() for m in ro.main_sub + ro.main_prod}
    for rid in sorted(one_sided):
        r, ro = pw.reactions[rid], roles[rid]
        sub_empty = not ro.main_sub and not ro.side_sub
        prod_empty = not ro.main_prod and not ro.side_prod
        present = ro.main_prod if sub_empty else ro.main_sub
        reg = None
        if r.modifiers and present and all(uses.get(m, 0) == 1 for m in present):
            reg = min(r.modifiers, key=lambda m: (m not in main_somewhere, KIND_RANK.get(r.modifiers[m], 9), m))
        if sub_empty:
            roles[rid] = Roles([reg] if reg else [], ro.main_prod, [EMPTY], ro.side_prod)
        elif prod_empty:
            roles[rid] = Roles(ro.main_sub, [reg] if reg else [], ro.side_sub, [EMPTY])


def compute_layout(pw: Pathway, P: Params | None = None, previous: Previous | None = None) -> Layout:
    """Classify, orient, split into connected components, lay out each one, pack them, finalize."""
    P = P or Params()
    prev = previous or Previous()
    orig_pw = pw          # written to SBML as given (modifiers included, even when not drawn)
    if not P.modifiers and any(r.modifiers for r in pw.reactions.values()):
        pw = replace(pw, reactions={rid: replace(r, modifiers={}) for rid, r in pw.reactions.items()},
                     metabolites=dict(pw.metabolites))
    if P.density != 1.0:      # density scales every distance; node and font sizes stay the same
        k_ = 1.0 / max(0.3, P.density)
        P = replace(P, spacing=P.spacing * k_, cycle_spacing=P.cycle_spacing * k_,
                    min_cycle_radius=P.min_cycle_radius * k_, parallel_step=P.parallel_step * k_,
                    density=1.0)
        P_meta = replace(P, density=1.0 / k_)
    else:
        P_meta = P
    roles = classify(pw, P.hub_degree, P.hub_ratio)          # hub degrees from the whole pathway
    # pseudo-reactions (biomass objective): not drawn as a node with lines across the whole map
    drains: dict[str, tuple] = {}
    if P.pseudo != "draw":
        pseudo = sorted(rid for rid, r in pw.reactions.items()
                        if r.pseudo or looks_pseudo(rid, r.name)
                        or len(roles[rid].main_sub) + len(roles[rid].main_prod) >= P.pseudo_min_partners)
        if pseudo and len(pseudo) < len(pw.reactions):
            if P.pseudo == "drain":
                drains = {rid: (list(roles[rid].main_sub), list(roles[rid].main_prod)) for rid in pseudo}
            rest = replace(pw, reactions={k: v for k, v in pw.reactions.items() if k not in pseudo})
            roles = classify(rest, P.hub_degree, P.hub_ratio)    # degrees without the pseudo-reaction
    forced = canonical_orientation(pw, roles) if P.cycles == "known" else {}
    if P.main_path != ["off"]:
        for rid, d in path_orientation(pw, roles, P.main_path or None).items():
            forced.setdefault(rid, d)
    flipped = orient_reversible(pw, roles, forced)
    for rid in flipped:   # draw reversible reactions in the direction of the overall flow
        ro = roles[rid]
        roles[rid] = Roles(ro.main_prod, ro.main_sub, ro.side_prod, ro.side_sub)
    add_regulators(pw, roles)

    # connected components of the main-compound network
    parent = {r: r for r in roles}

    def find(r):
        while parent[r] != r:
            parent[r] = parent[parent[r]]
            r = parent[r]
        return r
    first_rxn: dict[str, str] = {}
    for rid in sorted(roles):
        for m in roles[rid].main_sub + roles[rid].main_prod:
            if m in first_rxn:
                parent[find(rid)] = find(first_rxn[m])
            else:
                first_rxn[m] = rid
    comps: dict[str, list[str]] = {}
    for rid in sorted(roles):
        comps.setdefault(find(rid), []).append(rid)
    comp_list = sorted(comps.values(), key=lambda c: (-len(c), c[0]))

    P = replace(P, base_spacing=P.spacing)
    pinned_ids = set(prev.met) | set(prev.rxn) | set(pw.pins)
    participant = {m for ro in roles.values() for m in ro.main_sub + ro.main_prod}
    comp_of = {rid: i for i, rids in enumerate(comp_list) for rid in rids}
    owner: dict[str, int] = {}
    for rid in sorted(roles):
        for m in sorted(pw.reactions[rid].modifiers):
            if m not in participant:
                owner.setdefault(m, comp_of[rid])
    parts = []
    for ci_, rids in enumerate(comp_list):
        here = {m for r in rids for m in roles[r].main_sub + roles[r].main_prod}
        foreign = {m for r in rids for m in pw.reactions[r].modifiers
                   if (m in participant and m not in here) or (m not in participant and owner.get(m) != ci_)}
        mets = {m for r in rids for m in list(pw.reactions[r].substrates) + list(pw.reactions[r].products)
                + list(pw.reactions[r].modifiers)}
        if EMPTY in pw.metabolites:
            mets.add(EMPTY)
        sub_pw = replace(pw, reactions={r: pw.reactions[r] for r in rids},
                         metabolites={m: pw.metabolites[m] for m in sorted(mets)})
        sub_roles = lambda: {r: Roles(list(roles[r].main_sub), list(roles[r].main_prod),
                                      list(roles[r].side_sub), list(roles[r].side_prod)) for r in rids}
        Lc = _layout_component(sub_pw, P, prev, sub_roles(), flipped & set(rids), foreign, drains)
        pinned = bool(pinned_ids & (set(Lc.met_pos) | set(Lc.rxn_pos)))
        # too dense (many labels had to move away from their node): give it more room
        Pc = P
        for _ in range(2):
            leaders = sum(lb.leader for lb in Lc.labels.values())
            if pinned or not Lc.labels or leaders <= 0.08 * max(1.0, P_meta.density) * len(Lc.labels):
                break
            Pc = replace(Pc, spacing=Pc.spacing * 1.25, cycle_spacing=Pc.cycle_spacing * 1.25,
                         parallel_step=Pc.parallel_step * 1.25, grid=Pc.grid / 1.25)   # same grid cell
            Lc = _layout_component(sub_pw, Pc, prev, sub_roles(), flipped & set(rids), foreign, drains)
        parts.append((Lc, pinned))

    L = Layout(orig_pw, P, roles)
    L.flipped = flipped
    g_abs = P.grid * P.spacing
    _pack(L, parts, P, g_abs)
    finalize(L, keep_frame=bool(previous and (previous.met or previous.rxn)), parent=prev.layout_hash, grid=g_abs)
    return L


def _merge(L: Layout, Lc: Layout):
    off = len(L.cycles)
    L.met_pos.update(Lc.met_pos)
    L.rxn_pos.update(Lc.rxn_pos)
    L.rxn_axis.update(Lc.rxn_axis)
    L.rxn_side.update(Lc.rxn_side)
    L.clones += Lc.clones
    L.edges += Lc.edges
    L.labels.update(Lc.labels)
    L.cycles += Lc.cycles
    L.ring_met.update({k: (i + off, t) for k, (i, t) in Lc.ring_met.items()})
    L.ring_rxn.update({k: (i + off, t) for k, (i, t) in Lc.ring_rxn.items()})
    L.backbone += Lc.backbone if not L.backbone else []
    L.far_links.update(Lc.far_links)


def _pack(L: Layout, parts, P: Params, g: float = 0.0):
    """Largest component first; parts sharing a metabolite with it go next to that metabolite, the
    others are packed around it as tightly as possible. Components with pinned nodes keep their
    coordinates."""
    gap = 0.8 * P.spacing

    def gshift(v: float) -> float:    # shifts are whole grid cells, so all parts share one grid
        return float(math.ceil(v / g) * g) if g > 0 else float(math.ceil(v))
    fixed = [Lc for Lc, pinned in parts if pinned]
    free = [Lc for Lc, pinned in parts if not pinned]
    for Lc in fixed:
        _merge(L, Lc)
    if not free:
        return
    if not fixed:
        main = free.pop(0)
        x0, y0, _, _ = _bbox(main)
        _transform(main, (gshift(-x0) - g if g > 0 else -math.floor(x0), gshift(-y0) - g if g > 0 else -math.floor(y0)))
        _merge(L, main)
    placed = [_bbox(L)]                               # boxes of the parts drawn so far
    # 1. components that share a metabolite with what is already placed go next to that metabolite
    rest = []
    for Lc in free:
        t = _anchor_translation(L, Lc, P, g)
        if t is None:
            rest.append(Lc)
        else:
            _transform(Lc, t)
            placed.append(_bbox(Lc))
            _merge(L, Lc)
    free = rest
    if not free:
        return
    # 2. the remaining ones are packed around what is already drawn: largest first, each at the spot
    # (beside / above / below an already placed part, flush with its edges) that keeps the whole
    # drawing smallest, with a mild preference for a landscape-ish shape. Parts never overlap.
    gap = 0.6 * P.spacing

    def _area(Lc):
        x0, y0, x1, y1 = _bbox(Lc)
        return (x1 - x0) * (y1 - y0)
    for Lc in sorted(free, key=lambda Lc: (-_area(Lc), _bbox(Lc))):
        x0, y0, x1, y1 = _bbox(Lc)
        w, h = x1 - x0, y1 - y0
        # candidate corners: beside / flush with the edges of any placed part, combined across parts
        # (so a part can sit in the notch between two others)
        xs_ = sorted({v for (px0, py0, px1, py1) in placed for v in (px1 + gap, px0 - gap - w, px0, px1 - w)})
        ys_ = sorted({v for (px0, py0, px1, py1) in placed for v in (py1 + gap, py0 - gap - h, py0, py1 - h)})
        cands = [(cx, cy) for cx in xs_ for cy in ys_]
        ux0 = min(r[0] for r in placed)
        uy0 = min(r[1] for r in placed)
        ux1 = max(r[2] for r in placed)
        uy1 = max(r[3] for r in placed)
        best = None
        for cx, cy in cands:
            t = (gshift(cx - x0), gshift(cy - y0))
            r = (x0 + t[0], y0 + t[1], x1 + t[0], y1 + t[1])
            if any(r[0] < q[2] + gap * 0.5 and q[0] < r[2] + gap * 0.5 and
                   r[1] < q[3] + gap * 0.5 and q[1] < r[3] + gap * 0.5 for q in placed):
                continue
            W_ = max(ux1, r[2]) - min(ux0, r[0])
            H_ = max(uy1, r[3]) - min(uy0, r[1])
            score = W_ * H_ * (1.0 + 0.25 * abs(math.log(max(W_, 1.0) / max(H_, 1.0) / 1.4)))
            key = (round(score), abs(r[0] - ux0) + abs(r[1] - uy0))
            if best is None or key < best[0]:
                best = (key, t, r)
        if best is None:                               # cannot happen (outside candidates are free)
            t = (gshift(ux1 + gap - x0), gshift(uy0 - y0))
            r = (x0 + t[0], y0 + t[1], x1 + t[0], y1 + t[1])
        else:
            t, r = best[1], best[2]
        _transform(Lc, t)
        _merge(L, Lc)
        placed.append(r)


def _anchor_translation(L: Layout, Lc: Layout, P: Params, g: float = 0.0):
    """Find an integer shift that puts component Lc close to a metabolite it shares with L
    (e.g. a module whose link was drawn as a copy), without overlapping anything already drawn."""
    D = P.spacing
    pairs = []
    for c in Lc.clones:                         # shared species is a copy in Lc, a node in L
        if c.species in L.met_pos:
            pairs.append((c.species, L.met_pos[c.species], c.pos))
    for m, p in Lc.met_pos.items():             # shared species is a node in Lc, a copy in L
        for c in L.clones:
            if c.species == m:
                pairs.append((m, c.pos, p))
    if not pairs:
        return None
    occupied = list(L.met_pos.values()) + list(L.rxn_pos.values()) + [c.pos for c in L.clones]
    boxes = [(lb.x, lb.y, lb.x + lb.w, lb.y + lb.h) for lb in L.labels.values()]
    x0, y0, x1, y1 = _bbox(Lc)
    pad = 0.4 * D
    for _, A, B in sorted(pairs, key=lambda x: x[0]):
        for k in (1.5, 2.5, 3.5, 5.0, 7.0):
            for d in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                t = (round(A[0] + d[0] * k * D - B[0]), round(A[1] + d[1] * k * D - B[1]))
                if g > 0:
                    t = (round(t[0] / g) * g, round(t[1] / g) * g)
                bx0, by0, bx1, by1 = x0 + t[0] - pad, y0 + t[1] - pad, x1 + t[0] + pad, y1 + t[1] + pad
                if any(bx0 < q[0] < bx1 and by0 < q[1] < by1 for q in occupied):
                    continue
                if any(not (b[2] < bx0 or bx1 < b[0] or b[3] < by0 or by1 < b[1]) for b in boxes):
                    continue
                return (float(t[0]), float(t[1]))
    return None


def redox_ladder(rids: list[str], roles: dict[str, Roles], is_water=lambda m: False):
    """Chain of 2:2 reactions in which neighbours share one metabolite pair with swapped roles
    (X made by one, used by the next; Y used by one, made by the next). Returns (chain, cols) with
    cols: metabolite -> (column, top row?) - reactions sit on odd columns, pairs on even ones - or
    None if the reactions do not form one such chain."""
    quad = [r for r in rids if len(roles[r].main_sub) == 2 and len(roles[r].main_prod) == 2]
    if len(quad) < 2 or len(quad) != len(rids):
        return None
    nb: dict[str, list] = {r: [] for r in quad}
    for i, a in enumerate(quad):
        A_s, A_p = set(roles[a].main_sub), set(roles[a].main_prod)
        for b in quad[i + 1:]:
            B_s, B_p = set(roles[b].main_sub), set(roles[b].main_prod)
            x, y = A_p & B_s, A_s & B_p
            if len(x) == 1 and len(y) == 1 and not (A_s & B_s) and not (A_p & B_p):
                nb[a].append(b)
                nb[b].append(a)
    if any(len(v) > 2 for v in nb.values()) or sum(len(v) for v in nb.values()) != 2 * (len(quad) - 1):
        return None
    ends = sorted(r for r in quad if len(nb[r]) == 1)
    if len(ends) != 2:
        return None

    def free(r):                        # the pair of an end reaction not shared with its neighbour
        other = set(roles[nb[r][0]].main_sub + roles[nb[r][0]].main_prod)
        return ([m for m in roles[r].main_sub if m not in other], [m for m in roles[r].main_prod if m not in other])
    # reading direction: start where water is consumed (photosynthesis: H2O -> ... -> NADPH),
    # end where it is formed (respiratory chain: NADH -> ... -> H2O)
    ends.sort(key=lambda r: (-sum(1 for m in free(r)[0] if is_water(m)) + sum(1 for m in free(r)[1] if is_water(m)), r))
    chain = [ends[0]]
    while len(chain) < len(quad):
        nxt = [r for r in nb[chain[-1]] if r not in chain]
        if not nxt:
            return None
        chain.append(nxt[0])
    cols: dict[str, tuple] = {}
    prod_up = True                      # first reaction: substrates below, products above
    for k, r in enumerate(chain):
        left, right = 2 * k, 2 * k + 2
        prev_r = chain[k - 1] if k else None
        next_r = chain[k + 1] if k + 1 < len(chain) else None
        shared_next = set(roles[next_r].main_sub + roles[next_r].main_prod) if next_r else set()
        for m in roles[r].main_sub + roles[r].main_prod:
            up = (m in roles[r].main_prod) == prod_up
            shared_prev = prev_r is not None and m in set(roles[prev_r].main_sub + roles[prev_r].main_prod)
            if m in shared_next:
                col = right
            elif shared_prev:
                col = left
            else:                       # pair of the first (left end) or last reaction (right end)
                col = left if prev_r is None else right
            if m in cols and cols[m] != (col, up):
                return None                 # inconsistent: not a simple ladder
            cols[m] = (col, up)
        prod_up = not prod_up           # the next reaction uses this one's products: rows swap
    # two metabolites per column, one above and one below
    seen: dict[tuple, str] = {}
    for m, cu in cols.items():
        if cu in seen:
            return None
        seen[cu] = m
    return chain, cols


def _layout_component(pw: Pathway, P: Params, prev: Previous, roles: dict[str, Roles], flipped: set,
                      foreign: set | frozenset = frozenset(), drains: dict | None = None) -> Layout:
    L = Layout(pw, P, roles)
    L.flipped = flipped
    W = {rid: r.weight for rid, r in pw.reactions.items()}

    F: Pt = (0.0, 1.0) if P.orientation == "TB" else (1.0, 0.0)
    sigma = 1 if P.branch_side == "right" else -1
    S: Pt = mul((F[1], -F[0]), -sigma if P.orientation == "LR" else sigma)

    def side_rot(a: Pt) -> Pt:
        # rotation that maps the flow direction F onto -S (side compounds opposite to branches)
        r = (-a[1], a[0])
        return r if dot((-F[1], F[0]), neg(S)) > 0 else neg(r)

    dirc = 1 if P.cycle_direction == "cw" else -1
    D = P.spacing
    if P.orientation == "LR":  # horizontal: labels sit below nodes, so spacing must fit the widest name
        D = max(D, max((text_size(display_name(pw, m, P), 12)[0] for ro in roles.values()
                        for m in ro.main_sub + ro.main_prod), default=0) + 16)

    # boxes of labels fixed by a previous layout act as obstacles for new nodes
    fixed_boxes = []
    if prev.labels:
        for tg, (x, y) in prev.labels.items():
            if tg.startswith("TG_R_"):
                if not P.reaction_labels or tg[5:] not in pw.reactions:
                    continue
                w, h = text_size(short_text(pw.reactions[tg[5:]].name, P.label_max), 9)
            else:
                mid, _, rid = tg[5:].partition("__R_")
                if mid not in pw.metabolites:
                    continue
                w, h = text_size(display_name(pw, mid, P, clone=bool(rid)), 10 if rid else 12)
            fixed_boxes.append((x - 2, y - 2, w + 4, h + 4))

    def hits_fixed(p: Pt, pad=0.0) -> bool:
        return any(b[0] - pad < p[0] < b[0] + b[2] + pad and b[1] - pad < p[1] < b[1] + b[3] + pad
                   for b in fixed_boxes)

    # pins: explicit pins in input file override previous layout
    pin_met = dict(prev.met)
    pin_rxn = dict(prev.rxn)
    for k, v in pw.pins.items():
        if k in pw.reactions:
            pin_rxn[k] = tuple(v)
        else:
            pin_met[k] = tuple(v)

    # ---------------------------------------------------------------- skeleton
    G = nx.DiGraph()
    main_mets = sorted({m for ro in roles.values() for m in ro.main_sub + ro.main_prod})
    G.add_nodes_from(main_mets)
    for rid, ro in roles.items():
        for s in ro.main_sub:
            for p in ro.main_prod:
                if s == p:
                    continue
                if G.has_edge(s, p):
                    G[s][p]["rxns"].append(rid)
                    G[s][p]["weight"] = max(G[s][p]["weight"], W[rid])
                else:
                    G.add_edge(s, p, rxns=[rid], weight=W[rid])

    def best_rxn(rids):
        return sorted(rids, key=lambda r: (-W[r], r))[0]

    # ---------------------------------------------------------------- cycles
    raw_cycles = []
    for c in nx.simple_cycles(G, length_bound=P.max_cycle_length):
        if len(c) >= 3:
            i = c.index(min(c))
            raw_cycles.append(tuple(c[i:] + c[:i]))
    raw_cycles = sorted(set(raw_cycles), key=lambda c: (-len(c), c))
    cycles: list[tuple] = []
    cycle_info: list[dict] = []
    if P.cycles == "known":
        for c, info in select_known_cycles(raw_cycles, pw.metabolites):
            cycles.append(c)
            cycle_info.append(info)
    elif P.cycles == "all":
        used: set[str] = set()
        for c in raw_cycles:
            if not used & set(c):
                cycles.append(c)
                cycle_info.append({"name": "", "top": None, "member_of": {}})
                used |= set(c)

    sup = {m: m for m in main_mets}
    for k, c in enumerate(cycles):
        for m in c:
            sup[m] = f"@cyc{k}"
    H = nx.DiGraph()
    H.add_nodes_from(sorted(set(sup.values())))
    for u, v, d in sorted(G.edges(data=True), key=lambda e: (e[0], e[1])):
        su, sv = sup[u], sup[v]
        if su == sv:
            continue
        if H.has_edge(su, sv):
            H[su][sv]["weight"] = max(H[su][sv]["weight"], d["weight"])
        else:
            H.add_edge(su, sv, weight=d["weight"])
    while True:  # break remaining (short) cycles deterministically
        try:
            cyc = nx.find_cycle(H)
        except nx.NetworkXNoCycle:
            break
        u, v = max((e[0], e[1]) for e in cyc)
        H.remove_edge(u, v)

    # ---------------------------------------------------------------- backbone
    path: list[str] = []
    if H.number_of_nodes():
        topo = list(nx.lexicographical_topological_sort(H))
        path = nx.dag_longest_path(H, weight="weight", default_weight=1, topo_order=topo)
        # prefer the shortest route between any two backbone nodes over detours
        # (aldolase F16BP->GAP instead of via DHAP; G6P->...->GAP via glycolysis, not via the pentose phosphate pathway)
        changed = True
        while changed:
            changed = False
            for i in range(len(path)):
                lengths = nx.single_source_shortest_path_length(H, path[i])
                for j in range(len(path) - 1, i + 1, -1):
                    if lengths.get(path[j], j - i) < j - i:
                        short = nx.shortest_path(H, path[i], path[j])
                        path = path[:i] + short + path[j + 1:]
                        changed = True
                        break
                if changed:
                    break
    # a canonical main path (e.g. glycolysis D-glucose -> pyruvate) or a user-given one replaces the
    # longest chain; it is continued forward along the longest chain after its last member
    L.meta["main_path"] = None
    if H.number_of_nodes() and P.main_path != ["off"]:
        found = find_main_path(H, sup, pw.metabolites, P.main_path or None)
        if found:
            name_, chain = found
            # same compound just before the start in another compartment (e.g. glc__D_e -> glc__D)
            while True:
                pre = [u for u in H.predecessors(chain[0]) if u not in chain and not u.startswith("@cyc")
                       and pw.metabolites[u].base == pw.metabolites[chain[0]].base]
                if not pre:
                    break
                chain.insert(0, sorted(pre)[0])
            used_ = set(chain)
            best_tail: dict[str, list] = {}
            for n_ in reversed(list(nx.topological_sort(H))):
                if n_ in used_ and n_ != chain[-1]:
                    continue
                opts = [best_tail[v] for v in H.successors(n_) if v not in used_ and v in best_tail]
                best_tail[n_] = [n_] + (max(opts, key=lambda t: (len(t), t)) if opts else [])
            tail = best_tail.get(chain[-1], [chain[-1]])[1:]
            path = chain + tail
            L.meta["main_path"] = name_

    # a canonical cycle anchors the drawing: start the main line at it, so it keeps its conventional
    # rotation (citrate synthase at 12 o'clock); what leads into it is attached from outside as a branch
    def top_rxns_of(k: int) -> set:
        """reactions of the canonical top step (e.g. ornithine -> citrulline); a member may occur
        more than once in the cycle (mitochondrial and cytosolic ornithine), so all copies count"""
        info = cycle_info[k]
        if not info.get("top"):
            return set()
        a_set = [m for m, b in info["member_of"].items() if b == info["top"][0] and m in cycles[k]]
        b_set = [m for m, b in info["member_of"].items() if b == info["top"][1] and m in cycles[k]]
        return {r for a_ in a_set for b_ in b_set if G.has_edge(a_, b_) for r in G[a_][b_]["rxns"]}

    def enters_at_top(prev: str, k: int) -> bool:
        top_rxns = top_rxns_of(k)
        into = {r for m in cycles[k] if G.has_edge(prev, m) for r in G[prev][m]["rxns"]}
        return bool(into & top_rxns)

    for i_, n_ in enumerate(path):
        if n_.startswith("@cyc") and cycle_info[int(n_[4:])].get("name"):
            if i_ > 0 and not enters_at_top(path[i_ - 1], int(n_[4:])):
                # main path given (glycolysis ...): keep it and stop before the cycle, which is then
                # attached with its own conventional orientation; otherwise start at the cycle
                path = path[:i_] if L.meta.get("main_path") else path[i_:]
            break
    L.backbone = list(path)

    pos, rpos, raxis = L.met_pos, L.rxn_pos, L.rxn_axis
    ctx: dict[str, str] = {}
    attach: dict[str, Pt] = {}
    order: list[str] = []
    ring_side: dict[str, Pt] = {}

    par_side: dict[str, Pt] = {}    # parallel reactions: cofactors on their own outer side
    parent: dict[str, str] = {}      # branch metabolite -> the node it was attached to (for subtree moves)
    lab_pref: dict[str, Pt] = {}     # preferred label direction (rows of reaction partners: away from reaction)
    mid_of: dict[str, tuple] = {}    # reaction drawn halfway between two metabolites (+ offset)

    def place_rxn_between(u, v, rids):
        rids = sorted(rids, key=lambda r: (-W[r], r))
        n = len(rids)
        for j, rid in enumerate(rids):
            if rid in rpos:
                continue
            k_ = j - (n - 1) / 2
            off = mul(S, k_ * max(30.0, 2 * P.side_offset + 4) if n > 1 else 0.0)
            if n > 1 and abs(k_) > 1e-9:
                par_side[rid] = mul(S, 1 if k_ > 0 else -1)
            rpos[rid] = add(lerp(pos[u], pos[v], 0.5), off)
            mid_of[rid] = (u, v, off)
            ax = unit(sub(pos[v], pos[u]))
            raxis[rid] = ax if dot(ax, ax) > 0 else F

    cyc_of = {m: k for k, c in enumerate(cycles) for m in c}
    post_gap = [D]       # distance from a cycle's exit metabolite to the next backbone node
    placed_cycles: set[int] = set()

    def place_cycle(k: int, prev_node: str | None, F: Pt = F, gap: float = 1.0, exit_m: str | None = None):
        placed_cycles.add(k)
        ci = len(L.cycles)      # index of this circle in L.cycles
        cyc = cycles[k]
        n = len(cyc)
        slots = []
        for i, m in enumerate(cyc):
            slots.append(("M", m))
            slots.append(("R", best_rxn(G[m][cyc[(i + 1) % n]]["rxns"])))
        R = max(n * P.cycle_spacing / (2 * math.pi), P.min_cycle_radius)
        # room inside for stacks of parallel reactions
        max_par = max((len(G[cyc[i]][cyc[(i + 1) % n]]["rxns"]) - 1 for i in range(n)), default=0)
        R += P.parallel_step * max_par
        entry_idx, entry_pos, link = 0, (0.0, 0.0), None
        info = cycle_info[k]
        if prev_node is None and info.get("top"):
            # free-standing known cycle: conventional orientation (e.g. citrate synthase at 12 o'clock)
            tops_ = top_rxns_of(k)
            for j, (kind, id_) in enumerate(slots):
                if kind == "R" and id_ in tops_:
                    entry_idx = j
                    break
        if prev_node is not None:
            cands = sorted(((rid, m) for m in cyc for u, v in ((prev_node, m), (m, prev_node))
                            if G.has_edge(u, v) for rid in G[u][v]["rxns"]), key=lambda x: (-W[x[0]], x[0]))
            rid, m = cands[0]
            if ("R", rid) in slots:
                entry_idx = slots.index(("R", rid))
                entry_pos = add(pos[prev_node], mul(F, D * 0.8 * gap))
            else:
                entry_idx = slots.index(("M", m))
                entry_pos = add(pos[prev_node], mul(F, D * gap))
                link = (prev_node, m, rid) if G.has_edge(prev_node, m) else (m, prev_node, rid)
        center = add(entry_pos, mul(F, R))
        th0 = math.atan2(-F[1], -F[0])
        ns = len(slots)
        if exit_m is not None and link is not None and prev_node is not None and \
                link[1 if link[0] == prev_node else 0] == exit_m:
            # enters and leaves through the same metabolite (e.g. G1P for the Leloir loop):
            # that metabolite stays on the straight line, the circle hangs beside it
            def crowd(side):
                c_ = add(entry_pos, mul(side, R))
                return sum(1 for p in list(pos.values()) + list(rpos.values()) if dist(p, c_) < R + 0.6 * D)
            side = min((neg(S), S), key=crowd)
            # longer steps before/after the shared metabolite so the circle clears its neighbours
            entry_pos = add(pos[prev_node], mul(F, max(D, R + 0.45 * D)))
            post_gap[0] = max(D, R + 0.45 * D)
            center = add(entry_pos, mul(side, R))
            entry_idx = slots.index(("M", exit_m))
            th0 = math.atan2(-side[1], -side[0])
        elif exit_m is not None and exit_m in cyc and not info.get("name"):
            # (only for non-canonical circles; canonical ones keep their orientation)
            # the backbone continues below the circle: put the exit metabolite at the bottom
            entry_idx = slots.index(("M", exit_m))
            th0 = math.atan2(F[1], F[0])
        for j, (kind, id_) in enumerate(slots):
            th = wrap(th0 + (j - entry_idx) * 2 * math.pi / ns * dirc)
            p = (center[0] + R * math.cos(th), center[1] + R * math.sin(th))
            out = (math.cos(th), math.sin(th))
            if kind == "M":
                pos[id_] = p
                ctx[id_] = "ring"
                attach[id_] = out
                order.append(id_)
                L.ring_met[id_] = (ci, th)
            else:
                rpos[id_] = p
                raxis[id_] = mul((-math.sin(th), math.cos(th)), dirc)
                L.ring_rxn[id_] = (ci, th)
                ring_side[id_] = out
                # parallel reactions for the same step (isoenzymes, NAD/NADP variants, ...):
                # a small stack just inside the circle at the same angle
                jm = slots[j - 1][1]
                jn = slots[(j + 1) % ns][1]
                par = [r_ for r_ in sorted(G[jm][jn]["rxns"], key=lambda r_: (-W[r_], r_)) if r_ != id_]
                for kk, r_ in enumerate(par, start=1):
                    if r_ in rpos:
                        continue
                    rad = R - P.parallel_step * kk
                    rpos[r_] = (center[0] + rad * math.cos(th), center[1] + rad * math.sin(th))
                    raxis[r_] = raxis[id_]
                    ring_side[r_] = neg(out)
        L.cycles.append(Cycle(center, R, dirc, slots))
        if link:
            place_rxn_between(link[0], link[1], [link[2]])
        if info.get("name") and info.get("top"):
            # canonical cycle: an extra main partner of the top reaction (acetyl-CoA at citrate
            # synthase) sits directly outside it, like in textbooks
            kind_, top_r = slots[entry_idx]
            if kind_ == "R":
                th = L.ring_rxn[top_r][1]
                out = (math.cos(th), math.sin(th))
                extra = [m for m in roles[top_r].main_sub + roles[top_r].main_prod if m not in cyc and m not in pos]
                for j_, m in enumerate(extra):
                    pos[m] = add(rpos[top_r], add(mul(out, D), mul((-out[1], out[0]), (j_ - (len(extra) - 1) / 2) * 0.8 * D)))
                    ctx[m] = "anchor"
                    attach[m] = out
                    order.append(m)
                if extra:
                    # parallel reactions stacked inside the top reaction (ATP citrate synthase under
                    # citrate synthase) would sit on the line from acetyl-CoA: move them sideways,
                    # alternating left/right, so no reaction node lies on another reaction's line
                    stack = [r_ for r_ in rpos if r_ != top_r and r_ not in L.ring_rxn
                             and ring_side.get(r_) == neg(out) and abs(wrap(math.atan2(rpos[r_][1] - center[1],
                                                                                     rpos[r_][0] - center[0]) - th)) < 1e-6]
                    for kk, r_ in enumerate(sorted(stack, key=lambda r_: -dist(rpos[r_], center))):
                        rad = dist(rpos[r_], center)
                        dth = (1 if kk % 2 == 0 else -1) * (kk // 2 + 1) * 1.2 * P.parallel_step / max(rad, 1.0)
                        rpos[r_] = (center[0] + rad * math.cos(th + dth), center[1] + rad * math.sin(th + dth))

    # long backbones are folded into serpentine columns (aspect ratio / readability)
    n_line = sum(1 for n in path if not n.startswith("@cyc"))
    fold = P.fold or (math.ceil(n_line / math.ceil(n_line / P.auto_fold)) if n_line > P.auto_fold + 3 else 0)
    widest = max((text_size(display_name(pw, n, P), 12)[0] for n in path if not n.startswith("@cyc")), default=0)
    fold_gap_px = P.fold_gap * D if P.fold_gap else max(2.0 * D, widest + 2 * P.side_offset + 40)
    cur_F = F
    prev_node = None
    idx = 0
    after_cycle = D
    for i_node, node in enumerate(path):
        if node.startswith("@cyc"):
            k = int(node[4:])
            nxt = path[i_node + 1] if i_node + 1 < len(path) else None
            exit_m = None
            if nxt is not None:
                exits = sorted((-W[r], r, m) for m in cycles[k] if G.has_edge(m, nxt) for r in G[m][nxt]["rxns"])
                exit_m = exits[0][2] if exits else None
            post_gap[0] = D
            place_cycle(k, prev_node, cur_F, exit_m=exit_m)
            prev_node = exit_m
            after_cycle = post_gap[0]
            if exit_m is not None and exit_m in L.ring_met and dot(attach[exit_m], cur_F) < 0.7:
                # exit is not at the bottom (canonical orientation): leave the circle outward,
                # along the nearest horizontal/vertical direction
                exit_dir = attach[exit_m]
                cur_F = max((F, neg(F), S, neg(S)), key=lambda d: dot(d, exit_dir))
            continue
        if node in pos:           # already placed (e.g. acetyl-CoA above citrate synthase)
            prev_node = node
            continue
        if prev_node is None:
            pos[node] = (0.0, 0.0)
        else:
            if fold and idx % fold == 0:
                # U-turn into the next column; the gap is only as wide as labels and cofactors need
                pos[node] = add(pos[prev_node], mul(S, fold_gap_px))
                cur_F = neg(cur_F)
            else:
                pos[node] = add(pos[prev_node], mul(cur_F, after_cycle))
            after_cycle = D
            place_rxn_between(prev_node, node, G[prev_node][node]["rxns"]
                              if G.has_edge(prev_node, node) else [])
        ctx[node] = "backbone"
        attach[node] = cur_F
        order.append(node)
        prev_node = node
        idx += 1

    # ---------------------------------------------------------------- pins / previous layout
    deltas = [sub(pin_met[m], pos[m]) for m in pos if m in pin_met] + \
             [sub(pin_rxn[r], rpos[r]) for r in rpos if r in pin_rxn]
    if deltas:
        d = mean(deltas)
        for m in pos:
            pos[m] = add(pos[m], d)
        for r in rpos:
            rpos[r] = add(rpos[r], d)
        for c in L.cycles:
            c.center = add(c.center, d)
    for m in list(pos):
        if m in pin_met:
            if dist(pos[m], pin_met[m]) > 0.5:
                L.ring_met.pop(m, None)
            pos[m] = pin_met[m]
    for r in list(rpos):
        if r in pin_rxn:
            if dist(rpos[r], pin_rxn[r]) > 0.5:
                L.ring_rxn.pop(r, None)
            rpos[r] = pin_rxn[r]
    if pin_met or pin_rxn:
        # derived geometry follows the pinned positions: if a cycle's pinned members still form a
        # circle, re-derive centre, angles and direction from them; pinned reactions get their axis
        # from their pinned partners
        for ci, cy in enumerate(L.cycles):
            pts = [pos[i] if k == "M" else rpos[i] for k, i in cy.slots]
            c = mean(pts)
            rr = sum(dist(p, c) for p in pts) / len(pts)
            if max(abs(dist(p, c) - rr) for p in pts) > 1.0:
                continue
            ths = [math.atan2(p[1] - c[1], p[0] - c[0]) for p in pts]
            d_ = 1 if wrap(ths[1] - ths[0]) > 0 else -1
            cy.center, cy.radius, cy.direction = c, rr, d_
            for (k, i), th in zip(cy.slots, ths):
                if k == "M":
                    L.ring_met[i] = (ci, th)
                    attach[i] = (math.cos(th), math.sin(th))
                else:
                    L.ring_rxn[i] = (ci, th)
                    raxis[i] = mul((-math.sin(th), math.cos(th)), d_)
                    ring_side[i] = (math.cos(th), math.sin(th))
        for r in rpos:
            if r in pin_rxn and r not in L.ring_rxn:
                ms = [pos[m] for m in roles[r].main_sub if m in pos]
                mp = [pos[m] for m in roles[r].main_prod if m in pos]
                if ms and mp:
                    ax = unit(sub(mean(mp), mean(ms)))
                    if dot(ax, ax):
                        raxis[r] = ax

    # ---------------------------------------------------------------- branches
    def main_segments():
        segs = []
        for rid, rp in rpos.items():
            ro = roles[rid]
            for m in ro.main_sub + ro.main_prod:
                if m in pos:
                    segs.append((rp, pos[m]))
        return segs

    def violations(anchor: Pt, cand: Pt) -> int:
        v = 0
        v += sum(3 for p in pos.values() if dist(cand, p) < 0.7 * D)
        v += sum(3 for p in rpos.values() if dist(cand, p) < 0.45 * D)
        v += sum(5 for c in L.cycles if dist(cand, c.center) < c.radius + 0.4 * D)
        if fixed_boxes:
            v += sum(1 for t in range(2, 9) if hits_fixed(lerp(anchor, cand, t / 8), 2))
        for a, b in main_segments():
            if dist(a, anchor) < 1e-6 or dist(b, anchor) < 1e-6:
                continue
            if segs_cross(anchor, cand, a, b):
                v += 4
            elif seg_point_dist(cand, a, b) < 25:
                v += 2
        v += sum(2 for p in list(pos.values()) + list(rpos.values())
                 if dist(p, anchor) > 1e-6 and seg_point_dist(p, anchor, cand) < 20)
        return v

    def rot90(a: Pt, s: int) -> Pt:
        return (-a[1] * s, a[0] * s)

    def dir_prefs(anchor_id: str) -> list[Pt]:
        c = ctx.get(anchor_id, "branch")
        a = attach.get(anchor_id, F)
        if c == "ring":
            return [a, rot90(a, 1), rot90(a, -1)]
        if c == "backbone":
            return [S, neg(S), F, neg(F)]
        return [a, mul(neg(side_rot(a)), 1), side_rot(a), neg(a)]

    def step_len(anchor_id: str, d: Pt, new_id: str | None) -> float:
        """Distance for a new node: along label direction (perpendicular to labels) the plain spacing,
        across it enough room for both labels."""
        if new_id is None or abs(d[0]) < 0.7:
            return D
        w = (text_size(display_name(pw, anchor_id, P), 12)[0] + text_size(display_name(pw, new_id, P), 12)[0]) / 2
        return max(D, w + 30)

    def partners_of(m: str) -> list[str]:
        out = []
        for ro in roles.values():
            if m in ro.main_sub:
                out += ro.main_prod
            elif m in ro.main_prod:
                out += ro.main_sub
        return out

    def place_from(anchor_id: str, extra_len=(1.0, 1.5, 2.0, 3.0), new_id: str | None = None) -> tuple[Pt, Pt]:
        if ctx.get(anchor_id) == "ring" and extra_len[0] >= 1.0:
            extra_len = (1.6, 2.0, 2.6, 3.2)   # leave room for the ring's outward side compounds
        best = None
        for k in extra_len:
            for d in dir_prefs(anchor_id):
                cand = add(pos[anchor_id], mul(d, step_len(anchor_id, d, new_id) * k))
                # conflicts weigh more than length, but every extra step length costs:
                # no unnecessarily long edges, and never a disconnected node
                v = violations(pos[anchor_id], cand) * 3 + P.length_penalty * (k - extra_len[0])
                # other partners already drawn pull the node towards them (loops close compactly)
                v += P.partner_pull * sum(max(0.0, dist(cand, pos[q]) / D - 1.0) for q in partners_of(new_id)
                                          if q != anchor_id and q in pos) if new_id else 0.0
                if best is None or v < best[0]:
                    best = (v, d, cand)
        assert best is not None
        return best[1], best[2]   # least-bad candidate (direction, position)

    def place_siblings(rid: str):
        """Reaction template for several main substrates/products: the extra partners are placed in a
        row next to the one already drawn, at the same distance from the reaction (products on the
        product side, substrates on the substrate side), on the side opposite the cofactors."""
        if rid not in rpos:
            return
        ro, R, a = roles[rid], rpos[rid], raxis[rid]
        on_ring = rid in L.ring_rxn
        row0 = ring_side[rid] if on_ring else neg(side_rot(a))
        for lst, sgn in ((ro.main_prod, 1), (ro.main_sub, -1)):
            todo = [m for m in lst if m not in pos and m not in pin_met and cyc_of.get(m) not in
                    set(range(len(cycles))) - placed_cycles]
            if not todo:
                continue
            placed_side = [m for m in lst if m in pos]
            if on_ring:   # outside the ring, next to the reaction
                base = add(R, add(mul(row0, 0.0), mul(a, sgn * 0.35 * D)))
                row_dirs = [row0]
            else:
                base = pos[placed_side[0]] if placed_side else add(R, mul(a, sgn * D / 2))
                row_dirs = [row0, neg(row0)]
            ref = placed_side[0] if placed_side else None
            k = 0
            for m in todo:
                best = None
                for rd in row_dirs:
                    step = 0.8 * D        # partners of one reaction stay close together
                    for kk in range(k + 1, k + 4):
                        if step * kk > 2.2 * D:      # fork partners stay close; otherwise normal branch search
                            break
                        cand = add(base, mul(rd, step * kk if not on_ring else 0.9 * D + step * (kk - 1)))
                        v = violations(R, cand)
                        if best is None or v < best[0]:
                            best = (v, rd, cand, kk)
                        if v == 0:
                            break
                    if best is not None and best[0] == 0:
                        break
                if best is None:
                    continue
                _, rd, cand, kk = best
                if best[0] > 0:      # no free slot in the row (e.g. it would lie on the main line)
                    continue         # -> placed later by the normal branch search
                pos[m] = cand
                ctx[m], attach[m] = "branch", rd
                if ref:
                    parent[m] = ref
                lab_pref[m] = mul(a, sgn)
                order.append(m)
                k, ref = kk, m

    for rid in sorted(rpos):
        place_siblings(rid)

    unplaced = [m for m in main_mets if m not in pos]
    while unplaced:
        unplaced = [m for m in unplaced if m not in pos]
        if not unplaced:
            break
        rank = {m: i for i, m in enumerate(order)}
        cands = []
        for m in unplaced:
            for rid, ro in roles.items():
                if m in ro.main_sub:
                    partners = ro.main_prod
                elif m in ro.main_prod:
                    partners = ro.main_sub
                else:
                    continue
                for a in partners:
                    if a in pos:
                        cands.append(((0 if rid not in rpos else 1, rank.get(a, 1e9), rid, m), m, a, rid))
        if not cands:
            # disconnected remainder: start a new column right of the current drawing
            m = unplaced[0]
            xs = [p for p in list(pos.values()) + list(rpos.values())] or [(0.0, 0.0)]
            far = max(dot(p, S) for p in xs) + 2 * D
            start = min(dot(p, F) for p in xs)
            pos[m] = pin_met.get(m, add(mul(S, far), mul(F, start)))
            ctx[m], attach[m] = "backbone", F
            order.append(m)
            unplaced.remove(m)
            continue
        _, m, a, rid = min(cands, key=lambda c: c[0])
        k = cyc_of.get(m)
        if k is not None and k not in placed_cycles and not (set(cycles[k]) & set(pin_met)):
            # a cycle reached from a branch: draw it as a circle next to its anchor
            n_ = len(cycles[k])
            R_ = max(n_ * P.cycle_spacing / (2 * math.pi), P.min_cycle_radius)
            best_c = None
            for gap in (1.0, 1.6, 2.4):
                for d in dir_prefs(a):
                    c_ = add(pos[a], mul(d, D * gap + R_))
                    v = sum(3 for p in list(pos.values()) + list(rpos.values()) if dist(p, c_) < R_ + 0.6 * D)
                    v += sum(5 for c in L.cycles if dist(c.center, c_) < c.radius + R_ + 0.6 * D)
                    v += sum(2 for p1, p2 in main_segments() if seg_point_dist(c_, p1, p2) < R_ + 0.3 * D)
                    if best_c is None or v < best_c[0]:
                        best_c = (v, d, gap)
                if best_c is not None and best_c[0] == 0:
                    break
            assert best_c is not None
            place_cycle(k, a, best_c[1], best_c[2])
            for mm in cycles[k]:
                if mm in unplaced:
                    unplaced.remove(mm)
            for r_ in sorted(L.ring_rxn):
                place_siblings(r_)
            continue
        if m in pin_met:
            pos[m] = pin_met[m]
            d = unit(sub(pos[m], pos[a]))
        elif rid in rpos:
            # the reaction is already drawn: put the metabolite next to the reaction, not next to
            # some other partner of it (keeps the link short)
            R_ = rpos[rid]
            sgn_ = 1 if m in roles[rid].main_prod else -1
            best_r = None
            for kk in (1.0, 2.0, 3.0):
                for d_ in (mul(raxis[rid], sgn_), neg(side_rot(raxis[rid])), side_rot(raxis[rid]), mul(raxis[rid], -sgn_)):
                    cand = add(R_, mul(d_, kk * D / 2))
                    v = violations(R_, cand) * 3 + P.length_penalty * (kk - 1) / 2
                    if best_r is None or v < best_r[0]:
                        best_r = (v, d_, cand)
                if best_r is not None and best_r[0] == 0:
                    break
            assert best_r is not None
            _, d, p_r = best_r
            pos[m] = p_r
        else:
            d, p_r = place_from(a, new_id=m)
            pos[m] = p_r
        ctx[m], attach[m] = "branch", d
        parent[m] = a
        order.append(m)
        unplaced.remove(m)
        if rid not in rpos:
            ro = roles[rid]
            s_, p_ = (a, m) if a in ro.main_sub else (m, a)
            rpos[rid] = pin_rxn.get(rid, lerp(pos[s_], pos[p_], 0.5))
            mid_of[rid] = (s_, p_, (0.0, 0.0))
            ax = unit(sub(pos[p_], pos[s_]))
            raxis[rid] = ax if dot(ax, ax) else F
            place_siblings(rid)

    # remaining reactions (both ends already placed, one-sided, or without main participants)
    for rid, ro in roles.items():
        if rid in rpos:
            continue
        ms = [pos[m] for m in ro.main_sub if m in pos]
        mp = [pos[m] for m in ro.main_prod if m in pos]
        pairs = sorted((dist(pos[s_], pos[p_]), s_, p_) for s_ in ro.main_sub if s_ in pos
                       for p_ in ro.main_prod if p_ in pos)
        def _through_node(x):
            # the straight line between the pair runs over another metabolite (e.g. a GAPDH variant
            # GAP -> 3PG that skips 1,3-bisphosphoglycerate on the main line): never draw it there
            return any(k_ not in (x[1], x[2]) and seg_point_dist(q, pos[x[1]], pos[x[2]]) < 0.3 * D
                       for k_, q in pos.items())
        free_pairs = [x for x in pairs if all(dist(lerp(pos[x[1]], pos[x[2]], 0.5), q) > 0.3 * D
                                              for q in rpos.values()) and not _through_node(x)]
        pairs = free_pairs or pairs
        # prefer a pair whose connecting line does not run through other nodes
        others_pts = list(pos.items())
        pairs = sorted(pairs, key=lambda x: x[0] + 3 * D * sum(
            1 for k_, q in others_pts if k_ not in (x[1], x[2]) and seg_point_dist(q, pos[x[1]], pos[x[2]]) < 20))
        if rid in pin_rxn or pairs:
            # between the closest pair of partners that has a free midpoint; the rest become copies/forks
            if pairs:
                _, s_, p_ = pairs[0]
                off = (0.0, 0.0)
                if not free_pairs:
                    # a parallel reaction for the same pair is already there: put this one beside it
                    # (diamond), cofactors on its own outer side
                    mid = lerp(pos[s_], pos[p_], 0.5)
                    ax0 = unit(sub(pos[p_], pos[s_]))
                    perp = side_rot(ax0) if dot(ax0, ax0) else S
                    w_ = 2 * P.side_offset + 4
                    for k_ in (1, -1, 2, -2, 3, -3):
                        off = mul(perp, k_ * w_)
                        if all(dist(add(mid, off), q) > 0.3 * D for q in rpos.values()):
                            break
                    par_side[rid] = unit(off)
                rpos[rid] = pin_rxn.get(rid, add(lerp(pos[s_], pos[p_], 0.5), off))
                mid_of[rid] = (s_, p_, off)
                ax = unit(sub(pos[p_], pos[s_]))
            else:
                rpos[rid] = pin_rxn[rid]
                ax = unit(sub(mean(mp), mean(ms))) if ms and mp else F
            raxis[rid] = ax if dot(ax, ax) else F
        elif ms or mp:
            anchor = next(m for m in ro.main_sub + ro.main_prod if m in pos)
            d, cand = place_from(anchor, (0.5, 1.0, 1.5))
            rpos[rid] = pin_rxn.get(rid, cand)
            raxis[rid] = d if anchor in ro.main_sub else neg(d)
        else:
            xs = list(pos.values()) + list(rpos.values()) or [(0.0, 0.0)]
            rpos[rid] = pin_rxn.get(rid, add(mul(S, max(dot(p, S) for p in xs) + D),
                                             mul(F, max(dot(p, F) for p in xs) + D)))
            raxis[rid] = F

    # ---------------------------------------------------------------- shortcuts inside circles
    # a reaction linking two non-neighbouring members of a circle (e.g. isocitrate lyase:
    # isocitrate -> glyoxylate + succinate on the TCA cycle) is drawn inside the circle, like the
    # glyoxylate shunt in textbooks; branch metabolites that only connect to the circle (glyoxylate)
    # move inside with it. Nothing is duplicated.
    inside: dict[str, int] = {}
    inside_mets: set[str] = set()
    for rid in sorted(rpos):
        if rid in L.ring_rxn or rid in pin_rxn:
            continue
        ro = roles[rid]
        on_ring = [m for m in ro.main_sub + ro.main_prod if m in L.ring_met]
        if len(on_ring) < 2 or len({L.ring_met[m][0] for m in on_ring}) != 1:
            continue
        k_ = L.ring_met[on_ring[0]][0]
        cy = L.cycles[k_]
        step_ = 2 * math.pi / len(cy.slots)
        if len(on_ring) == 2 and abs(abs(wrap(L.ring_met[on_ring[0]][1] - L.ring_met[on_ring[1]][1])) - 2 * step_) < 1e-3:
            continue                     # neighbours on the circle: parallel reaction, handled elsewhere
        inside[rid] = k_
    for rid, k_ in inside.items():
        cy = L.cycles[k_]
        ro = roles[rid]
        ring_pts = [pos[m] for m in ro.main_sub + ro.main_prod if m in L.ring_met]
        rpos[rid] = lerp(cy.center, mean(ring_pts), 0.6)
        mid_of.pop(rid, None)
    for m in sorted(pos):
        if m in pin_met or ctx.get(m) != "branch":
            continue
        rids_m = [r for r, ro in roles.items() if m in ro.main_sub + ro.main_prod]
        ks = {inside.get(r) for r in rids_m}
        ring_partners = [q for r in rids_m for q in roles[r].main_sub + roles[r].main_prod if q in L.ring_met]
        if not rids_m or None in {inside.get(r) if inside.get(r) is not None else
                                  (L.ring_met[ring_partners[0]][0] if ring_partners else None) for r in rids_m[:1]}:
            continue
        if not any(r in inside for r in rids_m) or not ring_partners:
            continue
        cy = L.cycles[next(inside[r] for r in rids_m if r in inside)]
        pos[m] = lerp(cy.center, mean([pos[q] for q in ring_partners]), 0.3)
        inside_mets.add(m)
    for rid, k_ in inside.items():     # re-place the inside reactions between all their partners
        ro = roles[rid]
        rpos[rid] = mean([pos[m] for m in ro.main_sub + ro.main_prod if m in pos])
        ms_ = [pos[m] for m in ro.main_sub if m in pos]
        mp_ = [pos[m] for m in ro.main_prod if m in pos]
        if ms_ and mp_:
            ax = unit(sub(mean(mp_), mean(ms_)))
            raxis[rid] = ax if dot(ax, ax) else raxis[rid]
        ring_side[rid] = neg(side_rot(raxis[rid]))

    # ---------------------------------------------------------------- alignment
    # metabolites with nearly the same x (or y) are put on one vertical (horizontal) line; the main
    # line, circles and pinned nodes are anchors. Alignment matters more than equal edge lengths.
    movable = {m for m in pos if m not in pin_met and ctx.get(m) not in ("ring", "backbone")}
    tol = P.align_tol * D
    for axis_ in (0, 1):
        order_ = sorted(pos, key=lambda m: (pos[m][axis_], m))
        clusters, cur = [], []
        for m in order_:
            if cur and pos[m][axis_] - pos[cur[0]][axis_] > tol:
                clusters.append(cur)
                cur = []
            cur.append(m)
        if cur:
            clusters.append(cur)
        for cl in clusters:
            if len(cl) < 2 or not (set(cl) & movable):
                continue
            anchors = sorted(pos[m][axis_] for m in cl if m not in movable)
            vals = anchors or sorted(pos[m][axis_] for m in cl)
            target = vals[len(vals) // 2]
            for m in cl:
                if m not in movable or abs(pos[m][axis_] - target) < 1e-6:
                    continue
                new = (target, pos[m][1]) if axis_ == 0 else (pos[m][0], target)
                if all(dist(new, q) >= 0.6 * D for o, q in pos.items() if o != m) and \
                        all(dist(new, q) >= 0.35 * D for r_, q in rpos.items() if mid_of.get(r_, (None, None))[:2].count(m) == 0):
                    pos[m] = new
    # parallel reactions between the same two metabolites: always symmetric around the line
    groups: dict[tuple, list] = {}
    for rid in sorted(rpos):
        ro = roles[rid]
        if rid in L.ring_rxn or rid in pin_rxn or rid in ring_side or len(ro.main_sub) != 1 or len(ro.main_prod) != 1:
            continue
        key = tuple(sorted((ro.main_sub[0], ro.main_prod[0])))
        if all(k_ in pos for k_ in key):
            groups.setdefault(key, []).append(rid)
    for (u, v), rids_g in groups.items():
        if len(rids_g) < 2:
            continue
        ax0 = unit(sub(pos[v], pos[u]))
        perp = side_rot(ax0) if dot(ax0, ax0) else S
        n_ = len(rids_g)
        # the straight line u-v runs over another metabolite (reactions skipping a main-line
        # metabolite, e.g. GAPDH variants GAP -> 3PG past 1,3-bisphosphoglycerate): none on the line
        blocked = any(o not in (u, v) and seg_point_dist(q, pos[u], pos[v]) < 0.3 * D for o, q in pos.items())
        ks = [((i_ // 2) + 1) * (1 if i_ % 2 == 0 else -1) for i_ in range(n_)] if blocked else \
            [j - (n_ - 1) / 2 for j in range(n_)]
        for j, rid in enumerate(sorted(rids_g, key=lambda r: (-W[r], r))):
            k_ = ks[j]
            off = mul(perp, k_ * (2 * P.side_offset + 4))
            mid_of[rid] = (u, v, off)
            if abs(k_) > 1e-9:
                par_side[rid] = mul(perp, 1 if k_ > 0 else -1)
            else:
                par_side.pop(rid, None)

    multi = {rid for rid in rpos if rid not in L.ring_rxn and rid not in pin_rxn and rid not in ring_side
             and len(roles[rid].main_sub) + len(roles[rid].main_prod) >= 3
             and roles[rid].main_sub and roles[rid].main_prod}

    # ---------------------------------------------------------------- symmetric forks
    # reaction with a row of partners on one side and a single partner on the other: if that single
    # partner belongs to no other reaction, it moves to the middle of the row (Y shape: two substrates
    # side by side, reaction and product centred underneath)
    for rid in sorted(rpos):
        ro = roles[rid]
        if rid in L.ring_rxn or rid in pin_rxn:
            continue
        for one, row in ((ro.main_prod, ro.main_sub), (ro.main_sub, ro.main_prod)):
            if len(one) != 1 or len(row) < 2 or one[0] not in pos or one[0] in pin_met or one[0] in L.ring_met:
                continue
            m = one[0]
            if sum(1 for r_, ro_ in roles.items() if m in ro_.main_sub + ro_.main_prod) != 1:
                continue
            pts = [pos[q] for q in row if q in pos]
            if len(pts) < 2:
                continue
            c_ = mean(pts)
            # move along the row direction only
            rowdir = unit(sub(pts[-1], pts[0]))
            shift = dot(sub(c_, pos[m]), rowdir)
            new = add(pos[m], mul(rowdir, shift))
            if all(dist(new, q) >= 0.6 * D for o, q in pos.items() if o != m):
                pos[m] = new

    # ---------------------------------------------------------------- grid
    # metabolites snap to grid cells (nearest free cell if taken); circles and pinned nodes stay
    g = P.grid * D

    def snap(v: float, unit_: float) -> float:
        return round(v / unit_) * unit_ if unit_ > 0 else v

    if g > 0:
        for m in list(order) + sorted(set(pos) - set(order)):
            if m not in pos or m in pin_met or m in L.ring_met or ctx.get(m) == "anchor" or m in inside_mets:
                continue    # circles and nodes placed relative to a circle keep their exact position
            base = (snap(pos[m][0], g), snap(pos[m][1], g))
            best_ = None
            for rr in range(0, 6):
                for dx in range(-rr, rr + 1):
                    for dy in range(-rr, rr + 1):
                        if max(abs(dx), abs(dy)) != rr:
                            continue
                        q = (base[0] + dx * g, base[1] + dy * g)
                        if all(dist(q, p_) >= 0.6 * D for o, p_ in pos.items() if o != m):
                            cand = (dist(q, pos[m]), dx, dy, q)
                            if best_ is None or cand < best_:
                                best_ = cand
                if best_:
                    break
            if best_:
                pos[m] = best_[3]

    def snap_rxn(rid):
        if g > 0 and rid not in pin_rxn and rid not in L.ring_rxn and rid not in ring_side:
            rpos[rid] = (snap(rpos[rid][0], g / 2), snap(rpos[rid][1], g / 2))

    def recentre():
        for rid, (u, v, off) in mid_of.items():      # re-centre reactions between their metabolites
            if rid in pin_rxn or rid in multi or u not in pos or v not in pos:
                continue
            rpos[rid] = add(lerp(pos[u], pos[v], 0.5), off)
            ax = unit(sub(pos[v], pos[u]))
            raxis[rid] = ax if dot(ax, ax) else raxis[rid]
        for rid in multi:
            # several main partners: the reaction sits centred between the substrate row and the
            # product row (e.g. two substrates side by side, reaction centred underneath)
            ms_ = [pos[m] for m in roles[rid].main_sub if m in pos]
            mp_ = [pos[m] for m in roles[rid].main_prod if m in pos]
            if ms_ and mp_:
                rpos[rid] = lerp(mean(ms_), mean(mp_), 0.5)
                ax = unit(sub(mean(mp_), mean(ms_)))
                raxis[rid] = ax if dot(ax, ax) else raxis[rid]
        for rid in list(mid_of) + list(multi):           # reactions on the half grid
            snap_rxn(rid)
    recentre()

    # ---------------------------------------------------------------- regulons
    # a regulator with several targets that take part in nothing else (transcription factor ->
    # the genes it switches on): the targets sit evenly on an arc around it, on its free side, each
    # reaction on the straight ray between them - a tidy fan instead of a crowded star
    regulon: dict[str, list] = {}
    fans: list = []                                         # (regulator, radius): rigid for compaction
    for rid in sorted(rpos):
        ro = roles[rid]
        for m in ro.main_sub + ro.main_prod:
            if m in pw.reactions[rid].modifiers and m in pos and m not in pin_met:
                tgt = [x for x in ro.main_sub + ro.main_prod if x != m]
                if len(tgt) == 1 and tgt[0] in pos and tgt[0] not in pin_met and rid not in pin_rxn:
                    regulon.setdefault(m, []).append((rid, tgt[0]))
    for m, items in sorted(regulon.items()):
        if len(items) < 3:
            continue
        n_ = len(items)
        tg_set = {t for _, t in items} | {r for r, _ in items}
        others = [q for k_, q in pos.items() if k_ != m and k_ not in tg_set] + \
            [q for k_, q in rpos.items() if k_ not in tg_set]
        fan_step = 0.75 * D                                 # distance between neighbouring targets
        rad = max(1.8 * D, n_ * fan_step / math.pi)
        span = min(2 * math.pi * (n_ - 1) / n_, n_ * fan_step / rad)
        best_ = None
        for i_ in range(24):                                # centre direction of the fan
            th_c = 2 * math.pi * i_ / 24
            occ = 0.0
            for q in others:
                d_ = dist(q, pos[m])
                if d_ < 1e-6 or d_ > rad + 1.2 * D:
                    continue
                dth = abs(wrap(math.atan2(q[1] - pos[m][1], q[0] - pos[m][0]) - th_c))
                if dth <= span / 2 + 0.35:
                    occ += 1.0 + (rad + 1.2 * D - d_) / D
            key = (occ, abs(wrap(th_c - math.pi / 2)), i_)  # prefer below the regulator
            if best_ is None or key < best_[0]:
                best_ = (key, th_c)
        assert best_ is not None
        th_c = best_[1]
        fans.append((m, rad))
        ctx[m] = "regulon"            # the hub stays where its fan is (not moved by the crossing search)
        order = sorted(items, key=lambda it: (display_name(pw, it[1], P).lower(), it[0]))
        for j_, (rid, t) in enumerate(order):
            th = th_c - span / 2 + (span * j_ / (n_ - 1) if n_ > 1 else 0.0)
            pos[t] = add(pos[m], (rad * math.cos(th), rad * math.sin(th)))
            ctx[t] = "regulon"
            ro = roles[rid]
            u, v = (ro.main_sub[0], ro.main_prod[0])
            mid_of[rid] = (u, v, (0.0, 0.0))
        recentre()

    # ---------------------------------------------------------------- straight reactions
    # metabolite -> reaction -> metabolite is drawn as one straight line: a reaction with one main
    # substrate and one main product is centred between them (unless it belongs to a circle, a
    # parallel group or a pinned layout). Done before the crossing search, so metabolites are moved
    # to suit the straight line; only if that fails does the reaction return to its old spot.
    kinked0: dict[str, tuple] = {}
    pair_n: dict[tuple, int] = {}
    for rid_, (u_, v_, _o) in mid_of.items():
        pair_n[tuple(sorted((u_, v_)))] = pair_n.get(tuple(sorted((u_, v_))), 0) + 1
    for rid, (u, v, off) in sorted(mid_of.items()):
        # a single reaction put beside its line (its straight line ran over another reaction when it
        # was placed): straighten it too, unless the line runs over a metabolite (a real bypass)
        if (abs(off[0]) + abs(off[1]) < 1e-9 or rid in pin_rxn or rid in multi or rid in L.ring_rxn
                or pair_n[tuple(sorted((u, v)))] > 1 or u not in pos or v not in pos):
            continue
        if any(o not in (u, v) and seg_point_dist(q, pos[u], pos[v]) < 0.3 * D for o, q in pos.items()):
            continue
        kinked0[rid] = (rpos[rid], raxis[rid])
        mid_of[rid] = (u, v, (0.0, 0.0))
    for rid in sorted(rpos):
        ro = roles[rid]
        if (rid in mid_of or rid in pin_rxn or rid in L.ring_rxn or rid in ring_side or rid in multi
                or len(ro.main_sub) != 1 or len(ro.main_prod) != 1):
            continue
        u, v = ro.main_sub[0], ro.main_prod[0]
        if u not in pos or v not in pos:
            continue
        kinked0[rid] = (rpos[rid], raxis[rid])     # from now on it follows its metabolites
        mid_of[rid] = (u, v, (0.0, 0.0))
    if kinked0:
        recentre()

    # ---------------------------------------------------------------- no unnecessary crossings
    # local search over branch metabolites: single grid steps, moving or mirroring a whole branch,
    # and swaps. Cost (only for the lines that change): 10 x crossings + 10 x lines running through a
    # node + 3 x misaligned partners + length. Moves are horizontal/vertical, so alignment is kept.
    incident: dict[str, list] = {}
    for rid in rpos:
        for m in roles[rid].main_sub + roles[rid].main_prod:
            incident.setdefault(m, []).append(rid)

    def segs_of(rids):
        return [(rid, m, rpos[rid], pos[m]) for rid in rids for m in roles[rid].main_sub + roles[rid].main_prod
                if m in pos]

    def misaligned(rid):
        ro = roles[rid]
        if rid in L.ring_rxn:
            return 0
        bad = 0
        if len(ro.main_sub) == 1 and len(ro.main_prod) == 1:
            a_, b_ = pos.get(ro.main_sub[0]), pos.get(ro.main_prod[0])
            if a_ and b_ and abs(a_[0] - b_[0]) > 0.5 and abs(a_[1] - b_[1]) > 0.5:
                # a bent main line is worse than a bent branch
                bad += 3 if ctx.get(ro.main_sub[0]) == ctx.get(ro.main_prod[0]) == "backbone" else 1
        else:
            for row in (ro.main_sub, ro.main_prod):   # rows of partners should share x or y
                pts = [pos[m] for m in row if m in pos]
                for i in range(1, len(pts)):
                    if abs(pts[i][0] - pts[0][0]) > 0.5 and abs(pts[i][1] - pts[0][1]) > 0.5:
                        bad += 1
        return bad

    def local_cost(moved: set) -> float:
        rids = {r for m in moved for r in incident.get(m, [])}
        A = segs_of(rids)
        O = segs_of(set(rpos) - rids)
        cr = 0
        for i, (r1, m1, a1, b1) in enumerate(A):
            for (r2, m2, a2, b2) in O + A[i + 1:]:
                if r1 != r2 and m1 != m2 and segs_cross(a1, b1, a2, b2):
                    cr += 1
        hits = 0
        tol_m, tol_r = P.met_radius + 3, P.rxn_half + 3

        def near_seg(p_, a_, b_, tol):          # bounding-box test first (cheap), then exact
            return (min(a_[0], b_[0]) - tol <= p_[0] <= max(a_[0], b_[0]) + tol and
                    min(a_[1], b_[1]) - tol <= p_[1] <= max(a_[1], b_[1]) + tol and
                    seg_point_dist(p_, a_, b_) < tol)
        for i, (r1, m1, a1, b1) in enumerate(A):
            for q, pq in pos.items():
                if q != m1 and near_seg(pq, a1, b1, tol_m):
                    hits += 1
            for r_, pr in rpos.items():          # lines running through another reaction
                if r_ != r1 and near_seg(pr, a1, b1, tol_r):
                    hits += 1
            # lines lying on top of each other (collinear overlap) are as bad as crossings
            mid1 = lerp(a1, b1, 0.5)
            for (r2, m2, a2, b2) in O + A[i + 1:]:
                if r1 != r2 and near_seg(mid1, a2, b2, 5) and seg_point_dist(lerp(a2, b2, 0.5), a1, b1) < 5:
                    hits += 1
        for q in moved:
            for (r1, m1, a1, b1) in O:
                if q != m1 and near_seg(pos[q], a1, b1, tol_m):
                    hits += 1
        for r_ in rids:                          # moved reactions sitting on other lines
            for (r1, m1, a1, b1) in O:
                if r_ != r1 and near_seg(rpos[r_], a1, b1, tol_r):
                    hits += 1
        return 10.0 * (cr + hits) + 3.0 * sum(misaligned(r) for r in rids) + P.length_weight * sum(dist(a_, b_) for *_, a_, b_ in A) / D

    def global_bad() -> int:
        """crossings + 2 x reaction nodes on another reaction's line, for the whole component"""
        segs = segs_of(list(rpos))
        n = sum(1 for i, (r1, m1, a1, b1) in enumerate(segs) for (r2, m2, a2, b2) in segs[i + 1:]
                if r1 != r2 and m1 != m2 and segs_cross(a1, b1, a2, b2))
        return n + 2 * sum(1 for r1, m1, a1, b1 in segs for r_, q in rpos.items()
                           if r_ != r1 and seg_point_dist(q, a1, b1) < P.rxn_half + 4)

    def clear_all(moved: set) -> bool:
        for m in moved:
            q = pos[m]
            if any(dist(q, p_) < 0.7 * D for o, p_ in pos.items() if o != m):
                return False
            if any(dist(q, p_) < 0.4 * D for r_, p_ in rpos.items() if r_ not in incident.get(m, [])):
                return False
            if any(dist(q, cy.center) < cy.radius + 0.4 * D for cy in L.cycles):
                return False
        return True

    # branch metabolites, and main-line metabolites too (leaving the line costs via the alignment
    # term, but a crossing costs more); circles and pinned nodes stay
    movable2 = sorted(m for m in pos if m not in pin_met and ctx.get(m) in ("branch", "backbone"))
    children: dict[str, list] = {}
    for m, p_ in parent.items():
        if m in movable2:
            children.setdefault(p_, []).append(m)

    def subtree(m):
        out, stack = [], [m]
        while stack:
            x = stack.pop()
            if x in out or x not in movable2:
                continue
            out.append(x)
            stack += children.get(x, [])
        return out

    def try_move(new_pos: dict) -> bool:
        moved = set(new_pos)
        before = local_cost(moved)
        old = {m: pos[m] for m in moved}
        old_r = dict(rpos)
        old_ax = dict(raxis)
        pos.update(new_pos)
        if not clear_all(moved):
            pos.update(old)
            return False
        recentre()
        if local_cost(moved) < before - 1e-6:
            return True
        pos.update(old)
        rpos.update(old_r)
        raxis.update(old_ax)
        return False

    rigid_line = bool(L.meta.get("main_path"))   # known pathway (glycolysis) or --main-path
    gcell = g if g > 0 else D / 2          # all search moves are whole grid cells
    if movable2 and len(rpos) <= 150:
        st_ = 2 * gcell
        steps = ((st_, 0.0), (-st_, 0.0), (0.0, st_), (0.0, -st_))
        for _ in range(4):
            improved = False
            for m in movable2:
                for d in steps + tuple(mul(x, 2) for x in steps):
                    if rigid_line and ctx.get(m) == "backbone" and abs(dot(d, attach.get(m, F))) < 1e-9:
                        continue      # a known/user main path stays straight: only along it, never sideways
                    if try_move({m: add(pos[m], d)}):
                        improved = True
                        break
                sub_ = subtree(m)
                if len(sub_) > 1:
                    for d in steps:
                        if rigid_line and ctx.get(m) == "backbone" and abs(dot(d, attach.get(m, F))) < 1e-9:
                            continue
                        if try_move({x: add(pos[x], d) for x in sub_}):
                            improved = True
                            break
                    if m in parent and parent[m] in pos:
                        pp = pos[parent[m]]
                        ax_ = unit(sub(pos[m], pp))
                        if dot(ax_, ax_):
                            nrm = (-ax_[1], ax_[0])   # mirror the branch across the line parent->m
                            refl = {x: add(pp, sub(sub(pos[x], pp), mul(nrm, 2 * dot(sub(pos[x], pp), nrm))))
                                    for x in sub_}
                            if try_move(refl):
                                improved = True
                # a branch may move to the opposite side of the node it hangs from (point reflection
                # of the whole branch), e.g. (R)-lactoyl-CoA left instead of right of acryloyl-CoA
                # when its other link goes to the left
                if ctx.get(m) != "backbone" and m in parent and parent[m] in pos and \
                        all(x not in pin_met for x in sub_):
                    pp = pos[parent[m]]
                    out_ = unit(sub(pp, pos[m]))             # direction to the opposite side
                    side_ = (-out_[1], out_[0])
                    before_g = global_bad()
                    for extra in ((0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (0.0, -1.0), (1.0, 1.0), (1.0, -1.0)):
                        dv = add(mul(out_, extra[0] * gcell), mul(side_, extra[1] * gcell))
                        snap_ = (dict(pos), dict(rpos), dict(raxis))
                        opp = {x: add(sub(mul(pp, 2.0), pos[x]), dv) for x in sub_}
                        if try_move(opp):
                            if global_bad() >= before_g:     # only when it removes a crossing overall
                                pos.update(snap_[0]); rpos.update(snap_[1]); raxis.update(snap_[2])
                                recentre()
                                continue
                            improved = True
                            break
                # a leaf (one 1:1 reaction to one partner) may jump around its partner: right, left,
                # above or below, at 1-2x the spacing (frees the spot it blocked, e.g. a straight line)
                if ctx.get(m) != "backbone" and len(sub_) == 1 and len(incident.get(m, [])) == 1:
                    ro0 = roles[incident[m][0]]
                    others = [x for x in ro0.main_sub + ro0.main_prod if x != m and x in pos]
                    if len(others) == 1:
                        q0 = pos[others[0]]
                        done_ = False
                        for k_ in (1.0, 1.5, 2.0):
                            step_ = max(1, round(k_ * D / gcell)) * gcell
                            for dv in ((1.0, 0.0), (-1.0, 0.0), (0.0, 1.0), (0.0, -1.0)):
                                cand = add(q0, mul(dv, step_))
                                if dist(cand, pos[m]) > 1 and try_move({m: cand}):
                                    improved = done_ = True
                                    break
                            if done_:
                                break
            for i in range(len(movable2)):
                for j in range(i + 1, len(movable2)):
                    a_, b_ = movable2[i], movable2[j]
                    if ctx.get(a_) == "backbone" or ctx.get(b_) == "backbone":
                        continue          # never reorder the main line or put a branch onto it
                    if dist(pos[a_], pos[b_]) <= 3 * D and try_move({a_: pos[b_], b_: pos[a_]}):
                        improved = True
            if not improved:
                break
        recentre()

        # where crossings remain: simulated annealing over the metabolites involved (fixed seed,
        # so the result is reproducible)
        def crossing_nodes():
            segs = segs_of(list(rpos))
            bad = set()
            for i, (r1, m1, a1, b1) in enumerate(segs):
                for (r2, m2, a2, b2) in segs[i + 1:]:
                    if r1 != r2 and m1 != m2 and segs_cross(a1, b1, a2, b2):
                        bad |= {m1, m2} | set(roles[r1].main_sub + roles[r1].main_prod) | \
                            set(roles[r2].main_sub + roles[r2].main_prod)
            return sorted(m for m in bad if m in movable2 and (P.sa_backbone or ctx.get(m) != "backbone"))

        def n_crossings():
            segs = segs_of(list(rpos))
            n = sum(1 for i, (r1, m1, a1, b1) in enumerate(segs) for (r2, m2, a2, b2) in segs[i + 1:]
                    if r1 != r2 and m1 != m2 and segs_cross(a1, b1, a2, b2))
            return n

        def anneal(seed: int, focus: list):
            import random
            rng = random.Random(seed)
            iters = min(3000, 150 * len(focus))
            k_max = max(1, int(round(P.sa_max_step * D / gcell)))
            offsets = [x * gcell for x in range(-k_max, k_max + 1) if x]
            state_cost = 0.0
            best_state = ({m: pos[m] for m in pos}, dict(rpos), dict(raxis), 0.0)
            for it in range(iters):
                T = 8.0 * (1 - it / iters) + 0.05
                if rng.random() < 0.7:
                    m = focus[rng.randrange(len(focus))]
                    d_ = rng.choice(offsets)
                    new = {m: add(pos[m], (d_, 0.0) if rng.random() < 0.5 else (0.0, d_))}
                else:
                    if len(focus) < 2:
                        continue
                    a_, b_ = rng.sample(focus, 2)
                    new = {a_: pos[b_], b_: pos[a_]}
                moved = set(new)
                before = local_cost(moved)
                old = {k_: pos[k_] for k_ in moved}
                old_r, old_ax = dict(rpos), dict(raxis)
                pos.update(new)
                if not clear_all(moved):
                    pos.update(old)
                    continue
                recentre()
                delta = local_cost(moved) - before
                if delta < 0 or rng.random() < math.exp(-delta / T):
                    state_cost += delta
                    if state_cost < best_state[3] - 1e-6:
                        best_state = ({k_: pos[k_] for k_ in pos}, dict(rpos), dict(raxis), state_cost)
                else:
                    pos.update(old)
                    rpos.update(old_r)
                    raxis.update(old_ax)
            return best_state

        focus = crossing_nodes()
        if focus:
            start = ({m: pos[m] for m in pos}, dict(rpos), dict(raxis))
            seeds = (1729, 7, 4242) if len(focus) <= 15 or len(rpos) <= 40 else (1729,)
            results = []
            for seed in seeds:
                pos.update(start[0]); rpos.update(start[1]); raxis.update(start[2])
                st = anneal(seed, focus)
                pos.update(st[0]); rpos.update(st[1]); raxis.update(st[2])
                recentre()
                results.append((n_crossings(), st[3], seed, ({m: pos[m] for m in pos}, dict(rpos), dict(raxis))))
            best = min(results, key=lambda x: (x[0], x[1], x[2]))
            pos.update(best[3][0]); rpos.update(best[3][1]); raxis.update(best[3][2])
            recentre()

    # ---------------------------------------------------------------- side of offset reactions
    # reactions drawn beside the straight line between their metabolites (bypasses that skip a
    # main-line metabolite, parallel reactions) were put on a fixed side; flip a reaction (or its
    # whole parallel group) to the other side when that gives fewer crossings / lines through nodes
    def bad_total():
        """crossings + lines through a metabolite + reaction nodes lying on another reaction's line
        (that looks like one reaction with two substrates and two products, so it counts double)"""
        segs = segs_of(list(rpos))
        n = sum(1 for i, (r1, m1, a1, b1) in enumerate(segs) for (r2, m2, a2, b2) in segs[i + 1:]
                if r1 != r2 and m1 != m2 and segs_cross(a1, b1, a2, b2))
        for r_, m_, a_, b_ in segs:
            own = roles[r_].main_sub + roles[r_].main_prod
            n += sum(1 for o, q in pos.items() if o not in own and seg_point_dist(q, a_, b_) < 16)
            n += 2 * sum(1 for o, q in rpos.items() if o != r_ and seg_point_dist(q, a_, b_) < P.rxn_half + 4)
        return n

    # straightened reactions (see above) whose straight position still sits badly after the search
    # (on another line, through a node, causing a crossing): try the grid points around the midpoint
    # and take the one with the fewest problems, then the one closest to the straight line
    if kinked0:
        cur_bad_s = bad_total()
        hstep = (g if g > 0 else D / 2) / 2           # reactions sit on the half grid
        for rid in sorted(kinked0):
            if rid not in mid_of:
                continue
            u, v, _ = mid_of[rid]
            mid = lerp(pos[u], pos[v], 0.5)
            straight = (cur_bad_s, 0.0, None)
            best = straight
            keep = (mid_of.pop(rid), rpos[rid], raxis[rid])
            cands = [kinked0[rid][0]] + [add(mid, (i_ * hstep, j_ * hstep))
                                         for i_ in range(-6, 7) for j_ in range(-6, 7) if i_ or j_]
            for c_ in cands:
                if any(dist(c_, q) < 0.35 * D for q in pos.values()) or \
                        any(dist(c_, q) < 0.3 * D for o, q in rpos.items() if o != rid):
                    continue
                rpos[rid] = c_
                nb = bad_total()
                key = (nb, seg_point_dist(c_, pos[u], pos[v]), c_)
                if key[:2] < best[:2]:
                    best = key
            if best[2] is None:
                mid_of[rid], rpos[rid], raxis[rid] = keep
                recentre()
                continue
            rpos[rid] = best[2]
            raxis[rid] = unit(sub(pos[v], pos[u]))
            cur_bad_s = best[0]

    side_groups: dict[tuple, list] = {}
    for rid, (u, v, off) in sorted(mid_of.items()):
        if rid in pin_rxn or rid in multi or abs(off[0]) + abs(off[1]) < 1e-9:
            continue
        side_groups.setdefault(tuple(sorted((u, v))), []).append(rid)
    if side_groups:
        cur_bad = bad_total()
        for _ in range(2):
            changed_ = False
            for key in sorted(side_groups):
                if cur_bad == 0:
                    break
                grp = side_groups[key]
                old_ = {r: (mid_of[r], rpos[r], raxis[r], par_side.get(r)) for r in grp}
                for r in grp:
                    u, v, off = mid_of[r]
                    mid_of[r] = (u, v, neg(off))
                    if r in par_side:
                        par_side[r] = neg(par_side[r])
                recentre()
                nb = bad_total()
                if nb < cur_bad:
                    cur_bad = nb
                    changed_ = True
                    continue
                for r, (mo, rp, ra, ps) in old_.items():
                    mid_of[r], rpos[r], raxis[r] = mo, rp, ra
                    if ps is not None:
                        par_side[r] = ps
                recentre()
            if not changed_:
                break

    # ---------------------------------------------------------------- compaction
    # Empty strips (rows/columns without any node) left by grid snapping or by enlarging a crowded
    # component are removed: everything beyond such a strip moves back by whole grid cells. Every
    # link keeps at least the base spacing (1.5x with >= 3 parallel reactions), circles stay intact, and no crossing, line through a node or node contact is added.
    if P.compact and not pin_met and not pin_rxn and g > 0 and len(pos) > 2:
        D0 = P.base_spacing or D
        recentred = set(mid_of) | set(multi)
        n_par: dict[tuple, int] = {}
        for rid in rpos:
            ro = roles[rid]
            for s_ in ro.main_sub:
                for p_ in ro.main_prod:
                    key = tuple(sorted((s_, p_)))
                    n_par[key] = n_par.get(key, 0) + 1
        links = [k_ for k_ in n_par if all(x in pos for x in k_)]

        bad_count = bad_total            # crossings, lines through nodes, reactions on lines

        def snapshot():
            return dict(pos), dict(rpos), dict(raxis), [cy.center for cy in L.cycles]

        def restore(sn):
            pos.clear(); pos.update(sn[0]); rpos.clear(); rpos.update(sn[1])
            raxis.clear(); raxis.update(sn[2])
            for cy, c_ in zip(L.cycles, sn[3]):
                cy.center = c_

        # minimum length per link: base spacing; 1.5x with >= 3 parallel reactions (they fan out); a
        # link running past other main-line nodes (bypass, e.g. GAP -> 3PG past 1,3-BPG) 1.5x per step
        need_len: dict[tuple, float] = {}
        for x, y in links:
            inner = sum(1 for o, q in pos.items() if o not in (x, y) and seg_point_dist(q, pos[x], pos[y]) < 0.3 * D)
            f_ = 1.5 if (n_par[(x, y)] >= 3 or inner) else 1.0
            need_len[(x, y)] = D0 * f_ * (inner + 1)
        base_bad = bad_count()
        for axis_ in (1, 0):
            for _ in range(60):
                lv = sorted({round(q[axis_], 3) for q in pos.values()} |
                            {round(q[axis_], 3) for r_, q in rpos.items() if r_ not in recentred})
                blocks = [(cy.center[axis_] - cy.radius - 0.6 * D, cy.center[axis_] + cy.radius + 0.6 * D)
                          for cy in L.cycles] + \
                    [(pos[fm][axis_] - fr - 0.6 * D, pos[fm][axis_] + fr + 0.6 * D) for fm, fr in fans if fm in pos]
                done = False
                for a_, b_ in zip(lv, lv[1:]):
                    kmax = int((b_ - a_) / g + 1e-6) - 1
                    if kmax < 1 or any(lo < b_ and hi > a_ for lo, hi in blocks):
                        continue
                    for k_ in range(kmax, 0, -1):
                        sh = k_ * g
                        sn = snapshot()
                        mv = lambda q: (q[0] - sh, q[1]) if axis_ == 0 else (q[0], q[1] - sh)
                        for m_ in list(pos):
                            if pos[m_][axis_] >= b_ - 1e-6:
                                pos[m_] = mv(pos[m_])
                        for r_ in list(rpos):
                            if rpos[r_][axis_] >= b_ - 1e-6:
                                rpos[r_] = mv(rpos[r_])
                        for cy in L.cycles:
                            if cy.center[axis_] >= b_ - 1e-6:
                                cy.center = mv(cy.center)
                        recentre()
                        # links keep at least the base spacing (1.5x with >= 3 parallel reactions, which
                        # fan out sideways), or their current length if that is already shorter
                        ok = all(dist(pos[x], pos[y]) >= min(dist(sn[0][x], sn[0][y]), need_len[(x, y)]) - 1e-6
                                 for x, y in links)
                        ok = ok and all(dist(pos[x], pos[y]) >= 0.6 * D0 for i_, x in enumerate(sorted(pos))
                                        for y in sorted(pos)[i_ + 1:])
                        ok = ok and all(dist(q, p_) >= 0.3 * D0 for q in rpos.values() for p_ in pos.values())
                        ok = ok and bad_count() <= base_bad
                        if ok:
                            done = True
                            break
                        restore(sn)
                    if done:
                        break
                if not done:
                    break

    # ---------------------------------------------------------------- redox ladders
    # a chain of 2:2 reactions where neighbours hand a redox pair on (photosynthetic electron
    # transport: PSII -> plastoquinol -> cytochrome b6f -> plastocyanin -> PSI -> ferredoxin -> FNR):
    # reactions in a row, each pair in the gap between two reactions, one form above and one below
    ladder = redox_ladder(sorted(rpos), roles, lambda m: norm_name(pw.metabolites[m].name) in ("h2o", "water")) \
        if not pin_met and not pin_rxn else None
    if ladder is not None and set(ladder[0]) == set(rpos):
        chain, cols = ladder
        # compact: the pairs sit tight around the reaction row (independent of any automatic enlarging)
        D0 = P.base_spacing or D
        step, h_ = 0.6 * D0, 0.65 * D0           # column step, rows above / below the reactions
        for m_, (col, up) in cols.items():
            pos[m_] = (col * step, -h_ if up else h_)
            ctx[m_] = "ladder"
        recentre()

    # ---------------------------------------------------------------- modifiers
    # a modifier near its reaction gets a line; a modifier far away (or drawn in another part of the
    # map) gets a clone next to the reaction instead of a line across the drawing; an enzyme that
    # takes part in nothing else gets its own node next to the first reaction it acts on. Spots are
    # chosen across the reaction axis where the new line crosses nothing and the node touches nothing.
    mod_clone: dict[tuple, Pt] = {}
    mod_lines: list = []
    # regulators standing in for an empty side are laid out as main partners: always a line
    regulator_links = {(rid, m) for rid in rpos for m in roles[rid].main_sub + roles[rid].main_prod
                       if m in pw.reactions[rid].modifiers}
    all_segs_m = segs_of(list(rpos)) if any(pw.reactions[r].modifiers for r in rpos) else []

    def mod_spot(rid_, extra_rids=(), skip=None):
        R0 = rpos[rid_]
        ax_ = raxis.get(rid_, (1.0, 0.0))
        nrm = (-ax_[1], ax_[0])
        taken = [q for k_, q in pos.items() if k_ != skip] + [q for k_, q in mod_clone.items() if k_ != skip]
        skip_line = (rid_, skip) if skip else None
        best_ = None
        for k_ in (0.7, 1.0, 1.4, 1.9, 2.5):
            for sgn_ in (1, -1):
                for along in (0.0, 0.5, -0.5, 1.0, -1.0, 1.5, -1.5):
                    c_ = add(R0, add(mul(nrm, sgn_ * k_ * D), mul(ax_, along * D)))
                    if g > 0:
                        c_ = (round(c_[0] / (g / 2)) * (g / 2), round(c_[1] / (g / 2)) * (g / 2))
                    near = sum(1 for q in taken if dist(c_, q) < 0.45 * D) + \
                        sum(1 for q in rpos.values() if dist(c_, q) < 0.3 * D)
                    cross = sum(1 for r_ in (rid_,) + tuple(extra_rids)
                                for (r1, m1, a1, b1) in all_segs_m + mod_lines
                                if r1 != r_ and (r1, m1) != skip_line and segs_cross(c_, rpos[r_], a1, b1))
                    cross += sum(1 for (r1, m1, a1, b1) in all_segs_m + mod_lines if (r1, m1) != skip_line
                                 for r_ in (rid_,) if seg_point_dist(c_, a1, b1) < 2)
                    on_line = sum(1 for (r1, m1, a1, b1) in all_segs_m + mod_lines
                                  if (r1, m1) != skip_line and seg_point_dist(c_, a1, b1) < P.met_radius + 6)
                    thru = sum(1 for k2, q in pos.items() if k2 != skip and seg_point_dist(q, c_, R0) < P.met_radius + 4
                               and dist(q, c_) > 1) + \
                        sum(1 for r2, q in rpos.items() if r2 != rid_ and seg_point_dist(q, c_, R0) < P.rxn_half + 4)
                    # crossings cannot be repaired later, touching nodes can (overlap resolution)
                    cost = 150 * (cross + thru) + 60 * on_line + 40 * near + k_ + 0.5 * abs(along)
                    if best_ is None or cost < best_[0]:
                        best_ = (cost, c_)
        return best_[1]

    def line_ok(a_, rid_, m_):
        """a modifier line: short, crosses no other line, runs over no node"""
        b_ = rpos[rid_]
        if dist(a_, b_) > 2.5 * D:
            return False
        for (r1, m1, a1, b1) in all_segs_m + mod_lines:
            if r1 != rid_ and m1 != m_ and segs_cross(a_, b_, a1, b1):
                return False
        if any(k_ != m_ and seg_point_dist(q, a_, b_) < P.met_radius + 4 for k_, q in pos.items()):
            return False
        if any(r_ != rid_ and seg_point_dist(q, a_, b_) < P.rxn_half + 4 for r_, q in rpos.items()):
            return False
        return True

    for rid in sorted(rpos):
        for m in sorted(pw.reactions[rid].modifiers):
            if m not in pw.metabolites:
                continue
            if Layout.sg(m, rid) in prev.clones:         # kept from the previous layout
                px_, py_ = prev.clones[Layout.sg(m, rid)]
                mod_clone[(rid, m)] = (float(px_), float(py_))
                mod_lines.append((rid, m, mod_clone[(rid, m)], rpos[rid]))
                continue
            if m not in pos and m in pin_met and m not in foreign:
                px_, py_ = pin_met[m]
                pos[m] = (float(px_), float(py_))
                ctx[m] = "modifier"
            if m in pos and m not in foreign:
                if (rid, m) not in regulator_links and not line_ok(pos[m], rid, m):
                    mod_clone[(rid, m)] = mod_spot(rid)
                    mod_lines.append((rid, m, mod_clone[(rid, m)], rpos[rid]))
                else:
                    mod_lines.append((rid, m, pos[m], rpos[rid]))
            elif m in foreign:
                mod_clone[(rid, m)] = mod_spot(rid)
                mod_lines.append((rid, m, mod_clone[(rid, m)], rpos[rid]))
            else:
                pos[m] = mod_spot(rid)
                ctx[m] = "modifier"
                mod_lines.append((rid, m, pos[m], rpos[rid]))

    # second pass, now that every modifier line is known: a clone (or an enzyme node acting on one
    # reaction only) whose line crosses something moves to a better spot
    def crossing_mod(i_):
        r_, m_, a_, b_ = mod_lines[i_]
        return any(r1 != r_ and m1 != m_ and segs_cross(a_, b_, a1, b1)
                   for (r1, m1, a1, b1) in all_segs_m + mod_lines[:i_] + mod_lines[i_ + 1:])
    for _ in range(2):
        moved_any = False
        for i_ in range(len(mod_lines)):
            r_, m_, a_, b_ = mod_lines[i_]
            if not crossing_mod(i_):
                continue
            if Layout.sg(m_, r_) in prev.clones or m_ in pin_met:
                continue
            if (r_, m_) in mod_clone:
                mod_clone[(r_, m_)] = mod_spot(r_, skip=(r_, m_))
                mod_lines[i_] = (r_, m_, mod_clone[(r_, m_)], rpos[r_])
                moved_any = True
            elif ctx.get(m_) == "modifier" and sum(1 for x in mod_lines if x[1] == m_ and (x[0], m_) not in mod_clone) == 1:
                pos[m_] = mod_spot(r_, skip=m_)
                mod_lines[i_] = (r_, m_, pos[m_], rpos[r_])
                moved_any = True
        if not moved_any:
            break

    # main partners that ended up far away are not linked across the drawing: they become copies
    def far_of(rid):
        ro = roles[rid]
        far = set()
        for m in ro.main_sub + ro.main_prod:     # straight links through a circle (chords) -> copies
            if m not in pos or m in far:
                continue
            if rid in L.ring_rxn and m in L.ring_met and L.ring_rxn[rid][0] == L.ring_met[m][0]:
                cy = L.cycles[L.ring_rxn[rid][0]]
                if abs(abs(wrap(L.ring_met[m][1] - L.ring_rxn[rid][1])) - 2 * math.pi / len(cy.slots)) < 1e-3:
                    continue                      # neighbour on the ring: drawn as arc
            if P.split_chords and any(seg_point_dist(cy.center, rpos[rid], pos[m]) < 0.6 * cy.radius
                                      for cy in L.cycles):
                far.add(m)
        return far

    far_map = {rid: far_of(rid) for rid in sorted(rpos)}
    # very long, redundant lines: replace by a copy only if the reaction and the metabolite stay
    # connected through other drawn lines (never for a bridge such as a fold turn or a lone branch)
    if P.max_link > 0:
        Gd = nx.Graph()
        for rid in rpos:
            for m in roles[rid].main_sub + roles[rid].main_prod:
                if m in pos and m not in far_map[rid]:
                    Gd.add_edge(("R", rid), ("M", m))
        long_links = sorted(((dist(pos[m], rpos[rid]), rid, m) for rid in rpos
                             for m in roles[rid].main_sub + roles[rid].main_prod
                             if m in pos and m not in far_map[rid]), reverse=True)
        for d_, rid, m in long_links:
            if d_ <= P.max_link * D:
                break
            e = (("R", rid), ("M", m))
            Gd.remove_edge(*e)
            if nx.has_path(Gd, *e):
                far_map[rid].add(m)
            else:
                Gd.add_edge(*e)
    # a reaction whose partners on one side all became copies moves next to its remaining partner
    for rid in sorted(rpos):
        far = far_map[rid]
        if not far or rid in pin_rxn or rid in L.ring_rxn:
            continue
        ro = roles[rid]
        real = [m for m in ro.main_sub + ro.main_prod if m in pos and m not in far]
        if not real or (any(m in ro.main_sub for m in real) and any(m in ro.main_prod for m in real)):
            continue
        anchor = real[0]
        if dist(rpos[rid], pos[anchor]) <= 0.75 * D:
            continue
        best = None
        for d in (F, neg(F), S, neg(S)):
            cand = add(pos[anchor], mul(d, D / 2))
            v = sum(1 for q in list(pos.values()) + [p for r, p in rpos.items() if r != rid] if dist(cand, q) < 0.4 * D)
            v += sum(1 for a_, b_ in main_segments() if seg_point_dist(cand, a_, b_) < 15 and dist(a_, rpos[rid]) > 1e-6)
            if best is None or v < best[0]:
                best = (v, d, cand)
        rpos[rid] = best[2]
        raxis[rid] = best[1] if anchor in ro.main_sub else neg(best[1])

    for rid in sorted(rpos):
        far = far_map[rid]
        if far:
            ro = roles[rid]
            L.far_links[rid] = far
            roles[rid] = Roles([m for m in ro.main_sub if m not in far], [m for m in ro.main_prod if m not in far],
                               ro.side_sub + [m for m in ro.main_sub if m in far],
                               ro.side_prod + [m for m in ro.main_prod if m in far])

    # snap node positions to the export precision before deriving clones/edges from them,
    # so a layout re-read from SBML reproduces identical derived geometry
    for m in pos:
        pos[m] = (round(pos[m][0], 2), round(pos[m][1], 2))
    for r in rpos:
        rpos[r] = (round(rpos[r][0], 2), round(rpos[r][1], 2))

    # ---------------------------------------------------------------- side compounds (clones)
    base_segs = main_segments()

    def clone_slots(R: Pt, a: Pt, n: Pt, ro: Roles, far: set | frozenset = frozenset()):
        out = []
        for k, m in enumerate([m for m in ro.side_sub if m not in far]):
            out.append((m, "sidesubstrate",
                        add(R, add(mul(n, P.side_offset + 4 * k), mul(a, -(P.side_along + k * P.side_step))))))
        for k, m in enumerate([m for m in ro.side_prod if m not in far]):
            out.append((m, "sideproduct",
                        add(R, add(mul(n, P.side_offset + 4 * k), mul(a, P.side_along + k * P.side_step)))))
        return out

    def conflicts(rid, pts):
        score = 0
        others = list(pos.values()) + [p for r, p in rpos.items() if r != rid] + [c.pos for c in L.clones]
        for _, _, p in pts:
            score += 3 * hits_fixed(p, 4)
            score += sum(1 for q in others if dist(p, q) < 24)
            score += sum(1 for a, b in base_segs
                         if dist(a, rpos[rid]) > 1e-6 and seg_point_dist(p, a, b) < 14)
            # the hook from the reaction to the cofactor should not cross a main line
            score += 4 * sum(1 for a, b in base_segs
                             if dist(a, rpos[rid]) > 1e-6 and segs_cross(rpos[rid], p, a, b))
        return score

    for rid in sorted(rpos):
        ro, R, a = roles[rid], rpos[rid], raxis[rid]
        n0 = ring_side[rid] if rid in ring_side else par_side.get(rid, side_rot(a))
        far = L.far_links.get(rid, set())
        best = None
        for n in (n0, neg(n0)):
            pts = clone_slots(R, a, n, ro, far)
            old = [prev.clones[Layout.sg(m, rid)] for m, _, _ in pts if Layout.sg(m, rid) in prev.clones]
            if old:   # side used in the previous layout wins (consistency between versions)
                sc = sum(dist(o, p) for o, (_, _, p) in zip(old, pts))
            else:
                sc = conflicts(rid, pts)
            if best is None or sc < best[0]:
                best = (sc, n, pts)
        L.rxn_side[rid] = best[1]
        # a reaction kept where it was keeps its cofactor copies where they were, too (even if the
        # reaction's axis changed because new reactions made one of its partners a main compound)
        kept_here = rid in pin_rxn and dist(rpos[rid], pin_rxn[rid]) < 0.5
        for m, role, p in best[2]:
            g_id = Layout.sg(m, rid)
            if kept_here and g_id in prev.clones:
                p = (float(prev.clones[g_id][0]), float(prev.clones[g_id][1]))
            L.clones.append(Clone(g_id, m, rid, p, role, False))
        # copies of main partners (instead of long links): where the real node would be, i.e. on the
        # reaction axis at the normal distance, further copies in a row beside it
        n = best[1]
        for lst, role, sgn in ((ro.side_sub, "sidesubstrate", -1), (ro.side_prod, "sideproduct", 1)):
            j = 0
            for m in [m for m in lst if m in far]:
                offs = [0.0, -0.6, 0.6, -1.2, 1.2]
                p = None
                while j < len(offs):
                    cand = add(add(R, mul(a, sgn * D / 2)), mul(n, offs[j] * D))
                    j += 1
                    others = list(pos.values()) + list(rpos.values()) + [c.pos for c in L.clones]
                    if all(dist(cand, q) > 0.35 * D for q in others):
                        p = cand
                        break
                if p is None:
                    p = add(R, mul(a, sgn * D / 2))
                L.clones.append(Clone(Layout.sg(m, rid), m, rid, p, role, True))

    for (rid, m), p in sorted(mod_clone.items()):     # modifiers drawn as clones beside the reaction
        L.clones.append(Clone(Layout.sg(m, rid), m, rid, p, "modifier", True))
    # pseudo-reactions (biomass): a short arrow from each precursor to a small symbol right next to it
    # (into it for products), on the free side of the metabolite
    if drains:
        segs_d = segs_of(list(rpos))
        taken_d = list(pos.values()) + list(rpos.values()) + [c.pos for c in L.clones]
        for rid in sorted(drains):
            subs_d, prods_d = drains[rid]
            for m, role in [(x, "drain_in") for x in subs_d] + [(x, "drain_out") for x in prods_d]:
                if m not in pos:
                    continue
                g_id = Layout.sg(m, rid)
                if g_id in prev.clones:
                    px_, py_ = prev.clones[g_id]
                    p = (float(px_), float(py_))
                else:
                    best_ = None
                    for k_ in (0.55, 0.75, 1.0):
                        for i_ in range(16):
                            th = 2 * math.pi * i_ / 16
                            c_ = add(pos[m], (k_ * D * math.cos(th), k_ * D * math.sin(th)))
                            near = sum(1 for q in taken_d if dist(c_, q) < 0.35 * D)
                            cross = sum(1 for (r1, m1, a1, b1) in segs_d if m1 != m and segs_cross(pos[m], c_, a1, b1))
                            on_line = sum(1 for (r1, m1, a1, b1) in segs_d if seg_point_dist(c_, a1, b1) < 14)
                            # along a main line of the metabolite is where labels and cofactors go
                            along = sum(1 for (r1, m1, a1, b1) in segs_d if m1 == m and
                                        dot(unit(sub(c_, pos[m])), unit(sub(a1, b1) if dist(b1, pos[m]) < 1 else sub(b1, a1))) > 0.9)
                            cost = 100 * cross + 60 * near + 40 * on_line + 10 * along + k_ + 0.01 * i_
                            if best_ is None or cost < best_[0]:
                                best_ = (cost, c_)
                    assert best_ is not None
                    p = best_[1]
                taken_d.append(p)
                L.clones.append(Clone(g_id, m, rid, p, role, False))
    # ---------------------------------------------------------------- node overlaps (hard rule)
    resolve_node_overlaps(L, P, pos, rpos, raxis, ctx, pin_met, pin_rxn, side_rot, set(prev.clones),
                          grid_unit=P.grid * D)
    # ---------------------------------------------------------------- edges
    build_edges(L)
    # ---------------------------------------------------------------- labels
    place_labels(L, ctx, attach, S, side_rot, prev, lab_pref)
    return L


# ----------------------------------------------------------------------------- node overlaps
def resolve_node_overlaps(L: Layout, P: Params, pos, rpos, raxis, ctx, pin_met, pin_rxn, side_rot,
                          kept_clones: set | frozenset = frozenset(), gap: float = 4.0, grid_unit: float = 0.0):
    """No node may touch another. Overlaps are resolved by moving the most mobile of the two:
    cofactor copies first (to another free spot of their reaction's template), then copies of main
    metabolites, then reactions and branch metabolites (pushed apart). The main line, circles and
    pinned nodes only move if nothing else can."""
    D = P.spacing

    def items():
        out = [("M", m, P.met_radius) for m in sorted(pos)]
        out += [("R", r, P.rxn_half * 1.42) for r in sorted(rpos)]
        out += [("C", i, P.clone_radius + (1.5 if c.far else 0)) for i, c in enumerate(L.clones)]
        return out

    def getp(it):
        k, i, _ = it
        return pos[i] if k == "M" else rpos[i] if k == "R" else L.clones[i].pos

    def setp(it, p):
        k, i, _ = it
        p = (round(p[0], 2), round(p[1], 2))
        if k == "M":
            pos[i] = p
        elif k == "R":
            rpos[i] = p
        else:
            L.clones[i].pos = p

    def mobility(it):
        k, i, _ = it
        if k == "C":   # copies kept from a previous layout move last: new content gives way
            return 6 if L.clones[i].glyph in kept_clones else (0 if not L.clones[i].far else 1)
        if k == "R":
            return 6 if i in pin_rxn else 3 if i in L.ring_rxn else 2
        return 6 if i in pin_met else 5 if ctx.get(i) in ("ring", "backbone") else 4

    for _ in range(200):
        its = items()
        pts = [getp(t) for t in its]
        pair = None
        for i in range(len(its)):
            for j in range(i + 1, len(its)):
                if dist(pts[i], pts[j]) < its[i][2] + its[j][2] + gap:
                    pair = (its[i], its[j])
                    break
            if pair:
                break
        if pair is None:
            return
        a, b = pair
        ma, mb = mobility(a), mobility(b)
        # equally fixed (both from the previous layout): the small cofactor copy gives way, not the node
        if ma == mb and (a[0] == "C") != (b[0] == "C"):
            ma, mb = (ma - 0.5, mb) if a[0] == "C" else (ma, mb - 0.5)
        mover, other = (a, b) if ma <= mb else (b, a)

        def free(q, r):
            return all(dist(q, getp(t)) >= r + t[2] + gap for t in its if t != mover)

        if mover[0] == "C" and L.clones[mover[1]].role.startswith("drain"):
            c = L.clones[mover[1]]
            M_ = pos[c.species]
            cands = [add(M_, (k_ * P.spacing * math.cos(2 * math.pi * i_ / 16), k_ * P.spacing * math.sin(2 * math.pi * i_ / 16)))
                     for k_ in (0.55, 0.75, 1.0, 1.3) for i_ in range(16)]
            q = next((q for q in cands if free(q, mover[2])), None)
            if q is not None:
                setp(mover, q)
                continue
        if mover[0] == "C":
            c = L.clones[mover[1]]
            R_, a_ = rpos[c.rxn], raxis[c.rxn]
            n_: Pt = L.rxn_side[c.rxn] if c.rxn in L.rxn_side else side_rot(a_)
            sgn = -1 if c.role == "sidesubstrate" else 1
            cands = []
            if c.far:
                cands += [add(R_, add(mul(a_, sgn * D * f), mul(n_, o * D)))
                          for f in (0.5, 0.75, 1.0) for o in (0.0, -0.4, 0.4, -0.8, 0.8)]
            so, sa = P.side_offset, P.side_along
            for nn in (n_, neg(n_)):
                for off_n in (so, so + 10, so + 20, so + 32, so + 46):
                    for off_a in (sa, sa + 12, sa - 8, sa + 24, 0, sa + 36, sa + 48):
                        cands.append(add(R_, add(mul(nn, off_n), mul(a_, sgn * off_a))))
            q = next((q for q in cands if free(q, mover[2])), None)
            if q is not None:
                setp(mover, q)
                continue
        pa, pb = getp(mover), getp(other)
        v = unit(sub(pa, pb))
        if dot(v, v) == 0:
            v = (1.0, 0.0) if mover[0] != "R" else side_rot(raxis[mover[1]])
        if mover[0] in ("M", "R"):
            # keep alignment: move only horizontally or vertically, to the nearest spot free of all
            # nodes (searching all four directions avoids ping-pong between two neighbours)
            v0 = (math.copysign(1.0, v[0]), 0.0) if abs(v[0]) >= abs(v[1]) else (0.0, math.copysign(1.0, v[1]))
            dirs = [v0, (v0[1], v0[0]), (-v0[1], -v0[0]), neg(v0)]
            target = None
            unit_ = (grid_unit if mover[0] == "M" else grid_unit / 2) if grid_unit else 2.0
            for k in range(1, 60 if not grid_unit else 12):
                step = unit_ * k
                for d in dirs:
                    q = add(pa, mul(d, step))
                    if free(q, mover[2]):
                        target = q
                        break
                if target:
                    break
            setp(mover, target or add(pa, mul(v0, mover[2] + other[2] + gap)))
            continue
        need = mover[2] + other[2] + gap - dist(pa, pb) + 0.5
        setp(mover, add(pa, mul(v, need)))


# ----------------------------------------------------------------------------- edge geometry
SWAP = {"substrate": "product", "product": "substrate", "sidesubstrate": "sideproduct", "sideproduct": "sidesubstrate"}


def drain_label(name: str) -> str:
    return "biomass" if "biomass" in (name or "").lower() else short_text(name or "", 14)


MOD_ROLE = {"stimulation": "activator", "inhibition": "inhibitor", "catalysis": "modifier",
            "modulation": "modifier"}     # SBML Layout speciesReferenceGlyph roles


def build_edges(L: Layout):
    P, pw = L.params, L.pathway

    def true_role(rid, role, m):
        """Role as written in the model (layout roles are swapped for flipped reactions; far main
        partners are drawn as clones but keep their substrate/product role)."""
        if m in L.far_links.get(rid, ()):
            role = role.replace("side", "")
        return SWAP[role] if rid in L.flipped else role

    clone_at = {(c.rxn, c.species, c.role): c for c in L.clones}
    for c in L.clones:                         # pseudo-reaction drains: metabolite -> symbol (or back)
        if not c.role.startswith("drain") or c.species not in L.met_pos:
            continue
        mp = L.met_pos[c.species]
        tw = unit(sub(c.pos, mp))
        a_, b_ = add(mp, mul(tw, P.met_radius + 1)), add(c.pos, mul(tw, -(P.clone_radius + 3)))
        if c.role == "drain_in":
            L.edges.append(Edge(f"SRG_{c.rxn}_{c.species}_drain", c.rxn, c.species, Layout.sg(c.species), "substrate",
                                [("L", a_, b_)], 1.0, False, drawn="drain_in"))
        else:
            L.edges.append(Edge(f"SRG_{c.rxn}_{c.species}_drain", c.rxn, c.species, Layout.sg(c.species), "product",
                                [("L", b_, a_)], 1.0, False, drawn="drain_out"))
    for rid in sorted(L.rxn_pos):
        ro = L.roles[rid]
        r = pw.reactions[rid]
        R, a = L.rxn_pos[rid], L.rxn_axis[rid]
        p_in, p_out = add(R, mul(a, -P.rxn_half)), add(R, mul(a, P.rxn_half))
        for role, mets, port, sgn in (("substrate", ro.main_sub, p_in, -1), ("product", ro.main_prod, p_out, 1)):
            for m in mets:
                mp = L.met_pos[m]
                if m in r.modifiers:          # regulator standing in for an empty side: modulation arc
                    continue
                segs = None
                if rid in L.ring_rxn and m in L.ring_met and L.ring_rxn[rid][0] == L.ring_met[m][0]:
                    k, tr = L.ring_rxn[rid]
                    cyc = L.cycles[k]
                    tm = L.ring_met[m][1]
                    step = 2 * math.pi / len(cyc.slots)
                    if abs(abs(wrap(tm - tr)) - step) < 1e-3:
                        rr = cyc.radius
                        t0 = tr + sgn * cyc.direction * P.rxn_half / rr
                        t1 = tr + wrap(tm - tr) - sgn * cyc.direction * (P.met_radius + 1) / rr
                        segs = [arc_bezier(cyc.center, rr, t0, t0 + wrap(t1 - t0))]
                if segs is None:   # straight line from the reaction port to the metabolite
                    toward = unit(sub(mp, port))
                    if dot(toward, toward) == 0:
                        toward = mul(a, sgn)
                    segs = [("L", port, add(mp, mul(toward, -P.met_radius)))]
                tr = true_role(rid, role, m)
                L.edges.append(Edge(f"SRG_{rid}_{m}_{tr}", rid, m, Layout.sg(m), tr, segs,
                                    r.weight, r.reversible, drawn=role))
        # modifiers: a line from the regulator / enzyme to the reaction node, ending with the SBGN
        # symbol of its kind (catalysis, stimulation, inhibition, modulation)
        for m in sorted(r.modifiers):
            cl = clone_at.get((rid, m, "modifier"))
            if cl is not None:
                mp, rad, glyph = cl.pos, P.clone_radius + 2.5, cl.glyph
            elif m in L.met_pos:
                mp, rad, glyph = L.met_pos[m], P.met_radius + 1, Layout.sg(m)
            else:
                continue
            tw = unit(sub(R, mp))
            if dot(tw, tw) == 0:
                continue
            kind = r.modifiers[m]
            L.edges.append(Edge(f"SRG_{rid}_{m}_modifier", rid, m, glyph, MOD_ROLE.get(kind, "modifier"),
                                [("L", add(mp, mul(tw, rad)), add(R, mul(tw, -(P.rxn_half + 3))))],
                                1.0, False, drawn="mod:" + kind))
        n = L.rxn_side.get(rid, (0.0, 0.0))
        for role, mets, port, sgn in (("sidesubstrate", ro.side_sub, p_in, -1),
                                      ("sideproduct", ro.side_prod, p_out, 1)):
            for m in mets:
                c = clone_at[(rid, m, role)]
                if c.far:   # copy of a main partner: straight line like a normal main link
                    tw = unit(sub(c.pos, port))
                    tr = true_role(rid, role, m)
                    L.edges.append(Edge(f"SRG_{rid}_{m}_{tr}", rid, m, c.glyph, tr,
                                        [("L", port, add(c.pos, mul(tw, -(P.clone_radius + 1.5))))],
                                        r.weight, r.reversible, drawn=role, far=True))
                    continue
                c1 = add(port, mul(a, sgn * P.side_along * 0.8))
                back = unit(sub(add(c.pos, mul(n, -P.side_offset * 0.5)), c.pos))
                c2 = add(c.pos, mul(back, P.side_offset * 0.45))
                tr = true_role(rid, role, m)
                L.edges.append(Edge(f"SRG_{rid}_{m}_{tr}", rid, m, c.glyph, tr,
                                    [("C", port, c1, c2, add(c.pos, mul(back, P.clone_radius)))],
                                    r.weight, r.reversible, drawn=role, far=m in L.far_links.get(rid, ())))


# ----------------------------------------------------------------------------- labels
SLOT_DIRS = {"E": (1, 0), "W": (-1, 0), "N": (0, -1), "S": (0, 1),
             "NE": (0.71, -0.71), "NW": (-0.71, -0.71), "SE": (0.71, 0.71), "SW": (-0.71, 0.71)}
LABEL_SIZE = {"main": 12, "copy": 11, "side": 10, "reaction": 9}


def text_size(text: str, size: float) -> tuple[float, float]:
    # deliberately generous (Helvetica averages ~0.52 em): boxes never under-estimate the text
    return 0.58 * size * len(text) + 4, size + 4


def slot_box(slot: str, c: Pt, g: float, w: float, h: float) -> Pt:
    """top-left of a w x h box in direction `slot`, with gap g between node centre and box"""
    dx, dy = SLOT_DIRS[slot]
    x = c[0] + dx * g + (0 if dx > 0.1 else -w if dx < -0.1 else -w / 2)
    y = c[1] + dy * g + (0 if dy > 0.1 else -h if dy < -0.1 else -h / 2)
    return (x, y)


def slot_anchor(slot: str) -> str:
    dx = SLOT_DIRS[slot][0]
    return "start" if dx > 0.1 else "end" if dx < -0.1 else "middle"


def ordered_slots(pref: Pt) -> list[str]:
    return sorted(SLOT_DIRS, key=lambda s: (-dot(unit(SLOT_DIRS[s]), unit(pref)), s))


def place_labels(L: Layout, ctx, attach, S, side_rot, prev: Previous, lab_pref: dict | None = None):
    """Labels never overlap nodes or other labels (hard rule). Each label is searched ring by ring
    outward from its node; in the nearest ring with free positions the one crossing fewest edges wins
    (ties: preferred direction). Labels that end up away from their node get a leader line."""
    P, pw = L.params, L.pathway
    nodes = [(p, P.met_radius) for p in L.met_pos.values()]
    nodes += [(p, P.rxn_half * 1.42) for p in L.rxn_pos.values()]
    nodes += [(c.pos, P.clone_radius + (1.5 if c.far else 0)) for c in L.clones]
    edge_pts = [q for e in L.edges for s_ in e.segs for q in bezier_points(s_, 12)]
    placed: list[tuple] = []
    fixed_pending: dict[str, tuple] = {}
    PAD = 2.0

    def hits_node(box, own: Pt) -> bool:
        x0, y0, x1, y1 = box[0] - PAD, box[1] - PAD, box[0] + box[2] + PAD, box[1] + box[3] + PAD
        for p, r in nodes:
            if p is own or (abs(p[0] - own[0]) < 1e-6 and abs(p[1] - own[1]) < 1e-6):
                continue
            cx, cy = min(max(p[0], x0), x1), min(max(p[1], y0), y1)
            if math.hypot(cx - p[0], cy - p[1]) < r:
                return True
        return False

    def hits_label(box) -> bool:
        for o in placed + list(fixed_pending.values()):
            if not (box[0] + box[2] + PAD <= o[0] or o[0] + o[2] + PAD <= box[0] or
                    box[1] + box[3] + PAD <= o[1] or o[1] + o[3] + PAD <= box[1]):
                return True
        return False

    def edge_cost(box) -> int:
        return sum(1 for q in edge_pts if box[0] - 1 < q[0] < box[0] + box[2] + 1 and box[1] - 1 < q[1] < box[1] + box[3] + 1)

    def search(c: Pt, rr: float, w: float, h: float, pref: Pt, max_ring: int):
        """free positions only; cost = edge points under the label + distance (rings) + direction.
        A label one ring further out is preferred over one lying on a line."""
        order = ordered_slots(pref)
        best = None
        for ring in range(max_ring):
            g = rr + 3 + ring * 6
            for k, sl in enumerate(order):
                x0, y0 = slot_box(sl, c, g, w, h)
                box = (x0, y0, w, h)
                if hits_node(box, c) or hits_label(box):
                    continue
                cost = 2.0 * edge_cost(box) + 1.5 * ring + 0.2 * k
                if best is None or cost < best[0]:
                    best = (cost, x0, y0, sl, ring)
            if best is not None and best[0] < 1.5 * (ring + 1):   # nothing better can come further out
                break
        return best

    def put(tg, target, species, rxn, text, kind, c: Pt, rr, pref: Pt, fixed=None, shorter=None):
        size = LABEL_SIZE[kind]
        w, h = text_size(text, size)
        fixed_pending.pop(tg, None)
        anchor, leader = "start", False
        keep = tg in prev.labels and dist(add(prev.labels[tg], (w / 2, h / 2)), c) < rr + w / 2 + 30
        if keep:   # previous label still next to its object -> keep it exactly
            x, y = prev.labels[tg]
            # it had a leader line before if it sits clearly away from its node: keep that line too
            gx = max(x - c[0], 0.0, c[0] - (x + w))
            gy = max(y - c[1], 0.0, c[1] - (y + h))
            leader = math.hypot(gx, gy) > rr + 3 + 3 * 6 - 2       # same as ring >= 3 below
        elif isinstance(fixed, (list, tuple)):
            x, y = fixed
        elif isinstance(fixed, str) and fixed in SLOT_DIRS:
            x, y = slot_box(fixed, c, rr + 3, w, h)
            anchor = slot_anchor(fixed)
        else:
            best = search(c, rr, w, h, pref, 4)                 # close to the node
            if best is None and shorter and shorter != text:     # long enzyme name: try a shorter form
                w2, h2 = text_size(shorter, size)
                best2 = search(c, rr, w2, h2, pref, 4)
                if best2 is not None:
                    text, w, h, best = shorter, w2, h2, best2
            if best is None:
                best = search(c, rr, w, h, pref, 40)             # further out (up to ~250 px)
            if best is None:                                     # last resort: below everything
                ys = [o[1] + o[3] for o in placed] + [p[1] + r for p, r in nodes]
                best = (0, c[0] - w / 2, max(ys) + 6, "S", 99)
            _, x, y, sl, ring = best
            anchor = slot_anchor(sl)
            leader = ring >= 3
        placed.append((x, y, w, h))
        L.labels[tg] = Label(tg, target, species, rxn, text, x, y, w, h, size, kind, anchor, leader, c)

    # labels kept from a previous layout block their area for all new labels
    for m in L.met_pos:
        tg = f"TG_M_{m}"
        if tg in prev.labels:
            w, h = text_size(display_name(pw, m, P), LABEL_SIZE["main"])
            fixed_pending[tg] = (*prev.labels[tg], w, h)
    for c in L.clones:
        tg = "TG" + c.glyph[2:]
        if tg in prev.labels:
            w, h = text_size(display_name(pw, c.species, P, clone=not c.far), LABEL_SIZE["copy" if c.far else "side"])
            fixed_pending[tg] = (*prev.labels[tg], w, h)
    for rid in L.rxn_pos:
        tg = f"TG_R_{rid}"
        if tg in prev.labels and P.reaction_labels:
            w, h = text_size(short_text(pw.reactions[rid].name, P.label_max), LABEL_SIZE["reaction"])
            fixed_pending[tg] = (*prev.labels[tg], w, h)

    # order: main metabolites, reactions, copies of main metabolites, cofactors
    for m in sorted(L.met_pos, key=lambda m: (0 if ctx.get(m) != "branch" else 1, m)):
        p = L.met_pos[m]
        if lab_pref and m in lab_pref:      # partner rows: label away from the reaction
            pref = lab_pref[m]
        elif ctx.get(m) == "ring" and m in L.ring_met:
            pref = neg(attach[m])
        elif ctx.get(m) == "backbone":
            pref = S
        else:
            pref = neg(side_rot(attach.get(m, (0, 1))))
        put(f"TG_M_{m}", Layout.sg(m), m, None, display_name(pw, m, P), "main", p, P.met_radius, pref,
            pw.labels.get(m))
    if P.reaction_labels:
        for rid in sorted(L.rxn_pos):
            n = L.rxn_side.get(rid, (1, 0))
            name = pw.reactions[rid].name
            put(f"TG_R_{rid}", Layout.rg(rid), None, rid, short_text(name, P.label_max), "reaction",
                L.rxn_pos[rid], P.rxn_half + 1, neg(n), pw.labels.get(rid), shorter=short_text(name, 18))
    for c in sorted(L.clones, key=lambda c: (not c.far, c.glyph)):
        if c.species == EMPTY:          # source/sink symbol: no label
            continue
        if c.role.startswith("drain"):  # pseudo-reaction symbol: "biomass"
            r_ = pw.reactions.get(c.rxn)
            put("TG" + c.glyph[2:], c.glyph, None, c.rxn, drain_label(r_.name if r_ else c.rxn),
                "side", c.pos, P.clone_radius + 2, unit(sub(c.pos, L.met_pos[c.species])))
            continue
        n = L.rxn_side.get(c.rxn, (1, 0))
        put("TG" + c.glyph[2:], c.glyph, c.species, None, display_name(pw, c.species, P, clone=not c.far),
            "copy" if c.far else "side", c.pos, P.clone_radius + (1.5 if c.far else 0), n)


# ----------------------------------------------------------------------------- finalize
def _transform(L: Layout, d: Pt):
    for m in L.met_pos:
        L.met_pos[m] = add(L.met_pos[m], d)
    for r in L.rxn_pos:
        L.rxn_pos[r] = add(L.rxn_pos[r], d)
    for c in L.clones:
        c.pos = add(c.pos, d)
    for cy in L.cycles:
        cy.center = add(cy.center, d)
    for e in L.edges:
        e.segs = [(s[0],) + tuple(add(p, d) for p in s[1:]) for s in e.segs]
    for lb in L.labels.values():
        lb.x, lb.y = lb.x + d[0], lb.y + d[1]
        lb.node = add(lb.node, d)


def _bbox(L: Layout):
    P = L.params
    pts = []
    for p in L.met_pos.values():
        pts += [sub(p, (P.met_radius, P.met_radius)), add(p, (P.met_radius, P.met_radius))]
    pts += list(L.rxn_pos.values()) + [c.pos for c in L.clones]
    pts += [q for e in L.edges for s in e.segs for q in bezier_points(s, 6)]
    for lb in L.labels.values():
        pts += [(lb.x, lb.y), (lb.x + lb.w, lb.y + lb.h)]
    xs, ys = [p[0] for p in pts], [p[1] for p in pts]
    return min(xs), min(ys), max(xs), max(ys)


def finalize(L: Layout, keep_frame: bool, parent: str | None, grid: float = 0.0):
    P = L.params
    x0, y0, _, _ = _bbox(L)
    if keep_frame:   # keep previous coordinates; shift only if something left the canvas
        d = (max(0.0, P.margin - x0), max(0.0, P.margin - y0))
    else:
        d = (P.margin - x0, P.margin - y0)
    if grid > 0:     # whole grid cells: the drawing keeps its grid (positions stay multiples of it)
        d = (math.ceil(d[0] / grid) * grid, math.ceil(d[1] / grid) * grid)
    d = (float(math.ceil(d[0])) if not grid else d[0], float(math.ceil(d[1])) if not grid else d[1])
    if d != (0.0, 0.0):
        _transform(L, d)
    for m in L.met_pos:
        L.met_pos[m] = (round(L.met_pos[m][0], 2), round(L.met_pos[m][1], 2))
    for r in L.rxn_pos:
        L.rxn_pos[r] = (round(L.rxn_pos[r][0], 2), round(L.rxn_pos[r][1], 2))
    _, _, x1, y1 = _bbox(L)
    L.width, L.height = round(x1 + P.margin, 2), round(y1 + P.margin, 2)
    coords = {
        "m": {k: [round(v[0], 2), round(v[1], 2)] for k, v in sorted(L.met_pos.items())},
        "r": {k: [round(v[0], 2), round(v[1], 2)] for k, v in sorted(L.rxn_pos.items())},
        "l": {k: [round(v.x, 2), round(v.y, 2)] for k, v in sorted(L.labels.items())},
    }
    L.meta = {
        "algorithm": ALGORITHM_VERSION,
        "pathway_hash": L.pathway.content_hash(),
        "params_hash": sha(canonical(asdict(P))),
        "layout_hash": sha(canonical(coords)),
        "parent_layout_hash": parent or "",
        "shift": [round(d[0], 2), round(d[1], 2)],
    }
