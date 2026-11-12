"""Canonical cycles that are drawn as circles (Params.cycles = "known", the default).

A detected loop is drawn as a circle only if it contains EVERY member listed in CANONICAL_CYCLES
(so e.g. succinyl-CoA anchors the citric acid cycle). Members are recognised by BiGG id, KEGG
compound id or one of the listed names (case, stereo prefixes and charges ignored). If several
loops qualify, the one with the fewest extra metabolites is used. All other loops are laid out
like normal paths.
"""
from __future__ import annotations

import re

from .model import Metabolite


def norm_name(s: str) -> str:
    s = s.lower()
    s = re.sub(r"\[[^\]]*\]", "", s)                    # [compartment]
    s = re.sub(r"\((\d*[+\-−]|[+\-−]\d*)\)", "", s)      # charge "(2-)"
    s = re.sub(r"^(\((s|r|2s|2r|2r,3s|1r,2s)\)|[ld])-", "", s.strip())
    s = re.sub(r"ic acid\b", "ate", s)                   # "succinic acid" = "succinate"
    return re.sub(r"[^a-z0-9]", "", s)


def norm_id(s: str) -> str:
    """BiGG-like id without compartment ('argsuc[c]'), with '_L' the same as '__L'"""
    s = re.sub(r"\[[^\]]*\]$", "", s.lower())
    return re.sub(r"_+", "_", s)


# Canonical cycles drawn as circles: (BiGG id, KEGG compound, names). A loop is drawn as one of these
# only if it contains EVERY listed member.
CANONICAL_CYCLES = {
    "citric acid cycle": (
        ("oaa", "C00036", ("oxaloacetate", "oxalacetate")),
        ("cit", "C00158", ("citrate", "citric acid")),
        ("icit", "C00311", ("isocitrate", "d-threo-isocitrate", "isocitric acid")),
        ("akg", "C00026", ("2-oxoglutarate", "alpha-ketoglutarate", "2-oxoglutaric acid",
                           "akg", "oxoglutarate")),
        # An anchor, not a gap: every 2-oxoglutarate-dependent dioxygenase
        # turns akg into succinate, and with succinyl-CoA optional Recon3D's
        # ring closed through a prolyl hydroxylase instead of the TCA cycle.
        ("succoa", "C00091", ("succinyl-coa", "succinyl coenzyme a")),
        ("succ", "C00042", ("succinate", "succinic acid")),
        ("fum", "C00122", ("fumarate", "fumaric acid")),
        ("mal__L", "C00149", ("l-malate", "malate", "(s)-malate")),
    ),
    "urea cycle": (
        ("orn", "C00077", ("ornithine", "l-ornithine")),
        ("citr__L", "C00327", ("citrulline", "l-citrulline")),
        ("argsuc", "C03406", ("argininosuccinate", "arginosuccinate", "n(omega)-(l-arginino)succinate",
                              "l-argininosuccinate")),
        ("arg__L", "C00062", ("arginine", "l-arginine")),
    ),
    # Calvin-Benson cycle: carboxylation, reduction, regeneration. The regeneration part (aldolase,
    # FBPase/SBPase, transketolase, isomerase/epimerase) differs between models, so only the members
    # every version has are required; whatever lies between GAP and Ru5P is part of the ring.
    "calvin cycle": (
        ("rb15bp", "C01182", ("ribulose 1,5-bisphosphate", "ribulose-1,5-bisphosphate", "ribulose bisphosphate",
                              "ribulose 1,5-diphosphate")),
        ("3pg", "C00197", ("3-phospho-d-glycerate", "3-phosphoglycerate", "3-phospho-glycerate",
                           "glycerate 3-phosphate")),
        ("13dpg", "C00236", ("3-phospho-d-glyceroyl phosphate", "1,3-bisphosphoglycerate",
                             "1,3-bisphospho-d-glycerate", "3-phospho-d-glyceroyl-phosphate")),
        ("g3p", "C00118", ("glyceraldehyde 3-phosphate", "d-glyceraldehyde 3-phosphate",
                           "glyceraldehyde-3-phosphate")),
        ("ru5p__D", "C00199", ("ribulose 5-phosphate", "d-ribulose 5-phosphate", "ribulose-5-phosphate")),
    ),
    "methionine cycle": (
        ("met__L", "C00073", ("methionine", "l-methionine")),
        ("amet", "C00019", ("s-adenosyl-l-methionine", "s-adenosylmethionine")),
        ("ahcys", "C00021", ("s-adenosyl-l-homocysteine", "s-adenosylhomocysteine")),
        ("hcys__L", "C00155", ("homocysteine", "l-homocysteine")),
    ),
}
# Reaction between these two members is put at 12 o'clock when the cycle is drawn on its own.
CYCLE_TOP = {"citric acid cycle": ("oaa", "cit"), "urea cycle": ("orn", "citr__L"),
             "methionine cycle": ("met__L", "amet"), "calvin cycle": ("rb15bp", "3pg")}


# Canonical main paths: if a map contains (most of) such a chain, it becomes the main line (backbone)
# instead of the longest chain, e.g. glycolysis runs from D-glucose to pyruvate even when the galactose
# (Leloir) entry into G6P is longer. Members are matched like cycle members; at least
# MIN_PATH_MEMBERS of them must lie on one connected chain. Add your own pathways here.
CANONICAL_PATHS = {
    "glycolysis": (
        ("glc__D", "C00031", ("d-glucose", "glucose", "alpha-d-glucose", "beta-d-glucose")),
        ("g6p", "C00092", ("d-glucose 6-phosphate", "glucose 6-phosphate", "glucose-6-phosphate",
                           "alpha-d-glucose 6-phosphate", "beta-d-glucose 6-phosphate")),
        ("f6p", "C00085", ("d-fructose 6-phosphate", "fructose 6-phosphate", "fructose-6-phosphate",
                           "beta-d-fructose 6-phosphate")),
        ("fdp", "C00354", ("d-fructose 1,6-bisphosphate", "fructose 1,6-bisphosphate",
                           "beta-d-fructose 1,6-bisphosphate", "fructose-1,6-bisphosphate")),
        ("g3p", "C00118", ("glyceraldehyde 3-phosphate", "d-glyceraldehyde 3-phosphate",
                           "glyceraldehyde-3-phosphate")),
        ("13dpg", "C00236", ("3-phospho-d-glyceroyl phosphate", "1,3-bisphosphoglycerate",
                             "1,3-bisphospho-d-glycerate", "3-phospho-d-glyceroyl-phosphate")),
        ("3pg", "C00197", ("3-phospho-d-glycerate", "3-phosphoglycerate")),
        ("2pg", "C00631", ("2-phospho-d-glycerate", "2-phosphoglycerate")),
        ("pep", "C00074", ("phosphoenolpyruvate",)),
        ("pyr", "C00022", ("pyruvate", "pyruvic acid")),
    ),
}
MIN_PATH_MEMBERS = 5


def _member_keys(bigg: str, kegg: str, names) -> set[str]:
    return {f"id:{bigg.lower()}", f"id:{norm_id(bigg)}", f"kegg:{kegg}"} | {f"name:{norm_name(n)}" for n in names}


def met_keys(m: Metabolite) -> set[str]:
    keys = {f"id:{m.base.lower()}", f"id:{m.id.lower()}", f"name:{norm_name(m.name)}"}
    keys |= {f"id:{re.sub(r'_[a-z]$', '', x.lower())}" for x in (m.id, m.base)}   # BiGG compartment suffix
    keys |= {f"id:{norm_id(re.sub(r'_[a-z]$', '', x))}" for x in (m.id, m.base)} | \
        {f"id:{norm_id(x)}" for x in (m.id, m.base)}                               # "argsuc[c]", "citr_L"
    keys |= {f"chebi:{c}" for c in m.chebi} | {f"kegg:{k}" for k in m.kegg}
    return keys


def match_members(mets: dict[str, Metabolite], members) -> dict[str, str]:
    """metabolite id -> BiGG id of the canonical member it represents"""
    out = {}
    for mid, m in mets.items():
        ks = met_keys(m)
        for bigg, kegg, names in members:
            if ks & _member_keys(bigg, kegg, names):
                out[mid] = bigg
                break
    return out


def select_known_cycles(raw_cycles: list[tuple], mets: dict[str, Metabolite]) -> list[tuple[tuple, dict]]:
    """For each canonical cycle: among the detected loops that contain ALL its members, take the one
    with the fewest extra metabolites (e.g. cis-aconitate may sit between citrate and isocitrate)."""
    chosen: list[tuple[tuple, dict]] = []
    used: set[str] = set()
    for name, members in CANONICAL_CYCLES.items():
        member_of = match_members(mets, members)
        need = {m[0] for m in members}
        best = None
        for c in raw_cycles:
            if used & set(c):
                continue
            hit = {member_of[m] for m in c if m in member_of}
            if hit != need:
                continue
            extra = len(c) - sum(1 for m in c if m in member_of)
            if best is None or (extra, c) < best[0]:
                best = ((extra, c), c)
        if best:
            c = best[1]
            chosen.append((c, {"name": name, "top": CYCLE_TOP.get(name),
                               "member_of": {m: member_of[m] for m in c if m in member_of}}))
            used |= set(c)
    return chosen


def user_path_members(spec: list[str]):
    """--main-path A,B,C: members given by id or name (in this order)"""
    return tuple((x.strip(), "", (x.strip(),)) for x in spec if x.strip())


def find_main_path(H, sup: dict[str, str], mets: dict[str, Metabolite], user_spec=None):
    """Main line from a canonical (or user-given) chain of metabolites.

    H: contracted skeleton DAG (cycle members replaced by @cycK), sup: metabolite -> H node.
    Returns (name, chain) with chain a list of H nodes from the first to the last member found,
    connected through H; or None. Members that are missing or not reachable are skipped.
    """
    import networkx as nx
    tables = {"main path": user_path_members(user_spec)} if user_spec else CANONICAL_PATHS
    best = None
    for name, members in tables.items():
        member_of = match_members(mets, members)
        order = [m[0] for m in members]
        cands = {b: sorted({sup[m] for m, bb in member_of.items() if bb == b and m in sup
                            and not sup[m].startswith("@cyc")}) for b in order}
        present = [b for b in order if cands[b]]
        need = 2 if user_spec else MIN_PATH_MEMBERS
        if len(present) < need:
            continue
        starts = [(si, c) for si, b in enumerate(present) for c in cands[b]]
        for si, start in starts:              # a member may be missing or not lead on (e.g. no hexokinase)
            chain, hits, cur = [start], 1, start
            for b in present[si + 1:]:
                dist = nx.single_source_shortest_path_length(H, cur)
                reach = [(dist[c], c) for c in cands[b] if c in dist and c not in chain]
                if not reach:
                    continue
                _, c = min(reach)
                chain += nx.shortest_path(H, cur, c)[1:]
                cur, hits = c, hits + 1
            key = (hits, -si, len(chain))      # most members, then the earliest start
            if hits >= need and (best is None or key > best[0]):
                best = (key, name, chain)
    if best is None:
        return None
    return best[1], best[2]
