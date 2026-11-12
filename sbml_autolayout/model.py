"""Pathway data model, input loading and role classification (main vs. side compounds)."""
from __future__ import annotations

import hashlib
import re
import json
from dataclasses import dataclass, field

# Currency / cofactor metabolites that are drawn as per-reaction clones by default.
DEFAULT_CURRENCY = {
    "atp", "adp", "amp", "gtp", "gdp", "gmp", "utp", "udp", "ctp", "cdp",
    "nad", "nadh", "nadp", "nadph", "fad", "fadh2", "q", "qh2",
    "h2o", "h", "co2", "o2", "pi", "ppi", "nh4", "nh3", "coa", "hco3", "h2o2",
}


# Name-based fallback so cofactors are recognised in models with arbitrary ids (e.g. "M_s_0434").
CURRENCY_NAMES = {
    "atp", "adp", "amp", "gtp", "gdp", "utp", "udp", "ctp", "cdp", "nad", "nad+", "nadh", "nadp", "nadp+",
    "nadph", "fad", "fadh2", "ubiquinone", "ubiquinol", "h2o", "water", "h+", "h", "proton", "co2",
    "carbon dioxide", "o2", "oxygen", "phosphate", "orthophosphate", "pi", "diphosphate", "pyrophosphate",
    "ppi", "ammonium", "nh4+", "nh3", "ammonia", "coenzyme a", "coa", "hydrogen peroxide", "bicarbonate",
}


# ChEBI ids of cofactors / currency metabolites (charged and neutral forms).
CURRENCY_CHEBI = {
    "30616", "15422", "456216", "16761", "456215", "16027", "37565", "15996", "58189", "17552",
    "46398", "15713", "58223", "17659", "37563", "17677", "58069", "17239", "57540", "15846",
    "57945", "16908", "58349", "18009", "57783", "16474", "57692", "16238", "58307", "17877",
    "15377", "15378", "24636", "16526", "15379", "43474", "18367", "26078", "33019", "29888",
    "18361", "28938", "16134", "57287", "15346", "17544", "16240", "16389", "17976",
}
# Short display names for cofactor clones.
SHORT_NAMES = {
    "atp": "ATP", "adp": "ADP", "amp": "AMP", "gtp": "GTP", "gdp": "GDP", "utp": "UTP", "udp": "UDP",
    "ctp": "CTP", "cdp": "CDP", "nad": "NAD+", "nadh": "NADH", "nadp": "NADP+", "nadph": "NADPH",
    "fad": "FAD", "fadh2": "FADH2", "h2o": "H2O", "h": "H+", "co2": "CO2", "o2": "O2", "pi": "Pi",
    "ppi": "PPi", "nh4": "NH4+", "nh3": "NH3", "coa": "CoA", "hco3": "HCO3-", "h2o2": "H2O2",
    "q": "Q", "qh2": "QH2",
}
CHEBI_SHORT = {
    "30616": "ATP", "456216": "ADP", "456215": "AMP", "37565": "GTP", "58189": "GDP", "46398": "UTP",
    "58223": "UDP", "57540": "NAD+", "57945": "NADH", "58349": "NADP+", "57783": "NADPH", "57692": "FAD",
    "58307": "FADH2", "15377": "H2O", "15378": "H+", "16526": "CO2", "15379": "O2", "43474": "Pi",
    "33019": "PPi", "28938": "NH4+", "57287": "CoA", "17544": "HCO3-",
}


@dataclass
class Metabolite:
    id: str
    name: str
    compartment: str = "c"
    base: str = ""          # id without compartment suffix/prefix, used for cofactor matching
    chebi: tuple = ()       # ChEBI ids from annotations
    kegg: tuple = ()        # KEGG compound ids from annotations
    comp_label: str = ""    # short compartment tag shown in labels ("m", "e", ...)
    sbo: str = ""           # SBO term of the species (decides the SBGN glyph: simple chemical, protein, ...)
    sbml_id: str = ""       # id in the input SBML file (written back on export)

    def short(self) -> str:
        for c in self.chebi:
            if c in CHEBI_SHORT:
                return CHEBI_SHORT[c]
        return SHORT_NAMES.get(self.base.lower(), self.name)


@dataclass
class Reaction:
    id: str
    name: str
    substrates: dict[str, float]
    products: dict[str, float]
    reversible: bool = False
    weight: float = 1.0
    main: set[str] = field(default_factory=set)   # force these participants to be main
    side: set[str] = field(default_factory=set)   # force these participants to be side
    main_explicit: bool = False                   # main given by the source: everything else is side
    role_source: str = "rules"                    # where the main/side decision came from
    # species that influence the reaction without being consumed: id -> kind
    # ("catalysis", "stimulation", "inhibition", "modulation"; SBML modifiers / SBGN arcs)
    modifiers: dict[str, str] = field(default_factory=dict)
    # pseudo-reaction (biomass objective, SBO:0000629): not a biochemical step; drawn as small drains
    # next to its precursors instead of one node with lines across the whole map
    pseudo: bool = False
    sbml_id: str = ""             # id in the input SBML file (written back on export)
    # identical reactions of the same input file (same participants): drawn on top of this one, but
    # written to the output with their own ids: [{"id", "name", "substrates", "products", "reversible"}]
    aliases: list = field(default_factory=list)


@dataclass
class Pathway:
    id: str
    name: str
    metabolites: dict[str, Metabolite]
    reactions: dict[str, Reaction]
    currency: set[str]
    hidden: set[str]
    labels: dict = field(default_factory=dict)    # metabolite id -> slot ("N","E",..) or [x, y]
    pins: dict = field(default_factory=dict)      # metabolite/reaction id -> [x, y]
    raw: dict = field(default_factory=dict)

    def is_currency(self, mid: str) -> bool:
        m = self.metabolites[mid]
        if mid in self.currency or m.base.lower() in self.currency:
            return True
        if self.raw.get("currency_replace"):
            return False
        return m.name.strip().lower() in CURRENCY_NAMES or any(c in CURRENCY_CHEBI for c in m.chebi)

    def is_hidden(self, mid: str) -> bool:
        m = self.metabolites[mid]
        return mid in self.hidden or m.base.lower() in self.hidden or (
            "h" in self.hidden and ("15378" in m.chebi or "24636" in m.chebi))

    def content_hash(self) -> str:
        # the SBO terms only change how a species is drawn, not what the pathway is
        raw = dict(self.raw)
        if isinstance(raw.get("metabolites"), dict):
            raw["metabolites"] = {k: {a: b for a, b in v.items() if a not in ("sbo", "sbml_id")}
                                  if isinstance(v, dict) else v for k, v in raw["metabolites"].items()}
        if isinstance(raw.get("reactions"), list):     # original ids / duplicates: same pathway
            raw["reactions"] = [{a: b for a, b in r.items() if a not in ("sbml_id", "aliases")}
                                if isinstance(r, dict) else r for r in raw["reactions"]]
        return sha(canonical(raw))


def canonical(obj) -> str:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"))


def sha(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()[:16]


def load_pathway(path: str) -> Pathway:
    with open(path) as fh:
        return pathway_from_dict(json.load(fh))


def pathway_from_dict(raw: dict) -> Pathway:
    mets = {
        mid: Metabolite(mid, m.get("name", mid), m.get("compartment", "c"), m.get("base", mid),
                        chebi=tuple(m.get("chebi", [])), kegg=tuple(m.get("kegg", [])),
                        comp_label=m.get("comp_label", ""), sbo=m.get("sbo", ""),
                        sbml_id=m.get("sbml_id", ""))
        for mid, m in raw["metabolites"].items()
    }
    rxns = {}
    for r in raw["reactions"]:
        for mid in list(r["substrates"]) + list(r["products"]) + list(r.get("modifiers", {})):
            if mid not in mets:
                mets[mid] = Metabolite(mid, mid, "c", mid)
        rxns[r["id"]] = Reaction(
            id=r["id"], name=r.get("name", r["id"]),
            substrates=dict(r["substrates"]), products=dict(r["products"]),
            reversible=bool(r.get("reversible", False)),
            weight=float(r.get("weight", 1.0)),
            main=set(r.get("main", [])), side=set(r.get("side", [])),
            main_explicit=bool(r.get("main_explicit", False)),
            role_source=r.get("role_source", "explicit (JSON)" if r.get("main") or r.get("side") else "rules"),
            modifiers={m: k for m, k in r.get("modifiers", {}).items()
                       if m not in r["substrates"] and m not in r["products"]},
            pseudo=bool(r.get("pseudo", False)),
            sbml_id=r.get("sbml_id", ""), aliases=list(r.get("aliases", [])),
        )
    currency = set(raw.get("currency", [])) if raw.get("currency_replace") else (
        DEFAULT_CURRENCY | set(raw.get("currency", [])))
    return Pathway(
        id=raw.get("id", "pathway"), name=raw.get("name", ""), metabolites=mets,
        reactions=dict(sorted(rxns.items())), currency=currency,
        hidden=set(raw.get("hide", [])), labels=raw.get("labels", {}),
        pins=raw.get("pins", {}), raw=raw,
    )


@dataclass
class Roles:
    main_sub: list[str]
    main_prod: list[str]
    side_sub: list[str]
    side_prod: list[str]


def classify(pw: Pathway, hub_degree: int = 4, hub_ratio: float = 2.0) -> dict[str, Roles]:
    """Split each reaction's participants into main (backbone) and side (cloned) compounds.

    1. cofactors (id, name or ChEBI annotation) are side compounds
    2. hub demotion: if a reaction side has several main candidates, a candidate that takes part in
       >= hub_degree reactions and in >= hub_ratio x as many as the other candidates becomes a side
       compound in that reaction (e.g. acetyl-CoA released by every thiolase step of beta-oxidation).
       At least one main compound per side is kept.
    3. connectivity: a demotion is only applied if it does not cut the network apart, i.e. the reaction
       stays connected to the hub through other main links, or the hub would merely end up on its own
       (then it is drawn as copies everywhere). Example: G1P stays main in glycogen phosphorylase,
       because that link is what connects glycogen metabolism to the rest of the map.
    """
    import networkx as nx

    degree: dict[str, int] = {}
    for r in pw.reactions.values():
        for m in set(r.substrates) | set(r.products):
            if not pw.is_currency(m):
                degree[m] = degree.get(m, 0) + 1

    def is_main(r, m):
        if m in r.main:
            return True
        if m in r.side or r.main_explicit:
            return False
        return not pw.is_currency(m)

    # candidate demotions per reaction (as before), and the bipartite network of main links
    cands_all: list[tuple] = []
    G = nx.Graph()
    for rid, r in pw.reactions.items():
        G.add_node(("R", rid))
        for m in dict.fromkeys(list(r.substrates) + list(r.products)):
            if is_main(r, m):
                G.add_edge(("R", rid), ("M", m))
        if r.main_explicit:          # the source said which compounds are main: no hub rule
            continue
        for side in (r.substrates, r.products):
            cands = [m for m in side if is_main(r, m)]
            while len(cands) > 1:
                top = max(cands, key=lambda m: (degree.get(m, 0), m))
                rest = max(degree.get(m, 0) for m in cands if m != top)
                if top in r.main or degree.get(top, 0) < hub_degree or degree.get(top, 0) < hub_ratio * rest:
                    break
                cands_all.append((-degree.get(top, 0), rid, top))
                cands.remove(top)

    demoted: dict[str, set] = {}
    for _, rid, m in sorted(cands_all):
        e = (("R", rid), ("M", m))
        G.remove_edge(*e)
        if nx.has_path(G, *e) or len(nx.node_connected_component(G, ("M", m))) == 1:
            demoted.setdefault(rid, set()).add(m)
        else:
            G.add_edge(*e)      # would split the map: keep it as a main link

    out = {}
    for rid, r in pw.reactions.items():
        dm = demoted.get(rid, set())
        ms = [m for m in r.substrates if is_main(r, m) and m not in dm]
        mp = [m for m in r.products if is_main(r, m) and m not in dm]
        ss = [m for m in r.substrates if (not is_main(r, m) or m in dm) and not pw.is_hidden(m)]
        sp = [m for m in r.products if (not is_main(r, m) or m in dm) and not pw.is_hidden(m)]
        out[rid] = Roles(ms, mp, ss, sp)
    return out


# placeholder species for the empty side of a synthesis (-> X) or degradation (X ->) reaction; drawn as
# the SBGN source/sink symbol, never written to SBML/Escher as a species
EMPTY = "__empty__"

# SBO terms of modifier species references -> kind of influence (SBGN arc)
MODIFIER_SBO = {
    "SBO:0000013": "catalysis", "SBO:0000460": "catalysis", "SBO:0000172": "catalysis",
    "SBO:0000459": "stimulation", "SBO:0000461": "stimulation", "SBO:0000021": "stimulation",
    "SBO:0000170": "stimulation", "SBO:0000171": "stimulation", "SBO:0000636": "stimulation",
    "SBO:0000020": "inhibition", "SBO:0000206": "inhibition", "SBO:0000207": "inhibition",
    "SBO:0000169": "inhibition", "SBO:0000537": "inhibition", "SBO:0000638": "inhibition",
    "SBO:0000019": "modulation", "SBO:0000594": "modulation", "SBO:0000168": "modulation",
}

# placeholder species for the drain symbol of a pseudo-reaction (biomass); never written as a species
DRAIN = "__drain__"
PSEUDO_SBO = {"SBO:0000629", "SBO:0000630", "SBO:0000631", "SBO:0000632"}   # biomass production etc.


def looks_pseudo(rid: str, name: str) -> bool:
    return "biomass" in rid.lower() or "biomass" in (name or "").lower()


SID = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


def export_ids(pw: "Pathway") -> tuple[dict, dict]:
    """ids used in the output files: the original SBML ids where the input had them (so the result
    can be compared with the input as is), else M_<id> / R_<id>. Always valid and unique."""
    used: set[str] = set()

    def pick(orig: str, prefix: str, internal: str) -> str:
        cand = orig if orig and SID.match(orig) else prefix + re.sub(r"[^A-Za-z0-9_]", "_", internal)
        while cand in used:
            cand += "_"
        used.add(cand)
        return cand
    sp = {mid: pick(m.sbml_id, "M_", mid) for mid, m in sorted(pw.metabolites.items())}
    rx = {}
    for rid, r in pw.reactions.items():
        rx[rid] = pick(r.sbml_id, "R_", rid)
        for k, a in enumerate(r.aliases):
            rx[(rid, k)] = pick(a.get("id", ""), "R_", f"{rid}_alias{k + 1}")
    return sp, rx
