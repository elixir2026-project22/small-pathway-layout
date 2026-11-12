"""Escher map JSON (https://escher.github.io) import.

The reactions/metabolites of the map become a Pathway. The curator's drawing is kept as reference
information: which metabolite nodes are primary (main) per reaction, and their coordinates.
"""
from __future__ import annotations

import json
import re

from .model import Pathway, pathway_from_dict


def is_escher(path: str) -> bool:
    try:
        with open(path) as fh:
            d = json.load(fh)
    except Exception:
        return False
    return isinstance(d, list) and len(d) == 2 and isinstance(d[1], dict) and "reactions" in d[1]


def _split(bigg: str) -> tuple[str, str]:
    """BiGG id -> (base, compartment); both "glc__D_c" and "glc_D[c]" styles"""
    m = re.match(r"^(.*)\[([A-Za-z0-9_]{1,4})\]$", bigg)
    if m:
        return m.group(1), m.group(2)
    m = re.match(r"^(.*)_([a-z]{1,2}\d?)$", bigg)
    return (m.group(1), m.group(2)) if m else (bigg, "c")


def import_escher(path: str) -> tuple[Pathway, dict]:
    head, body = json.load(open(path))
    nodes, rxns = body["nodes"], body["reactions"]
    names = {n["bigg_id"]: n.get("name", n["bigg_id"]) for n in nodes.values() if n["node_type"] == "metabolite"}
    mets, reactions, primary, coords = {}, [], {}, {}
    for key in sorted(rxns, key=lambda k: rxns[k]["bigg_id"]):
        r = rxns[key]
        rid = re.sub(r"[^A-Za-z0-9_]", "_", r["bigg_id"])
        s, p = {}, {}
        for x in r["metabolites"]:
            mid = x["bigg_id"]
            base, comp = _split(mid)
            mets.setdefault(mid, {"name": names.get(mid, mid), "compartment": comp, "base": base})
            (s if x["coefficient"] < 0 else p)[mid] = abs(x["coefficient"])
        # nodes touched by this reaction's segments -> which participants the curator drew as primary
        touched = {n for seg in r["segments"].values() for n in (seg["from_node_id"], seg["to_node_id"])}
        primary[rid] = sorted({nodes[n]["bigg_id"] for n in touched
                               if nodes[n]["node_type"] == "metabolite" and nodes[n]["node_is_primary"]})
        for n in touched:
            if nodes[n]["node_type"] == "metabolite" and nodes[n]["node_is_primary"]:
                coords.setdefault(nodes[n]["bigg_id"], (nodes[n]["x"], nodes[n]["y"]))
        entry = {"id": rid, "name": r.get("name") or rid, "substrates": s, "products": p,
                 "reversible": bool(r.get("reversibility", False))}
        if primary[rid]:          # the curator's choice: primary nodes are main, the rest side
            entry.update(main=primary[rid], main_explicit=True, role_source="Escher node_is_primary")
        reactions.append(entry)
    comps = {m["compartment"] for m in mets.values()}
    counts = {}
    for v in mets.values():
        counts[v["compartment"]] = counts.get(v["compartment"], 0) + 1
    if len(counts) > 1:   # tag only metabolites outside the dominant compartment
        main_c = max(sorted(counts), key=lambda c: counts[c])
        for v in mets.values():
            if v["compartment"] != main_c:
                v["comp_label"] = v["compartment"]
    hide = sorted(k for k, v in mets.items() if v["base"] in ("h",))
    raw = {"id": re.sub(r"\W", "_", head.get("map_name", "escher_map")), "name": head.get("map_name", ""),
           "hide": hide, "metabolites": dict(sorted(mets.items())), "reactions": reactions,
           "source": {"escher": path.rsplit("/", 1)[-1]}}
    return pathway_from_dict(raw), {"primary": primary, "coords": coords}
