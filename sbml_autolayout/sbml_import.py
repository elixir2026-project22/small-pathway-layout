"""Import one or several SBML files (Level 2 or 3, any version) into a Pathway.

Merging across files
  Species are matched, per compartment, by (in this order) database annotations (ChEBI, KEGG,
  MetaNetX, BiGG, HMDB, InChIKey), by normalised id (prefix "M_" and compartment suffix removed),
  and optionally by name. Compartments are matched by id or name (cytosol/cytoplasm -> c, etc.).
  Reactions with the same stoichiometry (after species matching) are kept only once.

Filtering
  Boundary species are dropped; reactions left without substrates or products (exchange, demand,
  sink) are skipped unless keep_exchange=True. Subsystems are read from the SBML groups package
  and from "SUBSYSTEM:" notes; select them with subsystems=["TCA", ...] (case-insensitive substring).
"""
from __future__ import annotations

import os
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field

from .model import MODIFIER_SBO, PSEUDO_SBO, Pathway, looks_pseudo, pathway_from_dict

ID_PATTERNS = [
    ("chebi", re.compile(r"chebi[:/](?:CHEBI[:_])?(\d+)", re.I)),
    ("kegg", re.compile(r"kegg\.compound[:/]([CG]\d{5})", re.I)),
    ("mnx", re.compile(r"metanetx\.chemical[:/](MNXM\w+)", re.I)),
    ("bigg", re.compile(r"bigg\.metabolite[:/]([\w-]+)", re.I)),
    ("hmdb", re.compile(r"hmdb[:/](HMDB\d+)", re.I)),
    ("inchikey", re.compile(r"inchikey[:/]([A-Z]{14}-[A-Z]{10}-[A-Z])", re.I)),
]
COMPARTMENT_SYNONYMS = {
    "c": "c", "cytosol": "c", "cytoplasm": "c", "cyto": "c", "cytosolic": "c",
    "m": "m", "mitochondrion": "m", "mitochondria": "m", "mitochondrial matrix": "m",
    "e": "e", "extracellular": "e", "extracellular space": "e", "extracellular region": "e",
    "p": "p", "periplasm": "p", "n": "n", "nucleus": "n", "r": "r", "endoplasmic reticulum": "r",
    "x": "x", "peroxisome": "x", "g": "g", "golgi apparatus": "g", "v": "v", "vacuole": "v",
    "l": "l", "lysosome": "l", "h": "h", "chloroplast": "h",
}
PROTON = {"h", "h+", "proton", "hydron"}


def _local(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _kids(el, name):
    return [c for c in el if _local(c.tag) == name] if el is not None else []


def _kid(el, name):
    k = _kids(el, name)
    return k[0] if k else None


def _attr(el, name, default=None):
    for k, v in el.attrib.items():
        if _local(k) == name:
            return v
    return default


def _text(el) -> str:
    return " ".join(t.strip() for t in el.itertext() if t.strip()) if el is not None else ""


def _clean(s: str, lead: bool = False) -> str:
    """Make an SBML-safe id fragment. Internal ids get an "M_"/"R_" prefix on export, so a leading
    digit is fine unless lead=True (used for the model id)."""
    s = re.sub(r"[^A-Za-z0-9_]", "_", s) or "_"
    return "_" + s if lead and s[0].isdigit() else s


@dataclass
class ImportReport:
    files: list[str] = field(default_factory=list)
    merged_species: list[str] = field(default_factory=list)
    duplicate_reactions: list[str] = field(default_factory=list)
    skipped_reactions: list[str] = field(default_factory=list)
    subsystems: dict[str, int] = field(default_factory=dict)
    role_sources: dict[str, int] = field(default_factory=dict)
    layout_copies: list[str] = field(default_factory=list)
    added_main: list[str] = field(default_factory=list)     # main metabolites missing from their reaction

    def summary(self) -> str:
        lines = [f"read {len(self.files)} file(s): {', '.join(self.files)}"]
        if self.merged_species:
            lines.append(f"merged {len(self.merged_species)} species across files: "
                         + ", ".join(self.merged_species[:12]) + (" ..." if len(self.merged_species) > 12 else ""))
        if self.duplicate_reactions:
            lines.append("duplicate reactions kept once: " + ", ".join(self.duplicate_reactions))
        if self.skipped_reactions:
            lines.append(f"skipped {len(self.skipped_reactions)} exchange/boundary/filtered reactions: "
                         + ", ".join(self.skipped_reactions[:10]) + (" ..." if len(self.skipped_reactions) > 10 else ""))
        if self.added_main:
            lines.append("main metabolites named in the notes but missing from the reaction, added: "
                         + ", ".join(self.added_main))
        if self.role_sources:
            lines.append("main/side decided by: " + ", ".join(f"{k} ({v})" for k, v in sorted(
                self.role_sources.items(), key=lambda x: -x[1])))
        if self.layout_copies:
            lines.append("drawn several times in the input layout -> side compounds: " + ", ".join(self.layout_copies[:12]))
        return "\n".join(lines)


MAIN_NOTE = re.compile(r"main[_ ]metabolites?\s*:\s*([^|\n<]+)", re.I)


def _layout_info(model):
    """From an SBML Layout in the input: species glyph counts, and per reaction the species whose
    reference glyph has role side*/substrate/product."""
    glyph_species, counts, roles = {}, {}, {}
    for g in model.iter():
        if _local(g.tag) == "speciesGlyph":
            sp = _attr(g, "species")
            glyph_species[_attr(g, "id")] = sp
            counts[sp] = counts.get(sp, 0) + 1
    for rg in model.iter():
        if _local(rg.tag) != "reactionGlyph":
            continue
        rid = _attr(rg, "reaction")
        for srg in rg.iter():
            if _local(srg.tag) == "speciesReferenceGlyph":
                role = (_attr(srg, "role") or "").lower()
                sp = glyph_species.get(_attr(srg, "speciesGlyph"))
                if role and sp:
                    roles.setdefault(rid, {})[sp] = role
    return counts, roles


def _celldesigner_main(r) -> set[str]:
    """CellDesigner annotation: baseReactant / baseProduct are the main species of a reaction."""
    out = set()
    for e in r.iter():
        if _local(e.tag) in ("baseReactant", "baseProduct"):
            sp = _attr(e, "species")
            if sp:
                out.add(sp)
    return out


def list_subsystems(paths: list[str]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for p in paths:
        model = _kid(ET.parse(p).getroot(), "model")
        for rid, subs in _subsystems(model).items():
            for s in subs:
                counts[s] = counts.get(s, 0) + 1
    return dict(sorted(counts.items()))


def _subsystems(model) -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    for g in model.iter():
        if _local(g.tag) != "group":
            continue
        gname = _attr(g, "name") or _attr(g, "id") or "group"
        for mem in g.iter():
            if _local(mem.tag) == "member":
                ref = _attr(mem, "idRef")
                if ref:
                    out.setdefault(ref, set()).add(gname)
    for r in _kids(_kid(model, "listOfReactions"), "reaction"):
        m = re.search(r"SUBSYSTEM\s*:\s*([^\n<]+?)(?:\s{2,}|$|[A-Z_]+\s*:)", _text(_kid(r, "notes")))
        if m:
            out.setdefault(r.get("id"), set()).add(m.group(1).strip())
    return out


def import_sbml(paths: list[str], subsystems: list[str] | None = None, keep_exchange: bool = False,
                match_names: bool = True, show_protons: bool = False,
                keep_regulated: bool = True) -> tuple[Pathway, ImportReport]:
    rep = ImportReport()
    key_to_id: dict[tuple, str] = {}
    mets: dict[str, dict] = {}
    rxns: list[dict] = []
    rxn_ids: set[str] = set()
    signatures: dict[tuple, str] = {}
    names = []
    filt = [s.lower() for s in subsystems or []]

    for path in paths:
        rep.files.append(os.path.basename(path))
        root = ET.parse(path).getroot()
        level = int(root.get("level", "3"))
        model = _kid(root, "model")
        if model is None:
            raise ValueError(f"{path}: no <model> element - not an SBML model?")
        names.append(model.get("name") or model.get("id") or os.path.basename(path))

        comp = {}
        for c in _kids(_kid(model, "listOfCompartments"), "compartment"):
            cid = c.get("id")
            label = (c.get("name") or cid).strip().lower()
            comp[cid] = COMPARTMENT_SYNONYMS.get(label, COMPARTMENT_SYNONYMS.get(cid.lower()))
            if comp[cid] is None:   # unknown compartment: initials of its name ("endoplasmic reticulum lumen" -> erl)
                comp[cid] = "".join(w[0] for w in re.split(r"[\s_-]+", label) if w)[:4] or _clean(cid.lower())
        comp_names = {(c.get("name") or "").strip().lower() for c in _kids(_kid(model, "listOfCompartments"), "compartment")}

        local_map: dict[str, str | None] = {}
        for sp in _kids(_kid(model, "listOfSpecies"), "species"):
            sid = sp.get("id")
            if sp.get("boundaryCondition", "false") == "true":
                local_map[sid] = None
                continue
            cid = comp.get(sp.get("compartment"), "c")
            name = re.sub(r"\s+", " ", (sp.get("name") or sid)).strip()
            m_ = re.match(r"^(.*?)\s*\[([^\]]+)\]$", name)      # Reactome: "ATP [mitochondrial matrix]"
            if m_ and m_.group(2).strip().lower() in comp_names:
                name = m_.group(1)
            base = re.sub(r"^(M_|S_)", "", sid)
            for suffix in {sp.get("compartment", ""), cid}:
                if suffix and base.endswith("_" + suffix):
                    base = base[: -len(suffix) - 1]
                    break
            # identity annotations only (bqbiol:is); hasPart/hasVersion would e.g. make a complex containing ATP look like ATP
            is_q = [e for e in sp.iter() if _local(e.tag) == "is"]
            src = is_q if is_q else [sp]
            uris = " ".join(v for q in src for e in q.iter() for k, v in e.attrib.items() if _local(k) == "resource")
            keys = [(db, m.group(1).upper(), cid) for db, rx in ID_PATTERNS for m in rx.finditer(uris)]
            if not any(k[0] == "chebi" for k in keys):   # Reactome style: ...searchId.do?chebiId=CHEBI:15377
                keys += [("chebi", c, cid) for c in re.findall(r"CHEBI[:_](\d+)", uris)]
            keys.append(("id", base.lower(), cid))
            if match_names:
                keys.append(("name", name.lower(), cid))
            hit = next((key_to_id[k] for k in keys if k in key_to_id), None)
            if hit is None:
                new = _clean(base if cid == "c" else f"{base}_{cid}")
                while new in mets:
                    new += "_"
                chebi = sorted({k[1] for k in keys if k[0] == "chebi"})
                kegg = sorted({k[1] for k in keys if k[0] == "kegg"})
                mets[new] = {"name": name, "compartment": cid, "base": base, "chebi": chebi, "kegg": kegg}
                hit = new
            elif len(rep.files) > 1 and hit not in rep.merged_species and sid != hit:
                rep.merged_species.append(f"{sid}->{hit}" if sid != hit else hit)
            for k in keys:
                key_to_id.setdefault(k, hit)
            local_map[sid] = hit

        subs_of = _subsystems(model)
        # enzyme names from EC classification groups (BRENDA): "EC 4.2.1.3 aconitate hydratase"
        ec_name = {}
        for g in model.iter():
            if _local(g.tag) == "group":
                m_ = re.match(r"EC\s+([\d.A-Za-z-]+)\s+(.+)", _attr(g, "name") or "")
                if m_:
                    ec_name[m_.group(1)] = m_.group(2).strip()
        glyph_counts, layout_roles = _layout_info(model)
        sbo_of = {sp.get("id"): sp.get("sboTerm") for sp in _kids(_kid(model, "listOfSpecies"), "species")}
        # species drawn >= 3 times in the input layout (and not proteins/complexes): side compounds there
        layout_side = {sp for sp, c in glyph_counts.items() if c >= 3 and sbo_of.get(sp) in (None, "SBO:0000247")
                       and local_map.get(sp)}
        for sp in sorted(layout_side):
            mid_ = local_map.get(sp)
            if not mid_:
                continue
            nm = mets[mid_]["name"]
            if nm not in rep.layout_copies:
                rep.layout_copies.append(nm)
        for r in _kids(_kid(model, "listOfReactions"), "reaction"):
            rid0 = r.get("id")
            rsubs = sorted(subs_of.get(rid0, set()))
            if filt and not any(f in s.lower() for s in rsubs for f in filt):
                continue
            rev = r.get("reversible", "true" if level < 3 else "false") == "true"
            sides = []
            for lst in ("listOfReactants", "listOfProducts"):
                d = {}
                for ref in _kids(_kid(r, lst), "speciesReference"):
                    m = local_map.get(ref.get("species"))
                    if m is None:
                        continue
                    d[m] = d.get(m, 0) + float(ref.get("stoichiometry", "1"))
                sides.append(d)
            s, p = sides
            # modifiers (catalysts, transcription factors, inhibitors): species -> kind
            mods = {}
            for ref in _kids(_kid(r, "listOfModifiers"), "modifierSpeciesReference"):
                m = local_map.get(ref.get("species"))
                if m is not None and m not in s and m not in p:
                    mods.setdefault(m, MODIFIER_SBO.get(ref.get("sboTerm") or "", "modulation"))
            # one side empty: exchange/demand/sink reactions are skipped; a regulated synthesis or
            # degradation (transcription "-> protein" with a transcription factor) is kept when
            # modifiers are drawn (keep_regulated)
            if not keep_exchange and (not s or not p) and not (mods and keep_regulated):
                rep.skipped_reactions.append(rid0)
                continue
            sig = (tuple(sorted(s.items())), tuple(sorted(p.items())))
            alt = (sig[1], sig[0])
            if sig in signatures or (rev and alt in signatures):
                kept = signatures.get(sig) or signatures.get(alt)
                rep.duplicate_reactions.append(f"{rid0}={kept}")
                for e_ in rxns:                    # same reaction, other regulators: add them
                    if e_["id"] == kept:
                        for m, k in mods.items():
                            e_.setdefault("modifiers", {}).setdefault(m, k)
                continue
            rid = _clean(re.sub(r"^R_", "", rid0))
            while rid in rxn_ids:
                rid += "_2"
            rxn_ids.add(rid)
            signatures[sig] = rid
            for sub in rsubs:
                rep.subsystems[sub] = rep.subsystems.get(sub, 0) + 1
            notes = _text(_kid(r, "notes"))
            name = (r.get("name") or rid).strip()
            ec = re.search(r"preferred_ec\s*:\s*([\d.A-Za-z-]+)", notes)
            if ec and (re.fullmatch(r"(BR|BS)\d+(-\w+)?|NODE_\d+|rx_\w+", name) or name == rid0):
                name = ec_name.get(ec.group(1), f"EC {ec.group(1)}")   # BRENDA ids -> enzyme name
            entry = {"id": rid, "name": name, "substrates": s, "products": p,
                     "reversible": rev, "subsystem": rsubs}
            if mods:
                entry["modifiers"] = mods
            if r.get("sboTerm") in PSEUDO_SBO or looks_pseudo(rid0, name):
                entry["pseudo"] = True
            # main/side information carried by the input, in order of priority
            mm = MAIN_NOTE.search(notes)
            cd_main = _celldesigner_main(r)
            lroles = layout_roles.get(rid0, {})
            src, main_set = None, None
            if mm:
                toks = [t.strip() for t in re.split(r"[,;\s]+", mm.group(1)) if t.strip()]
                main_set = {local_map.get(t) or local_map.get("M_" + t) or local_map.get("S_" + t) for t in toks}
                src = "main_metabolites note"
            elif cd_main:
                main_set, src = {local_map.get(x) for x in cd_main}, "CellDesigner base reactants"
            elif any(v.startswith("side") for v in lroles.values()):
                main_set = {local_map.get(x) for x, v in lroles.items() if v in ("substrate", "product")}
                src = "SBML Layout roles"
            if mm and main_set:
                # BRENDA: a named main metabolite missing from the reaction (incomplete equation, e.g.
                # the polymer made by mannuronan synthase): add it on the other side of its partner,
                # otherwise the pathway would fall apart at this reaction
                for m in sorted(x for x in main_set if x and x not in s and x not in p):
                    if any(x in s for x in main_set if x):
                        p[m] = 1.0
                    else:
                        s[m] = 1.0
                    rep.added_main.append(f"{mets[m]['name']} ({'product' if m in p else 'substrate'} of {name})")
            main_set = {m for m in (main_set or set()) if m and (m in s or m in p)}
            if main_set:
                entry.update(main=sorted(main_set), main_explicit=True, role_source=src)
            elif layout_side:
                side = sorted(x_ for x in layout_side if (x_ := local_map.get(x)) and (x_ in s or x_ in p))
                if side:
                    entry.update(side=side, role_source="copies in input layout + rules")
            rep.role_sources[entry.get("role_source", "rules")] = rep.role_sources.get(entry.get("role_source", "rules"), 0) + 1
            rxns.append(entry)

    used = {m for r in rxns for m in list(r["substrates"]) + list(r["products"]) + list(r.get("modifiers", {}))}
    mets = {k: v for k, v in mets.items() if k in used}
    counts = {}
    for v in mets.values():
        counts[v["compartment"]] = counts.get(v["compartment"], 0) + 1
    if len(counts) > 1:   # tag only metabolites outside the dominant compartment
        main_c = max(sorted(counts), key=lambda c: counts[c])
        for v in mets.values():
            if v["compartment"] != main_c:
                v["comp_label"] = v["compartment"]
    hide = [] if show_protons else sorted(
        k for k, v in mets.items() if v["base"].lower() in PROTON or v["name"].lower() in PROTON
        or "15378" in v.get("chebi", []))
    raw = {
        "id": _clean("_".join(os.path.splitext(os.path.basename(p))[0] for p in paths), lead=True)[:60],
        "name": " + ".join(names),
        "hide": hide,
        "metabolites": dict(sorted(mets.items())),
        "reactions": rxns,
        "source": {"files": rep.files, "subsystems": filt},
    }
    return pathway_from_dict(raw), rep
