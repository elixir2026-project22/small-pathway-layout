"""SBML Level 3 Version 1 core + Layout package (v1) writer, and a reader for previous layouts.

Clones are expressed the standard way: several speciesGlyphs referencing the same species.
Layout metadata (algorithm, hashes, parent layout) is stored in an annotation on the layout element.
"""
from __future__ import annotations

import xml.etree.ElementTree as ET

from .layout import Layout, Previous
from .model import EMPTY

CORE = "http://www.sbml.org/sbml/level3/version1/core"
LAY = "http://www.sbml.org/sbml/level3/version1/layout/version1"
XSI = "http://www.w3.org/2001/XMLSchema-instance"
META = "urn:sbml-autolayout:layout-metadata:1"
META_OLD = "urn:metlayout:layout-metadata:1"   # files written under the old name

ET.register_namespace("", CORE)
ET.register_namespace("layout", LAY)
ET.register_namespace("xsi", XSI)
ET.register_namespace("autolayout", META)


def _c(tag):
    return f"{{{CORE}}}{tag}"


def _l(tag):
    return f"{{{LAY}}}{tag}"


def _f(v: float) -> str:
    return f"{round(v, 2):g}"


def _sid(prefix: str, x: str) -> str:
    return prefix + "".join(ch if ch.isalnum() or ch == "_" else "_" for ch in x)


def _bbox(parent, bid, x, y, w, h):
    bb = ET.SubElement(parent, _l("boundingBox"), {_l("id"): bid})
    ET.SubElement(bb, _l("position"), {_l("x"): _f(x), _l("y"): _f(y)})
    ET.SubElement(bb, _l("dimensions"), {_l("width"): _f(w), _l("height"): _f(h)})


def _curve(parent, segs):
    cv = ET.SubElement(parent, _l("curve"))
    lst = ET.SubElement(cv, _l("listOfCurveSegments"))
    for s in segs:
        if s[0] == "L":
            el = ET.SubElement(lst, _l("curveSegment"), {f"{{{XSI}}}type": "LineSegment"})
            pts = [("start", s[1]), ("end", s[2])]
        else:
            el = ET.SubElement(lst, _l("curveSegment"), {f"{{{XSI}}}type": "CubicBezier"})
            pts = [("start", s[1]), ("end", s[4]), ("basePoint1", s[2]), ("basePoint2", s[3])]
        for tag, p in pts:
            ET.SubElement(el, _l(tag), {_l("x"): _f(p[0]), _l("y"): _f(p[1])})


MOD_SBO = {"catalysis": "SBO:0000013", "stimulation": "SBO:0000459", "inhibition": "SBO:0000020",
           "modulation": "SBO:0000019"}


def write_sbml(L: Layout, path: str):
    pw, P = L.pathway, L.params
    root = ET.Element(_c("sbml"), {"level": "3", "version": "1", _l("required"): "false"})
    model = ET.SubElement(root, _c("model"), {"id": _sid("", pw.id) or "pathway", "name": pw.name})

    comps = sorted({m.compartment for k, m in pw.metabolites.items() if k != EMPTY})
    lc = ET.SubElement(model, _c("listOfCompartments"))
    for c in comps:
        ET.SubElement(lc, _c("compartment"), {"id": c, "constant": "true", "spatialDimensions": "3", "size": "1"})
    ls = ET.SubElement(model, _c("listOfSpecies"))
    used = {m for r in pw.reactions.values() for m in list(r.substrates) + list(r.products) + list(r.modifiers)}
    for mid in sorted(used):
        m = pw.metabolites[mid]
        ET.SubElement(ls, _c("species"), {
            "id": _sid("M_", mid), "name": m.name, "compartment": m.compartment,
            "hasOnlySubstanceUnits": "false", "boundaryCondition": "false", "constant": "false"})
    lr = ET.SubElement(model, _c("listOfReactions"))
    for rid, r in pw.reactions.items():
        re_ = ET.SubElement(lr, _c("reaction"), {"id": _sid("R_", rid), "name": r.name,
                                                 "reversible": str(r.reversible).lower(), "fast": "false"})
        for tag, parts, suf in (("listOfReactants", r.substrates, "s"), ("listOfProducts", r.products, "p")):
            if parts:
                lst = ET.SubElement(re_, _c(tag))
                for mid, st in parts.items():
                    ET.SubElement(lst, _c("speciesReference"), {
                        "id": _sid("SR_", f"{rid}_{mid}_{suf}"), "species": _sid("M_", mid),
                        "stoichiometry": _f(st), "constant": "true"})
        if r.modifiers:      # catalysts / regulators, with the SBO term of their kind
            lst = ET.SubElement(re_, _c("listOfModifiers"))
            for mid in sorted(r.modifiers):
                ET.SubElement(lst, _c("modifierSpeciesReference"), {
                    "id": _sid("SR_", f"{rid}_{mid}_m"), "species": _sid("M_", mid),
                    "sboTerm": MOD_SBO.get(r.modifiers[mid], "SBO:0000019")})

    lol = ET.SubElement(model, _l("listOfLayouts"))
    lay = ET.SubElement(lol, _l("layout"), {_l("id"): _sid("layout_", pw.id), _l("name"): "sbml-autolayout"})
    ann = ET.SubElement(lay, _c("annotation"))
    ET.SubElement(ann, f"{{{META}}}info", {f"{{{META}}}{k}": str(v) for k, v in L.meta.items()} |
                  {f"{{{META}}}orientation": P.orientation, f"{{{META}}}cycleDirection": P.cycle_direction})
    ET.SubElement(lay, _l("dimensions"), {_l("width"): _f(L.width), _l("height"): _f(L.height)})

    # compartment glyph spanning the drawing (most editors expect one)
    lcg = ET.SubElement(lay, _l("listOfCompartmentGlyphs"))
    for c in comps:
        cg = ET.SubElement(lcg, _l("compartmentGlyph"), {_l("id"): f"CG_{c}", _l("compartment"): c})
        _bbox(cg, f"CG_{c}_bb", 0, 0, L.width, L.height)

    lsg = ET.SubElement(lay, _l("listOfSpeciesGlyphs"))
    r_ = P.met_radius
    for mid in sorted(L.met_pos):
        x, y = L.met_pos[mid]
        g = ET.SubElement(lsg, _l("speciesGlyph"), {_l("id"): _sid("", L.sg(mid)), _l("species"): _sid("M_", mid)})
        _bbox(g, _sid("", L.sg(mid)) + "_bb", x - r_, y - r_, 2 * r_, 2 * r_)
    cr = P.clone_radius
    for c in sorted(L.clones, key=lambda c: c.glyph):
        if c.species == EMPTY or c.role.startswith("drain"):   # symbols, not species glyphs
            continue
        g = ET.SubElement(lsg, _l("speciesGlyph"), {_l("id"): _sid("", c.glyph), _l("species"): _sid("M_", c.species)})
        _bbox(g, _sid("", c.glyph) + "_bb", c.pos[0] - cr, c.pos[1] - cr, 2 * cr, 2 * cr)

    lrg = ET.SubElement(lay, _l("listOfReactionGlyphs"))
    edges_by_rxn: dict[str, list] = {}
    for e in L.edges:
        edges_by_rxn.setdefault(e.rxn, []).append(e)
    h = P.rxn_half
    for rid in sorted(L.rxn_pos):
        (x, y), (ax, ay) = L.rxn_pos[rid], L.rxn_axis[rid]
        g = ET.SubElement(lrg, _l("reactionGlyph"), {_l("id"): _sid("", L.rg(rid)), _l("reaction"): _sid("R_", rid)})
        _bbox(g, _sid("", L.rg(rid)) + "_bb", x - h, y - h, 2 * h, 2 * h)
        _curve(g, [("L", (x - ax * h, y - ay * h), (x + ax * h, y + ay * h))])
        lsrg = ET.SubElement(g, _l("listOfSpeciesReferenceGlyphs"))
        for e in edges_by_rxn.get(rid, []):
            if e.species == EMPTY:
                continue
            suf = "m" if (e.drawn or "").startswith("mod:") else ("s" if e.role in ("substrate", "sidesubstrate") else "p")
            el = ET.SubElement(lsrg, _l("speciesReferenceGlyph"), {
                _l("id"): _sid("", e.id), _l("speciesGlyph"): _sid("", e.glyph),
                _l("speciesReference"): _sid("SR_", f"{rid}_{e.species}_{suf}"), _l("role"): e.role})
            _curve(el, e.segs)

    # pseudo-reaction (biomass): one small reaction glyph per drained metabolite (SBML Layout allows
    # several glyphs for one reaction), each linked to that metabolite only
    for c in sorted((c for c in L.clones if c.role.startswith("drain")), key=lambda c: c.glyph):
        x, y = c.pos
        g = ET.SubElement(lrg, _l("reactionGlyph"), {_l("id"): _sid("RG_", c.glyph[3:]), _l("reaction"): _sid("R_", c.rxn)})
        _bbox(g, _sid("RG_", c.glyph[3:]) + "_bb", x - h, y - h, 2 * h, 2 * h)
        lsrg = ET.SubElement(g, _l("listOfSpeciesReferenceGlyphs"))
        for e in edges_by_rxn.get(c.rxn, []):
            if e.species != c.species:
                continue
            suf = "s" if e.role == "substrate" else "p"
            el = ET.SubElement(lsrg, _l("speciesReferenceGlyph"), {
                _l("id"): _sid("", e.id), _l("speciesGlyph"): _sid("", e.glyph),
                _l("speciesReference"): _sid("SR_", f"{c.rxn}_{c.species}_{suf}"), _l("role"): e.role})
            _curve(el, e.segs)

    ltg = ET.SubElement(lay, _l("listOfTextGlyphs"))
    drain_glyphs = {c.glyph for c in L.clones if c.role.startswith("drain")}
    for tg in sorted(L.labels):
        lb = L.labels[tg]
        if lb.target in drain_glyphs:          # label of a drain symbol -> its reaction glyph
            el = ET.SubElement(ltg, _l("textGlyph"), {_l("id"): _sid("", tg), _l("graphicalObject"): _sid("RG_", lb.target[3:]),
                                                      _l("text"): lb.text})
            _bbox(el, _sid("", tg) + "_bb", lb.x, lb.y, lb.w, lb.h)
            continue
        attrs = {_l("id"): _sid("", tg), _l("graphicalObject"): _sid("", lb.target)}
        full = pw.metabolites[lb.species].name if lb.species else None
        if lb.species and lb.text != full:      # shortened / cofactor short name: store the shown text
            attrs[_l("text")] = lb.text
        elif lb.species:
            attrs[_l("originOfText")] = _sid("M_", lb.species)
        elif lb.rxn:
            attrs[_l("originOfText")] = _sid("R_", lb.rxn)
        el = ET.SubElement(ltg, _l("textGlyph"), attrs)
        _bbox(el, _sid("", tg) + "_bb", lb.x, lb.y, lb.w, lb.h)

    ET.indent(root, "  ")
    ET.ElementTree(root).write(path, encoding="UTF-8", xml_declaration=True)


def read_previous(path: str) -> Previous:
    """Read glyph positions from an SBML layout written by sbml-autolayout (or edited in another tool).

    Main species glyphs (SG_M_<id>), reaction glyphs (RG_R_<id>) and text glyphs (TG_*) are matched by id;
    their positions become hard pins for the next layout run.
    """
    tree = ET.parse(path)
    prev = Previous()

    def num(el, name: str) -> float:
        v = el.get(_l(name)) if el is not None else None
        return float(v) if v is not None else 0.0

    def bb_center(g) -> tuple[float, float, float, float]:
        bb = g.find(_l("boundingBox"))
        if bb is None:
            return 0.0, 0.0, 0.0, 0.0
        p, d = bb.find(_l("position")), bb.find(_l("dimensions"))
        return num(p, "x"), num(p, "y"), num(d, "width"), num(d, "height")

    for g in tree.iter(_l("speciesGlyph")):
        gid = g.get(_l("id")) or ""
        x, y, w, h = bb_center(g)
        if gid.startswith("SG_M_") and "__R_" not in gid:
            prev.met[gid[5:]] = (round(x + w / 2, 2), round(y + h / 2, 2))
        elif gid.startswith("SG_M_"):
            prev.clones[gid] = (x + w / 2, y + h / 2)
    for g in tree.iter(_l("reactionGlyph")):
        gid = g.get(_l("id")) or ""
        cv = g.find(_l("curve"))
        if g.find(_l("boundingBox")) is not None and gid.startswith("RG_M_") and "__R_" in gid:
            x, y, w, h = bb_center(g)             # pseudo-reaction drain symbol (kept like a clone)
            prev.clones["SG_" + gid[3:]] = (x + w / 2, y + h / 2)
        elif g.find(_l("boundingBox")) is not None and gid.startswith("RG_R_"):
            x, y, w, h = bb_center(g)
            prev.rxn[gid[5:]] = (round(x + w / 2, 2), round(y + h / 2, 2))
        elif cv is not None and gid.startswith("RG_R_"):
            segs_ = cv.find(_l("listOfCurveSegments"))
            if segs_ is None or len(segs_) == 0:
                continue
            s, e = segs_[0].find(_l("start")), segs_[0].find(_l("end"))
            prev.rxn[gid[5:]] = ((num(s, "x") + num(e, "x")) / 2, (num(s, "y") + num(e, "y")) / 2)
    for g in tree.iter(_l("textGlyph")):
        x, y, _, _ = bb_center(g)
        prev.labels[g.get(_l("id")) or ""] = (x, y)
    for ns in (META, META_OLD):
        info = tree.find(f".//{{{ns}}}info")
        if info is not None:
            prev.layout_hash = info.get(f"{{{ns}}}layout_hash")
            break
    return prev
