"""SVG rendering of a computed layout (same geometry as the SBML export)."""
from __future__ import annotations

import math
from html import escape

from . import sbgn
from .layout import Layout
from .model import EMPTY

STYLE = """
  .main{fill:none;stroke:#3a4a5c;stroke-linecap:round}
  .side{fill:none;stroke:#8a96a3;stroke-width:1.1}
  .mod{fill:none;stroke:#5b6f8a;stroke-width:1.3}
  .mod-inhibition{stroke:#a34d4d}
  .sink{fill:#fff;stroke:#3a4a5c;stroke-width:1.2}
  .drain{fill:none;stroke:#7b8796;stroke-width:1.2;stroke-dasharray:3 2}
  .pheno{fill:#f3efe2;stroke:#7b8796;stroke-width:1.2}
  .far{fill:none;stroke:#3a4a5c;stroke-width:1.4}
  .met{fill:#fdf6e3;stroke:#c58b1c;stroke-width:1.6}
  .clone{fill:#eef2f6;stroke:#8a96a3;stroke-width:1}
  .clonebar{fill:#8a96a3}
  .rxn{fill:#fff;stroke:#3a4a5c;stroke-width:1.4}
  .cycle-hint{fill:none;stroke:#e9edf2;stroke-width:24}
  .lm{font:600 12px 'Helvetica Neue',Arial,sans-serif;fill:#1f2933}
  .lk{font:600 11px 'Helvetica Neue',Arial,sans-serif;fill:#1f2933}
  .lc{font:400 10px 'Helvetica Neue',Arial,sans-serif;fill:#8a96a3}
  .leader{stroke:#c3cad3;stroke-width:0.8;stroke-dasharray:2 2}
  .lr{font:italic 9px 'Helvetica Neue',Arial,sans-serif;fill:#9aa5b1}
  .meta{font:9px monospace;fill:#b0b8c1}
  .g-sc{fill:#fdf6e3;stroke:#a87613;stroke-width:1.3}
  .g-mm{fill:#e9f1fa;stroke:#3f6f9f;stroke-width:1.3}
  .g-cx{fill:#edf3e6;stroke:#5e7f45;stroke-width:1.3}
  .g-na{fill:#f5ecf5;stroke:#835a83;stroke-width:1.3}
  .g-ue{fill:#f2f3f5;stroke:#6b7785;stroke-width:1.3}
  .g-ph{fill:#f3efe2;stroke:#7b8796;stroke-width:1.3}
  .g-cof{fill:#f4f6f8;stroke:#8a96a3;stroke-width:1}
  .g-bar{fill:#b9c2cc}
  .port{stroke:#3a4a5c;stroke-width:1.4}
  .tm{font:400 10.5px 'Helvetica Neue',Arial,sans-serif;fill:#1f2933}
  .tk{font:400 10px 'Helvetica Neue',Arial,sans-serif;fill:#1f2933}
  .tc{font:400 9px 'Helvetica Neue',Arial,sans-serif;fill:#52606d}
"""

GLYPH_CSS = {"simple chemical": "g-sc", "macromolecule": "g-mm", "complex": "g-cx",
             "nucleic acid feature": "g-na", "unspecified entity": "g-ue", "phenotype": "g-ph"}


def _d(segs) -> str:
    out = []
    for i, s in enumerate(segs):
        p0 = s[1]
        if i == 0:
            out.append(f"M{p0[0]:.1f},{p0[1]:.1f}")
        if s[0] == "L":
            out.append(f"L{s[2][0]:.1f},{s[2][1]:.1f}")
        else:
            out.append(f"C{s[2][0]:.1f},{s[2][1]:.1f} {s[3][0]:.1f},{s[3][1]:.1f} {s[4][0]:.1f},{s[4][1]:.1f}")
    return " ".join(out)


def write_svg(L: Layout, path: str, show_meta: bool = True):
    P = L.params
    W, H = L.width, L.height + (16 if show_meta else 0)
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W:.0f}" height="{H:.0f}" '
         f'viewBox="0 0 {W:.1f} {H:.1f}">',
         f"<style>{STYLE}</style>",
         "<defs>",
         '<marker id="am" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" '
         'orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="#3a4a5c"/></marker>',
         '<marker id="as" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" '
         'orient="auto-start-reverse"><path d="M0,1 L10,5 L0,9 z" fill="#8a96a3"/></marker>',
         # SBGN arc ends: catalysis (open circle), stimulation (open arrow), inhibition (bar),
         # modulation (open diamond)
         '<marker id="m-catalysis" viewBox="0 0 12 12" refX="11" refY="6" markerWidth="9" markerHeight="9" '
         'markerUnits="userSpaceOnUse" orient="auto"><circle cx="6" cy="6" r="5" fill="#fff" stroke="#5b6f8a" stroke-width="1.4"/></marker>',
         '<marker id="m-stimulation" viewBox="0 0 12 12" refX="11" refY="6" markerWidth="10" markerHeight="10" '
         'markerUnits="userSpaceOnUse" orient="auto"><path d="M1,1 L11,6 L1,11 z" fill="#fff" stroke="#5b6f8a" stroke-width="1.4"/></marker>',
         '<marker id="m-inhibition" viewBox="0 0 12 12" refX="10" refY="6" markerWidth="10" markerHeight="10" '
         'markerUnits="userSpaceOnUse" orient="auto"><path d="M10,0 L10,12" stroke="#a34d4d" stroke-width="2.2"/></marker>',
         '<marker id="m-modulation" viewBox="0 0 14 12" refX="13" refY="6" markerWidth="11" markerHeight="10" '
         'markerUnits="userSpaceOnUse" orient="auto"><path d="M1,6 L7,1 L13,6 L7,11 z" fill="#fff" stroke="#5b6f8a" stroke-width="1.4"/></marker>',
         "</defs>",
         f'<rect width="{W:.1f}" height="{H:.1f}" fill="#ffffff"/>']
    for c in L.cycles:
        o.append(f'<circle class="cycle-hint" cx="{c.center[0]:.1f}" cy="{c.center[1]:.1f}" r="{c.radius:.1f}"/>')
    for e in L.edges:
        drawn = e.drawn or e.role
        if drawn.startswith("drain"):    # pseudo-reaction (biomass): short dashed arrow
            o.append(f'<path class="drain" d="{_d(e.segs)}" marker-end="url(#as)"><title>{escape(e.rxn)}</title></path>')
            continue
        if drawn.startswith("mod:"):
            kind = drawn[4:]
            o.append(f'<path class="mod mod-{kind}" d="{_d(e.segs)}" marker-end="url(#m-{kind})">'
                     f'<title>{escape(kind)}: {escape(L.pathway.metabolites[e.species].name)}</title></path>')
            continue
        if drawn in ("substrate", "product") and not e.far:
            sw = 1.8 * max(0.5, min(e.weight, 4.0)) ** 0.5
            end = ' marker-end="url(#am)"' if drawn == "product" or e.reversible else ""
            o.append(f'<path class="main" style="stroke-width:{sw:.2f}" d="{_d(e.segs)}"{end}><title>{escape(e.id)}</title></path>')
        else:
            mk = "am" if e.far else "as"
            end = f' marker-end="url(#{mk})"' if drawn in ("sideproduct", "product") or e.reversible else ""
            o.append(f'<path class="{"far" if e.far else "side"}" d="{_d(e.segs)}"{end}/>')
    if L.sbgn:
        _sbgn_nodes(L, o)
    else:
        h = P.rxn_half
        for rid, (x, y) in sorted(L.rxn_pos.items()):
            ax, ay = L.rxn_axis[rid]
            ang = math.degrees(math.atan2(ay, ax))
            o.append(f'<rect class="rxn" x="{x-h:.1f}" y="{y-h:.1f}" width="{2*h:.1f}" height="{2*h:.1f}" '
                     f'transform="rotate({ang:.1f} {x:.1f} {y:.1f})"><title>{escape(L.pathway.reactions[rid].name)}</title></rect>')
        cr = P.clone_radius
        for c in L.clones:
            x, y = c.pos
            if c.role.startswith("drain"):   # SBGN phenotype (hexagon): biomass / pseudo-reaction
                w_, h_ = cr * 1.6, cr * 1.0
                pts = [(x - w_, y), (x - w_ * 0.6, y - h_), (x + w_ * 0.6, y - h_), (x + w_, y), (x + w_ * 0.6, y + h_), (x - w_ * 0.6, y + h_)]
                o.append('<path class="pheno" d="M' + " L".join(f"{a:.1f},{b:.1f}" for a, b in pts) + ' Z">'
                         f'<title>{escape(L.pathway.reactions[c.rxn].name)}</title></path>')
                continue
            if c.species == EMPTY:     # SBGN source/sink: circle with a diagonal stroke
                o.append(f'<circle class="sink" cx="{x:.1f}" cy="{y:.1f}" r="{cr:.1f}"/>')
                o.append(f'<path class="sink" d="M{x-cr*0.95:.1f},{y+cr*0.95:.1f} L{x+cr*0.95:.1f},{y-cr*0.95:.1f}"/>')
                continue
            if c.far:   # main metabolite repeated here instead of a long link: main colours, clone marker
                o.append(f'<circle class="met" cx="{x:.1f}" cy="{y:.1f}" r="{cr+1.5:.1f}"/>')
                o.append(f'<path class="clonebar" d="M{x-(cr+1.5)*0.87:.1f},{y+(cr+1.5)*0.5:.1f} A{cr+1.5:.1f},{cr+1.5:.1f} 0 0 0 '
                         f'{x+(cr+1.5)*0.87:.1f},{y+(cr+1.5)*0.5:.1f} Z"/>')
                continue
            o.append(f'<circle class="clone" cx="{x:.1f}" cy="{y:.1f}" r="{cr:.1f}"/>')
            # SBGN-style clone marker: filled lower band
            o.append(f'<path class="clonebar" d="M{x-cr*0.87:.1f},{y+cr*0.5:.1f} A{cr:.1f},{cr:.1f} 0 0 0 '
                     f'{x+cr*0.87:.1f},{y+cr*0.5:.1f} Z"/>')
        r = P.met_radius
        for mid, (x, y) in sorted(L.met_pos.items()):
            o.append(f'<circle class="met" cx="{x:.1f}" cy="{y:.1f}" r="{r:.1f}"><title>{escape(mid)}</title></circle>')
        for lb in sorted(L.labels.values(), key=lambda l: l.glyph):
            if lb.leader:   # thin line from the node to a label that had to be placed further away
                cx, cy = min(max(lb.node[0], lb.x), lb.x + lb.w), min(max(lb.node[1], lb.y), lb.y + lb.h)
                o.append(f'<line class="leader" x1="{lb.node[0]:.1f}" y1="{lb.node[1]:.1f}" x2="{cx:.1f}" y2="{cy:.1f}"/>')
            cls = {"main": "lm", "copy": "lk", "side": "lc", "reaction": "lr"}[lb.kind]
            tx = lb.x + 1 if lb.anchor == "start" else lb.x + lb.w - 1 if lb.anchor == "end" else lb.x + lb.w / 2
            o.append(f'<text class="{cls}" text-anchor="{lb.anchor}" x="{tx:.1f}" y="{lb.y + lb.h - 4.5:.1f}">'
                     f'{escape(lb.text)}</text>')
    if show_meta:
        m = L.meta
        o.append(f'<text class="meta" x="8" y="{H-5:.1f}">{escape(m["algorithm"])} · layout {m["layout_hash"]}'
                 f' · pathway {m["pathway_hash"]}' + (f' · parent {m["parent_layout_hash"]}' if m["parent_layout_hash"] else "")
                 + "</text>")
    o.append("</svg>")
    with open(path, "w") as fh:
        fh.write("\n".join(o))


def _sbgn_nodes(L: Layout, o: list):
    """SBGN Process Description: process nodes with ports, entity pools as shapes with the name
    inside, cofactor copies with clone markers"""
    P = L.params
    h = P.rxn_half
    for rid, (x, y) in sorted(L.rxn_pos.items()):
        ax, ay = L.rxn_axis[rid]
        for sg_ in (-1, 1):
            a0 = (x + sg_ * ax * h, y + sg_ * ay * h)
            a1 = (x + sg_ * ax * (h + P.port), y + sg_ * ay * (h + P.port))
            o.append(f'<line class="port" x1="{a0[0]:.1f}" y1="{a0[1]:.1f}" x2="{a1[0]:.1f}" y2="{a1[1]:.1f}"/>')
        ang = math.degrees(math.atan2(ay, ax))
        o.append(f'<rect class="rxn" x="{x-h:.1f}" y="{y-h:.1f}" width="{2*h:.1f}" height="{2*h:.1f}" '
                 f'transform="rotate({ang:.1f} {x:.1f} {y:.1f})"><title>{escape(L.pathway.reactions[rid].name)}</title></rect>')
    cr = P.clone_radius
    for c in L.clones:
        x, y = c.pos
        if c.role.startswith("drain"):   # SBGN phenotype (hexagon): biomass / pseudo-reaction
            o.append(sbgn.outline_svg("phenotype", x, y, cr * 3.2, cr * 2.0, "pheno",
                                      escape(L.pathway.reactions[c.rxn].name)))
            continue
        if c.species == EMPTY:     # SBGN source/sink: circle with a diagonal stroke
            o.append(f'<circle class="sink" cx="{x:.1f}" cy="{y:.1f}" r="{cr:.1f}"/>')
            o.append(f'<path class="sink" d="M{x-cr*0.95:.1f},{y+cr*0.95:.1f} L{x+cr*0.95:.1f},{y-cr*0.95:.1f}"/>')
            continue
        w, hh, _, cls = L.clone_dims(c)
        css = GLYPH_CSS.get(cls, "g-sc") if (c.far or c.role == "modifier") else "g-cof"
        o.append(sbgn.outline_svg(cls, x, y, w, hh, css, escape(L.pathway.metabolites[c.species].name)))
        o.append(sbgn.clone_marker_svg(cls, x, y, w, hh, "g-bar"))
    for mid, (x, y) in sorted(L.met_pos.items()):
        w, hh, _, cls = L.met_dims(mid)
        o.append(sbgn.outline_svg(cls, x, y, w, hh, GLYPH_CSS.get(cls, "g-sc"), escape(mid)))
    for lb in sorted(L.labels.values(), key=lambda l: l.glyph):
        if not lb.inside:              # reaction names (and "biomass") stay beside their node
            if lb.leader:
                cx, cy = min(max(lb.node[0], lb.x), lb.x + lb.w), min(max(lb.node[1], lb.y), lb.y + lb.h)
                o.append(f'<line class="leader" x1="{lb.node[0]:.1f}" y1="{lb.node[1]:.1f}" x2="{cx:.1f}" y2="{cy:.1f}"/>')
            cls = {"main": "lm", "copy": "lk", "side": "lc", "reaction": "lr"}[lb.kind]
            tx = lb.x + 1 if lb.anchor == "start" else lb.x + lb.w - 1 if lb.anchor == "end" else lb.x + lb.w / 2
            o.append(f'<text class="{cls}" text-anchor="{lb.anchor}" x="{tx:.1f}" y="{lb.y + lb.h - 4.5:.1f}">'
                     f'{escape(lb.text)}</text>')
            continue
        cls = {"main": "tm", "copy": "tk", "side": "tc"}.get(lb.kind, "tm")
        lh = sbgn.LINE.get(lb.kind, 12.0)
        cx, cy = lb.x + lb.w / 2, lb.y + lb.h / 2
        y0 = cy - len(lb.lines) * lh / 2 + lh * 0.5 + lb.size * 0.35
        spans = "".join(f'<tspan x="{cx:.1f}" y="{y0 + i * lh:.1f}">{escape(t)}</tspan>' for i, t in enumerate(lb.lines))
        o.append(f'<text class="{cls}" text-anchor="middle"><title>{escape(lb.text)}</title>{spans}</text>')
