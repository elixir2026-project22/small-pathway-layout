"""Export a layout as an Escher map (https://escher.github.io), loadable via Map > Load map JSON.

Structure follows Escher's map schema 1-0-0: [header, {reactions, nodes, text_labels, canvas}].
Each reaction gets a midmarker and two multimarkers (substrate side / product side); every drawn
metabolite glyph becomes a metabolite node (main metabolites: node_is_primary = true, cofactor
copies: false). Escher draws larger than this layout, so coordinates are scaled (default 3x).
"""
from __future__ import annotations

import json

from .layout import Layout, add, mul, unit
from .model import EMPTY, export_ids


def write_escher(L: Layout, path: str, scale: float = 3.0, marker: float = 20.0):
    pw = L.pathway
    SP, RX = export_ids(pw)          # original SBML ids where known (bigg_id = id in the input)
    nodes: dict[str, dict] = {}
    reactions: dict[str, dict] = {}
    counter = [0]

    def new_id() -> str:
        counter[0] += 1
        return str(counter[0])

    S = lambda p: {"x": round(p[0] * scale, 2), "y": round(p[1] * scale, 2)}

    # Escher draws text left-aligned from (label_x, label_y = baseline), in its own font sizes. Labels
    # are therefore placed with Escher's convention, independent of text length (ids or names):
    # right of the node -> start just right of it; above/below -> start at the node's left edge;
    # labels that were left of a node go below (or above) it.
    R_NODE = {"main": 20.0, "copy": 20.0, "side": 12.0, "reaction": 10.0}
    FONT = {"main": 20.0, "copy": 20.0, "side": 20.0, "reaction": 30.0}

    def label_xy(tg, node_pos, kind):
        x, y = node_pos[0] * scale, node_pos[1] * scale
        r_, fs = R_NODE[kind], FONT[kind]
        lb = L.labels.get(tg)
        if lb is None:
            return round(x + r_ + 6, 2), round(y + fs * 0.35, 2)
        cx, cy = lb.x + lb.w / 2 - node_pos[0], lb.y + lb.h / 2 - node_pos[1]
        if lb.anchor == "start" and cx > 0 and abs(cy) <= lb.h:          # right of the node
            return round(x + r_ + 6, 2), round(y + fs * 0.35, 2)
        if cy < 0 and abs(cy) >= abs(cx) * 0.3:                          # above
            return round(x - r_, 2), round(y - r_ - 8, 2)
        return round(x - r_, 2), round(y + r_ + fs + 2, 2)               # below

    # metabolite nodes: one per drawn glyph
    # Escher has no modifiers: species that only act as modifiers, modifier clones, modifier arcs and
    # the source/sink symbol are left out
    participants = {m for r in pw.reactions.values() for m in list(r.substrates) + list(r.products)}
    glyph_node: dict[str, str] = {}
    for m, p in sorted(L.met_pos.items()):
        if m not in participants:
            continue
        nid = new_id()
        lx, ly = label_xy(f"TG_M_{m}", p, "main")
        nodes[nid] = {"node_type": "metabolite", **S(p), "bigg_id": SP[m], "name": pw.metabolites[m].name,
                      "label_x": lx, "label_y": ly, "node_is_primary": True}
        glyph_node[Layout.sg(m)] = nid
    for c in sorted(L.clones, key=lambda c: c.glyph):
        if c.species == EMPTY or c.role == "modifier" or c.role.startswith("drain"):
            continue
        nid = new_id()
        lx, ly = label_xy("TG" + c.glyph[2:], c.pos, "copy" if c.far else "side")
        nodes[nid] = {"node_type": "metabolite", **S(c.pos), "bigg_id": SP[c.species],
                      "name": pw.metabolites[c.species].name, "label_x": lx, "label_y": ly,
                      "node_is_primary": bool(c.far)}
        glyph_node[c.glyph] = nid

    edges_by_rxn: dict[str, list] = {}
    for e in L.edges:
        edges_by_rxn.setdefault(e.rxn, []).append(e)

    for rid in sorted(L.rxn_pos):
        r = pw.reactions[rid]
        R, a = L.rxn_pos[rid], unit(L.rxn_axis[rid])
        # the reaction itself, then identical reactions of the input (aliases) on top of it
        versions = [(RX[rid], r.name, r.reversible, r.substrates, r.products)]
        versions += [(RX[(rid, k)], al.get("name", ""), al.get("reversible", r.reversible),
                      al.get("substrates", r.substrates), al.get("products", r.products))
                     for k, al in enumerate(r.aliases)]
        for bigg, rname, rrev, subs, prods in versions:
            mid, m_in, m_out = new_id(), new_id(), new_id()
            nodes[mid] = {"node_type": "midmarker", **S(R)}
            nodes[m_in] = {"node_type": "multimarker", **S(add(R, mul(a, -marker / scale)))}
            nodes[m_out] = {"node_type": "multimarker", **S(add(R, mul(a, marker / scale)))}
            segments = {new_id(): {"from_node_id": m_in, "to_node_id": mid, "b1": None, "b2": None},
                        new_id(): {"from_node_id": mid, "to_node_id": m_out, "b1": None, "b2": None}}
            for e in edges_by_rxn.get(rid, []):
                nid = glyph_node.get(e.glyph)
                if nid is None or (e.drawn or "").startswith("mod:"):
                    continue
                drawn = e.drawn or e.role
                sub_side = drawn in ("substrate", "sidesubstrate")
                seg = e.segs[0]
                b1 = b2 = None
                if seg[0] == "C" and len(e.segs) == 1:      # curve from reaction port to metabolite
                    b1, b2 = S(seg[2]), S(seg[3])
                if sub_side:     # metabolite -> substrate-side multimarker (control points reversed)
                    segments[new_id()] = {"from_node_id": nid, "to_node_id": m_in, "b1": b2, "b2": b1}
                else:            # product-side multimarker -> metabolite
                    segments[new_id()] = {"from_node_id": m_out, "to_node_id": nid, "b1": b1, "b2": b2}
            lx, ly = label_xy(f"TG_R_{rid}", R, "reaction")
            mets = [{"bigg_id": SP[m], "coefficient": -float(v)} for m, v in subs.items()]
            mets += [{"bigg_id": SP[m], "coefficient": float(v)} for m, v in prods.items()]
            reactions[new_id()] = {"name": rname, "bigg_id": bigg, "reversibility": bool(rrev),
                                   "label_x": lx, "label_y": ly, "gene_reaction_rule": "", "genes": [],
                                   "metabolites": mets, "segments": segments}

    _place_labels_escher(nodes, reactions, R_NODE, FONT)

    xs = [n["x"] for n in nodes.values()]
    ys = [n["y"] for n in nodes.values()]
    pad = 200
    canvas = {"x": min(xs) - pad, "y": min(ys) - pad,
              "width": max(xs) - min(xs) + 2 * pad, "height": max(ys) - min(ys) + 2 * pad}
    header = {"map_name": pw.name or pw.id, "map_id": L.meta.get("layout_hash", "sbml-autolayout"),
              "map_description": f"Generated by {L.meta.get('algorithm', 'sbml-autolayout')} "
                                 f"(layout {L.meta.get('layout_hash', '')})",
              "homepage": "https://escher.github.io",
              "schema": "https://escher.github.io/escher/jsonschema/1-0-0#"}
    with open(path, "w") as fh:
        json.dump([header, {"reactions": reactions, "nodes": nodes, "text_labels": {}, "canvas": canvas}],
                  fh, indent=1)


def _place_labels_escher(nodes: dict, reactions: dict, R_NODE: dict, FONT: dict):
    """Collision-aware label placement in Escher's own geometry.

    Escher draws each label left-aligned from (label_x, baseline label_y). For every label a few
    candidate spots around its node are scored against all drawn lines, node circles and labels that
    are already placed; the cheapest wins. Text width is estimated from the bigg_id (what Escher shows
    by default); right-hand spots stay valid when names are switched on, so they are preferred.
    """
    # polylines of all segments (Bezier curves sampled)
    lines = []
    for r in reactions.values():
        for sg in r["segments"].values():
            a, b = nodes[sg["from_node_id"]], nodes[sg["to_node_id"]]
            P0, P3 = (a["x"], a["y"]), (b["x"], b["y"])
            if sg["b1"] and sg["b2"]:
                P1, P2 = (sg["b1"]["x"], sg["b1"]["y"]), (sg["b2"]["x"], sg["b2"]["y"])
                pts = []
                for i in range(11):
                    t = i / 10
                    u = 1 - t
                    pts.append((u**3 * P0[0] + 3*u*u*t * P1[0] + 3*u*t*t * P2[0] + t**3 * P3[0],
                                u**3 * P0[1] + 3*u*u*t * P1[1] + 3*u*t*t * P2[1] + t**3 * P3[1]))
            else:
                pts = [P0, P3]
            lines += list(zip(pts, pts[1:]))

    def seg_hits_box(p, q, bx):
        x0, y0, x1, y1 = bx
        # Liang-Barsky clipping
        dx, dy = q[0] - p[0], q[1] - p[1]
        t0, t1 = 0.0, 1.0
        for pp, qq in ((-dx, p[0] - x0), (dx, x1 - p[0]), (-dy, p[1] - y0), (dy, y1 - p[1])):
            if pp == 0:
                if qq < 0:
                    return False
            else:
                t = qq / pp
                if pp < 0:
                    t0 = max(t0, t)
                else:
                    t1 = min(t1, t)
                if t0 > t1:
                    return False
        return True

    circles = []
    for n in nodes.values():
        if n["node_type"] == "metabolite":
            circles.append((n["x"], n["y"], R_NODE["main"] if n["node_is_primary"] else R_NODE["side"]))

    def box_circle(bx, c):
        x, y, r = c
        nx, ny = min(max(x, bx[0]), bx[2]), min(max(y, bx[1]), bx[3])
        return (nx - x) ** 2 + (ny - y) ** 2 < r * r

    placed = []

    def overlap(a, b):
        w = min(a[2], b[2]) - max(a[0], b[0])
        h = min(a[3], b[3]) - max(a[1], b[1])
        return w * h if w > 0 and h > 0 else 0.0

    # items: (priority, holder dict, x, y, r, fs, text)
    items = []
    for n in nodes.values():
        if n["node_type"] != "metabolite":
            continue
        kind = "main" if n["node_is_primary"] else "side"
        items.append((0 if kind == "main" else 2, n, n["x"], n["y"], R_NODE[kind], FONT[kind], n["bigg_id"]))
    for r in reactions.values():
        mid = next(nodes[x] for sg in r["segments"].values() for x in (sg["from_node_id"], sg["to_node_id"])
                   if nodes[x]["node_type"] == "midmarker")
        items.append((1, r, mid["x"], mid["y"], R_NODE["reaction"], FONT["reaction"], r["bigg_id"]))
    items.sort(key=lambda it: (it[0], it[2], it[3]))

    for _, holder, x, y, r, fs, text in items:
        w = 0.55 * fs * max(len(text), 3)
        g = 6
        cands = [  # (start x, baseline y, preference penalty)
            (x + r + g, y + 0.35 * fs, 0.0),                     # right
            (x + r * 0.7, y - r - g, 0.5),                       # above right
            (x + r * 0.7, y + r + fs, 0.5),                      # below right
            (x - w / 2, y - r - g, 1.0),                         # above centred
            (x - w / 2, y + r + fs + 2, 1.0),                    # below centred
            (x - r - g - w, y + 0.35 * fs, 1.5),                 # left
            (x - r * 0.7 - w, y - r - g, 2.0),                   # above left
            (x - r * 0.7 - w, y + r + fs, 2.0),                  # below left
            (x + r + g, y - r - g - fs * 0.3, 3.0),              # further out, diagonal
            (x + r + g, y + r + g + fs, 3.0),
        ]
        for k in (2.2, 3.4):                                     # second ring, further out
            d = r * k + g
            cands += [(x + d, y + 0.35 * fs, 4.0 * k), (x - d - w, y + 0.35 * fs, 4.5 * k),
                      (x - w / 2, y - d, 4.0 * k), (x - w / 2, y + d + fs, 4.0 * k),
                      (x + d * 0.7, y - d * 0.7, 4.5 * k), (x + d * 0.7, y + d * 0.7 + fs, 4.5 * k)]
        old = (holder.get("label_x"), holder.get("label_y"))
        best = None
        for sx, by, pref in cands:
            bx = (sx - 2, by - 0.8 * fs, sx + w + 2, by + 0.25 * fs)
            cost = pref
            cost += 40 * sum(1 for p, q in lines if seg_hits_box(p, q, bx))
            cost += 60 * sum(1 for c in circles if box_circle(bx, c))
            cost += sum(150 + overlap(bx, o) / 20 for o in placed if overlap(bx, o) > 0)
            if old[0] is not None and abs(sx - old[0]) < 1 and abs(by - old[1]) < 1:
                cost -= 0.3                                  # tie-break: keep the layout's choice
            if best is None or cost < best[0]:
                best = (cost, sx, by, bx)
        assert best is not None
        _, sx, by, bx = best
        holder["label_x"], holder["label_y"] = round(sx, 2), round(by, 2)
        placed.append(bx)
