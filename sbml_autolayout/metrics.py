"""Layout quality metrics (used to compare layouts and to check changes against a curated map)."""
from __future__ import annotations

import math

from .layout import Layout, bezier_points, segs_cross


def _poly(segs):
    pts = []
    for s in segs:
        pts += bezier_points(s, 8)[(1 if pts else 0):]
    return pts


def _crossings(polys):
    """polys: list of (owner_set, polyline). Counts pairs of edges that properly cross."""
    n = 0
    boxes = []
    for owners, pl in polys:
        xs, ys = [p[0] for p in pl], [p[1] for p in pl]
        boxes.append((min(xs), min(ys), max(xs), max(ys)))
    for i in range(len(polys)):
        for j in range(i + 1, len(polys)):
            if polys[i][0] & polys[j][0]:
                continue
            a, b = boxes[i], boxes[j]
            if a[2] < b[0] or b[2] < a[0] or a[3] < b[1] or b[3] < a[1]:
                continue
            p, q = polys[i][1], polys[j][1]
            if any(segs_cross(p[k], p[k + 1], q[m], q[m + 1]) for k in range(len(p) - 1) for m in range(len(q) - 1)):
                n += 1
    return n


def _ortho(p0, p1, tol_deg=5.0):
    ang = math.degrees(math.atan2(abs(p1[1] - p0[1]), abs(p1[0] - p0[0])))
    return ang < tol_deg or ang > 90 - tol_deg


def layout_metrics(L: Layout, curator_primary: dict | None = None) -> dict:
    P = L.params
    main = [e for e in L.edges if e.role in ("substrate", "product") and not e.glyph.count("__R_")]
    polys = [({e.rxn, e.species}, _poly(e.segs)) for e in main]
    lens = [math.dist(e.segs[0][1], e.segs[-1][-1]) for e in main]
    labels = [lb for lb in L.labels.values() if not lb.inside]     # sbgn: names inside nodes are part of them
    lab_ov = sum(1 for i in range(len(labels)) for j in range(i + 1, len(labels))
                 if not (labels[i].x + labels[i].w <= labels[j].x or labels[j].x + labels[j].w <= labels[i].x or
                         labels[i].y + labels[i].h <= labels[j].y or labels[j].y + labels[j].h <= labels[i].y))
    # straight main edges outside cycles: are they horizontal/vertical?
    straight = [e for e in main if len(e.segs) == 1 and e.segs[0][0] == "L"]
    out = {
        "main_nodes": len(L.met_pos),
        "reactions": len(L.rxn_pos),
        "clones": len(L.clones),
        "main_edge_crossings": _crossings(polys),
        # modifier arcs (catalysts, regulators) crossing main lines or each other
        "modifier_crossings": _crossings(polys + [({e.rxn, e.species}, _poly(e.segs)) for e in L.edges
                                                  if (e.drawn or "").startswith("mod:")]) - _crossings(polys),
        "mean_main_edge_len/D": round(sum(lens) / len(lens) / P.spacing, 2) if lens else 0,
        "max_main_edge_len/D": round(max(lens) / P.spacing, 2) if lens else 0,
        "orthogonal_straight_edges_%": round(100 * sum(_ortho(e.segs[0][1], e.segs[0][2]) for e in straight)
                                             / max(1, len(straight)), 1),
        "label_overlaps": lab_ov,
        "width": round(L.width), "height": round(L.height),
        "aspect(w/h)": round(L.width / L.height, 2),
    }
    if curator_primary:
        tp = fp = fn = 0
        for rid, prim in curator_primary.items():
            ro = L.roles.get(rid)
            if not ro:
                continue
            ours = set(ro.main_sub + ro.main_prod) | L.far_links.get(rid, set())
            cur = set(prim)
            tp += len(ours & cur)
            fp += len(ours - cur)
            fn += len(cur - ours)
        out["main/side vs curator: precision"] = round(tp / max(1, tp + fp), 2)
        out["main/side vs curator: recall"] = round(tp / max(1, tp + fn), 2)
    return out


def escher_reference_metrics(path: str) -> dict:
    import json
    body = json.load(open(path))[1]
    N, R = body["nodes"], body["reactions"]
    polys, lens, straight = [], [], []
    for rid, r in R.items():
        for s in r["segments"].values():
            a, b = N[s["from_node_id"]], N[s["to_node_id"]]
            if any(n["node_type"] == "metabolite" and not n["node_is_primary"] for n in (a, b)):
                continue
            pa, pb = (a["x"], a["y"]), (b["x"], b["y"])
            if s.get("b1"):
                seg = ("C", pa, (s["b1"]["x"], s["b1"]["y"]), (s["b2"]["x"], s["b2"]["y"]), pb)
            else:
                seg = ("L", pa, pb)
                straight.append(_ortho(pa, pb))
            owners = {rid} | {s["from_node_id"], s["to_node_id"]}
            polys.append((owners, bezier_points(seg, 8)))
    xs = [n["x"] for n in N.values()]
    ys = [n["y"] for n in N.values()]
    w, h = max(xs) - min(xs), max(ys) - min(ys)
    return {"main_nodes": sum(1 for n in N.values() if n["node_type"] == "metabolite" and n["node_is_primary"]),
            "reactions": len(R),
            "clones (secondary nodes)": sum(1 for n in N.values() if n["node_type"] == "metabolite" and not n["node_is_primary"]),
            "main_segment_crossings": _crossings(polys),
            "orthogonal_straight_segments_%": round(100 * sum(straight) / max(1, len(straight)), 1),
            "aspect(w/h)": round(w / h, 2)}


def overlap_report(L: Layout) -> dict:
    """Hard overlaps that must be zero: nodes on nodes, labels on labels, labels on nodes."""
    P = L.params
    # (name, centre, half width, half height, round?)
    nodes = [(f"M:{k}", p, L.met_dims(k)[0] / 2, L.met_dims(k)[1] / 2, not L.sbgn) for k, p in L.met_pos.items()]
    nodes += [(f"R:{k}", p, P.rxn_half * 1.2, P.rxn_half * 1.2, True) for k, p in L.rxn_pos.items()]
    nodes += [(f"C:{c.glyph}", c.pos, L.clone_dims(c)[0] / 2, L.clone_dims(c)[1] / 2,
               L.clone_dims(c)[3] == "circle") for c in L.clones]

    def touch(a, b) -> bool:
        if a[4] and b[4]:
            return math.dist(a[1], b[1]) < a[2] + b[2] + 2
        return abs(a[1][0] - b[1][0]) < a[2] + b[2] + 1 and abs(a[1][1] - b[1][1]) < a[3] + b[3] + 1
    nn = [(a[0], b[0]) for i, a in enumerate(nodes) for b in nodes[i + 1:] if touch(a, b)]
    labels = [lb for lb in L.labels.values() if not lb.inside]     # sbgn: names inside nodes are part of them
    box = lambda l: (l.x, l.y, l.x + l.w, l.y + l.h)
    ll = [(a.glyph, b.glyph) for i, a in enumerate(labels) for b in labels[i + 1:]
          if not (box(a)[2] <= box(b)[0] or box(b)[2] <= box(a)[0] or box(a)[3] <= box(b)[1] or box(b)[3] <= box(a)[1])]
    ln = []
    for l in labels:
        x0, y0, x1, y1 = box(l)
        for name, p, hx, hy, round_ in nodes:
            if round_:
                cx, cy = min(max(p[0], x0), x1), min(max(p[1], y0), y1)
                if math.dist((cx, cy), p) < hx:
                    ln.append((l.glyph, name))
            elif x0 < p[0] + hx and p[0] - hx < x1 and y0 < p[1] + hy and p[1] - hy < y1:
                ln.append((l.glyph, name))
    return {"node_node": nn, "label_label": ll, "label_node": ln}
