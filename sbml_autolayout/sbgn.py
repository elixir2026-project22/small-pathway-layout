"""SBGN Process Description notation: entity pools are drawn as shapes with the name inside.

Shapes follow the SBGN PD glyphs, chosen by the species' SBO term:
  simple chemical (default) .......... stadium (rounded ends)
  macromolecule (proteins, enzymes) .. rounded rectangle
  complex ............................ rectangle with cut corners
  nucleic acid feature (genes, RNA) .. rectangle with rounded bottom corners
  unspecified entity ................. ellipse
  phenotype .......................... hexagon
Cofactor copies are small stadiums with the SBGN clone marker (filled lower band).
Box widths are limited (Params.box_max); long names wrap onto up to Params.box_lines lines and are
shortened with "…" after that.
"""
from __future__ import annotations

import math
import re
from functools import lru_cache

# SBO term -> SBGN entity pool class
SBO_CLASS = {
    "SBO:0000247": "simple chemical", "SBO:0000327": "simple chemical", "SBO:0000328": "simple chemical",
    "SBO:0000299": "simple chemical", "SBO:0000456": "simple chemical",
    "SBO:0000245": "macromolecule", "SBO:0000246": "macromolecule", "SBO:0000248": "macromolecule",
    "SBO:0000249": "macromolecule", "SBO:0000252": "macromolecule", "SBO:0000460": "macromolecule",
    "SBO:0000253": "complex", "SBO:0000296": "complex", "SBO:0000297": "complex",
    "SBO:0000250": "nucleic acid feature", "SBO:0000251": "nucleic acid feature",
    "SBO:0000243": "nucleic acid feature", "SBO:0000278": "nucleic acid feature",
    "SBO:0000354": "nucleic acid feature", "SBO:0000334": "nucleic acid feature",
    "SBO:0000335": "nucleic acid feature",
    "SBO:0000285": "unspecified entity",
    "SBO:0000358": "phenotype",
}

FONT = {"main": 10.5, "copy": 10.0, "side": 9.0}
LINE = {"main": 12.5, "copy": 12.0, "side": 10.5}
PAD_X, PAD_Y = 5.0, 3.5


def entity_class(sbo: str) -> str:
    return SBO_CLASS.get((sbo or "").strip(), "simple chemical")


def char_w(size: float) -> float:
    return 0.5 * size         # Helvetica/Arial lower case averages ~0.5 em


def _tokens(text: str) -> list[str]:
    """split after spaces and hyphens (chemical names break there; not inside "1,6")"""
    return [t for t in re.findall(r"[^\s\-]*[\s\-]*", text) if t]


def _splits(n: int, k: int):
    """all ways to cut n tokens into k non-empty consecutive groups (as cut positions)"""
    if k == 1:
        yield ()
        return
    for i in range(1, n - k + 2):
        for rest in _splits(n - i, k - 1):
            yield (i,) + tuple(i + r for r in rest)


def wrap(text: str, size: float, max_w: float, max_lines: int) -> list[str]:
    """Balanced line breaks: the longest line as short as possible, breaking at spaces rather than
    hyphens. A single word may stick out a little (15 %) instead of being cut; what still does not
    fit in max_lines is shortened with "…"."""
    limit = max(4, int(max_w / char_w(size)))
    text = text.strip()
    if len(text) <= limit or (len(text) <= limit * 1.15 and len(_tokens(text)) <= 2):
        return [text]                 # fits (a slightly long single word is not split)
    toks = _tokens(text)
    best = None
    if len(toks) <= 40:
        for k in range(2, max_lines + 1):
            if k > len(toks):
                break
            for cuts in _splits(len(toks), k):
                b = (0,) + cuts + (len(toks),)
                lines = ["".join(toks[b[j]:b[j + 1]]).strip() for j in range(k)]
                longest = max(len(x) for x in lines)
                if longest > limit * 1.15:
                    continue
                pen = sum(2 for x in lines[:-1] if x.endswith("-"))      # hyphen breaks: slightly worse
                pen += sum(12 for x in lines if len(x) <= 3)               # no "3-" / "L-" alone on a line
                cost = (max(longest - limit, 0) * 3 + longest + pen + 3 * (k - 2), )
                if best is None or cost < best[0]:
                    best = (cost, lines)
            if best is not None:
                break
    if best is not None:
        return best[1]
    # too long for max_lines: fill the lines greedily (very long words are split with "-" where the
    # line is full) and shorten the last one
    lines, cur = [], ""
    for t in toks:
        if len((cur + t).rstrip()) <= limit:
            cur += t
            continue
        if len(t.rstrip()) > limit or limit - len(cur) >= 6:     # split the word to use the line
            while len((cur + t).rstrip()) > limit:
                room = limit - len(cur) - 1
                if room < 3:
                    lines.append(cur.rstrip())
                    cur = ""
                    continue
                lines.append(cur + t[:room] + "-")
                cur, t = "", t[room:]
            cur += t
        else:
            lines.append(cur.rstrip())
            cur = t
    if cur:
        lines.append(cur.rstrip())
    if len(lines) > max_lines:
        last = lines[max_lines - 1]
        lines = lines[:max_lines - 1] + [last[: max(1, limit - 1)].rstrip(" ,-") + "…"]
    return lines


@lru_cache(maxsize=None)
def box(text: str, kind: str, cls: str, max_w: float, max_lines: int) -> tuple[float, float, tuple]:
    """(width, height, lines) of the node for this name. kind: main / copy (a far partner drawn
    again) / side (cofactor copy, one line, smaller)"""
    size = FONT[kind]
    if kind == "side":
        max_w, max_lines = min(max_w, 72.0), 1
    elif kind == "copy":
        max_w = min(max_w, 90.0)
    ends = (LINE[kind] + 2 * PAD_Y) * 0.35 if cls == "simple chemical" else 0.0
    lines = wrap(text, size, max_w - 2 * PAD_X - ends, max_lines)
    tw = max(len(s) for s in lines) * char_w(size)
    h = len(lines) * LINE[kind] + 2 * PAD_Y
    w = tw + 2 * PAD_X
    if cls == "simple chemical":
        w += h * 0.35                    # room for the round ends
    w = max(w, h * (1.0 if cls == "simple chemical" else 1.6), 22.0 if kind != "side" else 18.0)
    w = min(w, max_w * 1.12)             # a long single word may stick out a little, not more
    return round(w, 1), round(h, 1), tuple(lines)


def border(cls: str, hx: float, hy: float, u: tuple[float, float]) -> float:
    """distance from the centre to the outline in direction u (unit vector)"""
    ux, uy = abs(u[0]), abs(u[1])
    if ux < 1e-9 and uy < 1e-9:
        return 0.0
    if cls == "unspecified entity":
        return 1.0 / math.sqrt((ux / hx) ** 2 + (uy / hy) ** 2)
    rect = min(hx / ux if ux > 1e-9 else math.inf, hy / uy if uy > 1e-9 else math.inf)
    if cls != "simple chemical":
        return rect
    # stadium: straight middle part, half circles of radius hy at both ends
    a = max(0.0, hx - hy)
    if ux * rect <= a + 1e-9:
        return rect
    # ray from the origin hits the circle centred at (a, 0) with radius hy
    b = ux * a
    disc = b * b - (a * a - hy * hy)
    return b + math.sqrt(max(0.0, disc))


def outline_svg(cls: str, x: float, y: float, w: float, h: float, css: str, title: str = "") -> str:
    """SVG element(s) for an entity pool centred at (x, y)"""
    x0, y0 = x - w / 2, y - h / 2
    t = f"<title>{title}</title>" if title else ""
    if cls == "simple chemical":
        r = h / 2
        return f'<rect class="{css}" x="{x0:.1f}" y="{y0:.1f}" width="{w:.1f}" height="{h:.1f}" rx="{r:.1f}">{t}</rect>'
    if cls == "macromolecule":
        r = min(6.0, h / 3)
        return f'<rect class="{css}" x="{x0:.1f}" y="{y0:.1f}" width="{w:.1f}" height="{h:.1f}" rx="{r:.1f}">{t}</rect>'
    if cls == "complex":
        c = min(6.0, h / 3)
        pts = [(x0 + c, y0), (x0 + w - c, y0), (x0 + w, y0 + c), (x0 + w, y0 + h - c), (x0 + w - c, y0 + h),
               (x0 + c, y0 + h), (x0, y0 + h - c), (x0, y0 + c)]
        return f'<path class="{css}" d="M' + " L".join(f"{a:.1f},{b:.1f}" for a, b in pts) + f' Z">{t}</path>'
    if cls == "nucleic acid feature":
        r = min(7.0, h / 2.5)
        d = (f"M{x0:.1f},{y0:.1f} H{x0 + w:.1f} V{y0 + h - r:.1f} A{r:.1f},{r:.1f} 0 0 1 {x0 + w - r:.1f},{y0 + h:.1f} "
             f"H{x0 + r:.1f} A{r:.1f},{r:.1f} 0 0 1 {x0:.1f},{y0 + h - r:.1f} Z")
        return f'<path class="{css}" d="{d}">{t}</path>'
    if cls == "unspecified entity":
        return f'<ellipse class="{css}" cx="{x:.1f}" cy="{y:.1f}" rx="{w / 2:.1f}" ry="{h / 2:.1f}">{t}</ellipse>'
    if cls == "phenotype":
        c = h / 2
        pts = [(x0, y), (x0 + c, y0), (x0 + w - c, y0), (x0 + w, y), (x0 + w - c, y0 + h), (x0 + c, y0 + h)]
        return f'<path class="{css}" d="M' + " L".join(f"{a:.1f},{b:.1f}" for a, b in pts) + f' Z">{t}</path>'
    return f'<rect class="{css}" x="{x0:.1f}" y="{y0:.1f}" width="{w:.1f}" height="{h:.1f}">{t}</rect>'


def clone_marker_svg(cls: str, x: float, y: float, w: float, h: float, css: str) -> str:
    """SBGN clone marker: the lower part of the glyph filled (here: the lower 30 %)"""
    x0, y0 = x - w / 2, y - h / 2
    yb = y0 + h * 0.70
    clip = f"clip-{abs(hash((round(x, 1), round(y, 1)))) % 10**9}"
    shape = outline_svg(cls, x, y, w, h, "").replace(' class=""', "")
    return (f'<clipPath id="{clip}">{shape}</clipPath>'
            f'<rect class="{css}" x="{x0:.1f}" y="{yb:.1f}" width="{w:.1f}" height="{h * 0.30 + 0.5:.1f}" '
            f'clip-path="url(#{clip})"/>')
