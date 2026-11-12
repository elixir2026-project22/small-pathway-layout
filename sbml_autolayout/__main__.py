"""Command line: python -m sbml_autolayout INPUT [INPUT ...] -o out/name [--previous old.sbml] [options]

INPUT is a pathway JSON file or one or more SBML files (.xml/.sbml), which are merged."""
from __future__ import annotations

import argparse
import re
import os

from .layout import Params, compute_layout
from .model import load_pathway
from .sbml import read_previous, write_sbml
from .sbml_import import import_sbml, list_subsystems
from .escher import import_escher, is_escher
from .escher_export import write_escher
from .svg import write_svg


def add_layout_options(ap: argparse.ArgumentParser):
    """Options shared by the single-map command and the batch runner."""
    ap.add_argument("--subsystem", action="append", default=[],
                    help="only reactions whose subsystem/group contains this text (repeatable)")
    ap.add_argument("--keep-exchange", action="store_true",
                    help="also keep exchange/demand/sink reactions (one side empty and no modifier); off by "
                         "default. Regulated synthesis/degradation (e.g. transcription '-> protein' with a "
                         "transcription factor as modifier) is always kept")
    ap.add_argument("--no-name-match", action="store_true", help="merge species by annotation/id only")
    ap.add_argument("--pseudo", choices=("drain", "hide", "draw"), default="drain",
                    help="biomass / pseudo-reactions (SBO:0000629, 'biomass' in id or name, or >= 8 main "
                         "partners): 'drain' = a short arrow from each precursor to a small biomass symbol "
                         "next to it (default), 'hide', or 'draw' as a normal reaction node")
    ap.add_argument("--modifiers", action="store_true",
                    help="draw modifiers (catalysts, transcription factors, inhibitors) as SBGN arcs and keep "
                         "regulated synthesis/degradation ('-> protein'); off by default (modifiers stay in "
                         "the SBML output but are not drawn). Escher has no modifiers")
    ap.add_argument("--no-modifiers", action="store_true", help=argparse.SUPPRESS)   # old switch, now the default
    ap.add_argument("--show-protons", action="store_true",
                    help="draw H+ as a cofactor; off by default (protons are hidden)")
    ap.add_argument("--max-reactions", type=int, default=150,
                    help="refuse larger inputs (use --subsystem); this layout is meant for pathway-sized maps")
    ap.add_argument("--orientation", choices=["TB", "LR"], default="TB")
    ap.add_argument("--cycle-direction", choices=["cw", "ccw"], default="cw")
    ap.add_argument("--branch-side", choices=["right", "left"], default="right")
    ap.add_argument("--spacing", type=float, default=110.0)
    ap.add_argument("--no-reaction-labels", action="store_true")
    ap.add_argument("--main-path", default="",
                    help='metabolites of the main line, separated by ">" (ids or names), e.g. '
                         '"glc__D > pyr" or "D-glucose > pyruvate". Missing steps are filled in along the '
                         'network. Default: known pathways (glycolysis), else the longest chain. "off" = '
                         'always the longest chain')
    ap.add_argument("--fold", type=int, default=0, help="max metabolites per straight run before folding (0 = auto)")
    ap.add_argument("--labels", choices=["name", "id"], default="name", help="metabolite labels: names or ids")
    ap.add_argument("--hub-degree", type=int, default=4, help="metabolites in >= N reactions may become clones")
    ap.add_argument("--metrics", action="store_true", help="print layout quality metrics")
    ap.add_argument("--density", type=float, default=1.0,
                    help="how compact: >1 less white space (e.g. 1.3), <1 more air (e.g. 0.8)")
    ap.add_argument("--grid", type=float, default=0.5,
                    help="grid cell for metabolites in units of the spacing (reactions use half); 0 = off")
    ap.add_argument("--notation", choices=["default", "sbgn"], default="default",
                    help="default: small circles with the names beside them; sbgn: SBGN Process Description "
                         "(names inside shapes chosen by SBO term, process nodes with ports, clone markers)")
    ap.add_argument("--box-max", type=float, default=110.0,
                    help="sbgn: largest width of a name box in px; longer names wrap onto a second line "
                         "and are shortened with '…' after that")
    ap.add_argument("--box-lines", type=int, default=2, help="sbgn: max lines of text in a box")
    ap.add_argument("--cycles", choices=["known", "all", "none"], default="known",
                    help="draw as circles: curated cycles only (TCA, urea, glyoxylate, methionine, Calvin), "
                         "every loop found, or none")


def params_from_args(a) -> Params:
    return Params(orientation=a.orientation, cycle_direction=a.cycle_direction, branch_side=a.branch_side,
               spacing=a.spacing, reaction_labels=not a.no_reaction_labels, fold=a.fold,
               modifiers=a.modifiers and not a.no_modifiers, pseudo=a.pseudo,
               main_path=[x.strip() for x in re.split(r"[>;]", a.main_path) if x.strip()],
               label_mode=a.labels, hub_degree=a.hub_degree, cycles=a.cycles,
               density=a.density, grid=a.grid, notation=a.notation, box_max=a.box_max, box_lines=a.box_lines)


def run(inputs: list[str], out: str, a, verbose: bool = True):
    """Lay out one map (one JSON file, or one or more SBML files merged) and write all outputs.
    Raises ValueError with a readable message on bad input. Returns the Layout."""
    P = params_from_args(a)
    sbml_in = [p for p in inputs if not p.lower().endswith(".json")]
    if sbml_in:
        if len(sbml_in) != len(inputs):
            raise ValueError("give either one JSON file or SBML files, not both")
        pw, rep = import_sbml(sbml_in, a.subsystem, a.keep_exchange, not a.no_name_match, a.show_protons,
                              keep_regulated=a.modifiers and not a.no_modifiers)
        if verbose:
            print(rep.summary())
    else:
        if len(inputs) > 1:
            raise ValueError("only one JSON pathway at a time")
        pw = import_escher(inputs[0])[0] if is_escher(inputs[0]) else load_pathway(inputs[0])
    if not pw.reactions:
        raise ValueError("no reactions left after filtering" + (" (check --list-subsystems)" if a.subsystem else ""))
    if len(pw.reactions) > a.max_reactions:
        raise ValueError(f"{len(pw.reactions)} reactions: too large for a single pathway map; select part of "
                         "the model with --subsystem (see --list-subsystems) or raise --max-reactions")
    prev = read_previous(a.previous) if getattr(a, "previous", None) else None
    L = compute_layout(pw, P, prev)

    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    write_sbml(L, out + ".sbml")
    write_svg(L, out + ".svg")
    write_escher(L, out + ".escher.json")
    if verbose and a.metrics:
        from .metrics import layout_metrics
        for k, v in layout_metrics(L).items():
            print(f"  {k}: {v}")
    if verbose:
        print(f"{out}: {len(L.met_pos)} metabolites, {len(L.rxn_pos)} reactions, {len(L.clones)} clones, "
          f"{len(L.cycles)} cycle(s), backbone={L.backbone}, layout_hash={L.meta['layout_hash']}")
    return L


def main(argv=None):
    ap = argparse.ArgumentParser(prog="sbml-autolayout", description=__doc__)
    ap.add_argument("inputs", nargs="+", help="pathway JSON, or one or more SBML files to merge")
    ap.add_argument("-o", "--out", help="output prefix (writes .sbml, .svg, .escher.json)")
    ap.add_argument("--list-subsystems", action="store_true", help="list subsystems in the SBML input(s) and exit")
    ap.add_argument("--previous", help="previous SBML layout; its positions are kept fixed")
    add_layout_options(ap)
    a = ap.parse_args(argv)
    if a.list_subsystems:
        for name, n in list_subsystems([p for p in a.inputs if not p.lower().endswith(".json")]).items():
            print(f"{n:5d}  {name}")
        return
    if not a.out:
        ap.error("-o/--out is required")
    try:
        run(a.inputs, a.out, a)
    except ValueError as e:
        ap.error(str(e))
    except Exception as e:  # unreadable / malformed input
        ap.error(f"could not process {', '.join(a.inputs)}: {type(e).__name__}: {e}")


if __name__ == "__main__":
    main()
