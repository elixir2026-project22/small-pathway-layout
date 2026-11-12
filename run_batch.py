"""Lay out every pathway in a folder at once.

EASIEST: put your SBML (.xml/.sbml) or Escher/JSON files into the `examples` folder - or point
SETTINGS["inputs"] below to a folder anywhere on disk - then run

    uv run run_batch.py

Results go to `output/<folder name>/`: per map .svg, .sbml (SBML Layout), .escher.json, plus
`index.html` (all maps side by side - open it in a browser) and `summary.csv`.

Command line use (overrides SETTINGS), e.g.
    uv run sbml-autolayout-batch ~/data/my_examples                 # -> output/my_examples/
    uv run sbml-autolayout-batch ~/data/a ~/data/b one.xml -o ~/results --density 1.3 -j 4
    uv run sbml-autolayout-batch --help                             # all options
A subfolder of an input folder (e.g. examples/my_pathway/a.xml + b.xml) is merged into ONE map.
"""
import sys

# ----------------------------------------------------------------------------- SETTINGS
SETTINGS = {
    # files and/or folders, anywhere on disk, e.g. ["/home/anja/data/brenda_examples", "~/Desktop/tca.xml"]
    "inputs": ["examples"],
    "outdir": None,               # None = ./output/<input folder name>; or e.g. "/home/anja/results"
    "density": 1.0,               # >1 = more compact (e.g. 1.3), <1 = more white space (e.g. 0.8)
    "grid": 0.5,                  # grid cell in units of the spacing; 0 = no grid
    "orientation": "TB",          # "TB" top->bottom or "LR" left->right
    "main_path": "",              # main line, e.g. "D-glucose > pyruvate" (ids or names, ">" separated;
                                  # missing steps are filled in). "" = known pathways (glycolysis, see
                                  # sbml_autolayout/known.py CANONICAL_PATHS), else the longest chain; "off" = longest
    "fold": 0,                    # max metabolites in one straight run before the chain turns into the
                                  # next column (U-turn); 0 = auto (chains > 15 are split evenly into runs
                                  # of <= 12). Small value = wider/flatter map, large value = long and narrow
    "hub_degree": 4,              # hub rule, only for reactions WITHOUT main/side info in the source:
                                  # a metabolite in >= N reactions may become a cofactor clone. 999 = off
    "cycle_direction": "cw",      # "cw" or "ccw"
    "cycles": "known",            # "known" (curated cycles only), "all" or "none"
    "labels": "name",             # metabolite labels: "name" or "id"
    "reaction_labels": True,
    "modifiers": False,           # draw catalysts / transcription factors / inhibitors (SBGN arcs)
    "pseudo": "drain",            # biomass/pseudo-reactions: "drain", "hide" or "draw"
    "keep_exchange": False,       # also keep exchange/demand/sink reactions (regulated ones are always kept)
    "show_protons": False,        # draw H+ as a cofactor
    "merge_subfolders": True,     # each subfolder = one merged map
    "jobs": 1,                    # maps computed in parallel
}
# --------------------------------------------------------------------------------------


def settings_to_argv(s):
    argv = list(s["inputs"]) + (["-o", s["outdir"]] if s["outdir"] else []) + ["--density", str(s["density"]), "--grid", str(s["grid"]),
                                "--orientation", s["orientation"], "--fold", str(s["fold"]),
                                "--main-path", s["main_path"],
                                "--hub-degree", str(s["hub_degree"]), "--cycle-direction", s["cycle_direction"],
                                "--cycles", s["cycles"], "--labels", s["labels"], "-j", str(s["jobs"])]
    if not s["reaction_labels"]:
        argv.append("--no-reaction-labels")
    if s["modifiers"]:
        argv.append("--modifiers")
    argv += ["--pseudo", s["pseudo"]]
    if s["keep_exchange"]:
        argv.append("--keep-exchange")
    if s["show_protons"]:
        argv.append("--show-protons")
    if s["merge_subfolders"]:
        argv.append("--merge-subfolders")
    return argv


if __name__ == "__main__":
    import os
    here = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, here)
    from sbml_autolayout.batch import main
    if len(sys.argv) > 1:
        sys.exit(main(sys.argv[1:]))
    os.chdir(here)                 # SETTINGS paths are relative to this folder
    sys.exit(main(settings_to_argv(SETTINGS)))
