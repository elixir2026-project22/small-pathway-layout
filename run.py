"""Lay out ONE pathway (same as `uv run sbml-autolayout ...`) (one JSON/Escher file, or several SBML files merged into one map).

    uv run sbml-autolayout examples/brenda_citric_acid_cycle.xml -o output/tca
    uv run sbml-autolayout part1.xml part2.xml -o output/merged          # merge several SBML files
    uv run sbml-autolayout model.xml --list-subsystems                    # big models: pick a part ...
    uv run sbml-autolayout model.xml --subsystem "Glycolysis" -o output/glyc
    uv run sbml-autolayout examples/x.xml -o output/x_v2 --previous output/x.sbml   # keep an earlier layout fixed
    uv run sbml-autolayout --help                                         # all options
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sbml_autolayout.__main__ import main

if __name__ == "__main__":
    main()
