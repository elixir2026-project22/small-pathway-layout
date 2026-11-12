# sbml-autolayout – automatic layout for metabolic pathway maps

## Setup with uv (once)

Install uv first if needed: https://docs.astral.sh/uv/ . Then, in this folder:

    uv sync

This creates `.venv/` with Python and networkx (the only dependency). It also installs two commands:
- `sbml-autolayout`: one map;
- `sbml-autolayout-batch`: many maps or whole folders.

Prefix every command with `uv run`. No manual activation is needed.

## Many maps at once – data from any folder

    uv run sbml-autolayout-batch ~/data/my_examples

- Every `.xml`, `.sbml` (SBML) and `.json` (Escher) file in that folder becomes its own map.
- Results go to `output/my_examples/` in this project. Your data folder is never written to.
- Open `output/my_examples/index.html` in a browser to see all maps side by side with their statistics. `summary.csv` has the same numbers.

More examples:

    uv run sbml-autolayout-batch ~/data/brenda ~/data/reactome extra.xml -o ~/results   # several sources, own output folder
    uv run sbml-autolayout-batch ~/data/my_examples --density 1.3 -j 4                  # more compact, 4 maps in parallel
    uv run sbml-autolayout-batch ~/data/models --merge-subfolders                       # each subfolder = one merged map
    uv run sbml-autolayout-batch --help

Windows paths work too, e.g. `uv run sbml-autolayout-batch "C:\Users\anja\data\examples"`.

**Without typing paths:** set your data folder once in the `SETTINGS` block at the top of `run_batch.py`, e.g.

    "inputs": ["/home/anja/data/brenda_examples"]

The same block holds density, grid, orientation and the other options. Then just run:

    uv run run_batch.py

The bundled `examples/` folder (your four test files) is the default. Its results are already in `output/examples/`.

Output per map:

| File | Content |
|---|---|
| `NAME.svg` | picture |
| `NAME.sbml` | SBML with SBML Layout (CellDesigner, Newt, …) |
| `NAME.escher.json` | Escher map (Map → Load map JSON) |

A file that fails is listed as FAILED with the reason; all other maps are still produced.

## One map

    uv run sbml-autolayout ~/data/tca.xml -o output/tca
    uv run sbml-autolayout part1.xml part2.xml -o output/merged          # merge SBML files into one map
    uv run sbml-autolayout model.xml --list-subsystems                    # large models: choose a part
    uv run sbml-autolayout model.xml --subsystem "Glycolysis" -o output/glyc

## Adjusting a map by hand and keeping it

1. Open `NAME.sbml` in an SBML-Layout editor and move nodes or labels.
2. Re-run with that file as `--previous`. Everything already placed stays exactly where it is; only new parts are added:

       uv run sbml-autolayout ~/data/NAME.xml -o output/NAME_v2 --previous output/my_examples/NAME.sbml

## Main options (both commands)

| Option | Effect |
|---|---|
| `--density 1.3` / `0.8` | more compact / more white space |
| `--grid 0.5` (`0` = off) | grid cell for metabolites; reactions use half a cell |
| `--orientation TB\|LR` | main direction top→bottom or left→right |
| `--main-path "A > B"` | main line (ids or names, `>` separated, missing steps are filled in), e.g. `"D-glucose > pyruvate"`. Default: a known pathway if present (glycolysis), else the longest chain. `off` = always the longest chain |
| `--fold N` (`0` = auto) | max metabolites in one straight run before the chain makes a U-turn into the next column; auto splits chains longer than 15 evenly into runs of ≤ 12. Smaller N = wider, flatter map |
| `--cycle-direction cw\|ccw` | direction of drawn cycles |
| `--cycles known\|all\|none` | draw only curated cycles (TCA, urea, methionine, Calvin) as circles, every loop, or none |
| `--modifiers` | draw catalysts / transcription factors / inhibitors (SBML modifiers); off by default |
| `--notation sbgn` | SBGN Process Description: names inside shapes (default: small circles, names beside them) |
| `--box-max 110` / `--box-lines 2` | SBGN only: largest width of a name box (px) and number of text lines |
| `--pseudo drain\|hide\|draw` | biomass / pseudo-reactions: small drains next to the precursors (default), hidden, or a normal node |
| `--keep-exchange` | also keep exchange/demand/sink reactions (one side empty, no modifier); off by default |
| `--show-protons` | draw H⁺ as a cofactor; off by default |
| `--labels name\|id` | metabolite labels |
| `--no-reaction-labels` | hide enzyme/reaction names |
| `--hub-degree 4` | only for reactions where the source gives no main/side info: a metabolite in ≥ N reactions (and ≥ 2× as many as its co-substrate) may become a cofactor clone, never if that disconnects the map. `999` = off |
| `--keep-exchange`, `--show-protons` | keep exchange reactions / H⁺ |
| `-j 4` (batch only) | compute several maps in parallel |

## Which chain is the main line?

The main line is drawn as one straight run (folded only by `--fold`); everything else hangs off it.

1. **`--main-path`**, if you give one.
2. **A known pathway** from `CANONICAL_PATHS` in `sbml_autolayout/known.py`, if at least 5 of its members lie on one chain in the map. Currently that is glycolysis (D-glucose → … → pyruvate).
   - The chain that covers the most members wins. In the carbohydrate map there is no hexokinase, so the main line starts at glucose-6-phosphate; galactose and glycogen join as branches.
   - After its last member, the line continues along the longest chain (e.g. mitochondrial pyruvate, or into the TCA cycle at citrate synthase).
3. **Otherwise the longest chain** in the map.

For a known pathway or `--main-path`:
- **Direction:** reversible reactions along it are drawn in pathway order (D-glucose → pyruvate), even when the file writes them backwards. BRENDA does this, e.g. F6P → G6P or pyruvate → PEP. This applies only when at least 5 consecutive members are directly linked, so a pentose phosphate map containing glucose, G6P and F6P is not treated as glycolysis.
- **Straight line:** the line stays straight; crossing removal may only change its spacing. A longest-chain main line may still step aside to avoid a crossing.
- **Skip-over reactions:** reactions that skip a main-line metabolite (e.g. non-phosphorylating GAPDH, GAP → 3PG past 1,3-bisphosphoglycerate) are drawn beside the line, never over it.

## Main metabolites vs. cofactors

The source decides first:
- BRENDA `main_metabolites` notes;
- Escher `node_is_primary`;
- SBML Layout roles;
- CellDesigner main reactants.

Only where the source says nothing do built-in rules apply: currency metabolites by name/ChEBI, and hub demotion that never disconnects the map.

## Ids in the output (for comparisons and benchmarks)

The SBML and Escher outputs use the ids from the input SBML file:
- **SBML:** species and reactions have their original ids, and the glyphs point to them
  (`layout:species`, `layout:reaction`, `originOfText`).
- **Escher:** `bigg_id` is the original SBML id, both for metabolite nodes and for reactions and their
  metabolite lists. The internal glyph ids (`SG_M_…`, `RG_R_…`) are only used for `--previous`.
- **Species of one file are never merged**, even if they have the same name or annotation. For example,
  a protein and its active form keep two separate nodes. Species are only merged across files when
  several files are combined into one map, and then the id from the first file is used.
- **Identical reactions in one file** (same participants) are drawn once. Each one still gets its own
  reaction with its original id in the SBML (with its own reaction glyph on top of the drawn one) and in
  the Escher map.
- **Not drawn, so no glyph:** exchange/demand reactions (unless `--keep-exchange`), protons (unless
  `--show-protons`) and modifiers (unless `--modifiers`). With all three switches, every species and
  reaction of the input appears in the output. JSON or Escher inputs keep their own ids.

## Notation: default or SBGN

`--notation sbgn` (or `"notation": "sbgn"` in `run_batch.py`) draws the map in SBGN Process
Description style:
- **Species are shapes with the name inside**, chosen by the species' SBO term in the SBML file:
  - simple chemical (metabolites, the default when there is no SBO term): rounded ends;
  - macromolecule (SBO:0000245, 0000252, …): rounded rectangle;
  - complex (SBO:0000253, 0000296, 0000297): cut corners;
  - nucleic acid feature (genes, RNA): rounded bottom corners;
  - unspecified entity: ellipse;
  - phenotype: hexagon.
- **Long names wrap.** A box is at most `--box-max` px wide (default 110; a single long word may stick
  out by up to 12 %). Names break at spaces, then hyphens, onto at most `--box-lines` lines (default 2).
  What still doesn't fit is shortened with "…". The full name is kept as the hover text and in the SBML
  text glyph.
- **Cofactor copies** are small one-line boxes (at most 72 px) with the SBGN clone marker (shaded lower
  band).
- **Process nodes** are squares with two ports. Consumption arcs have no arrowhead, production arcs
  have a filled one, and reversible reactions have arrowheads on both sides.
- **Spacing:** the layout gives the boxes 30 % more room. Lines end at the box outline, and boxes never
  overlap each other or the reaction names.
- **SBML:** each species glyph's bounding box is the box. Shapes and colours are not stored, since that
  would need the SBML Render package.
- **Escher:** Escher has its own look, so the Escher map uses the same positions but Escher's usual
  circles.

You can switch the notation and still use `--previous`: positions are kept, and labels that were inside
the boxes are placed beside the nodes again.

## Modifiers (catalysts, transcription factors, inhibitors)

Off by default. Without `--modifiers` (or `"modifiers": True` in `run_batch.py`):
- modifiers are not drawn, but stay in the SBML output as `modifierSpeciesReference`;
- regulated synthesis/degradation reactions ("→ protein") are skipped like exchange reactions.

With `--modifiers`, SBML modifiers are drawn as SBGN arcs ending at the reaction node:

| Kind | SBO | Arc end |
|---|---|---|
| catalysis | 0000013 | open circle |
| stimulation | 0000459 | open arrow |
| inhibition | 0000020 | red bar |
| modulation | 0000019 / none | open diamond |

How they are placed:
- **Short arc or clone:** a modifier that is already in the map gets a line if that line is short and crosses nothing. Otherwise it gets a clone next to the reaction, never a long line across the map.
- **Enzyme-only species:** a species that only acts as a modifier (an enzyme) gets its own node next to the first reaction it acts on.
- **Regulated synthesis/degradation:** reactions like transcription ("→ protein") or degradation ("protein →") that have a regulator are kept, and the empty side is drawn as the SBGN source/sink symbol ∅.
- **Regulon fan:** if the protein made takes part in nothing else, it is placed next to its regulator. A transcription factor's target genes form a fan around it.

**Outputs:**
- **SBML:** modifiers are written as `modifierSpeciesReference` with their SBO term. In the Layout their roles are `activator`, `inhibitor` or `modifier`.
- **Escher:** Escher has no modifiers, so they are left out there.

## Biomass and other pseudo-reactions

A biomass objective function (SBO:0000629) is not a biochemical step. It drains a dozen precursors in fixed amounts. Drawn as a normal reaction it would have one node with lines across the whole map; in e_coli_core that caused 73 instead of 8 crossings.

A reaction counts as a pseudo-reaction if it has that SBO term, or "biomass" in its id or name, or at least 8 main partners. Instead of a node, each precursor gets a short dashed arrow to a small hexagon labelled "biomass" right next to it (SBGN phenotype); products get an arrow from it.
- **SBML:** the reaction stays unchanged in the model. In the SBML Layout every drain is its own glyph of that reaction.
- **Escher:** drains are left out.
- **Option:** `--pseudo hide` leaves the reaction out of the drawing, and `--pseudo draw` draws it as a normal reaction.

## Known cycles

The curated cycles are the citric acid cycle, the urea cycle, the methionine cycle and the Calvin cycle (RuBisCO at the top).
- **Recognition:** a loop is drawn as one of these only if it contains every required member. Matching uses ids, KEGG, ChEBI and names, and copies in several compartments may all be on the ring (e.g. mitochondrial and cytosolic ornithine).
- **Calvin cycle:** RuBP, 3-PG, 1,3-BPG, GAP and Ru5P are required. The regeneration steps between GAP and Ru5P vary between models and simply join the ring.
- **Adding cycles:** add your own in `sbml_autolayout/known.py` (`CANONICAL_CYCLES`, `CYCLE_TOP`).

## Incomplete reactions in the input

If a BRENDA `main_metabolites` note names a metabolite that is missing from the reaction, it is added on the other side, and the import report says so. Example: the polymer made by mannuronan synthase. Otherwise the map would fall apart at that reaction.

## Without uv

    pip install -e .        # then: sbml-autolayout ... / sbml-autolayout-batch ...
    # or, without installing:  python run_batch.py ...  /  python run.py ...

## Folder contents

- `pyproject.toml`: the uv/pip project, which defines the two commands.
- `run_batch.py`: batch run with the editable `SETTINGS` block.
- `run.py`: one map.
- `sbml_autolayout/`: the layout code.
  - `layout.py`: the algorithm.
  - `known.py`: the curated cycles (`CANONICAL_CYCLES`) and main pathways (`CANONICAL_PATHS`); extend them there.
  - `model.py`: the cofactor lists `CURRENCY_NAMES` / `CURRENCY_CHEBI`.
- `examples/`: your four test files.
- `output/examples/`: results for them.

The rules the layout follows, from most to least important, are described in `LAYOUT_RULES.md`.
