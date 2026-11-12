# sbml-autolayout – layout rules

These are the rules the layout algorithm (`sbml_autolayout/layout.py`, version sbml-autolayout-0.10) uses, ordered
from most to least important. When two rules conflict, the higher one wins. Each rule names the
option or code section where it can be changed.

---

## A. Hard rules

1. **Reproducible.** The same input with the same settings always gives the same picture. There is no
   randomness; even the annealing step uses a fixed seed. Every result carries `pathway_hash`,
   `params_hash` and `layout_hash` in the SBML layout annotation.
2. **Previous layouts and pins are respected.** With `--previous old.sbml`, every glyph that already
   existed (metabolites, reactions, clone sides, labels) stays exactly where it was; only new
   elements are placed. Pins given in the input file override even the previous layout.
3. **No node touches another node.** If two nodes overlap, the more movable one moves, in this order:
   cofactor copies, then copies of main metabolites, then reactions and branch metabolites. The main
   line, circles and pinned nodes move only as a last resort (`resolve_node_overlaps`).
4. **No label overlaps a node or another label.** Each label is tried in rings around its node: first
   N/E/S/W/…, then further out. In the nearest ring that has free spots, the spot crossing the
   fewest lines wins. A label that has to move away from its node gets a leader line
   (`place_labels`).
5. **Nothing is left unconnected that belongs together.** Whether a compound is main or side never
   splits the network: a compound stays main if demoting it would cut the map apart (see B2).

## B. What gets drawn as the structure (main vs. side compounds)

1. **Main vs. side compounds.** In each reaction, cofactors and currency metabolites (ATP, NAD(P)H,
   H2O, CO2, phosphate, …, found by id, name or ChEBI) are *side compounds*. They are drawn as
   small copies on short hooks next to the reaction. Everything else is *main* and forms the
   network. If the source file says which compounds are main (a `main_metabolites` note), that is
   used as given.
2. **Hub demotion.** If one side of a reaction has several main candidates, one that takes part in
   ≥ `--hub-degree` (4) reactions and in ≥ 2× as many as the others becomes a side compound *in that
   reaction*. Example: acetyl-CoA released by every β-oxidation step. At least one main compound
   per side is kept. A demotion is undone if it would disconnect the map.
3. **Pseudo-reactions (biomass).** Reactions with SBO:0000629–0000632, "biomass" in the id or name, or
   ≥ 8 main partners are not drawn as one node with lines across the map. Each precursor instead
   gets a short arrow to a small hexagon next to it (`--pseudo drain`).
4. **Exchange/proton reactions** are hidden unless `--keep-exchange` / `--show-protons` is given.
5. **Modifiers** (catalysts, inhibitors, stimulators) are off by default. With `--modifiers` they are
   drawn in SBGN style (see G2).

## C. Main line and cycles (the skeleton)

1. **Known cycles are drawn as circles in their textbook orientation.** These are TCA (citrate
   synthase at 12 o'clock), urea cycle, methionine cycle and Calvin cycle (RuBisCO at the top). A
   metabolite that occurs twice in a cycle (e.g. in two compartments) is allowed. `--cycles known`
   (default), `all` (every loop) or `none`.
2. **One main line (backbone).** The longest chain of the map, or a canonical path such as
   glycolysis (D-glucose → pyruvate), or your own `--main-path A,B,C`. It is drawn as one straight
   line in the flow direction (`--orientation TB` top→bottom, or `LR`). Between two backbone nodes
   the shortest route is preferred over detours: aldolase F16BP → GAP, not via DHAP.
3. **Reversible reactions follow the overall flow.** Their drawing direction is chosen so that the
   whole chain points the same way, including BRENDA files with mixed directions.
4. **Long main lines fold** into serpentine columns: above `--fold` metabolites, or automatically
   above 12. The U-turn gap is only as wide as labels and cofactors need.
5. **A cycle on the main line** is entered and left at the right places: the backbone continues below
   it. A loop that enters and leaves through the same metabolite hangs beside the line.
6. **Inside circles:**
   - Parallel reactions for the same step (isoenzymes, NAD/NADP variants) form a small stack
     inside the circle.
   - A shortcut between two non-neighbouring members (e.g. isocitrate lyase, the glyoxylate
     shunt) is drawn inside the circle.
   - An extra main partner of the top reaction (acetyl-CoA) sits directly outside it.

## D. Readability rules (weighted in the optimisation)

These are soft rules. The crossing search minimises

`10 × crossings + 10 × lines through a node + 3 × misaligned partners + 3 × length (per spacing)`,

and a bent main line costs more than a bent branch.

1. **No reaction node on another line.** A reaction square must never sit on a line that belongs to
   another reaction; it would look like a 2:2 reaction. This counts double in the final check.
2. **No line through a metabolite.** A bypass that skips a main-line metabolite (GAPDH variant
   GAP → 3PG past 1,3-BPG) is drawn beside the line, never over the metabolite. Its side is chosen
   by fewest problems.
3. **No unnecessary crossings.** Lines lying on top of each other count as crossings. The search
   uses these moves:
   - single grid steps;
   - moving or mirroring a whole branch;
   - reflecting a branch to the opposite side of its anchor (e.g. (R)-lactoyl-CoA left of
     acryloyl-CoA);
   - letting a leaf jump around its partner;
   - swaps (never involving the main line).

   Then simulated annealing (fixed seed) runs where crossings remain.
4. **Straight metabolite → reaction → metabolite.** A reaction with one main substrate and one main
   product sits exactly midway between them on a straight line. Metabolites are moved to allow
   this; only if that fails does the reaction move to the best nearby grid point.
5. **Alignment.** Metabolites with nearly the same x (or y), within 0.3 × spacing, are put on one
   line. Alignment matters more than equal edge lengths. Main line, circles and pins are the
   anchors.
6. **Short links.** Every extra step of length costs. Other partners that are already drawn pull a
   node towards them, so loops close compactly. A partner that ended up far away is shown as a copy
   next to the reaction instead of a line across the map, but only if the network stays visibly
   connected.
7. **Branches go to one side, cofactors to the other.** `--branch-side right` puts branches right of
   the main line; side compounds hook out on the opposite side. Cofactor hooks should not cross a
   main line.

## E. Symmetry patterns

1. **Parallel reactions** between the same two metabolites are placed symmetrically around the
   straight line (diamond shape).
2. **Y-shaped forks.** With two partners on one side and one on the other, the single partner moves
   to the middle of the row if it belongs to nothing else, so the reaction sits centred underneath.
   A reaction with ≥ 3 main partners sits centred between the substrate row and the product row.
3. **Regulons (fans).** A regulator with several targets that do nothing else (transcription factor →
   genes) puts the targets evenly on an arc on its free side, each reaction on its ray.
4. **Redox ladder.** A chain of 2:2 reactions handing a redox pair on (photosynthetic electron
   transport: PSII → plastoquinone → cyt b6f → plastocyanin → PSI → ferredoxin → FNR) is drawn as a
   zig-zag:
   - reactions in one row;
   - each pair in the gap between two reactions, one form above and one below;
   - water starts the chain;
   - drawn compactly (0.6 × spacing between columns, 0.65 × spacing above and below), independent of
     automatic enlarging.

## F. Compactness

1. **Grid.** Metabolites snap to a grid (`--grid`, 0.5 × spacing; reactions use half a cell). Circles
   and pins stay where they are.
2. **Compaction.** Empty rows/columns left over are removed in whole grid cells. Every link keeps at
   least the base spacing (1.5× for ≥ 3 parallel reactions and for bypasses). Circles and fans move as
   rigid blocks. No crossing, line through a node or node contact is added.
3. **Automatic enlarging.** If more than 8 % of a part's labels had to move away from their node
   (leader lines), that part gets 25 % more room (at most twice).
4. **Density.** `--density` scales all distances (> 1 = tighter); node and font sizes stay the same.

## G. Placement of separate parts

1. **Packing of unconnected components.** The largest part goes first. Parts sharing a metabolite
   with what is already drawn go next to that metabolite. All others are packed largest first into
   the free spot beside, above or below (flush with the edges of) an already placed part, including
   notches between two parts. The spot chosen keeps the whole map smallest, with a mild preference
   for a landscape-ish shape (about 1.4 : 1). Parts never overlap, and shifts are whole grid cells.
2. **Modifiers (if `--modifiers`).** A nearby modifier gets a line with an SBGN arrowhead
   (catalysis ○, stimulation, inhibition ⊣, modulation ◇). A far one gets a copy next to the
   reaction. An enzyme acting on nothing else gets its own small node. Spots are chosen where the
   new line crosses nothing. A regulator of a reaction with an empty side stands in for that side.
   Otherwise the empty side is drawn as the SBGN source/sink symbol (∅).

## H. Output

The SBML L3 file uses the Layout package: a copy is an extra speciesGlyph of the same species, one
reaction glyph is written per drain, and modifiers are written with SBO terms. There is also an
SVG, and an Escher JSON that can be loaded with Map → Load map JSON. Labels in the Escher map follow
Escher's own conventions and avoid collisions.
