# Standard Bias Studies

Four complementary checks of the discrete-profiling (`RooMultiPdf` envelope)
background model. All work on any 1D or 2D datacard whose per-category
background is a `RooMultiPdf` with a `pdfindex_*` discrete nuisance.

1. **Truth-vs-fit function matrix (per-category split datacard)** — for every
   pair of (truth function used to generate toys, function used to fit them)
   in the envelope, and for a few injected `r` values, measure the pull on
   `r`. Reveals which specific function pairings bias the fit. Uses a
   hand-split single-category datacard, which drops any nuisance whose name
   references another category (a known approximation — see caveat below).
2. **Masked envelope bias check (recommended)** — the same idea, but on the
   **full combined workspace** with the other categories switched off via
   `--channel-masks` instead of a split datacard, so every correlated
   systematic stays intact. Supports both "truth vs. fixed fit function"
   (comparable to #1) and, more importantly, **"truth vs. the envelope left
   floating"** — i.e. toys are thrown from one fixed truth function but fit
   exactly like combine would fit real data (envelope profiled, not frozen).
   This second mode is the flashggFinalFit convention and the one that
   actually matches analysis behavior, since real fits never freeze the
   discrete index.
3. **Full-workspace envelope bias check (no masking)** — same question as #2
   ("what if category J's truth function is i, i != 0"), but with every
   *other* category left completely normal: no `--channel-masks`, generation
   only ever perturbs category J's pdfindex (every other category is
   generated at its default index, 0, implicitly — nothing else needs to be
   set), and the envelope for every OTHER category is left floating in the
   fit rather than frozen. Two fit variants, both with the full envelope
   active everywhere except possibly category J:
   - `--fit-mode envelope`: fit with the envelope floating in **every**
     category, including J — the true "fit exactly like real data, only the
     toy's truth in J was function i" check.
   - `--fit-mode wrongfit`: freeze **only** category J's pdfindex to some
     function k != i (looping over every other function in J's envelope),
     while every other category's envelope still floats — tests how bad the
     bias gets if category J's discrete profiling ends up stuck on the wrong
     function while the rest of the fit behaves normally.
   Run over all `i` (every function in J's envelope) x all `J` (every
   category) by looping `--cat` and letting the script's own `TRUTH_FUNCTIONS`
   loop cover every `i` per category.
4. **Fixed-vs-envelope significance / spurious signal** — compare the
   background-only significance distribution when the envelope index is
   frozen to its best fit vs. left floating (profiled). If floating the
   envelope systematically inflates significance (Type-1 error rate) above
   what freezing it or a Gaussian would give, that's the classic "spurious
   signal from discrete profiling" bias.

## Changes for resonant HH (X → HH)

These changes were made to run the checks on the **resonant X → HH → bbγγ**
datacards; they were first used on the M900_v4 resolved card. Each one is
tagged **[resonant HH]** where it appears below. All new options are off by
default, so without them every script behaves exactly as before. The M900_v4
study itself is written up in [`M900_v4_res_bias_checks.md`](M900_v4_res_bias_checks.md).

Why they were needed: at high X mass the resolved categories have very little
background under the m_γγ peak. With the original fixed `--rMin -5 --rMax 5`
many background-only toys fail or end at the lower edge of the fit range, and
some categories' s+b pdf goes negative at small negative r (−0.55 for
cat0Resolved at M900). The new options make the fit range, the treatment of
those toys and the summary number match the existing resonant bias study,
`flashggFinalFit/submitBiasTests.py`.

**Fixes for these cards**
- `resolve_pdfindex_legs.py`, `build_dropfunc_workspace.py`: the category
  token (e.g. `cat2Resolved`) was compared against lower-cased pdf names, so
  bin names with capital letters never matched. Now case-insensitive.
- `run_dropfunc_envelope_bias.sh`: accepts absolute paths on `shapes` lines.

**New options** (in `run_dropfunc_envelope_bias.sh` and `run_full_envelope_bias.sh`)
- `--r-expected X`: fit each injected r in `[r_inj − X, r_inj + X]` instead
  of `--rmin/--rmax`; default injections become `0,X`. X is meant to be the
  median expected limit, as in `submitBiasTests.py --expRange 1`. Outputs get
  an `_rInjPmX` suffix.
- `--fixed-pull P`: also reports `median_pull_fixed`, the median pull with
  failed, edge, lopsided and |pull| ≥ 5 toys counted at ∓P instead of dropped.
  Same rules as `submitBiasTests.py --failedFitR -P --keepEdgeFits P
  --keepLopsidedFits P --clipOutlierPulls`, and it reproduces that script's
  medians exactly. Plots go to a `…_fixedPullP` folder.
- `--r-min-floor auto|R`: never let rMin go below the floor (and never above
  r_inj), as `submitBiasTests.py --rMinFloor`. `auto` = 0.9 × the least
  negative r at which some category's s+b density reaches zero. Outputs get
  a `_rMinPos` (`_rMinFloorR`) suffix.
- `run_full_envelope_bias.sh` only: `--r-values`, `--parallel` /
  `--max-jobs N`, and `--fit-mode samefit` (the fit = truth diagonal that
  `wrongfit` skips).

**Plotting** (`plot_bias_pull_FitDiagnostics.py`)
- `--fixed-pull P` (with `--r-min/--r-max`) computes `median_pull_fixed`.
- The r = 0 core-fit window is now also reused for non-integer injected r
  (e.g. r_exp = 2.3672), and plot names keep the decimal (`_r2.3672.png`).

**New helper scripts**
- `zero_density_r.py <ws.root>`: the rMin floor used by `--r-min-floor auto`
  (port of `submitBiasTests.py`'s `zeroDensityR()`).
- `make_median_tables.py`: truth × fit tables of `median_pull_fixed` for the
  full-envelope (envelope / wrongfit / samefit) and dropfunc runs, as Markdown
  plus one PNG per category and r, labelled with each category's best-fit function.
- `organize_plots_by_r.py <run_dir>`: moves each `bias_plots_*` folder's
  plots into one `r<value>/` subfolder per injected r.

## Step-by-step (new datacard, 1D or 2D)

```bash
# 0. Bootstrap a run directory (splits the datacard by category automatically)
./setup_new_run.sh /path/to/Datacard_newAnalysis.txt myRunName
cd myRunName

# 1. Build the per-category multipdf workspaces and print the RooMultiPdf map
#    (function name <-> discrete index), needed by parser_legs_matrix.py below.
./generate_workspaces_first.sh
```

### 1. Truth-vs-fit function matrix

`2b_standard_bias_studies_2D_1D_tryFitDiagnostic.sh` loops over all
(truth-function, fit-function, r) combinations for a category, using
`parser_legs_matrix.py` to look up each function's `pdfindex` value in the
workspace built above. Without `--run`, it only *prints* the commands
(dry-run) so you can sanity-check them before actually launching combine.

```bash
# Dry run first (prints the combine commands it would execute):
./2b_standard_bias_studies_2D_1D_tryFitDiagnostic.sh

# Actually run it for one category:
./2b_standard_bias_studies_2D_1D_tryFitDiagnostic.sh --run cat2

# If the category has independent HGG/HBB legs (combined datacards) and you
# only want to scan the HGG leg to save time:
./2b_standard_bias_studies_2D_1D_tryFitDiagnostic.sh --run cat2 --only-hgg
```

Edit the `TRUTH_FUNCTIONS` / `FIT_FUNCTIONS` / `R_TRUTH_VALUES` arrays at the
top of the script to match the function families available in your envelope
and the `r` values you want to scan. Each combination runs
`combine -M GenerateOnly --saveToys`, then `combine -M FitDiagnostics
--toysFile ... --toysFrequentist --robustFit 1`, then plots the pull with
`plot_bias_pull_FitDiagnostics.py` (and the r-only variant
`plot_bias_pull_FitDiagnostics_onlyrfit.py`) into
`./bias_plots_FitDiagnostics_allCombinations/`, appending each result to
`pull_summary.json` in that folder.

Already ran the generation+fit and just want to redo the plots quickly
(e.g. after tweaking a plotting cut)? Use the plotting-only variant instead
of re-running combine:

```bash
./2b_standard_bias_studies_2D_1D_tryFitDiagnostic_onlyplottingoption.sh --run cat2
```

`2b_standard_bias_studies_2D_1D_tryFitDiagnostic_vs_envelope.sh` is the same
matrix but additionally compares against a frozen-envelope fit per
combination — an optional, more expensive cross-check.

To see, per category, how often the envelope chose a non-truth function and
what the residual pull looks like once those "wrong pick" toys are filtered
out:

```bash
python3 analyze_envelope_choice_bias.py \
    -i fitDiagnostics.fit_cat2_both_truth_exp1_fit_bern1_r0.root \
    -v pdfindex_hgg_cat2_13TeV --truth-idx 1 --r-truth 0
```

Once you have enough entries in `pull_summary.json`, build the summary
tables (median pull ± error per truth/fit pair, per category and `r`):

```bash
python3 create_table.py
```

### 2. Masked envelope bias check (recommended)

`run_envelope_bias.sh` runs on the **combined** workspace (not a split
per-category datacard) and uses `--channel-masks` to switch off every
category except the one under study, so all cross-category correlated
nuisances stay exactly as in the real analysis. It supports the same
truth-vs-fixed-fit matrix as script #1 *and* the truth-vs-floating-envelope
check, from one script:

```bash
# Build the combined, channel-masked workspace once:
text2workspace.py Datacard_newAnalysis.txt --channel-masks -m 125 -o Datacard_newAnalysis.root

# Recommended: throw toys from one fixed truth function, fit with the
# envelope left floating (matches what combine does on real data):
./run_envelope_bias.sh --ws Datacard_newAnalysis.root --cat cat2 --fit-mode envelope

# Diagnostic: truth vs. a fixed (possibly different) fit function, same as
# script #1's matrix but on the combined+masked workspace:
./run_envelope_bias.sh --ws Datacard_newAnalysis.root --cat cat2 --fit-mode fixed

# Dry-run by default (prints the combine commands); add --run to execute:
./run_envelope_bias.sh --ws Datacard_newAnalysis.root --cat cat2 --fit-mode envelope --run
```

Under the hood it:
- Discovers `mask_*` parameters in the workspace and sets `mask_cat2=0`
  (active) / all others `=1` (masked out).
- Uses `resolve_pdfindex_legs.py <ws> <cat_number>` to map each envelope
  function name (`exp1`, `bern1`, `bern2`, `pow1`, `lau1`, `lau2`) to its
  `pdfindex` value **for that specific category only** — unlike
  `parser_legs_matrix.py`, which assumes a single-category workspace and can
  mis-associate function indices across categories if reused directly on a
  combined workspace.
- Freezes every *other* category's `pdfindex` (masked channels don't
  contribute to the likelihood, so their discrete index is both irrelevant
  and, left floating, a source of minimizer instability) and, on the fit
  step, adds the same `--X-rtd MINIMIZER_*` stability flags used by
  `RunBiasInSignificance.py` below.
- Reuses `plot_bias_pull_FitDiagnostics.py` for plotting, writing to
  `./bias_plots_masked_<cat>_<fit-mode>/`.

`--ntoys`, `--seed`, `--rmin`/`--rmax` and `--only-hgg` (for categories with
both an HGG and HBB leg) are all overridable flags — see `--help`.

### 3. Full-workspace envelope bias check (no masking)

```bash
# Build the combined workspace once -- plain, NO --channel-masks needed:
text2workspace.py Datacard_newAnalysis.txt -m 125 -o Datacard_newAnalysis.root

# Fit with the full envelope floating everywhere (recommended: matches real data):
./run_full_envelope_bias.sh --ws Datacard_newAnalysis.root --cat cat2 --fit-mode envelope

# Freeze ONLY cat2 to a wrong function (k != truth i), envelope floats elsewhere:
./run_full_envelope_bias.sh --ws Datacard_newAnalysis.root --cat cat2 --fit-mode wrongfit

# [resonant HH] Freeze ONLY cat2 to the truth function itself (k == i), the
# diagonal wrongfit skips; envelope floats elsewhere:
./run_full_envelope_bias.sh --ws Datacard_newAnalysis.root --cat cat2 --fit-mode samefit

# Dry-run by default (prints the combine commands); add --run to execute:
./run_full_envelope_bias.sh --ws Datacard_newAnalysis.root --cat cat2 --fit-mode envelope --run
```

Under the hood, for every truth function `i` present in cat2's envelope
(discovered via `resolve_pdfindex_legs.py`, same as `run_envelope_bias.sh`):
- **Generation** (identical in both fit modes): `--setParameters` sets ONLY
  cat2's pdfindex to `i`; `--freezeParameters` freezes it there, plus every
  other category's pdfindex with no explicit value (so they stay at the
  workspace default, index 0). No `--channel-masks`, no `mask_*` parameters
  at all — every category's data goes into the toy exactly as it would for
  a real combined dataset.
- **Fit, `--fit-mode envelope`**: only `MH` is frozen. Every category's
  pdfindex (including cat2's) is left floating, so combine discrete-profiles
  the envelope everywhere, exactly like it would on real data.
- **Fit, `--fit-mode wrongfit`**: `MH` plus cat2's pdfindex are frozen, set
  to some function `k != i` (loops over every other function in cat2's
  envelope); every other category's pdfindex is still left floating.

Because the envelope floats in every *other* category in both modes here
(unlike `run_envelope_bias.sh`, which freezes them since they're masked
out), this is a much bigger discrete-combinatorial fit if you have several
multi-function categories — expect it to run slower. The `--X-rtd
MINIMIZER_*` stability flags are therefore ON by default (they were
deliberately left off in `run_envelope_bias.sh` only because they crash
`--channel-masks` + `FitDiagnostics` specifically, per that script's
comments — that combination doesn't apply here). Use
`--no-stability-flags` if you run into trouble.

Reuses `plot_bias_pull_FitDiagnostics.py` for plotting, writing to
`./bias_plots_full_<cat>_<fit-mode>/`. Same `--ntoys`, `--seed`,
`--rmin`/`--rmax`, `--only-hgg` flags as `run_envelope_bias.sh` — see
`--help`.

**[resonant HH]** `--r-values`, `--r-expected X` (fit window `r_inj -/+ X`,
default injections `0,X`), `--fixed-pull P` (extra `median_pull_fixed` in
`pull_summary.json` with failed/edge/lopsided/outlier toys counted at -/+P,
same rules as `submitBiasTests.py --failedFitR -P --keepEdgeFits P
--keepLopsidedFits P --clipOutlierPulls`) and `--r-min-floor auto|R` work as
in the dropfunc check. `--parallel` runs every combination's GEN+FIT as a
background job (plots afterward, sequentially); `--max-jobs N` caps how many
run at once. `--fit-mode samefit` uses the same generation and seed as the
other modes, so its toys are identical and the three modes fill one
truth × fit table (see `make_median_tables.py`).

### 4. Fixed-vs-envelope significance (spurious signal)

`RunBiasInSignificance.py` runs the classic combine discrete-profiling bias
recipe against the **full/combined workspace** (`text2workspace.py` output
for the whole datacard, not a single category):

```bash
text2workspace.py Datacard_newAnalysis.txt -o Datacard.root -m 125

python3 RunBiasInSignificance.py --inputWSFile Datacard.root --mode setup \
    --initial-fit-param <low_impact_nuisance_name>
python3 RunBiasInSignificance.py --inputWSFile Datacard.root --mode generate --nToys 2000
python3 RunBiasInSignificance.py --inputWSFile Datacard.root --mode fixed    --nToys 2000
python3 RunBiasInSignificance.py --inputWSFile Datacard.root --mode envelope --nToys 2000

python3 SummaryBiasSignificance.py
```

`--initial-fit-param` must name any single low-impact nuisance already in the
datacard (e.g. a luminosity term) — combine needs at least one floating
parameter to seed the initial fit. `SummaryBiasSignificance.py` writes
`plots/summary_significance.{png,pdf}`: the Type-1 error rate vs. Z-score for
`fixed` vs `envelope` vs the Gaussian expectation. Envelope tracking above the
`fixed`/Gaussian curves by more than ~10-20% indicates the discrete profiling
is inflating significance (i.e. injecting bias / look-elsewhere-like effect).

### 5. Follow-up checks used for the v7 final results

These scripts were added for the v7 final datacard
(`v7_SnTV7F1E500T0821thTest_resample_DDonly_min`) and assume its category
layout (cat0 boosted with a single fixed function, cat1-4 with a
`RooMultiPdf` envelope). Run them from the bootstrapped run directory.

- **cat0 fixed model.** `run_cat0_fixedmodel_bias.sh` /
  `run_cat0_fixedmodel_bias_default0.sh` test truth = exp1/pow1;
  `run_cat0_fixedmodel_bias_extra.sh` / `run_cat0_fixedmodel_bias_default0_extra.sh`
  add the missing truth = bern1/lau1 and append to the same
  `pull_summary.json`. Other categories are frozen to their best-fit index
  (`bestfit`) or to index 0 (`default0`). Summary plots:
  `plot_mean_sigma_vs_pdfindex_cat0_bestfit.py`,
  `plot_mean_sigma_vs_pdfindex_cat0_default0.py`.
- **dropfunc summaries.** `plot_mean_sigma_vs_pdfindex_dropfunc.py` and
  `plot_mean_sigma_vs_pdfindex_dropfunc_bestfit.py` now include cat0 (read
  from the cat0 fixed-model `pull_summary.json`) as the first panel when it
  exists. The `_bestfit` variant reads the published `pull_summary.json`
  files under `~/public_html`.
- **cat1 checks.** `run_cat1_standalone_check.sh` (cat1 workspace alone),
  `run_cat1_altws_gen_nofreezeothers.sh` + `run_cat1_altws_nofreezeothers.sh`
  (alternative workspace, other categories not frozen; set `ALTWS=` to point
  to it) and the pull summaries `make_pull_summary_cat1_envelope.py` /
  `make_pull_summary_cat1_reproducibility_check.py` (pull =
  (r - r_inj) / ((|rLoErr| + |rHiErr|) / 2), fixed Gaussian range [-2, 2]).

`plot_bias_pull_FitDiagnostics.py` draws the fitted Gaussian on a clone over
the full window (the previous version clipped the curve at the core-fit
range) and, for r_truth != 0, reuses the core-fit window of the r_truth = 0
plot of the same truth/fit combination (cached in
`<outdir>/.core_range_cache.json`), so plots at different injected r are
directly comparable. Run the r = 0 case first.

### 6. Final summary (dropfunc + cat0 fixed model)

This is how the v7 final summary
(`~/public_html/BiasSensitivityStudies_v7_final_summary_rMin-5rMax5_1000toys/StandardStudies_summary/`)
was produced. Everything runs from the bootstrapped run directory, on the
combined workspace **without channel masks**: one category J is perturbed at a
time, the fit uses all categories, so `r` is the combined signal strength.
How the dropfunc check works step by step: section 7.

```bash
# 0. Combined workspace + best-fit pdfindex of every category on data
#    (writes higgsCombine_initial.MultiDimFit.mH125.root, then cached by
#    get_bestfit_pdfindices.py into bestfit_pdfindices.sh on first use).
text2workspace.py Datacard.txt -m 125 -o Datacard.root
python3 RunBiasInSignificance.py --inputWSFile Datacard.root --mode setup \
    --initial-fit-param <low_impact_nuisance_name>

# 1. dropfunc, cat1-4: truth function i generated in cat J, then REMOVED from
#    the RooMultiPdf of cat J before the fit (truth not in the envelope).
#    Other categories: set + frozen to their data best-fit index.
#    Defaults: 1000 toys, r_truth = 0,2, rMin/rMax = -5/5. Dry run without --run.
for c in cat1 cat2 cat3 cat4; do
    ./run_dropfunc_envelope_bias.sh --ws Datacard.root --datacard Datacard.txt \
        --cat $c --run --parallel --publish
done
#    -> bias_plots_dropfunc_<cat>/pull_summary.json + pull_*.png

# 2. cat0 (single fixed function, nothing to drop): toys from its standard
#    4-function envelope, fit with cat0 fixed (see section 5).
./run_cat0_fixedmodel_bias.sh          && ./run_cat0_fixedmodel_bias_extra.sh
./run_cat0_fixedmodel_bias_default0.sh && ./run_cat0_fixedmodel_bias_default0_extra.sh

# 3. Summary plots: pull mean/sigma vs truth pdfindex, cat0 first, then cat1-4.
python3 plot_mean_sigma_vs_pdfindex_dropfunc_bestfit.py   # reads the published dropfunc/cat0 pull_summary.json
python3 plot_mean_sigma_vs_pdfindex_dropfunc.py           # reads local bias_plots_dropfunc_cat*/ (default0 runs)
python3 plot_mean_sigma_vs_pdfindex_cat0_bestfit.py
python3 plot_mean_sigma_vs_pdfindex_cat0_default0.py
```

Notes:
- The `default0` dropfunc results (other categories at pdfindex 0, published
  in `..._dropfunc_cat*_default0_...`) come from an earlier pass of
  `run_dropfunc_envelope_bias.sh`; the current script always freezes the
  other categories to their best-fit index.
- The AN table `bias_tests_dropfunc_cat0.tex` in the summary page was written
  by hand from the `pull_summary.json` files; no script in this repo produces
  it (`create_table.py` makes the truth-vs-fit matrix table of section 1).
- Still open: `run_cat0_wrongfit_bias*.sh` (truth = cat0 fixed function, fit
  forced on a wrong one) was never completed or published.

### 7. The dropfunc check in detail

**Question it answers.** With discrete profiling, the fit can only choose
among the functions in the envelope. If the true background shape is *not*
one of them, how much is `r` biased? The dropfunc check makes this happen
on purpose: toys are generated from function `i` in category J, and the fit
is done on a workspace where `i` has been **removed** from category J's
`RooMultiPdf`, so the fit can never pick the truth.

Compare with the other checks of this README:

| Check | Truth in the fit envelope? | Fit pdfindex of cat J |
|---|---|---|
| envelope (sections 2, 3) | yes | floating |
| wrongfit (section 3) | yes, but not reachable | frozen to a wrong function k != i |
| **dropfunc** | **no** | floating among the remaining functions |

wrongfit only tests a fit stuck on the wrong function; dropfunc tests an
envelope that is an incomplete set of shapes, which is the realistic worry
on data.

**What `run_dropfunc_envelope_bias.sh --cat catJ` does**, for each truth
function `i` in {exp1, bern1, bern2, bern3, pow1, lau1, lau2} that cat J
actually has (others are skipped; `bern1` also matches a `flat0` pdf):

1. *Find cat J's envelope.* Reads the datacard line
   `shapes bkg_mass catJ <file> multipdf:<name>` to get the raw per-category
   multipdf file and the pdf name (categories may point to different
   directories, e.g. `background` vs `background_boosted`, so nothing is
   hard-coded). `resolve_pdfindex_legs.py` maps function names to pdfindex
   values in the combined workspace.
2. *Remove function i.* `build_dropfunc_workspace.py` rebuilds that
   per-category file with `i` taken out of the `RooMultiPdf`, keeping every
   object name (workspace `multipdf`, pdf, `pdfindex`, `_norm`) unchanged.
   Output: `dropfunc_workspaces/<file>_drop_<i>.root`.
   `visualize_envelope_drop.py` plots the envelope before/after the drop
   (fit to data) if you want to look at it.
3. *Patch the datacard.* `patch_datacard_dropfunc.py` copies the datacard and
   changes only the file path on cat J's `bkg_mass` shapes line (its
   `data_obs` line and all other categories are untouched):
   `<datacard>_catJ_drop_<i>.txt`, then `text2workspace.py` ->
   `<datacard>_catJ_drop_<i>.root`. Steps 2-3 are cached: existing files
   are reused.
4. *Generate* on the **original** combined workspace (`--ws`), no channel
   masks: cat J pdfindex set + frozen to `i`; every other category set +
   frozen to its best-fit index on data (from `bestfit_pdfindices.sh`, see
   section 6 step 0). `GenerateOnly`, `--expectSignal r_truth`,
   `--saveToys`.
5. *Fit* the same toys on the **reduced** workspace with `FitDiagnostics`,
   only `MH` frozen: every pdfindex floats, including cat J among its
   remaining functions. `--robustFit 1 --cminDefaultMinimizerStrategy 0`.
6. *Pull.* `plot_bias_pull_FitDiagnostics.py` with
   `--fit-func envelope_minus_<i>`:
   - toys kept if `fit_status` is 0 or 300 and both `rLoErr`, `rHiErr` > 0.05;
   - pull = (r - r_truth) / sigma, with sigma = `rLoErr` if r > r_truth,
     `rHiErr` otherwise, or their average if the two differ by more than a
     factor 3;
   - Gaussian fit in an iterative core window (mean +- 2 sigma; for
     r_truth != 0 the window of the r_truth = 0 plot is reused, so run r = 0
     first). The fit is drawn over the full [-5, 5] range.
   - Each (i, r_truth) appends one entry to
     `bias_plots_dropfunc_catJ/pull_summary.json` (`mean_pull`,
     `median_pull`, `rms_pull`, `fit_mean`, `fit_sigma`, `n_toys`).

Defaults: 1000 toys, seed 123456, r_truth = 0 and 2, `--rMin -5 --rMax 5`.
Useful options: dry run without `--run` (still builds the reduced
workspaces so you can inspect them), `--parallel` (all GEN+FIT jobs at
once; plotting stays sequential because all entries go to the same JSON),
`--truth-funcs` / `--r-values` (fill in missing combinations),
**[resonant HH]** `--r-expected X` (fit each injection in the window
`[r_inj - X, r_inj + X]` instead of `--rmin/--rmax`, default injections become
`0,X`; outputs get an `_rInjPmX` suffix so they never mix with fixed-range runs),
**[resonant HH]** `--fixed-pull P` (extra `median_pull_fixed` with bad toys at
-/+P, see section 3), **[resonant HH]** `--r-min-floor auto|R` (rMin never below
the floor where the s+b pdf stays positive, `_rMinPos` suffix),
`--publish` (copies plots + JSON to
`~/public_html/BiasSensitivityStudies_v7_final_dropfunc_<cat>_rMin<X>rMax<Y>_<N>toys/StandardStudies_dropfunc/`).
`--stability-flags` is off by default because
`MINIMIZER_freezeDisassociatedParams` makes FitDiagnostics segfault on this
analysis. `--leg hbb` (2D categories) has never been tested on a real
two-leg workspace.

**Reading the result.** `fit_mean` of the pull is the bias in units of the
statistical uncertainty on `r`; `fit_sigma` close to 1 means the
uncertainty is correctly estimated. The summary plots (section 6, step 3)
show `fit_mean` +- `fit_sigma` vs the pdfindex of the removed function,
per category, for r_truth = 0 and 2. In the v7 final AN table the bias is
quoted as max and mean |`fit_mean`| per category in %.

**cat0** has a single fixed function, so there is nothing to drop: it is
covered by the fixed-model check (section 5), and shown as the first panel
of the dropfunc summary plots.

Note: the dropfunc runner does not pass `--r-min/--r-max` to the plotting
script. This has no effect: toys with `r` stuck at rMin/rMax (~2%) always
have a near-zero error on that side and are already removed by the
`rLoErr, rHiErr > 0.05` cut (checked on the v7 final outputs: identical
`fit_mean`, `fit_sigma` and toy counts with and without the boundary cut).
