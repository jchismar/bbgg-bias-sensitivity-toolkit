# M900_v4 resolved: toolkit bias checks vs `submitBiasTests.py`

Launched 2026-11-19 12:14 PST on uaf-4 via `run_M900_v4_res_all.sh`.

- **Card:** `Combine/Datacard_2016-2025_cats_M900_v4/Datacard_2016-2025_allcats.txt`. Copied to `M900_v4_dropfunc/Datacard.txt` with absolute shape paths.
- **Signal:** `xhh` (run2 + run3).
- **Bins:** `cat0Resolved` … `cat4Resolved`.

The previous study being compared against is
`flashggFinalFit/biasTests/M900_v4_res_percatInComb_1000toys_rngExp1`:

```
submitBiasTests.py --perCategoryInComb --nToys 1000 --nJobs 10 --expRange 1 --rinj 0,1
                   --failedFitR -5 --keepEdgeFits -5 --keepLopsidedFits 5 --clipOutlierPulls
```

## Envelopes in this card

| Category | Functions (pdfindex order) | Best fit on data |
|---|---|---|
| cat0Resolved | bern1 | bern1 |
| cat1Resolved | bern1 | bern1 |
| cat2Resolved | bern1, exp1, pow1, lau1 | bern1 |
| cat3Resolved | bern1, exp1, pow1, lau1 | bern1 |
| cat4Resolved | bern1, exp1, pow1, lau1 | pow1 |

The best fit comes from `RunBiasInSignificance.py --mode setup --MH 125 --initial-fit-param lumi_13TeV_corr_151617`.

These best-fit indices are identical to the workspace defaults that `submitBiasTests.py` generates the
other categories with (0, 0, 0, 0, 2). So for this card, the other categories are generated the
same way in both setups.

## What is running

All runs use 1000 toys and inject r = 0 and r = r_exp = 2.3672. That is the median expected limit on this card, the same number
your run used. Each fit uses the window r_inj ± r_exp:

| r_inj | Fit range |
|---|---|
| 0 | [−2.3672, 2.3672] |
| 2.3672 | [0, 4.7344] |

| Test | README | Categories | Combinations | Output directory |
|---|---|---|---|---|
| dropfunc | §7 | cat2–4Resolved | 4 truth × 2 r = 8 per category | `M900_v4_dropfunc/bias_plots_dropfunc_<cat>_rInjPm2.3672_fixedPull5/` |
| full envelope, `envelope` mode | §3 | all 5 | cat0/1: 2, cat2–4: 8 | `M900_v4_fullenv/bias_plots_full_<cat>_envelope_rInjPm2.3672_fixedPull5/` |
| full envelope, `wrongfit` mode | §3 | cat2–4Resolved | 4 truth × 3 wrong fits × 2 r = 24 per category | `M900_v4_fullenv/bias_plots_full_<cat>_wrongfit_rInjPm2.3672_fixedPull5/` |

There's no dropfunc or wrongfit for cat0/cat1Resolved: with a single function there is nothing to remove or swap.

Each output directory has one pull plot per combination and a `pull_summary.json`. Progress logs:

- Overall: `M900_v4_res_all.log`, one line per finished category/test
- Per category: `M900_v4_dropfunc/dropfunc_<cat>.log`, `M900_v4_fullenv/fullenv_<mode>_<cat>.log`
- Per combine job: `combine_gen_*.log` / `combine_fit_*.log` next to them

### What each test asks

All three tests use the full combined workspace with no channel masks. Only category J's
background is perturbed; every other category is generated at its best-fit function and fitted
with its envelope floating.

- **dropfunc.** Generate cat J from function *i*, then fit on a workspace where *i* has been
  **removed** from cat J's envelope. The fit can never pick the truth. This tests an
  incomplete envelope. There is no counterpart in `submitBiasTests.py`.
- **full envelope, `envelope` mode.** Generate cat J from *i*, fit with the **full** envelope
  floating (truth available). This corresponds to your **gen-*i* / fit-envelope** entries.
- **full envelope, `wrongfit` mode.** Generate cat J from *i*, fit with cat J's pdfindex
  **frozen** to a function *k ≠ i*. This corresponds to your **gen-*i* / fit-*k*** off-diagonal entries.

## How these differ from the `submitBiasTests.py` runs

| | `submitBiasTests.py --perCategoryInComb` | Toolkit runs (here) |
|---|---|---|
| Grid | Every gen model (envelope + each function) × every fit model, including the diagonal | `envelope`: gen *i* × fit envelope. `wrongfit`: gen *i* × fit *k ≠ i*. dropfunc: gen *i* × fit envelope − *i*. **No gen-envelope row, no diagonal (fit = gen).** |
| Truth in the fit envelope | Yes | `envelope`: yes. `wrongfit`: yes, but not reachable. dropfunc: **no** |
| Other categories, generation | Workspace default pdfindex | Best fit on data. **Same indices for this card.** |
| Toys | `--toysFrequentist --bypassFrequentistFit` (generation and fit) | **Combine default (not frequentist)**: nuisances randomized from their priors |
| Fit options | `FitDiagnostics --robustFit 1 --cminDefaultMinimizerStrategy 0 --skipBOnlyFit` + `--X-rtd MINIMIZER_multiMin_*` flags | `FitDiagnostics --robustFit 1 --cminDefaultMinimizerStrategy 0`, background-only fit **not** skipped, no `X-rtd` flags |
| Range used to generate toys | r_inj − 4σ to r_inj + 6σ | Same as the fit window, r_inj ± r_exp. Only the POI bounds differ. |
| Injected r / fit window | 0 and 1 × r_exp, r_inj ± 1 × r_exp | Same |
| Toy splitting / seeds | 10 jobs, seeds 123456 … 123465 | 1 job per combination, seed 123456. For a given truth function, the three tests use the same generation command and seed, so they should be throwing the same toys. |
| Error used in the pull | Side facing r_in, symmetric error if missing or > 3× lopsided | Same |
| Bad toys | Failed = −5, at rMin/rMax with no error = ∓5, lopsided clipped at the edge = ∓5, \|pull\| ≥ 5 → ±5. Median only. | Two numbers, both in `pull_summary.json`. **`median_pull_fixed`: same rules, same classification.** `fit_mean` / `fit_sigma` / `median_pull`: bad toys dropped (`fit_status` 0 or 300 and both errors > 0.05). |
| Headline number | Median pull (flag at \|median\| ≥ 0.2) | Compare with your median: `median_pull_fixed`. The toolkit's own: `fit_mean` (Gaussian core fit) |

`median_pull_fixed` was checked against your existing fits. Re-scoring
`cat0Resolved` gen-envelope / fit-envelope gave −0.116 ± 0.028 (N = 999) at r_inj = 0
and −0.065 ± 0.029 (N = 1000) at r_inj = r_exp. Both are identical to your `biasSummary` JSON,
including every count.

### Comparing the numbers

- For a given category and truth *i*, compare:
  - full-envelope `envelope` with your gen-*i* / fit-envelope,
  - `wrongfit` with your gen-*i* / fit-*k*,
  - each using `median_pull_fixed` against your `median`.
- What can still make them differ: frequentist vs. default toys, the skipped background-only fit, the `X-rtd`
  flags, and different toys (seeds and job splitting).
- For cat0/cat1Resolved, the toolkit's gen-bern1 / fit-envelope is the same model as your
  gen-envelope / fit-envelope (one function), so those should agree within toy statistics.
- dropfunc has no counterpart in your study. To read it, compare it with the full-envelope
  `envelope` result for the same truth *i*. The difference is the effect of the truth missing from the envelope.
- Don't read `fit_mean` as the bias at this mass. Your cat0Resolved run had 280 of 1000
  fits at rMin for r_inj = 0, and dropping the low toys pushes the Gaussian mean up: +0.34, against a
  median of −0.12. Check `n_toys` (kept) against `n_in_median_fixed` and the
  `n_failed` / `n_edge_noerr` / `n_at_rmin` counts.
- For r_inj = r_exp the window starts at r = 0, so toys that fluctuate low pile up at rMin. That
  is the same window as yours, so the comparison is like for like, but expect many edge toys.

## Changes made to the toolkit

All in `BiasStudiesToolkit/StandardStudies/`; the full diff is in `git diff`.

- `resolve_pdfindex_legs.py`, `build_dropfunc_workspace.py`: the category token (`cat2Resolved`)
  was compared against lower-cased pdf names, so bin names with capitals never matched. Now case-insensitive.
- `run_dropfunc_envelope_bias.sh`: accepts absolute paths on `shapes` lines.
- `run_dropfunc_envelope_bias.sh`, `run_full_envelope_bias.sh`: new `--r-expected X` (window
  r_inj ± X, default injections 0 and X, `_rInjPmX` output suffix) and `--fixed-pull P`
  (`median_pull_fixed`, `_fixedPullP` plot directory).
- `run_full_envelope_bias.sh`: new `--r-values`, `--parallel`, and `--max-jobs N`. Plots still run one at a time
  after the fits, because they all append to one `pull_summary.json`.
- `plot_bias_pull_FitDiagnostics.py`: `--fixed-pull P` (needs `--r-min/--r-max`), and the r = 0
  core-window reuse now works for non-integer r and suffixed names.

## Rerunning or extending

```bash
cd /home/users/jachisma/HH/finalfits/bbgg-bias-sensitivity-toolkit/BiasStudiesToolkit/StandardStudies
nohup ./run_M900_v4_res_all.sh > M900_v4_res_all.log 2>&1 &
```

The runners regenerate toys on every `--run`; only the dropfunc workspaces are cached. To rerun only
part of a test:
- dropfunc: `--truth-funcs` / `--r-values`.
- full envelope: `--r-values`.

Before rerunning a combination, remove its entries from `pull_summary.json`. New entries are appended,
not replaced.
