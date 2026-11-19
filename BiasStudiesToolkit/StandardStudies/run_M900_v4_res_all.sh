#!/bin/bash
# M900_v4 resolved, all categories: dropfunc (README section 7) and
# full-envelope (section 3, envelope + wrongfit modes) bias checks.
# 1000 toys, r_inj = 0 and r_exp, fit window r_inj -/+ r_exp, bad toys also
# summarised at fixed pull -/+5 (median_pull_fixed), like submitBiasTests.py.
#
# r_exp = 2.3672: median expected limit (AsymptoticLimits -t -1) on this card,
# from biasTests/M900_v4_res_percatInComb_1000toys_rngExp1.
#
# Prerequisites already in M900_v4_dropfunc/: Datacard.txt (absolute shape
# paths), Datacard.root, higgsCombine_initial.MultiDimFit.mH125.root and
# bestfit_pdfindices.sh (RunBiasInSignificance.py --mode setup --MH 125).
#
# Usage: nohup ./run_M900_v4_res_all.sh > M900_v4_res_all.log 2>&1 &

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DROP_DIR="$HERE/M900_v4_dropfunc"
FULL_DIR="$HERE/M900_v4_fullenv"
mkdir -p "$FULL_DIR"

source /home/users/jachisma/HH/finalfits/setup_ff.sh /home/users/jachisma/HH/finalfits/CMSSW_14_1_0_pre4 >/dev/null 2>&1 \
    || { echo "[ERROR] CMSSW setup failed"; exit 1; }

COMMON="--r-expected 2.3672 --fixed-pull 5 --ntoys 1000 --run --parallel"
BEST="--bestfit-ws $DROP_DIR/higgsCombine_initial.MultiDimFit.mH125.root --bestfit-cache $DROP_DIR/bestfit_pdfindices.sh"

echo "=== $(date) launching M900_v4 resolved bias checks ==="

# dropfunc: cat2-4Resolved (cat0/cat1Resolved have a single function, nothing to drop).
# 8 combinations per category, all at once.
for c in cat2Resolved cat3Resolved cat4Resolved; do
    ( cd "$DROP_DIR" && "$HERE/run_dropfunc_envelope_bias.sh" --ws Datacard.root --datacard Datacard.txt \
          --cat $c $COMMON $BEST > "dropfunc_$c.log" 2>&1
      echo "=== $(date) dropfunc $c finished (rc=$?) ===" ) &
done

# full-envelope, envelope mode: all 5 categories (cat0/cat1: 2 combinations, cat2-4: 8).
for c in cat0Resolved cat1Resolved cat2Resolved cat3Resolved cat4Resolved; do
    ( cd "$FULL_DIR" && "$HERE/run_full_envelope_bias.sh" --ws "$DROP_DIR/Datacard.root" --cat $c \
          --fit-mode envelope $COMMON --max-jobs 8 $BEST > "fullenv_envelope_$c.log" 2>&1
      echo "=== $(date) full-envelope envelope $c finished (rc=$?) ===" ) &
done

# full-envelope, wrongfit mode: cat2-4Resolved (24 combinations each, 12 at once).
for c in cat2Resolved cat3Resolved cat4Resolved; do
    ( cd "$FULL_DIR" && "$HERE/run_full_envelope_bias.sh" --ws "$DROP_DIR/Datacard.root" --cat $c \
          --fit-mode wrongfit $COMMON --max-jobs 12 $BEST > "fullenv_wrongfit_$c.log" 2>&1
      echo "=== $(date) full-envelope wrongfit $c finished (rc=$?) ===" ) &
done

wait
echo "=== $(date) all M900_v4 resolved bias checks finished ==="
