#!/bin/bash
# Envelope discrete-profiling bias check on the full COMBINED workspace, with
# NO channel masking -- every category contributes to the toy and to the fit
# exactly as it would for real data. Only category J's truth function is
# perturbed; every other category is generated (and, in --fit-mode envelope,
# fit) with its RooMultiPdf envelope handled normally.
#
# This is the un-masked sibling of run_envelope_bias.sh (see that script /
# README.md check #2 for the --channel-masks version, which isolates one
# category by switching every other one off). Use THIS script when you want
# to know: "if category J's background actually came from truth function i,
# while every other category's background is exactly what the envelope
# already assumes (index 0), how does the fit behave -- (a) with the full
# discrete-profiling envelope active everywhere, exactly like on real data,
# vs (b) if category J's envelope got stuck on/forced onto the WRONG
# function, with every other category still profiled normally?"
#
#   --fit-mode envelope  (fit with the full envelope, everywhere)
#       Every category's pdfindex is left FLOATING (discrete-profiled), same
#       as combine does on real data -- including category J. Loops over
#       TRUTH_FUNCTIONS x R_TRUTH_VALUES.
#
#   --fit-mode wrongfit  (fit with a PDF != i only in cat J, envelope
#                         elsewhere)
#       Category J's pdfindex is FROZEN to some function k != i (looping
#       over every OTHER function present in J's envelope); every other
#       category's pdfindex is left floating (discrete-profiled), same as
#       'envelope' mode. Loops over TRUTH_FUNCTIONS x FIT_FUNCTIONS (k != i)
#       x R_TRUTH_VALUES.
#
# In BOTH modes, generation touches category J AND every other category:
# --setParameters sets J's pdfindex to i (--freezeParameters freezes it),
# and every OTHER category's pdfindex is explicitly set + frozen to its
# BEST-FIT function as determined from real data (a MultiDimFit to data
# with the envelope floating everywhere, e.g. via RunBiasInSignificance.py
# --mode setup), read via RooMultiPdf::getCurrentPdf()->GetName() -- see
# get_bestfit_pdfindices.py. This replaces the old behavior of leaving
# other categories at the workspace's structural default index (0, i.e.
# the first function in each envelope, which need not match what the data
# actually prefers): now "every other category's background is exactly
# what the discrete profiling already assumes" means the data's own
# best-fit choice, not just index 0. See --bestfit-ws / --bestfit-cache.
#
# Because every other category's envelope is left floating in BOTH fit
# modes here (unlike run_envelope_bias.sh, which freezes them since they're
# masked out), a full run over many categories x many functions is a much
# bigger discrete-combinatorial fit (every category's envelope is profiled
# simultaneously), so expect it to be slower than run_envelope_bias.sh. The
# --X-rtd MINIMIZER_* stability flags below are OFF by default, though:
# MINIMIZER_freezeDisassociatedParams segfaults FitDiagnostics on this
# analysis even WITHOUT --channel-masks (confirmed: exit code 139 on the
# plain, unmasked full-envelope workspace too -- not just the masked case
# run_envelope_bias.sh already avoids them for), and the remaining X-rtd
# flags gave no measurable speedup here. Opt in with --stability-flags if
# you want to re-test them on your own workspace.
#
# Prerequisite (once per datacard, NO --channel-masks needed):
#   text2workspace.py Datacard_xxx.txt -m 125 -o Datacard_xxx.root
#
# Usage:
#   ./run_full_envelope_bias.sh --ws Datacard_xxx.root --cat cat2 --fit-mode envelope [--run] [--only-hgg]
#   ./run_full_envelope_bias.sh --ws Datacard_xxx.root --cat cat2 --fit-mode wrongfit [--run] [--only-hgg]
#
# Without --run, only prints the combine commands (dry-run).

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

show_help() {
    echo "Usage: $0 --ws <workspace.root> --cat <catN> [--fit-mode envelope|wrongfit] [--run] [--only-hgg] [--ntoys N] [--seed S] [--rmin X] [--rmax Y] [--stability-flags]"
    echo ""
    echo "  --ws <file>       Combined workspace (plain text2workspace.py output, no --channel-masks needed) (required)"
    echo "  --cat <catN>      Category J whose truth function is perturbed, e.g. cat2 (required)"
    echo "  --fit-mode        'envelope' (default): fit with the full envelope floating everywhere."
    echo "                    'wrongfit': freeze ONLY cat J's pdfindex to a function != truth, envelope floats everywhere else."
    echo "  --run             Actually execute combine; otherwise dry-run (prints commands only)"
    echo "  --parallel        Run every (truth x fit x r) combination's GEN+FIT as a background job instead of one at"
    echo "                    a time; plotting runs afterward, sequentially (all plots append to one pull_summary.json)"
    echo "  --max-jobs <N>    With --parallel: at most N combinations running at once (default 0 = all at once)"
    echo "  --only-hgg        Ignore an HBB leg on this category, if present"
    echo "  --ntoys           Toys per combine call (default 1000)"
    echo "  --seed            combine -s seed (default 123456)"
    echo "  --rmin / --rmax   POI range passed to combine (default -5 / 5)"
    echo "  --stability-flags      Add the --X-rtd MINIMIZER_* flags (OFF by default: MINIMIZER_freezeDisassociatedParams"
    echo "                         segfaults FitDiagnostics on this analysis, confirmed even on the untouched full-envelope"
    echo "                         workspace -- exit code 139 -- and the other X-rtd flags gave no measurable speedup here)"
    echo "  --bestfit-ws <file>    Fitted-to-data workspace to read other categories' best-fit pdfindex from, if the cache"
    echo "                         doesn't exist yet (default: \$SCRIPT_DIR/higgsCombine_initial.MultiDimFit.mH125.root,"
    echo "                         i.e. the output of RunBiasInSignificance.py --mode setup)"
    echo "  --bestfit-cache <file> Cached best-fit pdfindex bash snippet (default: \$SCRIPT_DIR/bestfit_pdfindices.sh);"
    echo "                         computed once from --bestfit-ws via get_bestfit_pdfindices.py and reused after that"
    echo "  --r-values <list>      Comma-separated r_truth values to run (default: 0,2; 0,X with --r-expected X)"
    echo "  --r-expected <X>       Per-injection POI window: each combination is generated + fit with"
    echo "                         --rMin r_inj-X --rMax r_inj+X instead of --rmin/--rmax (which are then ignored)."
    echo "                         Default --r-values become 0,X. Typically X = the expected upper limit on r."
    echo "                         Outputs get an _rInjPm<X> suffix."
    echo "  --fixed-pull <P>       Also report median_pull_fixed in pull_summary.json: bad toys counted at a fixed pull"
    echo "                         instead of dropped, like submitBiasTests.py --failedFitR -P --keepEdgeFits P"
    echo "                         --keepLopsidedFits P --clipOutlierPulls (failed = -P, at/near rMin = -P, at/near"
    echo "                         rMax = +P, |pull|>=5 = +-P). Median only; histogram and Gaussian fit unchanged."
    echo "                         Plots go to bias_plots_full_<cat>_<mode>..._fixedPull<P>/ (fit outputs are shared)."
    exit 0
}

WS=""; CAT=""; FIT_MODE="envelope"; DO_RUN=false; ONLY_HGG=false; PARALLEL=false; MAX_JOBS=0
PAR_NAMES=(); PAR_PLOTS=()
NTOYS=1000; SEED=123456; RMIN=-5; RMAX=5; STABILITY=false
BESTFIT_WS=""; BESTFIT_CACHE=""
R_VALUES_ARG=""; R_EXPECTED=""; FIXED_PULL=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --help) show_help ;;
        --ws) WS="$2"; shift 2 ;;
        --cat) CAT="$2"; shift 2 ;;
        --fit-mode) FIT_MODE="$2"; shift 2 ;;
        --run) DO_RUN=true; shift ;;
        --only-hgg) ONLY_HGG=true; shift ;;
        --ntoys) NTOYS="$2"; shift 2 ;;
        --seed) SEED="$2"; shift 2 ;;
        --rmin) RMIN="$2"; shift 2 ;;
        --rmax) RMAX="$2"; shift 2 ;;
        --stability-flags) STABILITY=true; shift ;;
        --parallel) PARALLEL=true; shift ;;
        --max-jobs) MAX_JOBS="$2"; shift 2 ;;
        --bestfit-ws) BESTFIT_WS="$2"; shift 2 ;;
        --bestfit-cache) BESTFIT_CACHE="$2"; shift 2 ;;
        --r-values) R_VALUES_ARG="$2"; shift 2 ;;
        --r-expected) R_EXPECTED="$2"; shift 2 ;;
        --fixed-pull) FIXED_PULL="$2"; shift 2 ;;
        *) echo "[ERROR] Unknown option: $1"; show_help ;;
    esac
done

[ -z "$WS" ] && { echo "[ERROR] --ws <workspace.root> is required"; exit 1; }
[ -z "$CAT" ] && { echo "[ERROR] --cat <catN> is required"; exit 1; }
[ ! -f "$WS" ] && { echo "[ERROR] Workspace $WS not found. Build it first with:"; echo "  text2workspace.py <datacard>.txt -m 125 -o $WS"; exit 1; }
if [ "$FIT_MODE" != "wrongfit" ] && [ "$FIT_MODE" != "envelope" ]; then
    echo "[ERROR] --fit-mode must be 'envelope' or 'wrongfit'"; exit 1
fi

CAT_NUM="${CAT#cat}"
TRUTH_FUNCTIONS=("exp1" "bern1" "bern2" "bern3" "pow1" "lau1" "lau2")
FIT_FUNCTIONS=("exp1" "bern1" "bern2" "bern3" "pow1" "lau1" "lau2")
[[ "$MAX_JOBS" =~ ^[0-9]+$ ]] || { echo "[ERROR] --max-jobs must be a non-negative integer, got '$MAX_JOBS'"; exit 1; }
for _opt in R_EXPECTED FIXED_PULL; do
    if [ -n "${!_opt}" ]; then
        awk -v x="${!_opt}" 'BEGIN{exit !(x+0 > 0)}' || { echo "[ERROR] --$(echo $_opt | tr 'A-Z_' 'a-z-') must be a positive number, got '${!_opt}'"; exit 1; }
    fi
done
if [ -n "$R_VALUES_ARG" ]; then
    IFS=',' read -r -a R_TRUTH_VALUES <<< "$R_VALUES_ARG"
elif [ -n "$R_EXPECTED" ]; then
    R_TRUTH_VALUES=(0 "$R_EXPECTED")
else
    R_TRUTH_VALUES=(0 2)
fi

# POI window for one injected r: [r_inj - R_EXPECTED, r_inj + R_EXPECTED] with
# --r-expected, else the global --rmin/--rmax. Prints "<rmin> <rmax>".
r_window() {
    if [ -n "$R_EXPECTED" ]; then
        awk -v r="$1" -v x="$R_EXPECTED" 'BEGIN{printf "%g %g\n", r - x, r + x}'
    else
        echo "$RMIN $RMAX"
    fi
}
# Extra plotting args: the fit window (boundary cut + fixed-pull edge
# classification) and --fixed-pull. Nothing without either option, so the
# default behaviour is unchanged.
plot_range_args() {
    [ -z "$R_EXPECTED" ] && [ -z "$FIXED_PULL" ] && return 0
    local lo hi; read -r lo hi <<< "$(r_window "$1")"
    echo " --r-min $lo --r-max $hi${FIXED_PULL:+ --fixed-pull $FIXED_PULL}"
}
NAME_TAG="${R_EXPECTED:+_rInjPm${R_EXPECTED}}"

echo "=== Resolving pdfindex legs and function map for $CAT (mode: $FIT_MODE) ==="
eval "$(python3 "$SCRIPT_DIR/resolve_pdfindex_legs.py" "$WS" "$CAT_NUM")"

HAS_HGG=false; HAS_HBB=false
[ -n "$HGG_VAR" ] && HAS_HGG=true
if [ -n "$HBB_VAR" ] && [ "$ONLY_HGG" = false ]; then HAS_HBB=true; fi

if [ "$HAS_HGG" = false ] && [ "$HAS_HBB" = false ]; then
    echo "[ERROR] No pdfindex variable found for $CAT in $WS."
    echo "        Check that $CAT is a real bin name in this workspace."
    exit 1
fi

TARGET_LIST=""
[ "$HAS_HGG" = true ] && TARGET_LIST="$HGG_VAR"
[ "$HAS_HBB" = true ] && TARGET_LIST="${TARGET_LIST:+$TARGET_LIST }$HBB_VAR"
TARGET_LIST_CSV="$(echo "$TARGET_LIST" | tr ' ' ',')"

# --- Best-fit pdfindex for every OTHER category, from real data. Computed
#     once (via get_bestfit_pdfindices.py, reading RooMultiPdf::
#     getCurrentPdf()->GetName() off a fitted-to-data workspace snapshot)
#     and cached, since it doesn't depend on CAT/FIT_MODE/t_func/etc. ---
[ -z "$BESTFIT_CACHE" ] && BESTFIT_CACHE="$SCRIPT_DIR/bestfit_pdfindices.sh"
[ -z "$BESTFIT_WS" ] && BESTFIT_WS="$SCRIPT_DIR/higgsCombine_initial.MultiDimFit.mH125.root"
if [ ! -f "$BESTFIT_CACHE" ]; then
    if [ ! -f "$BESTFIT_WS" ]; then
        echo "[ERROR] No best-fit pdfindex cache ($BESTFIT_CACHE) and no fitted workspace ($BESTFIT_WS)."
        echo "        Fit the envelope to real data first, e.g.:"
        echo "          python3 RunBiasInSignificance.py --inputWSFile $WS --mode setup --initial-fit-param <low_impact_nuisance>"
        echo "        then re-run this script (or pass --bestfit-ws <file> / --bestfit-cache <file>)."
        exit 1
    fi
    echo "=== Computing best-fit pdfindex per category from $BESTFIT_WS (caching to $BESTFIT_CACHE) ==="
    python3 "$SCRIPT_DIR/get_bestfit_pdfindices.py" "$BESTFIT_WS" > "$BESTFIT_CACHE"
fi
source "$BESTFIT_CACHE"

OTHER_SET=""; OTHER_FREEZE=""
for v in $ALL_PDFINDEX_VARS; do
    [ "$v" == "$HGG_VAR" ] && continue
    [ "$v" == "$HBB_VAR" ] && continue
    bidx_var="BESTIDX_${v}"
    bidx_val="${!bidx_var}"
    if [ -z "$bidx_val" ]; then
        echo "[WARN] No cached best-fit index for $v -- leaving it at the workspace's structural default for generation."
        continue
    fi
    OTHER_SET="${OTHER_SET:+$OTHER_SET,}${bidx_val}"
    OTHER_FREEZE="${OTHER_FREEZE:+$OTHER_FREEZE,}${v}"
done

# --- Generation touches category J (freeze + set its pdfindex to the truth
#     index i) AND every other category (freeze + set its pdfindex to its
#     best-fit function from data, OTHER_SET/OTHER_FREEZE above). No channel
#     masking: every category's data is generated and enters the toy. ---
FREEZE_GEN="MH${TARGET_LIST_CSV:+,$TARGET_LIST_CSV}${OTHER_FREEZE:+,$OTHER_FREEZE}"

# --- Fit: 'envelope' floats every category's pdfindex, including J.
#     'wrongfit' freezes ONLY J's pdfindex (to a function != truth); every
#     other category's pdfindex is left floating in both modes. ---
if [ "$FIT_MODE" == "wrongfit" ]; then
    FREEZE_FIT="MH${TARGET_LIST_CSV:+,$TARGET_LIST_CSV}"
else
    FREEZE_FIT="MH"
fi

STABILITY_FLAGS=""
if [ "$STABILITY" = true ]; then
    STABILITY_FLAGS="--X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2"
fi

echo ""
if [ "$DO_RUN" = false ]; then
    echo "DRY-RUN MODE: only printing commands. Add --run to execute."
fi

OUTDIR="./bias_plots_full_${CAT}_${FIT_MODE}${NAME_TAG}${FIXED_PULL:+_fixedPull${FIXED_PULL}}"

for t_func in "${TRUTH_FUNCTIONS[@]}"; do
    if [ "$FIT_MODE" == "envelope" ]; then
        FIT_LOOP=("envelope")
    else
        FIT_LOOP=("${FIT_FUNCTIONS[@]}")
    fi

    for f_func in "${FIT_LOOP[@]}"; do
        if [ "$FIT_MODE" == "wrongfit" ] && [ "$t_func" == "$f_func" ]; then continue; fi

        for r_truth in "${R_TRUTH_VALUES[@]}"; do
            SET_GEN=""; SET_FIT=""
            SKIP=false

            if [ "$HAS_HGG" = true ]; then
                idx_gen="IDX_${HGG_VAR}_${t_func}"
                [ -z "${!idx_gen}" ] && SKIP=true
                SET_GEN="${!idx_gen}"
                if [ "$FIT_MODE" == "wrongfit" ]; then
                    idx_fit="IDX_${HGG_VAR}_${f_func}"
                    [ -z "${!idx_fit}" ] && SKIP=true
                    SET_FIT="${!idx_fit}"
                fi
            fi
            if [ "$HAS_HBB" = true ]; then
                idx_gen_b="IDX_${HBB_VAR}_${t_func}"
                [ -z "${!idx_gen_b}" ] && SKIP=true
                SET_GEN="${SET_GEN:+$SET_GEN,}${!idx_gen_b}"
                if [ "$FIT_MODE" == "wrongfit" ]; then
                    idx_fit_b="IDX_${HBB_VAR}_${f_func}"
                    [ -z "${!idx_fit_b}" ] && SKIP=true
                    SET_FIT="${SET_FIT:+$SET_FIT,}${!idx_fit_b}"
                fi
            fi
            [ "$SKIP" = true ] && continue

            SET_GEN="${SET_GEN:+$SET_GEN,}${OTHER_SET}"

            NAME="${CAT}_full_${FIT_MODE}_truth_${t_func}_fit_${f_func}_r${r_truth}${NAME_TAG}"
            read -r R_LO R_HI <<< "$(r_window "$r_truth")"
            echo "-------------------------------------------------------------------------"
            echo "   [full-${FIT_MODE}] $CAT | Truth [${t_func^^}] vs Fit [${f_func^^}] | r=${r_truth} | r range [${R_LO}, ${R_HI}]"
            echo "-------------------------------------------------------------------------"

            CMD_GEN="combine $WS -M GenerateOnly -m 125 -t $NTOYS -s $SEED -n .gen_${NAME} --expectSignal ${r_truth} --setParameters ${SET_GEN} --freezeParameters ${FREEZE_GEN} --cminDefaultMinimizerStrategy 0 --rMin $R_LO --rMax $R_HI --saveToys"
            CMD_FIT="combine -M FitDiagnostics $WS -m 125 -t $NTOYS -s $SEED -n .fit_${NAME} --expectSignal ${r_truth} --toysFile higgsCombine.gen_${NAME}.GenerateOnly.mH125.${SEED}.root ${SET_FIT:+--setParameters ${SET_FIT}} --freezeParameters ${FREEZE_FIT} --cminDefaultMinimizerStrategy 0 --robustFit 1 --rMin $R_LO --rMax $R_HI $STABILITY_FLAGS"
            CMD_PLOT="python3 $SCRIPT_DIR/plot_bias_pull_FitDiagnostics.py -i fitDiagnostics.fit_${NAME}.root --r-truth ${r_truth} --truth-func ${t_func} --fit-func ${f_func} --cat ${NAME} -o $OUTDIR$(plot_range_args "$r_truth")"

            echo "  [GEN]  $CMD_GEN"
            echo "  [FIT]  $CMD_FIT"
            echo "  [PLOT] $CMD_PLOT"

            if [ "$DO_RUN" = true ] && [ "$PARALLEL" = true ]; then
                # Throttle to --max-jobs (0 = no limit), then launch GEN+FIT
                # in the background; plotting is deferred (shared JSON).
                if [ "$MAX_JOBS" -gt 0 ]; then
                    while [ "$(jobs -rp | wc -l)" -ge "$MAX_JOBS" ]; do wait -n; done
                fi
                (
                    eval $CMD_GEN > combine_gen_${NAME}.log 2>&1 || { echo "[ERROR] Generation failed for ${NAME}. See combine_gen_${NAME}.log"; exit 1; }
                    eval $CMD_FIT > combine_fit_${NAME}.log 2>&1 || { echo "[ERROR] Fit failed for ${NAME}. See combine_fit_${NAME}.log"; exit 1; }
                    echo "[DONE] ${NAME}"
                ) &
                PAR_NAMES+=("$NAME"); PAR_PLOTS+=("$CMD_PLOT")
            elif [ "$DO_RUN" = true ]; then
                eval $CMD_GEN > combine_gen_${NAME}.log 2>&1
                if [ $? -ne 0 ]; then echo "[ERROR] Generation failed for ${NAME}. See combine_gen_${NAME}.log"; continue; fi
                eval $CMD_FIT > combine_fit_${NAME}.log 2>&1
                if [ $? -eq 0 ]; then
                    mkdir -p "$OUTDIR"
                    eval $CMD_PLOT
                else
                    echo "[ERROR] Fit failed for ${NAME}. See combine_fit_${NAME}.log"
                fi
            fi
        done
    done
done

if [ "$DO_RUN" = true ] && [ "$PARALLEL" = true ]; then
    echo "=== [parallel] ${#PAR_NAMES[@]} GEN+FIT jobs launched (max at once: ${MAX_JOBS/#0/unlimited}); waiting ==="
    wait
    echo "=== [parallel] all GEN+FIT jobs finished -- plotting sequentially ==="
    mkdir -p "$OUTDIR"
    for i in "${!PAR_NAMES[@]}"; do
        if [ -f "fitDiagnostics.fit_${PAR_NAMES[$i]}.root" ]; then
            eval ${PAR_PLOTS[$i]}
        else
            echo "[SKIP] ${PAR_NAMES[$i]}: no fitDiagnostics output (see combine_gen/fit_${PAR_NAMES[$i]}.log)"
        fi
    done
fi
