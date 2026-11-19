#!/bin/bash
# Envelope discrete-profiling bias check where the truth function is
# GENUINELY ABSENT from the fit's envelope, not just frozen away from.
#
# This is the sibling of run_full_envelope_bias.sh's --fit-mode wrongfit,
# but strictly stronger: wrongfit freezes cat J's pdfindex to a function
# k != i while i is still a candidate in the RooMultiPdf (the fit is simply
# forced off it); here, function i is surgically removed from cat J's
# RooMultiPdf itself, so the fit's discrete profiling can never recover the
# truth at all. That isolates bias coming from the envelope being an
# incomplete truncation of function space, as opposed to discrete
# profiling merely landing on the wrong candidate.
#
#   Generation: on the ORIGINAL (full-envelope) combined workspace --
#   --setParameters sets cat J's pdfindex to i and --freezeParameters
#   freezes it there, plus every other category's pdfindex set + frozen to
#   its BEST-FIT function from real data (same --bestfit-ws/--bestfit-cache
#   mechanism as run_full_envelope_bias.sh, via get_bestfit_pdfindices.py).
#   No --channel-masks.
#
#   Fit: on a SEPARATE combined workspace built from a patched copy of the
#   datacard, whose only difference is that cat J's `bkg_mass` shapes line
#   points at a rebuilt per-category multipdf file with function i removed
#   (see build_dropfunc_workspace.py). Every other category's envelope
#   floats normally; only MH is frozen -- fit exactly like real data, on a
#   workspace that structurally cannot pick i for cat J.
#
# Prerequisites:
#   - The ORIGINAL combined workspace already built:
#       text2workspace.py Datacard_xxx.txt -m 125 -o Datacard_xxx.root
#   - The datacard's own `shapes bkg_mass <cat> ...` line resolves to a
#     real per-category multipdf file (read directly from the datacard --
#     categories may point at different subdirectories, e.g. background vs
#     background_alt vs background_boosted, each with a different function
#     count, so this is never assumed/hardcoded), and its "multipdf:<name>"
#     token gives the exact top-level pdf object name passed through to
#     build_dropfunc_workspace.py.
#
# 2D categories (hgg+hbb legs): pass --leg hgg|hbb to say which leg's
# envelope to drop a function from; the OTHER leg is generated and fit with
# its full envelope, untouched -- deliberately NOT a (i_hgg, i_hbb) cross
# product (see build_dropfunc_workspace.py's docstring for why: one leg at
# a time keeps this linear in the number of categories x functions instead
# of combinatorial, and answers "which leg's mismodeling biases r" directly
# rather than conflating the two). NOTE: the --leg hbb path has not been
# exercised against a real two-leg workspace -- this toolkit's current
# datacards are all 1D/hgg-only -- so verify against a real file first.
#
# Usage:
#   ./run_dropfunc_envelope_bias.sh --ws Datacard_xxx.root --datacard Datacard_xxx.txt --cat cat2 [--leg hgg|hbb] [--run] [--parallel]
#
# Without --run, only prints the commands it would execute (dry-run), and
# still performs the (cheap) workspace-surgery + text2workspace build steps
# so you can inspect the resulting files before committing to toys.
#
# --parallel launches every (truth function x r-truth) combination's
# GenerateOnly+FitDiagnostics pair as a background job at once instead of
# looping through them one at a time -- each combination is fully
# independent (different truth function and/or different injected r), so
# on a multi-core machine this takes roughly as long as ONE combination
# instead of N x that. Plotting still runs afterward, one combination at a
# time -- every combination's plot appends to the SAME pull_summary.json
# in the output directory, so plotting concurrently would race on that
# file; only the (much cheaper) toy generation+fitting is parallelized.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

show_help() {
    echo "Usage: $0 --ws <workspace.root> --datacard <datacard.txt> --cat <catN> [--leg hgg|hbb] [--run] [--ntoys N] [--seed S] [--rmin X] [--rmax Y] [--stability-flags]"
    echo ""
    echo "  --ws <file>        Original, full-envelope combined workspace (text2workspace.py output) (required)"
    echo "  --datacard <file>  The datacard --ws was built from -- its 'shapes bkg_mass <cat>' line is"
    echo "                     read to find cat J's real per-category multipdf file (required)"
    echo "  --cat <catN>       Category J to drop a function from, e.g. cat2 (required)"
    echo "  --leg hgg|hbb      Which leg's envelope to drop a function from (default hgg). Only matters"
    echo "                     for a 2D category with both legs; the OTHER leg's envelope is left intact."
    echo "  --run              Actually execute combine; otherwise dry-run (still builds the dropped-function workspaces)"
    echo "  --parallel         Run every (truth function x r) combination's GEN+FIT concurrently instead of one at a time"
    echo "                     (plotting still happens afterward, sequentially -- see header comment)"
    echo "  --publish          After plotting, copy the pull plots + pull_summary.json to"
    echo "                     ~/public_html/BiasSensitivityStudies_v7_final_dropfunc_<cat>_rMin<X>rMax<Y>_<N>toys/StandardStudies_dropfunc/"
    echo "                     (reuses whatever gallery index.php already exists under ~/public_html, if any)"
    echo "  --ntoys            Toys per combine call (default 1000)"
    echo "  --seed             combine -s seed (default 123456)"
    echo "  --rmin / --rmax    POI range passed to combine (default -5 / 5)"
    echo "  --stability-flags      Add the --X-rtd MINIMIZER_* flags (OFF by default: MINIMIZER_freezeDisassociatedParams"
    echo "                         segfaults FitDiagnostics on this analysis, confirmed even on the untouched full-envelope"
    echo "                         workspace -- exit code 139 -- and the other X-rtd flags gave no measurable speedup here)"
    echo "  --bestfit-ws <file>    Fitted-to-data workspace to read other categories' best-fit pdfindex from, if the cache"
    echo "                         doesn't exist yet (default: \$SCRIPT_DIR/higgsCombine_initial.MultiDimFit.mH125.root,"
    echo "                         i.e. the output of RunBiasInSignificance.py --mode setup)"
    echo "  --bestfit-cache <file> Cached best-fit pdfindex bash snippet (default: \$SCRIPT_DIR/bestfit_pdfindices.sh);"
    echo "                         computed once from --bestfit-ws via get_bestfit_pdfindices.py and reused after that"
    echo "  --truth-funcs <list>   Comma-separated subset of truth functions to run (default: all 7) -- use this to"
    echo "                         resume a run and only fill in combinations missing from a previous pass"
    echo "  --r-values <list>      Comma-separated subset of r_truth values to run (default: 0,2)"
    echo "  --r-expected <X>       Per-injection POI window: each (truth, r_inj) combination is generated + fit with"
    echo "                         --rMin r_inj-X --rMax r_inj+X instead of the global --rmin/--rmax (which are then"
    echo "                         ignored). Also changes the default --r-values to 0,X. Typically X = the expected"
    echo "                         upper limit on r. The window is also passed to the plotting script's boundary cut."
    echo "  --fixed-pull <P>       Also report median_pull_fixed in pull_summary.json: bad toys counted at a fixed pull"
    echo "                         instead of dropped, like submitBiasTests.py --failedFitR -P --keepEdgeFits P"
    echo "                         --keepLopsidedFits P --clipOutlierPulls (failed = -P, at/near rMin = -P, at/near"
    echo "                         rMax = +P, |pull|>=5 = +-P). Median only; histogram and Gaussian fit unchanged."
    echo "                         Plots go to bias_plots_dropfunc_<cat>..._fixedPull<P>/ (fit outputs are shared)."
    exit 0
}

WS=""; DATACARD=""; CAT=""; LEG="hgg"; DO_RUN=false; PARALLEL=false; PUBLISH=false
NTOYS=1000; SEED=123456; RMIN=-5; RMAX=5; STABILITY=false
BESTFIT_WS=""; BESTFIT_CACHE=""
TRUTH_FUNCS_ARG=""; R_VALUES_ARG=""; R_EXPECTED=""; FIXED_PULL=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --help) show_help ;;
        --ws) WS="$2"; shift 2 ;;
        --datacard) DATACARD="$2"; shift 2 ;;
        --cat) CAT="$2"; shift 2 ;;
        --leg) LEG="$2"; shift 2 ;;
        --run) DO_RUN=true; shift ;;
        --parallel) PARALLEL=true; shift ;;
        --publish) PUBLISH=true; shift ;;
        --ntoys) NTOYS="$2"; shift 2 ;;
        --seed) SEED="$2"; shift 2 ;;
        --rmin) RMIN="$2"; shift 2 ;;
        --rmax) RMAX="$2"; shift 2 ;;
        --stability-flags) STABILITY=true; shift ;;
        --bestfit-ws) BESTFIT_WS="$2"; shift 2 ;;
        --bestfit-cache) BESTFIT_CACHE="$2"; shift 2 ;;
        --truth-funcs) TRUTH_FUNCS_ARG="$2"; shift 2 ;;
        --r-values) R_VALUES_ARG="$2"; shift 2 ;;
        --r-expected) R_EXPECTED="$2"; shift 2 ;;
        --fixed-pull) FIXED_PULL="$2"; shift 2 ;;
        *) echo "[ERROR] Unknown option: $1"; show_help ;;
    esac
done

[ -z "$WS" ] && { echo "[ERROR] --ws <workspace.root> is required"; exit 1; }
[ -z "$DATACARD" ] && { echo "[ERROR] --datacard <datacard.txt> is required"; exit 1; }
[ -z "$CAT" ] && { echo "[ERROR] --cat <catN> is required"; exit 1; }
[ ! -f "$WS" ] && { echo "[ERROR] Workspace $WS not found."; exit 1; }
[ ! -f "$DATACARD" ] && { echo "[ERROR] Datacard $DATACARD not found."; exit 1; }
[ "$LEG" != "hgg" ] && [ "$LEG" != "hbb" ] && { echo "[ERROR] --leg must be 'hgg' or 'hbb'"; exit 1; }

CAT_NUM="${CAT#cat}"
if [ -n "$TRUTH_FUNCS_ARG" ]; then
    IFS=',' read -r -a TRUTH_FUNCTIONS <<< "$TRUTH_FUNCS_ARG"
else
    TRUTH_FUNCTIONS=("exp1" "bern1" "bern2" "bern3" "pow1" "lau1" "lau2")
fi
if [ -n "$FIXED_PULL" ]; then
    awk -v x="$FIXED_PULL" 'BEGIN{exit !(x+0 > 0)}' || { echo "[ERROR] --fixed-pull must be a positive number, got '$FIXED_PULL'"; exit 1; }
fi
if [ -n "$R_EXPECTED" ]; then
    awk -v x="$R_EXPECTED" 'BEGIN{exit !(x+0 > 0)}' || { echo "[ERROR] --r-expected must be a positive number, got '$R_EXPECTED'"; exit 1; }
fi
if [ -n "$R_VALUES_ARG" ]; then
    IFS=',' read -r -a R_TRUTH_VALUES <<< "$R_VALUES_ARG"
elif [ -n "$R_EXPECTED" ]; then
    R_TRUTH_VALUES=(0 "$R_EXPECTED")
else
    R_TRUTH_VALUES=(0 2)
fi

# POI window for one injected r: [r_inj - R_EXPECTED, r_inj + R_EXPECTED] when
# --r-expected is given, else the global --rmin/--rmax. Prints "<rmin> <rmax>".
r_window() {
    if [ -n "$R_EXPECTED" ]; then
        awk -v r="$1" -v x="$R_EXPECTED" 'BEGIN{printf "%g %g\n", r - x, r + x}'
    else
        echo "$RMIN $RMAX"
    fi
}
# With --r-expected the window is also passed to the plotting script's
# boundary-pinned-toy cut (a narrower window makes pinned toys more common).
# Without it, nothing is passed, so default behaviour is unchanged.
plot_range_args() {
    [ -z "$R_EXPECTED" ] && [ -z "$FIXED_PULL" ] && return 0
    local lo hi; read -r lo hi <<< "$(r_window "$1")"
    echo " --r-min $lo --r-max $hi${FIXED_PULL:+ --fixed-pull $FIXED_PULL}"
}
# Tags so --r-expected runs never share output files / pull_summary.json
# with fixed-range runs (or with runs using a different --r-expected).
if [ -n "$R_EXPECTED" ]; then
    RANGE_TAG="rInjPm${R_EXPECTED}"
    NAME_TAG="_${RANGE_TAG}"
else
    RANGE_TAG="rMin${RMIN}rMax${RMAX}"
    NAME_TAG=""
fi

DATACARD_DIR="$(cd "$(dirname "$DATACARD")" && pwd)"
DATACARD_BASE="$(basename "$DATACARD" .txt)"

# Real per-category source file AND top-level pdf object name, read
# directly off the datacard -- never assumed, since categories can point
# at different subdirectories (and, for 2D categories, the object after
# "multipdf:" may be a wrapper combining both legs, not a RooMultiPdf
# itself).
SHAPES_LINE=$(awk -v c="$CAT" '$1=="shapes" && $2=="bkg_mass" && $3==c {print; exit}' "$DATACARD")
[ -z "$SHAPES_LINE" ] && { echo "[ERROR] No 'shapes bkg_mass $CAT ...' line found in $DATACARD"; exit 1; }
SRC_REL=$(echo "$SHAPES_LINE" | awk '{print $4}')
TOP_PDF_TOKEN=$(echo "$SHAPES_LINE" | awk '{print $5}')
TOP_PDF_NAME="${TOP_PDF_TOKEN#*:}"
[ -z "$TOP_PDF_NAME" ] && { echo "[ERROR] Could not parse top-level pdf name out of: $SHAPES_LINE"; exit 1; }
if [[ "$SRC_REL" = /* ]]; then SRC_ABS="$SRC_REL"; else SRC_ABS="$DATACARD_DIR/$SRC_REL"; fi
[ ! -f "$SRC_ABS" ] && { echo "[ERROR] $CAT's bkg_mass source file not found: $SRC_ABS"; exit 1; }
SRC_BASENAME="$(basename "$SRC_REL" .root)"

echo "=== Resolving pdfindex legs and function map for $CAT ==="
eval "$(python3 "$SCRIPT_DIR/resolve_pdfindex_legs.py" "$WS" "$CAT_NUM")"

if [ "$LEG" = "hbb" ]; then
    TARGET_VAR="$HBB_VAR"
else
    TARGET_VAR="$HGG_VAR"
fi
if [ -z "$TARGET_VAR" ]; then
    echo "[ERROR] $CAT has no '$LEG' pdfindex variable in $WS."
    exit 1
fi

# --- Best-fit pdfindex for every OTHER category, from real data. Same
#     mechanism as run_full_envelope_bias.sh: computed once (via
#     get_bestfit_pdfindices.py, reading RooMultiPdf::getCurrentPdf()->
#     GetName() off a fitted-to-data workspace snapshot) and cached, since
#     it doesn't depend on CAT/t_func/etc. ---
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

# Every category OTHER than J (both of J's own legs are excluded, even the
# untouched one if J is 2D -- that leg stays fully free, per the header
# comment above) is set + frozen to its best-fit function from data.
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

# Generation touches cat J (freeze+set ONLY its pdfindex to the truth i)
# AND every other category (freeze+set its pdfindex to its best-fit
# function from data, OTHER_SET/OTHER_FREEZE above) -- GenerateOnly doesn't
# fit, so this only affects what background shape the toy is drawn from in
# categories other than J, making them "exactly what the data already
# prefers" instead of the workspace's structural default (index 0).
FREEZE_GEN="MH,$TARGET_VAR${OTHER_FREEZE:+,$OTHER_FREEZE}"
FREEZE_FIT="MH"

STABILITY_FLAGS=""
if [ "$STABILITY" = true ]; then
    STABILITY_FLAGS="--X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2"
fi

DROPWS_DIR="$DATACARD_DIR/dropfunc_workspaces"
mkdir -p "$DROPWS_DIR"

echo ""
if [ "$DO_RUN" = false ]; then
    echo "DRY-RUN MODE: only printing combine commands. Add --run to execute."
fi

LEG_TAG=""
[ "$LEG" = "hbb" ] && LEG_TAG="_hbbleg"
OUTDIR="./bias_plots_dropfunc_${CAT}${LEG_TAG}${NAME_TAG}${FIXED_PULL:+_fixedPull${FIXED_PULL}}"

# --- Incremental publish target: set up once up front (not just at the end)
#     so run_one() can copy each combination's plot to public_html as soon
#     as it's ready, instead of waiting for the whole category to finish. ---
if [ "$PUBLISH" = true ]; then
    PUB_DIR="$HOME/public_html/BiasSensitivityStudies_v7_final_dropfunc_${CAT}${LEG_TAG}_${RANGE_TAG}_${NTOYS}toys"
    mkdir -p "$PUB_DIR/StandardStudies_dropfunc"
    TEMPLATE_INDEX=$(find "$HOME/public_html" -maxdepth 3 -iname "index.php" 2>/dev/null | head -1)
    if [ -n "$TEMPLATE_INDEX" ]; then
        cp "$TEMPLATE_INDEX" "$PUB_DIR/StandardStudies_dropfunc/index.php"
        cp "$TEMPLATE_INDEX" "$PUB_DIR/index.php"
    fi
fi

# --- Phase 1 (always sequential -- cheap): for every truth function cat J
#     actually has, drop it from the envelope, patch a copy of the
#     datacard, build the reduced combined workspace. Cached: skipped if
#     the file already exists from a previous run. ---
VALID_TFUNCS=()
declare -A DROP_COMBINED_WS_OF SET_GEN_OF
for t_func in "${TRUTH_FUNCTIONS[@]}"; do
    idx_gen="IDX_${TARGET_VAR}_${t_func}"
    [ -z "${!idx_gen}" ] && continue

    DROP_WS="$DROPWS_DIR/${SRC_BASENAME}_drop_${t_func}${LEG_TAG}.root"
    if [ ! -f "$DROP_WS" ]; then
        echo "[BUILD] Dropping '$t_func' from $CAT's $LEG-leg envelope -> $DROP_WS"
        python3 "$SCRIPT_DIR/build_dropfunc_workspace.py" "$SRC_ABS" "$CAT_NUM" "$TOP_PDF_NAME" "$t_func" "$DROP_WS" --leg "$LEG" || { echo "[ERROR] build_dropfunc_workspace.py failed for $t_func"; continue; }
    fi

    PATCHED_DATACARD="$DATACARD_DIR/${DATACARD_BASE}_${CAT}_drop_${t_func}${LEG_TAG}.txt"
    if [ ! -f "$PATCHED_DATACARD" ]; then
        python3 "$SCRIPT_DIR/patch_datacard_dropfunc.py" "$DATACARD" "$CAT" "dropfunc_workspaces/$(basename "$DROP_WS")" "$PATCHED_DATACARD" || { echo "[ERROR] patch_datacard_dropfunc.py failed for $t_func"; continue; }
    fi

    DROP_COMBINED_WS="$DATACARD_DIR/${DATACARD_BASE}_${CAT}_drop_${t_func}${LEG_TAG}.root"
    if [ ! -f "$DROP_COMBINED_WS" ]; then
        echo "[BUILD] text2workspace.py -> $DROP_COMBINED_WS"
        ( cd "$DATACARD_DIR" && text2workspace.py "$(basename "$PATCHED_DATACARD")" -m 125 -o "$(basename "$DROP_COMBINED_WS")" ) || { echo "[ERROR] text2workspace.py failed for $t_func"; continue; }
    fi

    VALID_TFUNCS+=("$t_func")
    DROP_COMBINED_WS_OF[$t_func]="$DROP_COMBINED_WS"
    SET_GEN_OF[$t_func]="${!idx_gen}${OTHER_SET:+,$OTHER_SET}"
done

run_one() {
    local t_func="$1" r_truth="$2" drop_ws="$3" set_gen="$4"
    local name="${CAT}_dropfunc${LEG_TAG}_truth_${t_func}_r${r_truth}${NAME_TAG}"
    local rmin rmax; read -r rmin rmax <<< "$(r_window "$r_truth")"
    local cmd_gen="combine $WS -M GenerateOnly -m 125 -t $NTOYS -s $SEED -n .gen_${name} --expectSignal ${r_truth} --setParameters ${set_gen} --freezeParameters ${FREEZE_GEN} --cminDefaultMinimizerStrategy 0 --rMin $rmin --rMax $rmax --saveToys"
    local cmd_fit="combine -M FitDiagnostics $drop_ws -m 125 -t $NTOYS -s $SEED -n .fit_${name} --expectSignal ${r_truth} --toysFile higgsCombine.gen_${name}.GenerateOnly.mH125.${SEED}.root --freezeParameters ${FREEZE_FIT} --cminDefaultMinimizerStrategy 0 --robustFit 1 --rMin $rmin --rMax $rmax $STABILITY_FLAGS"
    local cmd_plot="python3 $SCRIPT_DIR/plot_bias_pull_FitDiagnostics.py -i fitDiagnostics.fit_${name}.root --r-truth ${r_truth} --truth-func ${t_func} --fit-func envelope_minus_${t_func} --cat ${name} -o $OUTDIR$(plot_range_args "$r_truth")"

    echo "-------------------------------------------------------------------------"
    echo "   [dropfunc] $CAT ($LEG leg) | Truth [${t_func^^}] absent from fit envelope | r=${r_truth} | r range [${rmin}, ${rmax}]"
    echo "-------------------------------------------------------------------------"
    echo "  [GEN]  $cmd_gen"
    echo "  [FIT]  $cmd_fit"
    echo "  [PLOT] $cmd_plot"

    [ "$DO_RUN" = false ] && return 0

    eval $cmd_gen > "combine_gen_${name}.log" 2>&1
    if [ $? -ne 0 ]; then echo "[ERROR] Generation failed for ${name}. See combine_gen_${name}.log"; return 1; fi
    eval $cmd_fit > "combine_fit_${name}.log" 2>&1
    if [ $? -ne 0 ]; then echo "[ERROR] Fit failed for ${name}. See combine_fit_${name}.log"; return 1; fi
    if [ "$5" = "plot_now" ]; then
        mkdir -p "$OUTDIR"
        eval $cmd_plot
        # Incremental publish: copy just THIS combination's plot(s) + the
        # running pull_summary.json to public_html right away, instead of
        # waiting for every combination in the category to finish.
        if [ "$PUBLISH" = true ]; then
            cp "$OUTDIR"/*"${name}"*.png "$PUB_DIR/StandardStudies_dropfunc/" 2>/dev/null
            cp "$OUTDIR"/pull_summary.json "$PUB_DIR/StandardStudies_dropfunc/" 2>/dev/null
        fi
    fi
    return 0
}

if [ "$PARALLEL" = false ]; then
    # --- Sequential (default): GEN, FIT, PLOT, one combination at a time. ---
    for t_func in "${VALID_TFUNCS[@]}"; do
        for r_truth in "${R_TRUTH_VALUES[@]}"; do
            run_one "$t_func" "$r_truth" "${DROP_COMBINED_WS_OF[$t_func]}" "${SET_GEN_OF[$t_func]}" plot_now
        done
    done
else
    # --- Parallel: launch every combination's GEN+FIT at once, wait for
    #     all of them, THEN plot sequentially (shared pull_summary.json). ---
    echo "=== [parallel] launching $(( ${#VALID_TFUNCS[@]} * ${#R_TRUTH_VALUES[@]} )) GEN+FIT combinations at once ==="
    PIDS=()
    for t_func in "${VALID_TFUNCS[@]}"; do
        for r_truth in "${R_TRUTH_VALUES[@]}"; do
            run_one "$t_func" "$r_truth" "${DROP_COMBINED_WS_OF[$t_func]}" "${SET_GEN_OF[$t_func]}" no_plot &
            PIDS+=($!)
        done
    done
    if [ "$DO_RUN" = true ]; then
        echo "Launched PIDs: ${PIDS[@]}"
        for pid in "${PIDS[@]}"; do wait "$pid"; done
        echo "=== [parallel] all GEN+FIT jobs finished -- plotting sequentially ==="
        mkdir -p "$OUTDIR"
        for t_func in "${VALID_TFUNCS[@]}"; do
            for r_truth in "${R_TRUTH_VALUES[@]}"; do
                name="${CAT}_dropfunc${LEG_TAG}_truth_${t_func}_r${r_truth}${NAME_TAG}"
                if [ -f "fitDiagnostics.fit_${name}.root" ]; then
                    python3 "$SCRIPT_DIR/plot_bias_pull_FitDiagnostics.py" -i "fitDiagnostics.fit_${name}.root" --r-truth "${r_truth}" --truth-func "${t_func}" --fit-func "envelope_minus_${t_func}" --cat "${name}" -o "$OUTDIR" $(plot_range_args "$r_truth")
                else
                    echo "[SKIP] ${name}: no fitDiagnostics output (see combine_fit_${name}.log)"
                fi
            done
        done
    fi
fi

if [ "$DO_RUN" = true ] && [ "$PUBLISH" = true ]; then
    PUB_DIR="$HOME/public_html/BiasSensitivityStudies_v7_final_dropfunc_${CAT}${LEG_TAG}_${RANGE_TAG}_${NTOYS}toys"
    mkdir -p "$PUB_DIR/StandardStudies_dropfunc"
    cp "$OUTDIR"/*.png "$PUB_DIR/StandardStudies_dropfunc/" 2>/dev/null
    cp "$OUTDIR"/pull_summary.json "$PUB_DIR/StandardStudies_dropfunc/" 2>/dev/null
    TEMPLATE_INDEX=$(find "$HOME/public_html" -maxdepth 3 -iname "index.php" 2>/dev/null | head -1)
    if [ -n "$TEMPLATE_INDEX" ]; then
        cp "$TEMPLATE_INDEX" "$PUB_DIR/StandardStudies_dropfunc/index.php"
        cp "$TEMPLATE_INDEX" "$PUB_DIR/index.php"
    fi
    echo "=== [publish] copied plots + pull_summary.json to $PUB_DIR/StandardStudies_dropfunc/ ==="
fi
