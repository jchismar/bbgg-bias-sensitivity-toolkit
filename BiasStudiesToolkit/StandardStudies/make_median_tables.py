#!/usr/bin/env python3
"""Median-pull tables for the full-envelope (README section 3, envelope,
wrongfit and samefit modes) and dropfunc (section 7) checks, from their pull_summary.json
files. Writes one Markdown file with every table, plus one PNG per
(category, r_inj).

Each table: rows = truth function i generated in the category under test;
columns =
  envelope          full-envelope --fit-mode envelope (fit with the full envelope)
  fit <k>           full-envelope --fit-mode wrongfit (category frozen to k != i);
                    the k == i diagonal comes from --fit-mode samefit (shown as
                    '—' if that mode was not run)
  envelope - truth  dropfunc (truth i removed from the category's envelope)
Each cell is median_pull_fixed +- its error: the median with failed / edge /
lopsided / |pull|>=5 toys counted at -/+P (same rules as submitBiasTests.py).
Every run must have been made with --fixed-pull; a record without
median_pull_fixed is an error (the kept-toys-only median_pull is never used).

Usage:
  make_median_tables.py --fullenv-dir M900_v4_fullenv --dropfunc-dir M900_v4_dropfunc \
      [-o median_tables.md] [--png-dir median_tables_png] [--threshold 0.2]
"""
import argparse
import glob
import json
import os
import re

FUNC_ORDER = ["bern1", "bern2", "bern3", "exp1", "pow1", "lau1", "lau2"]


def cat_key(c):
    m = re.match(r"cat(\d+)", c)
    return (int(m.group(1)) if m else 99, c)


def func_key(f):
    return FUNC_ORDER.index(f) if f in FUNC_ORDER else len(FUNC_ORDER)


def load(pattern, regex):
    """{(cat, r_truth, truth, fit): record} from every pull_summary.json matching
    pattern; regex pulls the category out of the directory name. A later entry
    for the same key (a rerun appended to the json) wins."""
    out = {}
    for path in sorted(glob.glob(pattern)):
        m = re.search(regex, path)
        if not m:
            continue
        with open(path) as fh:
            for rec in json.load(fh):
                out[(m.group(1), round(rec["r_truth"], 4), rec["truth_func"], rec["fit_func"])] = rec
    return out


def load_bestfit(path):
    """{category: best-fit function tag} from a get_bestfit_pdfindices.py cache
    (bestfit_pdfindices.sh), e.g.
    BESTIDX_pdfindex_hgg_cat4Resolved_13TeV='...=2'  # env_pdf_..._pow1 (tag=pow1)"""
    best = {}
    if not path or not os.path.exists(path):
        return best
    with open(path) as fh:
        for line in fh:
            m = re.match(r"BESTIDX_pdfindex_(?:hgg_)?(.+?)_13TeV\S*=.*\(tag=([^)]+)\)", line)
            if m:
                best[m.group(1)] = m.group(2)
    return best


def value(rec):
    """(median_pull_fixed, its error, fixed pull P). Never falls back to median_pull."""
    if "median_pull_fixed" not in rec:
        raise SystemExit(f"[ERROR] {rec['category']}: no median_pull_fixed -- rerun its plots with --fixed-pull")
    return rec["median_pull_fixed"], rec["error_median_pull_fixed"], rec["fixed_pull"]


def build_tables(env, wrong, drop, thr):
    """List of tables: dict(cat, r, cols, rows=[(truth, [cell, ...])]); a cell is
    None (missing), '—' (diagonal), or (median_pull_fixed, error, fixed pull P, flagged)."""
    allkeys = list(env) + list(wrong) + list(drop)
    tables = []
    for cat in sorted({k[0] for k in allkeys}, key=cat_key):
        for r in sorted({k[1] for k in allkeys if k[0] == cat}):
            truths = sorted({k[2] for k in allkeys if k[0] == cat and k[1] == r}, key=func_key)
            fits = sorted({k[3] for k in wrong if k[0] == cat and k[1] == r}, key=func_key)
            cols = [("env", None, "envelope")] + [("wrong", f, f"fit {f}") for f in fits]
            if any(k[0] == cat and k[1] == r for k in drop):
                cols.append(("drop", None, "envelope − truth"))
            rows = []
            for t in truths:
                cells = []
                for kind, f, _ in cols:
                    if kind == "wrong" and f == t and (cat, r, t, f) not in wrong:
                        cells.append("—")  # samefit not run
                        continue
                    rec = {"env": env.get((cat, r, t, "envelope")),
                           "wrong": wrong.get((cat, r, t, f)),
                           "drop": drop.get((cat, r, t, f"envelope_minus_{t}"))}[kind]
                    if rec is None:
                        cells.append(None)
                    else:
                        med, err, vk = value(rec)
                        # flag on the displayed (2-decimal) value, so equal-looking cells agree
                        cells.append((med, err, vk, abs(round(med, 2)) >= thr))
                rows.append((t, cells))
            tables.append({"cat": cat, "r": r, "cols": [c[2] for c in cols], "rows": rows})
    return tables


def fmt(c, bold_md=False):
    if c is None:
        return ""
    if isinstance(c, str):
        return c
    med, err, _, flag = c
    s = f"{med:+.2f} ± {err:.2f}"
    return f"**{s}**" if (flag and bold_md) else s


def write_markdown(tables, path, src_full, src_drop, thr, pulls, best, extra=""):
    out = ["# Median-pull tables" + (f" ({extra})" if extra else ""), "",
           f"- full envelope: `{src_full}`",
           f"- dropfunc: `{src_drop}`", "",
           f"Each cell: `median_pull_fixed` ± error, i.e. the median pull with failed / edge / lopsided / "
           f"|pull| ≥ 5 toys counted at ∓{'/'.join(f'{x:g}' for x in sorted(pulls))} "
           "(same rules as `submitBiasTests.py --failedFitR -5 --keepEdgeFits 5 --keepLopsidedFits 5 "
           "--clipOutlierPulls`). Rows: the truth function generated in the category under test. "
           "Columns: **envelope** = fit with the full envelope floating (full envelope, `envelope` mode); "
           "**fit k** = that category frozen to function k (`wrongfit` mode for k ≠ truth, `samefit` mode for "
           "the k = truth diagonal); **envelope − truth** = the truth "
           "removed from that category's envelope (dropfunc). Other categories are generated at their best fit "
           f"on data and float in every fit. |median| ≥ {thr:g} in **bold**.", ""]
    cur = None
    for tb in tables:
        if tb["cat"] != cur:
            cur = tb["cat"]
            out += [f"## {cur}", ""]
        out += [f"#### {tb['cat']}, r_inj = {tb['r']:g} (best fit on data: {best.get(tb['cat'], '?')})", "",
                "| truth \\ fit | " + " | ".join(tb["cols"]) + " |",
                "|---|" + "---|" * len(tb["cols"])]
        for t, cells in tb["rows"]:
            out.append(f"| {t} | " + " | ".join(fmt(c, bold_md=True) for c in cells) + " |")
        out.append("")
    with open(path, "w") as fh:
        fh.write("\n".join(out) + "\n")


def write_pngs(tables, outdir, thr, label, best, extra=""):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(outdir, exist_ok=True)
    written = []
    for tb in tables:
        cols = ["truth \\ fit"] + tb["cols"]
        text = [[t] + [fmt(c) for c in cells] for t, cells in tb["rows"]]
        widths = [1.1] + [1.9 if c.startswith("envelope") else 1.5 for c in tb["cols"]]
        fig = plt.figure(figsize=(max(sum(widths), 7.0), 0.4 * (len(text) + 1) + 0.75), dpi=200)
        ax = fig.add_axes([0.0, 0.12, 1.0, 0.68])
        ax.axis("off")
        table = ax.table(cellText=text, colLabels=cols, cellLoc="center", bbox=[0, 0, 1, 1],
                         colWidths=[w / sum(widths) for w in widths])
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        for (row, col), c in table.get_celld().items():
            c.set_edgecolor("#cbd5e0")
            if row == 0 or col == 0:
                c.set_text_props(weight="bold", color="#1a365d")
                c.set_facecolor("#f7fafc")
                continue
            truth, cells = tb["rows"][row - 1]
            cellv = cells[col - 1]
            diag = tb["cols"][col - 1] == f"fit {truth}"
            if cellv == "—":
                c.set_facecolor("#edf2f7")
            elif isinstance(cellv, tuple) and cellv[3]:
                c.set_text_props(weight="bold", color="#9b2c2c")
                c.set_facecolor("#fed7d7")
            elif diag:
                c.set_facecolor("#ebf8ff")
            if diag and isinstance(cellv, tuple):
                c.set_edgecolor("#2b6cb0")
                c.set_linewidth(1.5)
        fig.text(0.5, 0.90, f"{tb['cat']} (best fit: {best.get(tb['cat'], '?')})  |  r_inj = {tb['r']:g}  |  {label}"
                 + (f"  |  {extra}" if extra else ""),
                 ha="center", va="center",
                 fontsize=11, weight="bold", color="#2b6cb0")
        fig.text(0.5, 0.04, f"median_pull_fixed ± error (bad toys at ∓fixed pull);  |median| ≥ {thr:g} highlighted;  "
                            "blue border: fit = truth", ha="center",
                 fontsize=8, color="#4a5568")
        fn = os.path.join(outdir, f"median_table_{tb['cat']}_r{tb['r']:g}.png")
        fig.savefig(fn, bbox_inches="tight")
        plt.close(fig)
        written.append(fn)
    return written


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fullenv-dir", required=True)
    p.add_argument("--dropfunc-dir", required=True)
    p.add_argument("-o", "--output", default="median_tables.md", help="Markdown output (default median_tables.md)")
    p.add_argument("--png-dir", default="median_tables_png", help="Directory for the PNG tables")
    p.add_argument("--bestfit-cache", default=None,
                   help="bestfit_pdfindices.sh with each category's best-fit function "
                        "(default: <dropfunc-dir>/bestfit_pdfindices.sh)")
    p.add_argument("--label", default="", help="Extra text for every table title, e.g. 'rMin floor -0.50'")
    p.add_argument("--threshold", type=float, default=0.2, help="|median| >= this is highlighted (default 0.2)")
    args = p.parse_args()

    env = load(os.path.join(args.fullenv_dir, "bias_plots_full_*_envelope*/pull_summary.json"),
               r"bias_plots_full_([^_/]+)_envelope")
    wrong = load(os.path.join(args.fullenv_dir, "bias_plots_full_*_wrongfit*/pull_summary.json"),
                 r"bias_plots_full_([^_/]+)_wrongfit")
    same = load(os.path.join(args.fullenv_dir, "bias_plots_full_*_samefit*/pull_summary.json"),
                r"bias_plots_full_([^_/]+)_samefit")
    # the fit = truth diagonal goes into the same "fit k" columns as wrongfit
    wrong.update({k: v for k, v in same.items() if k[2] == k[3]})
    drop = load(os.path.join(args.dropfunc_dir, "bias_plots_dropfunc_*/pull_summary.json"),
                r"bias_plots_dropfunc_([^_/]+)")
    if not (env or wrong or drop):
        raise SystemExit("No pull_summary.json found under the given directories")

    tables = build_tables(env, wrong, drop, args.threshold)
    best = load_bestfit(args.bestfit_cache or os.path.join(args.dropfunc_dir, "bestfit_pdfindices.sh"))
    missing = sorted({tb["cat"] for tb in tables} - set(best))
    if missing:
        print(f"[WARN] no best-fit function found for {missing}; labelled '?'")
    pulls = {c[2] for tb in tables for _, cells in tb["rows"] for c in cells if isinstance(c, tuple)}
    if len(pulls) > 1:
        raise SystemExit(f"[ERROR] runs used different --fixed-pull values {sorted(pulls)}; not comparable in one table")

    write_markdown(tables, args.output, os.path.abspath(args.fullenv_dir), os.path.abspath(args.dropfunc_dir),
                   args.threshold, pulls, best, args.label)
    pngs = write_pngs(tables, args.png_dir, args.threshold, f"median_pull_fixed (P = {pulls.pop():g})", best, args.label)
    print(f"Wrote {args.output} and {len(pngs)} PNGs in {args.png_dir}/")


if __name__ == "__main__":
    main()
