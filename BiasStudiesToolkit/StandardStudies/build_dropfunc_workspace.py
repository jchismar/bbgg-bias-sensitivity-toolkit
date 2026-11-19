#!/usr/bin/env python3
"""Rebuild a per-category RooMultiPdf workspace with ONE envelope function
permanently removed, for the "generate with truth i, fit with i genuinely
absent from the envelope" bias check.

This is a stronger test than freezing pdfindex away from i (the existing
--fit-mode wrongfit in run_full_envelope_bias.sh): there, the true function
is still in the RooMultiPdf's candidate list, just forced off; here, the
fit's discrete profiling can never recover the truth at all, because it was
never a candidate. That isolates the bias coming from the envelope itself
being an incomplete truncation of function space, as opposed to the
discrete-profiling algorithm merely landing on the wrong choice.

Input is the RAW per-category multipdf file a datacard's `shapes bkg_mass`
line points to (e.g. .../background/CMS-HGG_multipdf_cat2.root), NOT the
combined text2workspace.py output. Output keeps the same RooWorkspace name
("multipdf") and the same names for every object involved -- only the
number of candidate pdfs in one leg's RooMultiPdf shrinks by one -- so a
datacard needs to change nothing but the file path on its `bkg_mass` shapes
line for this category to pick this up.

1D categories (one leg): <top_pdf_name> IS the RooMultiPdf's own name (the
part after "multipdf:" on the datacard's shapes line) -- this is the
originally-tested path, unchanged.

2D categories (hgg+hbb legs): <top_pdf_name> is whatever object combines
both legs' RooMultiPdfs (e.g. a RooProdPdf) -- pass --leg hgg|hbb to say
which leg's envelope to drop a function from; the OTHER leg is carried
through completely untouched. This relies on RooFit::RecycleConflictNodes()
substituting-by-name: the reduced leg-multipdf is imported first under its
original name, then the original top-level object is imported on top of
it, so anywhere it references that name it picks up the reduced version
instead of re-cloning the original. NOTE: this 2D path has not been
exercised against a real two-leg workspace (this toolkit's current
datacards are all 1D/hgg-only) -- verify against a real file before
trusting its output.

Usage:
  build_dropfunc_workspace.py <in_ws.root> <cat_number> <top_pdf_name> <drop_tag> <out_ws.root> [--leg hgg|hbb]

<cat_number>     e.g. 2 (matches the "catN" token in the RooMultiPdf/pdfindex names)
<top_pdf_name>   exact object name after "multipdf:" on the datacard's shapes line
<drop_tag>       one of exp1, bern1, bern2, bern3, pow1, lau1, lau2 -- same tag
                 vocabulary as resolve_pdfindex_legs.py ("bern1" also matches a
                 bare "flat0" pdf, same aliasing as that script)
<--leg hgg|hbb>  only needed when more than one RooMultiPdf matches <cat_number>
                 (i.e. a 2D category) -- picks which leg's envelope to modify
"""
import argparse
import sys

import ROOT

TAG_MATCHERS = [
    ("bern1", lambda n: "bern1" in n or "flat0" in n),
    ("bern2", lambda n: "bern2" in n),
    ("bern3", lambda n: "bern3" in n),
    ("exp1", lambda n: "exp1" in n),
    ("pow1", lambda n: "pow1" in n),
    ("lau1", lambda n: "lau1" in n),
    ("lau2", lambda n: "lau2" in n),
]


def tag_of(pdf_name):
    name = pdf_name.lower()
    for tag, match in TAG_MATCHERS:
        if match(name):
            return tag
    return "unknown"


def find_multipdf_names(w, cat_token):
    names = []
    it = w.allPdfs().createIterator()
    obj = it.Next()
    while obj:
        if obj.ClassName() == "RooMultiPdf" and cat_token.lower() in obj.GetName().lower():
            names.append(obj.GetName())
        obj = it.Next()
    return names


def find_pdfindex_names(w, cat_token):
    names = []
    cit = w.allCats().createIterator()
    cobj = cit.Next()
    while cobj:
        name = cobj.GetName()
        if "pdfindex" in name and cat_token in name:
            names.append(name)
        cobj = cit.Next()
    return names


def filter_by_leg(names, leg):
    if leg is None:
        return names
    is_hbb = leg == "hbb"
    return [n for n in names if ("hbb" in n.lower()) == is_hbb]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("in_ws")
    p.add_argument("cat_number")
    p.add_argument("top_pdf_name")
    p.add_argument("drop_tag")
    p.add_argument("out_ws")
    p.add_argument("--leg", choices=["hgg", "hbb"], default=None)
    args = p.parse_args()

    cat_token = f"cat{args.cat_number}"

    ROOT.gErrorIgnoreLevel = ROOT.kError

    f = ROOT.TFile.Open(args.in_ws)
    if not f or f.IsZombie():
        sys.exit(f"[ERROR] Could not open {args.in_ws}")
    w = f.Get("multipdf")
    if not w:
        sys.exit(f"[ERROR] No RooWorkspace 'multipdf' in {args.in_ws}")

    mpdf_candidates = filter_by_leg(find_multipdf_names(w, cat_token), args.leg)
    if len(mpdf_candidates) == 0:
        sys.exit(
            f"[ERROR] No RooMultiPdf tagged '{cat_token}'"
            + (f" leg={args.leg}" if args.leg else "")
            + f" found in {args.in_ws}"
        )
    if len(mpdf_candidates) > 1:
        sys.exit(
            f"[ERROR] {len(mpdf_candidates)} RooMultiPdf objects match '{cat_token}': "
            f"{mpdf_candidates}. Disambiguate with --leg hgg|hbb."
        )
    mpdf_name = mpdf_candidates[0]
    mpdf_old = w.pdf(mpdf_name)
    n_pdfs = mpdf_old.getNumPdfs()
    all_names = [mpdf_old.getPdf(i).GetName() for i in range(n_pdfs)]

    kept = ROOT.RooArgList()
    kept_names = []
    dropped_name = None
    for i in range(n_pdfs):
        pdf_i = mpdf_old.getPdf(i)
        if tag_of(pdf_i.GetName()) == args.drop_tag:
            dropped_name = pdf_i.GetName()
            continue
        kept.add(pdf_i)
        kept_names.append(pdf_i.GetName())

    if dropped_name is None:
        sys.exit(
            f"[ERROR] Tag '{args.drop_tag}' not found among {mpdf_name}'s {n_pdfs} "
            f"functions ({all_names})"
        )
    if kept.getSize() < 2:
        sys.exit(
            f"[ERROR] Dropping '{args.drop_tag}' from {mpdf_name} would leave only "
            f"{kept.getSize()} function(s) left -- not a meaningful envelope anymore."
        )

    cat_candidates = filter_by_leg(find_pdfindex_names(w, cat_token), args.leg)
    if len(cat_candidates) != 1:
        sys.exit(
            f"[ERROR] Expected exactly 1 pdfindex category for '{cat_token}'"
            + (f" leg={args.leg}" if args.leg else "")
            + f", found {cat_candidates}"
        )
    old_cat_name = cat_candidates[0]

    # Normalization can live on the top-level object (1D: top_pdf_name IS
    # mpdf_name) or per-leg (2D: a shared top-level norm may not exist).
    norm_old = w.arg(f"{args.top_pdf_name}_norm")
    if norm_old is None and args.top_pdf_name != mpdf_name:
        norm_old = w.arg(f"{mpdf_name}_norm")
    if norm_old is None:
        sys.exit(
            f"[ERROR] No normalization found as '{args.top_pdf_name}_norm' or "
            f"'{mpdf_name}_norm' in {args.in_ws}"
        )

    # RooMultiPdf's constructor calls index.defineType("_pdfN", N) itself for
    # every entry in the pdf list it is given, so a freshly-constructed,
    # state-less RooCategory with the same name is all that's needed here.
    new_cat = ROOT.RooCategory(old_cat_name, old_cat_name)
    mpdf_new = ROOT.RooMultiPdf(mpdf_name, mpdf_old.GetTitle(), new_cat, kept)

    wout = ROOT.RooWorkspace("multipdf", "multipdf")
    getattr(wout, "import")(mpdf_new, ROOT.RooFit.RecycleConflictNodes())
    getattr(wout, "import")(norm_old)

    if args.top_pdf_name != mpdf_name:
        # 2D case: re-import the wrapper that combines both legs. Because
        # the reduced leg-multipdf was already imported above under its
        # original name, RecycleConflictNodes() splices it in wherever the
        # wrapper's original (full-envelope) copy of that name would
        # otherwise have been cloned -- the other leg is untouched.
        top_old = w.pdf(args.top_pdf_name)
        if top_old is None:
            sys.exit(f"[ERROR] Top-level pdf '{args.top_pdf_name}' not found in {args.in_ws}")
        getattr(wout, "import")(top_old, ROOT.RooFit.RecycleConflictNodes())

    fout = ROOT.TFile.Open(args.out_ws, "RECREATE")
    wout.Write()
    fout.Close()

    print(f"[OK] Wrote {args.out_ws}")
    print(f"     Modified leg-multipdf: {mpdf_name}" + (f" (leg={args.leg})" if args.leg else ""))
    print(f"     Dropped: {dropped_name}  (tag={args.drop_tag})")
    print(f"     Kept ({kept.getSize()}/{n_pdfs}): {kept_names}")
    if args.top_pdf_name != mpdf_name:
        print(f"     Top-level pdf re-imported, other leg unchanged: {args.top_pdf_name}")


if __name__ == "__main__":
    main()
