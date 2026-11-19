#!/usr/bin/env python3
"""Category-scoped RooMultiPdf function<->index resolver, for the combined
(multi-category) channel-masked workspace built with
`text2workspace.py --channel-masks`.

parser_legs_matrix.py was written for a single-category workspace (one hgg
and/or hbb RooMultiPdf) and associates a function tag to a RooMultiPdf with a
loose "is this an hgg-like name" heuristic. On a combined workspace, where
several categories' RooMultiPdfs all exist at once, that heuristic can
silently pick up the wrong category's function ordering. This script
restricts the association to a single requested category by requiring both
the pdfindex name and the RooMultiPdf name to contain that category's token
(e.g. "cat2"), which removes the ambiguity.

Usage: resolve_pdfindex_legs.py <workspace.root> <cat_number>

Prints bash-sourceable assignments (same IDX_<var>_<tag>, HGG_VAR, HBB_VAR
convention as parser_legs_matrix.py), plus ALL_PDFINDEX_VARS listing every
pdfindex-like category in the whole workspace, so the caller can freeze every
*other* category's discrete index (masked channels don't need it floating).
"""

import sys

import ROOT


def main():
    if len(sys.argv) != 3:
        print("[ERROR] Usage: resolve_pdfindex_legs.py <workspace.root> <cat_number>", file=sys.stderr)
        sys.exit(1)

    ws_file, cat_number = sys.argv[1], sys.argv[2]
    token = f"cat{cat_number}"

    f = ROOT.TFile.Open(ws_file)
    if not f or f.IsZombie():
        print(f"[ERROR] Could not open {ws_file}", file=sys.stderr)
        sys.exit(1)
    w = f.Get("w")
    if not w:
        print(f"[ERROR] No RooWorkspace 'w' in {ws_file}", file=sys.stderr)
        sys.exit(1)

    # 1. Function tag -> index, per RooMultiPdf (name kept, not yet filtered by category)
    pdf_map = {}
    it = w.allPdfs().createIterator()
    obj = it.Next()
    while obj:
        if obj.ClassName() == "RooMultiPdf":
            for i in range(obj.getNumPdfs()):
                pdf_name = obj.getPdf(i).GetName().lower()
                tag = "unknown"
                if "exp1" in pdf_name:
                    tag = "exp1"
                elif "bern1" in pdf_name or "flat0" in pdf_name:
                    tag = "bern1"
                elif "bern2" in pdf_name:
                    tag = "bern2"
                elif "bern3" in pdf_name:
                    tag = "bern3"
                elif "pow1" in pdf_name:
                    tag = "pow1"
                elif "lau1" in pdf_name:
                    tag = "lau1"
                elif "lau2" in pdf_name:
                    tag = "lau2"
                pdf_map[(obj.GetName(), tag)] = i
        obj = it.Next()

    # 2. Every pdfindex-like category in the whole workspace
    all_pdfindex = []
    cit = w.allCats().createIterator()
    cobj = cit.Next()
    while cobj:
        name = cobj.GetName()
        if "pdfindex" in name:
            all_pdfindex.append(name)
        cobj = cit.Next()

    # 3. Restrict to the requested category, split into hgg/hbb legs if both exist
    matched = [p for p in all_pdfindex if token in p]
    multi_leg = len({("hbb" in p.lower()) for p in matched}) > 1

    hgg_var, hbb_var = "", ""
    lines = []
    for p_name in matched:
        is_hbb = "hbb" in p_name.lower()
        if is_hbb:
            hbb_var = p_name
        else:
            hgg_var = p_name
        leg_kw = "hbb" if is_hbb else "hgg"

        for (mpdf_name, tag), idx in pdf_map.items():
            if tag == "unknown":
                continue
            mpdf_lower = mpdf_name.lower()
            if token.lower() not in mpdf_lower:
                continue
            if multi_leg and leg_kw not in mpdf_lower:
                continue
            lines.append(f"IDX_{p_name}_{tag}='{p_name}={idx}'")

    for line in lines:
        print(line)
    print(f"HGG_VAR='{hgg_var}'")
    print(f"HBB_VAR='{hbb_var}'")
    print(f"ALL_PDFINDEX_VARS='{' '.join(all_pdfindex)}'")


if __name__ == "__main__":
    main()
