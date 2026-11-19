#!/usr/bin/env python3
"""Least negative r at which the s+b density of some category first reaches
zero (nominal parameters), i.e. max over categories and observable values of
-b(x)/s(x). Port of submitBiasTests.py's zeroDensityR(), used for
--r-min-floor auto: the floor is FACTOR x that r (default 0.9, as there).

Usage: zero_density_r.py <combined_workspace.root> [--sig-proc xhh[,name2,...]] [--factor 0.9]

Per-category values go to stderr; stdout is one line, the floor, so a shell
can do FLOOR=$(python3 zero_density_r.py ws.root).
"""
import argparse
import sys

import numpy as np
import ROOT


def main():
    p = argparse.ArgumentParser()
    p.add_argument("ws")
    p.add_argument("--sig-proc", default="xhh",
                   help="Signal process name(s) or substring(s), comma-separated; a coefficient whose name "
                        "contains any of them counts as signal (default xhh)")
    p.add_argument("--factor", type=float, default=0.9, help="Floor = factor x zero-density r (default 0.9)")
    args = p.parse_args()

    ROOT.gErrorIgnoreLevel = ROOT.kError
    f = ROOT.TFile.Open(args.ws)
    w = f.Get("w") if f and not f.IsZombie() else None
    if not w:
        sys.exit(f"[ERROR] no RooWorkspace 'w' in {args.ws}")
    w.var("r").setVal(1.0)
    sig_keys = [x for x in args.sig_proc.split(",") if x]

    best, where = -1e30, None
    for t in w.cat("CMS_channel"):
        c = str(t.first)
        top = w.pdf(f"pdf_bin{c}_obsOnly") or w.pdf(f"pdf_bin{c}")
        obs = [o for o in top.getObservables(w.data("data_obs").get()) if o.GetName() != "CMS_channel"][0]
        nodes = ROOT.RooArgList()
        top.branchNodeServerList(nodes)
        adds = [nodes[i] for i in range(nodes.getSize()) if nodes[i].InheritsFrom("RooAddPdf")]
        add = max(adds, key=lambda a: a.pdfList().getSize())
        xs = np.linspace(obs.getMin(), obs.getMax(), 600)
        nset = ROOT.RooArgSet(obs)
        sig, bkg = np.zeros_like(xs), np.zeros_like(xs)
        for i in range(add.pdfList().getSize()):
            pdf, co = add.pdfList()[i], add.coefList()[i]
            y = []
            for x in xs:
                obs.setVal(x)
                y.append(pdf.getVal(nset))
            y = np.array(y) * co.getVal()
            if any(k in co.GetName() for k in sig_keys):
                sig += y
            else:
                bkg += y
        m = sig > 1e-9 * sig.max()
        rz = float(np.max(-bkg[m] / sig[m]))
        print(f"  {c:<20s} s+b density reaches 0 at r = {rz:.4f}", file=sys.stderr)
        if rz > best:
            best, where = rz, c
    f.Close()
    floor = args.factor * best
    print(f"rMin floor: {args.factor:g} x {best:.4f} ({where}) = {floor:.4f}", file=sys.stderr)
    print(f"{floor:.4f}")


if __name__ == "__main__":
    main()
