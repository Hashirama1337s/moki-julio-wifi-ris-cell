"""Lever test (2026-10-03): chip inductor Lc in series with each diode (same gap, so exact in the multiport algebra).
Loss R = w*Lc/Q (Q 30, conservative); tolerance box = MADP-000907 (C, L corners, D1/D2 +-5 %) x Lc +-0.1 nH.
Spacer G optimised analytically (sheet_scan transfer). usage: python lever_chip.py <mp.json> [<mp.json> ...]"""
import sys
import numpy as np
from sheet_scan import sheet, r_of

D = dict(Ron=5.2, C=0.025e-12, L=0.3e-9, Roff=60.0, Cb=(0.022e-12, 0.030e-12), Lb=(0.2e-9, 0.4e-9))
GG = np.arange(8.0, 12.01, 0.25); Q = 30.0


def worst(f, G, yo, zc, n, Lc):
    w = 2 * np.pi * f; out = []
    for C in D["Cb"] + (D["C"],):
        for L in D["Lb"] + (D["L"],):
            for dl in ((-0.1e-9, 0.0, 0.1e-9) if Lc > 0 else (0.0,)):
                for mm in (0.95, 1.0, 1.05):
                    lc = Lc + dl; zl = 1j * w * lc + w * lc / Q
                    zon = D["Ron"] + 1j * w * L + zl; zoff = D["Roff"] + 1 / (1j * w * C * mm) + 1j * w * L + zl
                    out.append(-20 * np.log10(abs(r_of(f, G, yo, zc, n, zon) - r_of(f, G, yo, zc, n, zoff)) / 2))
    return np.max(out, axis=0)


if __name__ == "__main__":
  for path in sys.argv[1:]:
    f, yo, zc, n = sheet(path, 10.0)
    print(path.split("/")[-1])
    for Lc in (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0):
        best = min(((worst(f, G, yo, zc, n, Lc * 1e-9), G) for G in GG), key=lambda t: t[0].max())
        print(f"  Lc {Lc:3.1f} nH: best G {best[1]:5.2f}  worst {best[0].max():.2f} dB  per-f {np.round(best[0], 2)}")
