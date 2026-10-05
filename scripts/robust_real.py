"""Realism screen for the best design (s2 pattern + series chip inductor), free algebra on the Palace sheet (2026-10-03, after
independent lever review): (1) spacer tolerance G +-0.2 mm (foam sheet); (2) mounting parasitic C_par across each diode junction;
(3) chip inductor with Q(f) and self-resonance (parallel C); (4) loss split: absorption vs co-pol leakage.
Worst case = diode box (C, L) x D1/D2 +-5 % x Lc +-0.1 nH x G +-0.2 mm."""
import numpy as np
from sheet_scan import sheet, r_of
import os
DATA = os.environ.get("RIS_DATA", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"))   # network JSON files

f, yo, zc, n = sheet(f"{DATA}/explore/omp3_2p000_0p400_1p275_0p644_0p300_10p000.json", 10.0); w = 2 * np.pi * f
D = dict(Ron=5.2, C=0.025e-12, L=0.3e-9, Roff=60.0, Cb=(0.022e-12, 0.030e-12), Lb=(0.2e-9, 0.4e-9))


def zl_chip(Lc, Q0=30.0, srf=None):
    """chip inductor: Q rising ~sqrt(f) from Q0 at 2.4 GHz (typical 0201/0402 curves); optional SRF via parallel C."""
    q = Q0 * np.sqrt(f / 2.4e9); z = 1j * w * Lc + w * Lc / q
    if srf:
        cp = 1 / ((2 * np.pi * srf) ** 2 * Lc); z = 1 / (1 / z + 1j * w * cp)
    return z


def states(G, Lc, C, L, mm, cpar, chip):
    zl = chip(Lc)
    def dio(z):
        return 1 / (1 / z + 1j * w * cpar) if cpar else z
    zon = dio(D["Ron"] + 1j * w * L) + zl; zoff = dio(D["Roff"] + 1 / (1j * w * C * mm) + 1j * w * L) + zl
    return r_of(f, G, yo, zc, n, zon), r_of(f, G, yo, zc, n, zoff)


def worst(G, Lc, cpar=0.0, chip=zl_chip, dG=0.2):
    out = []
    for g in (G - dG, G, G + dG):
        for C in D["Cb"] + (D["C"],):
            for L in D["Lb"] + (D["L"],):
                for dl in (-0.1e-9, 0.0, 0.1e-9):
                    for mm in (0.95, 1.0, 1.05):
                        a, b = states(g, Lc + dl, C, L, mm, cpar, chip)
                        out.append(-20 * np.log10(abs(a - b) / 2))
    return np.max(out, axis=0)


def best(cpar=0.0, chip=zl_chip, dG=0.2):
    return min(((worst(G, Lc * 1e-9, cpar, chip, dG), G, Lc) for G in np.arange(9.0, 13.01, 0.25)
                for Lc in np.arange(1.0, 3.01, 0.25)), key=lambda x: x[0].max())

if __name__ == "__main__":
    print("f GHz", np.round(f / 1e9, 3))
    for label, kw in (("ideal (Q 30 flat-ish), no G tol", dict(dG=0.0)), ("+ G +-0.2 mm", {}),
                      ("+ G tol + C_par 0.02 pF", dict(cpar=0.02e-12)), ("+ G tol + C_par 0.05 pF", dict(cpar=0.05e-12)),
                      ("+ G tol + C_par 0.05 + Q 20 + SRF 9 GHz", dict(cpar=0.05e-12, chip=lambda Lc: zl_chip(Lc, 20.0, 9e9)))):
        b = best(**kw)
        print(f"{label:42s}: worst {b[0].max():.2f} dB  G {b[1]:.2f}  Lc {b[2]:.2f} nH  per-f {np.round(b[0], 2)}")
    # loss split at the nominal point of the plain optimum
    b = best(); a, bb = states(b[1], b[2] * 1e-9, D["C"], D["L"], 1.0, 0.0, zl_chip)
    absorb = 1 - (abs(a) ** 2 + abs(bb) ** 2) / 2; copol = abs(a + bb) ** 2 / 4; cross = abs(a - bb) ** 2 / 4
    print("nominal split at G %.2f Lc %.2f: cross %s  absorbed %s  co-pol %s" % (b[1], b[2], np.round(cross, 3), np.round(absorb, 3), np.round(copol, 3)))
