"""Spacer-height transfer from ONE full-wave multiport: de-embed the diagonal-chain sheet as an exact per-frequency Moebius
map  y_sheet(Z) = y_open + eta0 / (zc + n*Z)  (3 complex numbers per frequency, from diode short / open / 10 ohm), then
re-attach a different air spacer G' analytically and score the MADP-000907 tolerance box (C, L corners, D1/D2 +-5 %).
Assumes the sheet does not depend on G (evanescent coupling to the ground ~exp(-2*pi*G/15) ~ 2 %): validate first.
usage: python sheet_scan.py <mp.json> <G0> [G1 G2 ...]          (validate: python sheet_scan.py A.json 10 9 vs B.json at 9)"""
import sys
import numpy as np
from mp_pol_analyze import load
from eig_extract import eig, y_down, ETA0

D = dict(Ron=5.2, C=0.025e-12, L=0.3e-9, Roff=60.0, Cb=(0.022e-12, 0.030e-12), Lb=(0.2e-9, 0.4e-9))


def sheet(path, G0):
    d, f, S = load(path); nf = len(f)
    SH, OP, R10 = np.zeros(nf, complex), np.full(nf, 1e12 + 0j), np.full(nf, 10.0 + 0j)
    a1, b1 = eig(S, f, SH, OP); a2, b2 = eig(S, f, SH, SH)
    pick = (lambda a, b: a) if np.max(abs(a1 - a2)) < np.max(abs(b1 - b2)) else (lambda a, b: b)
    yin = lambda r: (1 - r) / (1 + r)
    ys, yo, y10 = (yin(pick(*eig(S, f, z, OP))) - y_down(f, G0) for z in (SH, OP, R10))
    zc = ETA0 / (ys - yo); n = (ETA0 / (y10 - yo) - zc) / 10.0
    return f, yo, zc, n


def r_of(f, G, yo, zc, n, Z):
    y = yo + ETA0 / (zc + n * Z) + y_down(f, G)
    return (1 - y) / (1 + y)


def score(f, G, yo, zc, n, d=D):
    w = 2 * np.pi * f; worst = np.zeros(len(f)); nom = None
    for C in d["Cb"] + (d["C"],):
        for L in d["Lb"] + (d["L"],):
            for mm in (0.95, 1.0, 1.05):
                zon = d["Ron"] + 1j * w * L; zoff = d["Roff"] + 1 / (1j * w * C * mm) + 1j * w * L
                l = -20 * np.log10(abs(r_of(f, G, yo, zc, n, zon) - r_of(f, G, yo, zc, n, zoff)) / 2)
                worst = np.maximum(worst, l)
                if C == d["C"] and L == d["L"] and mm == 1.0:
                    nom = l
    return worst, nom


if __name__ == "__main__":
    path, G0 = sys.argv[1], float(sys.argv[2])
    f, yo, zc, n = sheet(path, G0)
    print(path, " f GHz", np.round(f / 1e9, 3))
    for G in [G0] + [float(x) for x in sys.argv[3:]]:
        worst, nom = score(f, G, yo, zc, n)
        print(f"  G {G:5.2f}: worst {worst.max():.2f} dB  per-f worst {np.round(worst, 2)}  nominal {np.round(nom, 2)}")
