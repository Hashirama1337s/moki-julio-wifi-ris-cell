"""Bias-network screen (2026-10-04): boss design with 4 corner bias posts whose top gap (ports 7-10) holds the bias choke.
Ports: 2-5 diodes (NE, SW, NW, SE), 6 centre via (shorted), 7-10 bias chokes (NE, SW, NW, SE corner pads).
Choke models (Palace exp(+jwt)): open (no bias line), ideal 47 nH, 1 kohm || 0.05 pF (resistive feed), 22 nH chip with SRF 4 GHz and Q 30,
short (worst case). Worst case over the MADP-000907 box x D1/D2 +-5 % x Lc +-0.1 nH, steerable loss (TE in).
usage: python bias_an.py <mp.json> [Lc_nH]"""
import sys
import numpy as np
from mp_pol_analyze import load, R0

ROT = [1, 0, 3, 2, 4, 6, 5, 8, 7]
D = dict(Ron=5.2, C=0.025e-12, L=0.3e-9, Roff=60.0, Cb=(0.022e-12, 0.030e-12), Lb=(0.2e-9, 0.4e-9)); Q = 30.0


def prep(path):
    d, f, S = load(path)
    if d.get("theta", 0) != 0:
        k = S[:, 0, 1:][:, ROT] / S[:, 2:, 0]; wt = abs(S[:, 2:, 0]) ** 2; km = (k * wt).sum(1) / wt.sum(1)
        S = S.copy(); S[:, 2:, 0] = S[:, 2:, 0] * (-km)[:, None]
    return d, f, S


def choke(f, kind):
    w = 2 * np.pi * f
    if kind == "open":
        return np.full(len(f), 1e12 + 0j)
    if kind == "short":
        return np.full(len(f), 1e-3 + 0j)
    if kind == "L47":
        return 1j * w * 47e-9 + 1.0
    if kind == "R1k":
        return 1 / (1 / 1000.0 + 1j * w * 0.05e-12)
    if kind == "L22srf4":
        L = 22e-9; cp = 1 / ((2 * np.pi * 4e9) ** 2 * L); zl = 1j * w * L + w * L / Q
        return 1 / (1 / zl + 1j * w * cp)
    raise ValueError(kind)


def te_in(S, i, Z):
    g = (Z - R0) / (Z + R0); Gm = np.diag(g)
    b = np.linalg.solve(np.eye(len(Z)) - S[i, 2:, 1:] @ Gm, S[i, 2:, 0]); return S[i, :2, 0] + S[i, :2, 1:] @ (Gm @ b)


def worst(f, S, Lc, kind):
    w = 2 * np.pi * f; out = []; zc = choke(f, kind)
    for C in D["Cb"] + (D["C"],):
        for L in D["Lb"] + (D["L"],):
            for dl in (-0.1e-9, 0.0, 0.1e-9):
                for mm in (0.95, 1.0, 1.05):
                    lc = Lc + dl; zl = 1j * w * lc + w * lc / Q
                    on = D["Ron"] + 1j * w * L + zl; o1 = D["Roff"] + 1 / (1j * w * C) + 1j * w * L + zl; o2 = D["Roff"] + 1 / (1j * w * C * mm) + 1j * w * L + zl
                    RA, RB = [], []
                    for i in range(len(f)):
                        zb = [zc[i]] * 4
                        RA.append(te_in(S, i, np.array([on[i], on[i], o2[i], o2[i], 1e-3] + zb)))
                        RB.append(te_in(S, i, np.array([o1[i], o1[i], on[i], on[i], 1e-3] + zb)))
                    RA, RB = np.array(RA), np.array(RB)
                    out.append(-20 * np.log10(np.linalg.norm((RA - RB) / 2, axis=1)))
    return np.max(out, axis=0)


if __name__ == "__main__":
    d, f, S = prep(sys.argv[1]); Lc = float(sys.argv[2]) * 1e-9 if len(sys.argv) > 2 else 2.5e-9
    print(sys.argv[1], "theta", d.get("theta"), "f GHz", np.round(f / 1e9, 3), f"Lc {Lc * 1e9:.2f} nH")
    for kind in ("open", "L47", "R1k", "L22srf4", "short"):
        print(f"  choke {kind:8s}: {np.round(worst(f, S, Lc, kind), 2)}")
