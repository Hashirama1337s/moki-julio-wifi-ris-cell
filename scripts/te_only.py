"""TE-in steerable loss at oblique incidence from the TE multiport alone (no TM column needed). usage: python te_only.py <mp.json> [Lc...]"""
import sys
import numpy as np


def wmean(k, a):
    """reciprocity factor averaged with weights |S|^2: ports not excited by symmetry (e.g. phi 45) carry ~0 and garbage ratios"""
    wt = abs(a) ** 2; return (k * wt).sum(1) / wt.sum(1)
from mp_pol_analyze import load
import angle_an as a

d, f, S = load(sys.argv[1]); kte = S[:, 0, 1:][:, a.ROT] / S[:, 2:, 0]
S6 = np.zeros((len(f), 6, 6), complex); S6[:, :2, 0] = S[:, :2, 0]; S6[:, :, 2:] = S[:, :, 1:]; S6[:, 2:, 0] = S[:, 2:, 0] * (-wmean(kte, S[:, 2:, 0]))[:, None]
S6[:, :2, 1] = np.nan
print(sys.argv[1], "theta", d.get("theta"), "f", np.round(f / 1e9, 3) if len(f) < 12 else f"{len(f)} pts {f[0]/1e9:.2f}-{f[-1]/1e9:.2f}")
D, Q = a.D, a.Q; w = 2 * np.pi * f
for Lc in [float(x) * 1e-9 for x in sys.argv[2:]] or [2.0e-9, 2.25e-9, 3.0e-9]:
    out, pe = [], []
    for C in D["Cb"] + (D["C"],):
        for L in D["Lb"] + (D["L"],):
            for dl in (-0.1e-9, 0.0, 0.1e-9):
                for mm in (0.95, 1.0, 1.05):
                    lc = Lc + dl; zl = 1j * w * lc + w * lc / Q
                    on = D["Ron"] + 1j * w * L + zl; o1 = D["Roff"] + 1 / (1j * w * C) + 1j * w * L + zl; o2 = D["Roff"] + 1 / (1j * w * C * mm) + 1j * w * L + zl
                    gA = [(z - 50) / (z + 50) for z in (on, on, o2, o2)]; gB = [(z - 50) / (z + 50) for z in (o1, o1, on, on)]
                    RA, RB = [], []
                    for gl, R in ((gA, RA), (gB, RB)):
                        for i in range(len(f)):
                            Gm = np.diag([g[i] for g in gl]); ap = Gm @ np.linalg.solve(np.eye(4) - S6[i, 2:, 2:] @ Gm, S6[i, 2:, 0])
                            R.append(S6[i, :2, 0] + S6[i, :2, 2:] @ ap)
                    RA, RB = np.array(RA), np.array(RB)
                    out.append(-20 * np.log10(np.linalg.norm((RA - RB) / 2, axis=1)))
                    pe.append(180 - abs(np.degrees(np.angle(RA[:, 1] / RB[:, 1]))))
    wv = np.max(out, axis=0)
    if len(f) <= 12:
        print(f"  Lc {Lc*1e9:.2f}: steerable loss {np.round(wv, 2)}  phase err {np.round(np.max(pe, axis=0), 2)}")
    else:
        m1 = f <= 2.49e9; m2 = (f >= 5.149e9) & (f <= 7.126e9); band = m2 if m2.any() else ~m1
        print(f"  Lc {Lc*1e9:.2f}: worst in 5.15-7.125 {wv[band].max():.2f} dB at {f[band][np.argmax(wv[band])]/1e9:.2f} GHz; curve:",
              " ".join(f"{x/1e9:.2f}:{v:.1f}" for x, v in zip(f[::2], wv[::2])))
