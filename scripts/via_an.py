"""Analysis with the centre-via port (port 6). S (nf, 2+NP, 1+NP): rows TE,TM,b2..b6; cols TE_in,a2..a6 (TE incidence).
Ports 2-5 diodes (NE,SW,NW,SE), port 6 = via. Via load: 'short' (metal via), 'open' (no via), or an inductance in nH (choke).
At oblique incidence the port<-Floquet entries are renormalised by reciprocity + C2 as in angle_an (weights |S|^2; the via maps
onto itself under C2). usage: python via_an.py <mp.json> [Lc_nH ...]"""
import sys
import numpy as np
from mp_pol_analyze import load, R0

ROT = [1, 0, 3, 2, 4]
D = dict(Ron=5.2, C=0.025e-12, L=0.3e-9, Roff=60.0, Cb=(0.022e-12, 0.030e-12), Lb=(0.2e-9, 0.4e-9)); Q = 30.0


def prep(path):
    d, f, S = load(path)
    if d.get("theta", 0) != 0:
        k = S[:, 0, 1:][:, ROT] / S[:, 2:, 0]; wt = abs(S[:, 2:, 0]) ** 2; km = (k * wt).sum(1) / wt.sum(1)
        S = S.copy(); S[:, 2:, 0] = S[:, 2:, 0] * (-km)[:, None]
    return d, f, S


def te_in(S, i, Z):
    g = (Z - R0) / (Z + R0); Gm = np.diag(g)
    b = np.linalg.solve(np.eye(len(Z)) - S[i, 2:, 1:] @ Gm, S[i, 2:, 0]); return S[i, :2, 0] + S[i, :2, 1:] @ (Gm @ b)


def worst(f, S, Lc, via):
    w = 2 * np.pi * f; out = []
    for C in D["Cb"] + (D["C"],):
        for L in D["Lb"] + (D["L"],):
            for dl in (-0.1e-9, 0.0, 0.1e-9):
                for mm in (0.95, 1.0, 1.05):
                    lc = Lc + dl; zl = 1j * w * lc + w * lc / Q
                    on = D["Ron"] + 1j * w * L + zl; o1 = D["Roff"] + 1 / (1j * w * C) + 1j * w * L + zl; o2 = D["Roff"] + 1 / (1j * w * C * mm) + 1j * w * L + zl
                    zv = np.zeros(len(f)) + 1e-3 if via == "short" else (np.full(len(f), 1e12) if via == "open" else 1j * w * float(via) * 1e-9 + 0.2)
                    RA = np.array([te_in(S, i, np.array([on[i], on[i], o2[i], o2[i], zv[i]])) for i in range(len(f))])
                    RB = np.array([te_in(S, i, np.array([o1[i], o1[i], on[i], on[i], zv[i]])) for i in range(len(f))])
                    out.append(-20 * np.log10(np.linalg.norm((RA - RB) / 2, axis=1)))
    return np.max(out, axis=0)


if __name__ == "__main__":
    d, f, S = prep(sys.argv[1]); band = (f >= 5.149e9) & (f <= 7.126e9)
    print(sys.argv[1], "theta", d.get("theta"), f"{len(f)} freqs")
    for Lc in [float(x) for x in sys.argv[2:]] or [2.0, 2.5, 3.0]:
        for via in ("open", "short", "2", "5"):
            wv = worst(f, S, Lc * 1e-9, via)
            if len(f) <= 12:
                print(f"  Lc {Lc:.2f} via {via:5s}: {np.round(wv, 2)}")
            else:
                print(f"  Lc {Lc:.2f} via {via:5s}: worst 5.15-7.125 {wv[band].max():.2f} dB @ {f[band][np.argmax(wv[band])]/1e9:.2f};", " ".join(f"{x/1e9:.1f}:{v:.1f}" for x, v in zip(f[::4], wv[::4])))
