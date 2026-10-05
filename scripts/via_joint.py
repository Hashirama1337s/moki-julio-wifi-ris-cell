"""Joint spacer G x chip-inductor Lc scan for the s2 pattern WITH a solid centre via (2026-10-03).
Normal incidence: dense_s2_g12_n (no via port; the via is invisible at normal incidence, shown to 0.01 dB) via net6.regap.
45 deg (phi 0): dense TE multiport via_d45_g12 + TM-in direct run via_d45_g12_tm (all 5 ports on 50 ohm) -> full 7x7
network (reciprocity + C2 renormalisation, weights |S|^2), oblique shunt-stub spacer transfer (! oblique transfer not yet
validated against a real run at another G), via port shorted. Worst case over the MADP-000907 box x D1/D2 +-5 % x Lc +-0.1 nH.
usage: python via_joint.py"""
import json
import numpy as np
from mp_pol_analyze import load, R0
import net6, one_scan
from angle_an import y_stub
import os
DATA = os.environ.get("RIS_DATA", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"))   # network JSON files

ROT = [1, 0, 3, 2, 4]
D = dict(Ron=5.2, C=0.025e-12, L=0.3e-9, Roff=60.0, Cb=(0.022e-12, 0.030e-12), Lb=(0.2e-9, 0.4e-9)); Q = 30.0
C0 = 299792458.0


def net45():
    d, f, S = load(f"{DATA}/via_d45_g12.json"); t = json.load(open(f"{DATA}/via_d45_g12_tm.json")); c = lambda x: np.array([complex(*z) for z in x])
    Vtm = np.stack([c(v) for v in t["V"]], 1) / np.sqrt(2 * R0); n = 2 + 5
    S7 = np.zeros((len(f), n, n), complex); S7[:, :2, 0] = S[:, :2, 0]; S7[:, :, 2:] = S[:, :, 1:]
    S7[:, 0, 1] = c(t["TE"]); S7[:, 1, 1] = c(t["TM"])
    wm = lambda k, a: (k * abs(a) ** 2).sum(1) / (abs(a) ** 2).sum(1)
    kte = S[:, 0, 1:][:, ROT] / S[:, 2:, 0]; ktm = S[:, 1, 1:][:, ROT] / Vtm
    S7[:, 2:, 0] = S[:, 2:, 0] * (-wm(kte, S[:, 2:, 0]))[:, None]; S7[:, 2:, 1] = Vtm * (-wm(ktm, Vtm))[:, None]
    return f, S7


def regap(S, f, th, G0, G1):
    n = S.shape[1]; k = 2 * np.pi * f / C0 * 1e-3 * np.cos(np.radians(th)); out = np.zeros_like(S); I = np.eye(n)
    for i in range(len(f)):
        Dm = np.diag([np.exp(1j * k[i] * 20.0)] * 2 + [1] * (n - 2)); s = Dm @ S[i] @ Dm
        y = (I - s) @ np.linalg.inv(I + s)
        for j, pol in enumerate(("TE", "TM")):
            y[j, j] += y_stub(f[i:i + 1], G1, th, pol)[0] - y_stub(f[i:i + 1], G0, th, pol)[0]
        out[i] = np.linalg.inv(Dm) @ ((I - y) @ np.linalg.inv(I + y)) @ np.linalg.inv(Dm)
    return out


def worst45(f, S, Lc):
    w = 2 * np.pi * f; out = []; zv = np.full(len(f), 1e-3 + 0j)
    for C in D["Cb"] + (D["C"],):
        for L in D["Lb"] + (D["L"],):
            for dl in (-0.1e-9, 0.0, 0.1e-9):
                for mm in (0.95, 1.0, 1.05):
                    lc = Lc + dl; zl = 1j * w * lc + w * lc / Q
                    on = D["Ron"] + 1j * w * L + zl; o1 = D["Roff"] + 1 / (1j * w * C) + 1j * w * L + zl; o2 = D["Roff"] + 1 / (1j * w * C * mm) + 1j * w * L + zl
                    RR = []
                    for Z in (np.stack([on, on, o2, o2, zv], 1), np.stack([o1, o1, on, on, zv], 1)):
                        g = (Z - R0) / (Z + R0); R = np.zeros((len(f), 2, 2), complex)
                        for i in range(len(f)):
                            Gm = np.diag(g[i]); R[i] = S[i, :2, :2] + S[i, :2, 2:] @ Gm @ np.linalg.solve(np.eye(5) - S[i, 2:, 2:] @ Gm, S[i, 2:, :2])
                        RR.append(R)
                    sw = (RR[0] - RR[1]) / 2
                    out.append(-20 * np.log10(np.maximum(np.linalg.norm(sw[:, :, 0], axis=1), 1e-9)))
                    out.append(-20 * np.log10(np.maximum(np.linalg.norm(sw[:, :, 1], axis=1), 1e-9)))
    return np.max(out, axis=0)


if __name__ == "__main__":
    _, fn, Sn, _ = net6.s6(f"{DATA}/dense_s2_g12_n.json"); f4, S4 = net45()
    bn = ((fn >= 2.399e9) & (fn <= 2.4836e9)) | ((fn >= 5.149e9) & (fn <= 7.126e9)); b4 = (f4 >= 5.149e9) & (f4 <= 7.126e9)
    rows = []
    for G in np.arange(8.0, 13.01, 0.5):
        Sng = net6.regap(Sn, fn, 12.0, G); S4g = regap(S4, f4, 45.0, 12.0, G)
        for Lc in np.arange(1.5, 3.51, 0.5):
            n0 = one_scan.worst(Sng, fn, Lc * 1e-9, "two"); n45 = worst45(f4, S4g, Lc * 1e-9)
            rows.append((max(n0[bn].max(), n45[b4].max()), G, Lc, n0[bn].max(), n45[b4].max(), f4[b4][np.argmax(n45[b4])]))
            print(f"G {G:4.1f} Lc {Lc:.1f}: normal {n0[bn].max():.2f}  45deg(5.15-7.125) {n45[b4].max():.2f} @ {rows[-1][5]/1e9:.2f}", flush=True)
    rows.sort(key=lambda r: r[0])
    print("BEST:", [f"G {r[1]} Lc {r[2]}: joint {r[0]:.2f} (normal {r[3]:.2f}, 45 {r[4]:.2f})" for r in rows[:5]])
