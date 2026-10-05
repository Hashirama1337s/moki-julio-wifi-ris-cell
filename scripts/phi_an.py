"""Full-network oblique analysis for ANY plane (incl. phi 45): TE multiport <te.json> + TM-in direct run with all ports on 50 ohm <tm.json>
(5 ports incl. the via). Reciprocity + C2 renormalisation (weights |S|^2). Steerable loss for TE and TM incidence = |R_sw e|.
usage: python phi_an.py <te.json> <tm.json> [Lc_nH ...]"""
import json, sys
import numpy as np
from mp_pol_analyze import load, R0
from via_joint import worst45

def net(te, tm):
    d, f, S = load(te); t = json.load(open(tm)); c = lambda x: np.array([complex(*z) for z in x]); ROT = [1, 0, 3, 2, 4]
    Vtm = np.stack([c(v) for v in t["V"]], 1) / np.sqrt(2 * R0); n = 7
    S7 = np.zeros((len(f), n, n), complex); S7[:, :2, 0] = S[:, :2, 0]; S7[:, :, 2:] = S[:, :, 1:]
    S7[:, 0, 1] = c(t["TE"]); S7[:, 1, 1] = c(t["TM"])
    wm = lambda k, a: (k * abs(a) ** 2).sum(1) / (abs(a) ** 2).sum(1)
    kte = S[:, 0, 1:][:, ROT] / S[:, 2:, 0]; ktm = S[:, 1, 1:][:, ROT] / Vtm
    S7[:, 2:, 0] = S[:, 2:, 0] * (-wm(kte, S[:, 2:, 0]))[:, None]; S7[:, 2:, 1] = Vtm * (-wm(ktm, Vtm))[:, None]
    return d, f, S7

if __name__ == "__main__":
    d, f, S7 = net(sys.argv[1], sys.argv[2]); b5 = (f >= 5.149e9) & (f <= 7.126e9)
    print(sys.argv[1], "theta", d.get("theta"), "phi", d.get("phi"), len(f), "pts")
    for Lc in [float(x) for x in sys.argv[3:]] or [2.0, 2.5, 3.0]:
        wv = worst45(f, S7, Lc * 1e-9)
        print(f"  Lc {Lc}: worst (TE and TM in) 5.15-7.125 {wv[b5].max():.2f} dB @ {f[b5][np.argmax(wv[b5])]/1e9:.2f};", " ".join(f"{x/1e9:.1f}:{y:.1f}" for x, y in zip(f[::5], wv[::5])))
