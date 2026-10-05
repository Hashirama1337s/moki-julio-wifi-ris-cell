"""Realistic-parts scoring of a saved multiport (closure with lumped loads at the diode sites).
Box: diode R_s list (default 5.2 and 7.0 ohm), C_T list `ct` (default 0.022 / 0.030 / 0.025 pF), package L_s 0.2-0.4 nH, +-5 % OFF-pair
mismatch, series chip inductor Lnom +-0.1 nH through the module-level `chip(w, L)` model (callers set it, e.g. real_all.chip_model(20e9)),
optional mounting capacitance `cpar` across each diode. Works on via-port files (5 ports, via shorted) and full 7-port oblique files.
`worst(..., phase=True)` also returns the worst equivalent bit-phase error 2 atan(|R_fix e| / |R_sw e|).
usage: python real_box.py <mp.json> [tm.json] [Lnom ...]"""
import sys
import numpy as np
from mp_pol_analyze import load, R0

ROT = [1, 0, 3, 2, 4]


def prep(path):
    d, f, S = load(path)
    if d.get("theta", 0) != 0:
        k = S[:, 0, 1:][:, ROT] / S[:, 2:, 0]; wt = abs(S[:, 2:, 0]) ** 2; km = (k * wt).sum(1) / wt.sum(1)
        S = S.copy(); S[:, 2:, 0] = S[:, 2:, 0] * (-km)[:, None]
    return d, f, S


def chip(w, L, Q=14.0, srf=10e9):
    zl = 1j * w * L + w * L / Q; cp = 1 / ((2 * np.pi * srf) ** 2 * L); return 1 / (1 / zl + 1j * w * cp)


def worst(f, S, Lnom, cpar=0.0, rs_list=(5.2, 7.0), full=False, phase=False, ct=(0.022e-12, 0.030e-12, 0.025e-12)):
    """Worst-case loss (dB) per frequency over the box. phase=True also returns the worst equivalent bit-phase error (deg),
    2 atan(|R_fix e| / |R_sw e|): the A-vs-B phase error of an equal-amplitude bit with the same unswitched (specular) leak. It does not
    depend on how the output polarisations are booked (for E along a cell diagonal at 0 deg it is close to the direct A-vs-B phase error;
    it also counts amplitude imbalance, so it is the stricter of the two)."""
    w = 2 * np.pi * f; out = []; perr = np.zeros(len(f)); zv = np.full(len(f), 1e-3 + 0j); npt = S.shape[2] - (2 if full else 1)
    for Rs in rs_list:
        for C in ct:
            for L in (0.2e-9, 0.4e-9, 0.3e-9):
                for dl in (-0.1e-9, 0.0, 0.1e-9):
                    for mm in (0.95, 1.0, 1.05):
                        zc = chip(w, Lnom + dl)
                        dio = lambda z: 1 / (1 / z + 1j * w * cpar) if cpar else z
                        on = dio(Rs + 1j * w * L) + zc; o1 = dio(60.0 + 1 / (1j * w * C) + 1j * w * L) + zc
                        o2 = dio(60.0 + 1 / (1j * w * C * mm) + 1j * w * L) + zc
                        RR = []
                        for Z in (np.stack([on, on, o2, o2, zv], 1), np.stack([o1, o1, on, on, zv], 1)):
                            g = (Z - R0) / (Z + R0); R = []
                            for i in range(len(f)):
                                Gm = np.diag(g[i])
                                if full:
                                    R.append(S[i, :2, :2] + S[i, :2, 2:] @ Gm @ np.linalg.solve(np.eye(5) - S[i, 2:, 2:] @ Gm, S[i, 2:, :2]))
                                else:
                                    b = np.linalg.solve(np.eye(5) - S[i, 2:, 1:] @ Gm, S[i, 2:, 0]); R.append((S[i, :2, 0] + S[i, :2, 1:] @ (Gm @ b))[:, None])
                            RR.append(np.array(R))
                        sw = (RR[0] - RR[1]) / 2
                        for j in range(sw.shape[2]):
                            out.append(-20 * np.log10(np.linalg.norm(sw[:, :, j], axis=1)))
                            if phase:
                                fx = np.linalg.norm((RR[0][:, :, j] + RR[1][:, :, j]) / 2, axis=1)
                                perr = np.maximum(perr, np.degrees(2 * np.arctan(fx / np.linalg.norm(sw[:, :, j], axis=1))))
    return (np.max(out, axis=0), perr) if phase else np.max(out, axis=0)


if __name__ == "__main__":
    args = sys.argv[1:]; full = len(args) > 1 and args[1].endswith(".json")
    if full:
        from phi_an import net
        d, f, S = net(args[0], args[1]); Ls = [float(x) for x in args[2:]] or [1.8, 2.0, 2.2, 2.5]
    else:
        d, f, S = prep(args[0]); Ls = [float(x) for x in args[1:]] or [1.8, 2.0, 2.2, 2.5]
    b24 = (f >= 2.399e9) & (f <= 2.4836e9); b5 = (f >= 5.149e9) & (f <= 7.126e9)
    print(args[0], "theta", d.get("theta"), "phi", d.get("phi"), len(f), "pts", "(both pols)" if full else "(TE in)")
    for Ln in Ls:
        for cpar in (0.0, 0.05e-12):
            wv = worst(f, S, Ln * 1e-9, cpar, full=full)
            m = lambda b: f"{wv[b].max():.2f}" if b.any() else "  - "
            print(f"  L {Ln:.1f} nH (0201, Q14, SRF10){' +0.05pF' if cpar else '         '}: 2.4 band {m(b24)} | 5.15-7.125 {m(b5)}  (Rs 5.2 & 7.0)")
