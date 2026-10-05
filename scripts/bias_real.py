"""Bias network with realistic parts (2026-10-04; packages the calculation behind the README bias numbers).
Boss design + 4 corner bias posts (ports 7-10 hold the choke). Series inductor: generic 0201 model (SRF 20 GHz, Q(f)); diode R_s 5.2 and
7.0 ohm, C_T 0.022-0.030 pF, L_s 0.2-0.4 nH, +-5 % pair mismatch, inductor +-0.1 nH. Chokes: open, ideal 47 nH, and two 0201 chips in
series, 33 nH (SRF 3.2 GHz) + 10 nH (SRF 8.5 GHz), Q(f) as the series inductor. Reference: the same cell without posts.
usage: python bias_real.py <bias.json> <no-post.json> [Lc_nH ...]"""
import sys
import numpy as np
import bias_an as ba
import real_box as rb
from real_all import chip_model

chip = chip_model(20e9)


def chip_srf(w, L, srf):
    q = np.minimum(14 * np.sqrt(w / (2 * np.pi * 0.5e9)), 40); zl = 1j * w * L + w * L / q
    return 1 / (1 / zl + 1j * w / ((2 * np.pi * srf) ** 2 * L))


def chokes(f):
    w = 2 * np.pi * f
    return {"open": np.full(len(f), 1e12 + 0j), "ideal 47 nH": 1j * w * 47e-9 + 1.0,
            "33 nH + 10 nH chips": chip_srf(w, 33e-9, 3.2e9) + chip_srf(w, 10e-9, 8.5e9)}


def worst(f, S, Lc, zch):
    w = 2 * np.pi * f; out = []
    for Rs in (5.2, 7.0):
        for C in (0.022e-12, 0.025e-12, 0.030e-12):
            for L in (0.2e-9, 0.3e-9, 0.4e-9):
                for dl in (-0.1e-9, 0.0, 0.1e-9):
                    for mm in (0.95, 1.0, 1.05):
                        zl = chip(w, Lc + dl); on = Rs + 1j * w * L + zl
                        o1 = 60 + 1 / (1j * w * C) + 1j * w * L + zl; o2 = 60 + 1 / (1j * w * C * mm) + 1j * w * L + zl
                        RA = np.array([ba.te_in(S, i, np.array([on[i], on[i], o2[i], o2[i], 1e-3] + [zch[i]] * 4)) for i in range(len(f))])
                        RB = np.array([ba.te_in(S, i, np.array([o1[i], o1[i], on[i], on[i], 1e-3] + [zch[i]] * 4)) for i in range(len(f))])
                        out.append(-20 * np.log10(np.linalg.norm((RA - RB) / 2, axis=1)))
    return np.max(out, axis=0)


if __name__ == "__main__":
    d, f, S = ba.prep(sys.argv[1]); d0, f0, S0 = rb.prep(sys.argv[2]); idx = [int(np.argmin(abs(f0 - x))) for x in f]
    rb.chip = chip; Ls = [float(x) * 1e-9 for x in sys.argv[3:]] or [2.0e-9, 2.2e-9]
    print(sys.argv[1], "theta", d.get("theta"), "f GHz", np.round(f / 1e9, 3))
    for Lc in Ls:
        ref = rb.worst(f0, S0, Lc, 0.0)[idx]
        ok = np.array([abs(f0[i] - x) < 0.03e9 for i, x in zip(idx, f)])          # reference only where the no-post file has that frequency
        print(f"  Lc {Lc * 1e9:.1f} nH  no posts (reference)   : " + str([round(float(r), 2) if o else 'n/a' for r, o in zip(ref, ok)]))
        for k, z in chokes(f).items():
            print(f"  Lc {Lc * 1e9:.1f} nH  posts, {k:20s}: {np.round(worst(f, S, Lc, z), 2)}")
