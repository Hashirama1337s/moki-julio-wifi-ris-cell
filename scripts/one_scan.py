"""One vs two diodes per diagonal with spacer G and series chip inductor Lc both optimised (net6 transfer).
Worst case over the MADP-000907 box x D1/D2 +-5 % x Lc +-0.1 nH. usage: python one_scan.py <mp.json at G 10> ..."""
import sys
import numpy as np
from net6 import s6, regap, close6

D = dict(Ron=5.2, C=0.025e-12, L=0.3e-9, Roff=60.0, Cb=(0.022e-12, 0.030e-12), Lb=(0.2e-9, 0.4e-9)); Q = 30.0


def worst(S6, f, Lc, mode):
    w = 2 * np.pi * f; out = []; SH = np.zeros(len(f), complex)
    for C in D["Cb"] + (D["C"],):
        for L in D["Lb"] + (D["L"],):
            for dl in ((-0.1e-9, 0.0, 0.1e-9) if Lc > 0 else (0.0,)):
                for mm in (0.95, 1.0, 1.05):
                    lc = Lc + dl; zl = 1j * w * lc + w * lc / Q
                    on = D["Ron"] + 1j * w * L + zl; off1 = D["Roff"] + 1 / (1j * w * C) + 1j * w * L + zl
                    off2 = D["Roff"] + 1 / (1j * w * C * mm) + 1j * w * L + zl
                    if mode == "two":
                        ZA, ZB = [on, on, off2, off2], [off1, off1, on, on]
                    else:
                        ZA, ZB = [on, SH, SH, off2], [off1, SH, SH, on]
                    A = close6(S6, np.stack(ZA, 1)); B = close6(S6, np.stack(ZB, 1))
                    out.append(-20 * np.log10(abs(A[:, 1, 0] - B[:, 1, 0]) / 2))      # switched cross-pol component
    return np.max(out, axis=0)


if __name__ == "__main__":
    for path in sys.argv[1:]:
        d, f, S6, _ = s6(path); print(path.split("/")[-1])
        for mode, lcs in (("two", np.arange(0, 3.51, 0.5)), ("one", np.arange(0, 8.01, 1.0))):
            best = None
            for G in np.arange(6.0, 14.01, 0.5):
                Sg = regap(S6, f, 10.0, G)
                for Lc in lcs:
                    wv = worst(Sg, f, Lc * 1e-9, mode)
                    if best is None or wv.max() < best[0].max():
                        best = (wv, G, Lc)
            print(f"  {mode} diode(s)/diag: worst {best[0].max():.2f} dB  G {best[1]:.1f}  Lc {best[2]:.1f} nH  per-f {np.round(best[0], 2)}")
