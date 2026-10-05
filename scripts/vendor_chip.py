"""Measured vendor series-inductor model: Coilcraft 0201DS 2.3 nH (02DS-2N3.s2p), Z(f) interpolated; optional scaling to a nearby
nominal value (Z * L_nom / 2.3 nH, a first-order stand-in for the sibling part)."""
import os, sys
import numpy as np
_csv = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "vendor", "coilcraft_02DS-2N3_Zseries.csv")
if not os.path.exists(_csv):
    sys.exit("vendor_chip: vendor/coilcraft_02DS-2N3_Zseries.csv missing; build it with scripts/coilcraft_to_csv.py from the Coilcraft .s2p")
_d = np.loadtxt(_csv, delimiter=",", comments="#")


def chip(w, L, Q=None, srf=None):
    f = w / (2 * np.pi) / 1e9
    z = np.interp(f, _d[:, 0], _d[:, 1]) + 1j * np.interp(f, _d[:, 0], _d[:, 2])
    return z * (L / 2.3e-9)


if __name__ == "__main__":
    # this shape scored with the measured Coilcraft 0201DS-2N3 (2.3 nH, +-0.1 nH by scaling) over the same box as reproduce.py
    import real_box as rb
    from phi_an import net
    rb.chip = chip; D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data") + "/"
    for name, te, tm in (("0 deg", "v4_d0", None), ("30 deg edge", "v4_d30", None), ("45 deg edge", "v4_d45", None),
                         ("30 deg diagonal", "v4_p45t30_te", "v4_p45t30_tm"), ("45 deg diagonal", "v4_p45_te", "v4_p45_tm")):
        d, f, S = net(D + te + ".json", D + tm + ".json") if tm else rb.prep(D + te + ".json"); lw = np.zeros(len(f))
        for cp in (0.013e-12, 0.015e-12, 0.017e-12):
            lw = np.maximum(lw, rb.worst(f, S, 2.3e-9, cp, full=bool(tm)))
        m = (f >= 5.149e9) & (f <= 7.126e9)
        print(f"{name:16s} measured 0201DS-2N3: worst loss 5.15-7.125 GHz {lw[m].max():.2f} dB")
