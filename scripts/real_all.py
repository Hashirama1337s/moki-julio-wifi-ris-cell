"""All angles with a REAL series chip inductor (self-resonance + Q) and diode R_s 5.2 / 7.0 ohm; optional mounting C across each diode.
usage: python real_all.py --srf 10 --L 1.0 1.2 --cp 0 0.05 [--qflat 14]
Q model: Q(f) = 14 * sqrt(f / 0.5 GHz) capped at 40 (datasheet minimum 14 at 500 MHz, rising), or flat with --qflat."""
import argparse
import numpy as np
import real_box as rb
from phi_an import net
import os
DATA = os.environ.get("RIS_DATA", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"))   # network JSON files

CASES = [("0 deg", f"{DATA}/boss_d0.json", None), ("30 edge", f"{DATA}/boss_d30.json", None), ("45 edge", f"{DATA}/boss_d45.json", None),
         ("45 diag", f"{DATA}/boss_p45_te.json", f"{DATA}/boss_p45_tm.json"), ("2.4@30e", f"{DATA}/boss24_t30.json", None),
         ("2.4@45e", f"{DATA}/boss24_t45.json", None), ("2.4@30d", f"{DATA}/boss24_p45t30_te.json", f"{DATA}/boss24_p45t30_tm.json"),
         ("2.4@45d", f"{DATA}/boss24_p45.json", f"{DATA}/boss24_p45_tm.json")]


def chip_model(srf, qflat=None):
    def chip(w, L, Q=14, srf_=None):
        q = qflat if qflat else np.minimum(14 * np.sqrt(w / (2 * np.pi * 0.5e9)), 40)
        zl = 1j * w * L + w * L / q
        return 1 / (1 / zl + 1j * w / ((2 * np.pi * srf) ** 2 * L))
    return chip


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--srf", type=float, default=10.0); ap.add_argument("--qflat", type=float, default=None)
    ap.add_argument("--L", type=float, nargs="+", default=[1.0, 1.2, 1.4]); ap.add_argument("--cp", type=float, nargs="+", default=[0.0])
    a = ap.parse_args(); rb.chip = chip_model(a.srf * 1e9, a.qflat)
    data = []
    for name, te, tm in CASES:
        if tm:
            d, f, S = net(te, tm); full = True
        else:
            d, f, S = rb.prep(te); full = False
        data.append((name, f, S, full))
    for cp in a.cp:
        for Ln in a.L:
            row = []
            for name, f, S, full in data:
                wv = rb.worst(f, S, Ln * 1e-9, cp * 1e-12, full=full)
                b = ((f >= 2.399e9) & (f <= 2.4836e9)) | ((f >= 5.149e9) & (f <= 7.126e9))
                row.append(f"{name}:{wv[b].max():.2f}")
            print(f"SRF {a.srf:g} GHz Q {'flat %g' % a.qflat if a.qflat else 'f'} Cpad {cp:.2f} pF  L {Ln:.1f}: " + " ".join(row), flush=True)


if __name__ == "__main__":
    main()
