"""Exact shape scorer (2026-10-04, replaces shape_tune's one-chain approximation): full 6-port network from a 0-deg 4-port run (net6.s6),
exact spacer transfer (net6.regap, validated at normal incidence), every diode on its own port (so the OFF-pair mismatch is exact),
mounting capacitance across each diode as a DESIGN value with a +-tolerance, realistic series chip L.
Score = worst equivalent bit-phase error over both bands and the box; loss must stay <= LMAX.
usage: python shape_exact.py <mp.json> [--G 9:13:0.25] [--L 1.0:2.6:0.1] [--pad 0.015 --padtol 0.005] [--lmax 1.5] [--full]"""
import argparse, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from net6 import s6, regap
from mp_pol_analyze import R0
from real_all import chip_model

chip = chip_model(20e9)


def close6v(S6, ZL):
    nf = len(ZL); g = (ZL - R0) / (ZL + R0); G = np.zeros((nf, 4, 4), complex); G[:, range(4), range(4)] = g
    X = np.linalg.solve(np.eye(4) - S6[:, 2:, 2:] @ G, S6[:, 2:, :2]); return S6[:, :2, :2] + S6[:, :2, 2:] @ G @ X


CT = (0.022e-12, 0.025e-12, 0.030e-12)


def box(full):
    if full:
        return [(Rs, C, Ls, dl, mm) for Rs in (5.2, 7.0) for C in CT for Ls in (0.2e-9, 0.3e-9, 0.4e-9)
                for dl in (-0.1e-9, 0.0, 0.1e-9) for mm in (0.95, 1.0, 1.05)]
    return [(Rs, C, Ls, dl, mm) for Rs in (5.2, 7.0) for C in (CT[0], CT[-1]) for Ls in (0.2e-9, 0.4e-9)
            for dl in (-0.1e-9, 0.1e-9) for mm in (0.95, 1.05)]


def score(S6, f, L, pads, corners):
    w = 2 * np.pi * f; pe = np.zeros(len(f)); ls = np.zeros(len(f))
    for cp in pads:
        dio = (lambda z: 1 / (1 / z + 1j * w * cp)) if cp else (lambda z: z)
        for Rs, C, Ls, dl, mm in corners:
            z = chip(w, L + dl); on = dio(Rs + 1j * w * Ls) + z
            o1 = dio(60 + 1 / (1j * w * C) + 1j * w * Ls) + z; o2 = dio(60 + 1 / (1j * w * C * mm) + 1j * w * Ls) + z
            A = close6v(S6, np.stack([on, on, o2, o2], 1)); B = close6v(S6, np.stack([o1, o1, on, on], 1))
            sw = np.linalg.norm((A - B)[:, :, 0] / 2, axis=1); fx = np.linalg.norm((A + B)[:, :, 0] / 2, axis=1)   # TE in (0 deg: any pol same loss)
            pe = np.maximum(pe, np.degrees(2 * np.arctan(fx / sw))); ls = np.maximum(ls, -20 * np.log10(sw))
    return pe, ls


def rng(s):
    a, b, c = (float(x) for x in s.split(":")); return np.arange(a, b + 1e-9, c)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("path"); ap.add_argument("--G", default="9:13:0.25"); ap.add_argument("--L", default="1.0:2.6:0.1")
    ap.add_argument("--pad", type=float, default=0.015); ap.add_argument("--padtol", type=float, default=0.005); ap.add_argument("--lmax", type=float, default=1.5)
    ap.add_argument("--full", action="store_true"); ap.add_argument("--G0", type=float, default=12.0)
    ap.add_argument("--ct", default="0.022,0.025,0.030", help="diode C_T box, pF")
    a = ap.parse_args(); CT = tuple(float(x) * 1e-12 for x in a.ct.split(","))
    d, f, S6, rec = s6(a.path)
    m = ((f >= 2.399e9) & (f <= 2.4836e9)) | ((f >= 5.149e9) & (f <= 7.126e9)); f, S6 = f[m], S6[m]; b24 = f < 3e9
    pads = sorted({max(0.0, a.pad - a.padtol), a.pad, a.pad + a.padtol}) if a.padtol else [a.pad]
    pads = [p * 1e-12 for p in pads]; corners = box(a.full); res = []
    for G in rng(a.G):
        SG = regap(S6, f, a.G0, G)
        for L in rng(a.L):
            pe, ls = score(SG, f, L * 1e-9, pads, corners)
            res.append((max(pe.max(), 0) + 100 * max(0, ls.max() - a.lmax), G, L, pe[b24].max(), pe[~b24].max(), ls[b24].max(), ls[~b24].max()))
    res.sort()
    print(f"{a.path} (pad {a.pad}+-{a.padtol} pF, {'full' if a.full else 'reduced'} box, reciprocity residual {rec:.1e})")
    for r in res[:3]:
        print(f"   G {r[1]:5.2f} mm  L {r[2]:.1f} nH -> eq phase {r[3]:.1f} / {r[4]:.1f} deg, loss {r[5]:.2f} / {r[6]:.2f} dB")
