"""Regenerate every headline number in README.md (v4) from the network files in data/ (no field solver needed).
usage:  python reproduce.py            (from the package top folder; needs numpy)
Design v4 (new shape: stub 2.3 mm, diode gap 1.7 mm, spacer 12.75 mm, boss + via): series inductor 2.3 nH (generic 0201 model, self-resonance
20 GHz, Q(f) = 14 sqrt(f / 0.5 GHz) <= 40); diode box R_s 5.2 / 7.0 ohm, C_T 0.022 / 0.025 / 0.030 pF, L_s 0.2-0.4 nH, +-5 % pair
mismatch, inductor +-0.1 nH; mounting capacitance across each diode 0.015 +- 0.002 pF (designed in).
Each line: worst-case loss (dB) and worst equivalent bit-phase error 2 atan(|R_fix e| / |R_sw e|) (deg), per band.
The earlier v3 shape (stub 0.3, gap 0.644, spacer 12 mm) is scored under the same rules for comparison, and the lower-capacitance diode
option (MACOM MA4AGFCP910: C_T 0.016 / 0.018 / 0.021 pF, R_s 5.2 / 6.0 ohm) on this shape, phase-tuned (1.9 nH) and loss-tuned (2.5 nH)."""
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "scripts")); os.chdir(HERE)
import real_box as rb
from real_all import chip_model
from phi_an import net
from mp_pol_analyze import load, close, R0

rb.chip = chip_model(20e9)
D = "data/"
PADS = (0.013e-12, 0.015e-12, 0.017e-12)
MADP = dict(ct=(0.022e-12, 0.030e-12, 0.025e-12), rs=(5.2, 7.0))
FCP910 = dict(ct=(0.016e-12, 0.021e-12, 0.018e-12), rs=(5.2, 6.0))
CASES = [("0 deg", "d0", None), ("30 deg edge", "d30", None), ("45 deg edge", "d45", None),
         ("30 deg diagonal", "p45t30_te", "p45t30_tm"), ("45 deg diagonal", "p45_te", "p45_tm"),
         ("2.4 GHz, 30 deg edge", "24_t30", None), ("2.4 GHz, 45 deg edge", "24_t45", None),
         ("2.4 GHz, 30 deg diagonal (dense)", "d24_p45t30_te", "d24_p45t30_tm"), ("2.4 GHz, 45 deg diagonal", "24_p45", "24_p45_tm")]
OLD = {"d0": "boss_d0", "d30": "boss_d30", "d45": "boss_d45", "p45t30_te": "boss_p45t30_te", "p45t30_tm": "boss_p45t30_tm", "p45_te": "boss_p45_te",
       "p45_tm": "boss_p45_tm", "24_t30": "boss24_t30", "24_t45": "boss24_t45", "d24_p45t30_te": "d24_p45t30_te", "d24_p45t30_tm": "d24_p45t30_tm",
       "24_p45": "boss24_p45", "24_p45_tm": "boss24_p45_tm"}
NEW = {k: "v4_" + k for k in OLD}


def band_max(f, wv, nd=1):
    out = []
    for a, b in ((2.399e9, 2.4836e9), (5.149e9, 7.126e9)):
        m = (f >= a) & (f <= b); out.append(f"{wv[m].max():5.{nd}f}" if m.any() else "  -  ")
    return out


def table(title, names, L, diode, pads):
    print(f"\n{title}\n{'case':34s} loss 2.4 / 5-7 GHz (dB)    eq. phase 2.4 / 5-7 GHz (deg)")
    for name, te, tm in CASES:
        pte = D + names[te] + ".json"; ptm = D + names[tm] + ".json" if tm else None
        if not os.path.exists(pte) or (ptm and not os.path.exists(ptm)):
            print(f"{name:34s} (data not present)"); continue
        d, f, S = net(pte, ptm) if tm else rb.prep(pte)
        lw = np.zeros(len(f)); pw = np.zeros(len(f))
        for cp in pads:
            wv, pe = rb.worst(f, S, L, cp, rs_list=diode["rs"], full=bool(tm), phase=True, ct=diode["ct"]); lw = np.maximum(lw, wv); pw = np.maximum(pw, pe)
        a = band_max(f, lw, 2); p = band_max(f, pw)
        print(f"{name:34s} {a[0]} / {a[1]}                 {p[0]} / {p[1]}     ({len(f)} frequencies, {'TE+TM' if tm else 'TE'} incidence)")


table("DESIGN v4: new shape, MADP-000907, 2.3 nH, mounting 0.015 +- 0.002 pF", NEW, 2.3e-9, MADP, PADS)
table("COMPARISON: earlier v3 shape, same rules (2.3 nH, mounting 0.015 +- 0.002 pF)", OLD, 2.3e-9, MADP, PADS)
table("COMPARISON: earlier v3 shape as designed in v3 (2.0 nH, no mounting capacitance)", OLD, 2.0e-9, MADP, (0.0,))
table("OPTION A (phase-tuned): this shape with MA4AGFCP910, 1.9 nH, mounting 0.015 +- 0.002 pF", NEW, 1.9e-9, FCP910, PADS)
table("OPTION B (loss-tuned): this shape with MA4AGFCP910, 2.5 nH, mounting 0.015 +- 0.002 pF", NEW, 2.5e-9, FCP910, PADS)

# typical parts, 0 deg
for lab, names, L, diode, C0 in (("v4 design", NEW, 2.3e-9, MADP, 0.025e-12), ("v3 shape", OLD, 2.3e-9, MADP, 0.025e-12),
                                  ("v4 + MA4AGFCP910, 1.9 nH", NEW, 1.9e-9, FCP910, 0.018e-12),
                                  ("v4 + MA4AGFCP910, 2.5 nH", NEW, 2.5e-9, FCP910, 0.018e-12)):
    p = D + names["d0"] + ".json"
    if not os.path.exists(p):
        continue
    d, f, S = rb.prep(p); w = 2 * np.pi * f; zv = np.full(len(f), 1e-3 + 0j); cp = 0.015e-12; z = rb.chip(w, L)
    dio = lambda x: 1 / (1 / x + 1j * w * cp)
    on = dio(5.2 + 1j * w * 0.3e-9) + z; off = dio(60 + 1 / (1j * w * C0) + 1j * w * 0.3e-9) + z; R = []
    for Z in (np.stack([on, on, off, off, zv], 1), np.stack([off, off, on, on, zv], 1)):
        g = (Z - R0) / (Z + R0)
        R.append(np.array([S[i, :2, 0] + S[i, :2, 1:] @ (np.diag(g[i]) @ np.linalg.solve(np.eye(5) - S[i, 2:, 1:] @ np.diag(g[i]), S[i, 2:, 0])) for i in range(len(f))]))
    A, B = R; e = np.degrees(2 * np.arctan(np.linalg.norm((A + B) / 2, axis=1) / np.linalg.norm((A - B) / 2, axis=1))); q = band_max(f, e)
    print(f"0 deg, typical parts, {lab:20s}: equivalent bit-phase error {q[0]} / {q[1]} deg")

print("\nvalidation: network model vs direct Palace run with passive loads (earlier design, 0 deg, 12 frequencies)")
d, f, S = load(D + "v12_s4.json"); dd, _, _ = load(D + "vdir_s4.json")
R = np.array([dd["direct"][str(k)] for k in (2, 3, 4, 5)], float)
te, tm = close(S, np.tile(R, (len(f), 1)).astype(complex))
TEd = np.array([complex(*z) for z in dd["TE"]]); TMd = np.array([complex(*z) for z in dd["TM"]])
print(f"  max |model - direct|: TE {np.max(abs(te - TEd)):.1e}, TM {np.max(abs(tm - TMd)):.1e}")
print("mesh convergence (v3 boss design, 0 deg, 4 frequencies, 2.0 nH, no mounting C): h0 0.25 vs 0.18 mm")
for p in ("boss_n", "bossfine_n"):
    d, f, S = rb.prep(D + p + ".json"); print(f"  {p:11s} ({d['tets']} elements):", np.round(rb.worst(f, S, 2.0e-9, 0.0), 2))
