"""Eigen-axis extraction for the X pol-rotation cell (design study, 2026-10-03).
By mirror symmetry r_u depends only on the NE/SW diodes and r_v only on NW/SE; co = (r_u + r_v)/2, cross = +-(r_u - r_v)/2.
From a saved multiport (palace_mp_pol.py) this: (1) identifies r_u and checks the symmetry claim (r_u must not move when the
OFF loads change); (2) de-embeds the port and the spacer to get the sheet admittance of one diagonal chain with the diode
SHORTED and OPEN (normalised to 1/eta0); (3) fits the chain X = w*La - 1/(w*Cp) and reports its series resonance;
(4) prints delta = arg(r_u/r_v) for real diodes (lossless conversion >= -2 dB needs |delta - 180| <= 75 deg).
usage: python eig_extract.py <mp.json> <G_mm> [--diode MADP000907]"""
import sys
import numpy as np
from mp_pol_analyze import load, close

ETA0 = 376.730313; C0 = 299792458.0; HUP = 20.0; T = 0.508; ER = 3.55
DIO = {"MADP000907": (5.2, 0.025e-12, 0.3e-9, 60.0), "MPP4203": (2.0, 0.09e-12, 0.2e-9, 20.0)}


def eig(S, f, z_u, z_v):
    """r at the sheet plane for both basis combinations; z_u on ports 2,3 (NE,SW), z_v on ports 4,5."""
    te, tm = close(S, np.stack([z_u, z_u, z_v, z_v], 1))
    k = 2 * np.pi * f / C0 * 1e-3                                  # rad/mm
    ph = np.exp(2j * k * HUP)                                      # Palace exp(+jwt): r_sheet = r_port * exp(+2jkh)
    return (te + tm) * ph, (te - tm) * ph


def y_down(f, G):
    """normalised admittance looking down from the sheet: RO4003C (T) over air G, shorted (exp(+jwt))."""
    k = 2 * np.pi * f / C0 * 1e-3; ks = k * np.sqrt(ER); zs = 1 / np.sqrt(ER)
    za = 1j * np.tan(k * G)
    zd = zs * (za + 1j * zs * np.tan(ks * T)) / (zs + 1j * za * np.tan(ks * T))
    return 1 / zd


if __name__ == "__main__":
    path, G = sys.argv[1], float(sys.argv[2])
    name = sys.argv[sys.argv.index("--diode") + 1] if "--diode" in sys.argv else "MADP000907"
    d, f, S = load(path); w = 2 * np.pi * f; nf = len(f)
    SH, OP = np.zeros(nf, complex), np.full(nf, 1e12 + 0j)
    a1, b1 = eig(S, f, SH, OP); a2, b2 = eig(S, f, SH, SH)          # change only the v-diagonal loads
    da, db = np.max(abs(a1 - a2)), np.max(abs(b1 - b2))
    ru_is_a = da < db
    print(f"{path}  G {G} mm   f GHz {np.round(f / 1e9, 3)}")
    print(f"  symmetry control: max change of (co+cross) {da:.2e}, of (co-cross) {db:.2e} when only v loads change"
          f" -> r_u = co{'+' if ru_is_a else '-'}cross")
    pick = (lambda a, b: a) if ru_is_a else (lambda a, b: b)
    r_sh = pick(*eig(S, f, SH, OP)); r_op = pick(*eig(S, f, OP, SH))
    print("  |r_u| diode shorted", np.round(abs(r_sh), 3), " open", np.round(abs(r_op), 3), "(lossless loads -> ~1)")
    yin = lambda r: (1 - r) / (1 + r)
    ys, yo = yin(r_sh) - y_down(f, G), yin(r_op) - y_down(f, G)
    print("  b_down (spacer)      ", np.round(y_down(f, G).imag, 3))
    print("  b_sheet diode SHORT  ", np.round(ys.imag, 3))
    print("  b_sheet diode OPEN   ", np.round(yo.imag, 3), " -> C0 eff pF", np.round(yo.imag / (w * ETA0) * 1e12, 4))
    zc = 1 / (ys - yo) * ETA0                                       # chain impedance (ohm-equivalent per square)
    X = zc.imag
    w9 = w / 1e9                                                    # scaled: [w, 1/w] in SI is ill-conditioned
    A = np.stack([w9, -1 / w9], 1); (a, b), *_ = np.linalg.lstsq(A, X, rcond=None)
    La, invCp = a * 1e-9, b * 1e9
    fr = 1 / (2 * np.pi * np.sqrt(La / invCp)) if La > 0 and invCp > 0 else float("nan")
    X = X; A = np.stack([w, -1 / w], 1)
    print("  chain X (ohm)        ", np.round(X, 1), "  R", np.round(zc.real, 2))
    print(f"  fit X = wLa - 1/(wCp): La {La * 1e9:.2f} nH, Cp {1 / invCp * 1e12:.3f} pF -> series resonance {fr / 1e9:.2f} GHz;"
          f" fit resid {np.round(X - A @ np.array([La, invCp]), 1)}")
    r10 = pick(*eig(S, f, np.full(nf, 10.0 + 0j), OP))              # diode -> chain transformer: dZ_chain = n * dZ_diode
    n = ((1 / (yin(r10) - y_down(f, G) - yo)) * ETA0 - zc) / 10.0
    print("  diode->chain ratio n (re)", np.round(n.real, 3), " im", np.round(n.imag, 3))
    Ron, C, L, Roff = DIO[name]
    zon = Ron + 1j * w * L; zoff = Roff + 1 / (1j * w * C) + 1j * w * L
    ru = pick(*eig(S, f, zon, zoff)); rv_from_B = pick(*eig(S, f, zoff, zon))   # r_v(A) = r_u(B) by C4
    delta = np.degrees(np.angle(ru / rv_from_B)) % 360
    print(f"  {name}: delta = arg(r_on/r_off) deg", np.round(delta, 1), " |r_on|", np.round(abs(ru), 3), " |r_off|", np.round(abs(rv_from_B), 3))
    print("  cross loss from delta+|r| dB", np.round(-20 * np.log10(abs(ru - rv_from_B) / 2), 2))
