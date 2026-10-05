"""Supercell analysis (2026-10-04): neighbour (opposite-code) coupling check for the boss design.
Codes per cell: A = NE+SW on, B = NW+SE on; vias shorted. From one supercell multiport:
  AA, BB  -> specular R00 must equal the single-cell periodic result (checks the supercell itself);
  AB      -> NOTE: the specular test R00 = R_fix is only meaningful when the (-1,0) order propagates; at normal incidence with a
             30 mm supercell it is evanescent, the switched power must return specularly (checked: |R00|^2 = 0.70-0.93) -> ignore there.
             Local periodicity predicts order (-1,0) = R_sw = (R_AA - R_BB)/2 (oblique, if
             propagating). The printed 'diff vs uniform prediction' uses |R_sw| as the reference, which ignores the 1-bit 2/pi factor and obliquity;
             use sup_energy.py for the fair comparison (energy accounting + physical-optics 1-bit reference).
Parts: realistic series inductor (2.0 nH, SRF 20 GHz, Q(f)), diode nominal and R_s 7 ohm / C_T corners (worst case over the box).
usage: python super_an.py <sup.json> [--single <single-cell.json>]"""
import sys
import numpy as np
from mp_pol_analyze import load, R0
from real_all import chip_model

chip = chip_model(20e9)


def prep(path):
    d, f, S = load(path)
    if d.get("theta", 0) != 0:                     # reciprocity + C2 renormalisation of the Floquet-driven port waves (TE column)
        # C2 about each cell centre maps the cell onto itself (NE<->SW, NW<->SE); about the supercell centre it swaps the cells.
        # Use the supercell-centre rotation: cell0 NE <-> cell1 SW etc.
        ROT = [5, 4, 7, 6, 1, 0, 3, 2, 9, 8]
        k = S[:, 0, 1:][:, ROT] / S[:, 4:, 0]; wt = abs(S[:, 4:, 0]) ** 2; km = (k * wt).sum(1) / wt.sum(1)
        print("  renorm factor", np.round(km, 3), " spread", np.round(np.max(abs(k - km[:, None]) * (wt > 0.05 * wt.max(1, keepdims=True))), 3))
        S = S.copy(); S[:, 4:, 0] = S[:, 4:, 0] * (-km)[:, None]
    return d, f, S


def close(S, i, Z):
    g = (Z - R0) / (Z + R0); Gm = np.diag(g)
    b = np.linalg.solve(np.eye(len(Z)) - S[i, 4:, 1:] @ Gm, S[i, 4:, 0]); out = S[i, :4, 0] + S[i, :4, 1:] @ (Gm @ b)
    return out                                       # [TE00, TM00, TE(-1,0), TM(-1,0)]


def loads(f, i, code, Rs, C, L, Lc):
    w = 2 * np.pi * f[i]; zs = chip(np.array([w]), Lc)[0]
    on = Rs + 1j * w * L + zs; off = 60.0 + 1 / (1j * w * C) + 1j * w * L + zs
    z = []
    for c in code:
        z += [on, on, off, off] if c == "A" else [off, off, on, on]
    return np.array(z + [1e-3, 1e-3], complex)


if __name__ == "__main__":
    d, f, S = prep(sys.argv[1]); Lc = 2.0e-9
    print(sys.argv[1], "theta", d.get("theta"), "f GHz", np.round(f / 1e9, 3), "orders", d.get("orders"))
    corners = [(Rs, C, L) for Rs in (5.2, 7.0) for C in (0.022e-12, 0.025e-12, 0.030e-12) for L in (0.2e-9, 0.3e-9, 0.4e-9)]
    for i in range(len(f)):
        worst_pred, worst_ab, err_spec, err_m1 = 0, 0, 0, 0
        for Rs, C, L in corners:
            R = {c: close(S, i, loads(f, i, c, Rs, C, L, Lc)) for c in ("AA", "BB", "AB")}
            Rsw = (R["AA"][:2] - R["BB"][:2]) / 2; Rfix = (R["AA"][:2] + R["BB"][:2]) / 2
            pred = -20 * np.log10(np.linalg.norm(Rsw))                       # single-cell-equivalent steerable loss
            spec_err = np.linalg.norm(R["AB"][:2] - Rfix) / max(np.linalg.norm(Rsw), 1e-9)
            worst_pred = max(worst_pred, pred); err_spec = max(err_spec, spec_err)
            if not np.isnan(R["AB"][2]):
                ab = -20 * np.log10(np.linalg.norm(R["AB"][2:]))              # steered (-1,0) order, both pols
                worst_ab = max(worst_ab, ab); err_m1 = max(err_m1, abs(ab - pred))
        line = f"  {f[i] / 1e9:.2f} GHz: uniform-code steerable loss (worst) {worst_pred:.2f} dB | AB specular error |R00 - R_fix|/|R_sw| <= {err_spec:.3f}"
        if worst_ab:
            line += f" | AB steered (-1,0) loss (worst) {worst_ab:.2f} dB, |diff vs uniform prediction| <= {err_m1:.2f} dB"
        print(line)
    if "--single" in sys.argv:
        import real_box as rb
        rb.chip = chip
        d1, f1, S1 = rb.prep(sys.argv[sys.argv.index("--single") + 1])
        idx = [int(np.argmin(abs(f1 - x))) for x in f]
        w1 = rb.worst(f1, S1, Lc, 0.0)
        print("  single-cell periodic (same parts, worst):", np.round(w1[idx], 2), "at", np.round(f1[idx] / 1e9, 3))
