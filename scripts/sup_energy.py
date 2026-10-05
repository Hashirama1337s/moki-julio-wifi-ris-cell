"""Energy accounting for the supercell (2026-10-04): where does the power go for AA, BB, AB codes?
Absorbed = 1 - |R00|^2 - |R(-1,0)|^2 (TE+TM). If AB absorbs about as much as AA/BB, neighbour coupling adds no loss;
the specular/steered split is then grating physics (identical for an ideal cell), not a cell defect.
Also prints the ideal physical-optics 1-bit reference (2/pi)^2 = -3.92 dB per first-order beam and the same scaled by this cell's
|R_sw|^2 (a reference, not a ceiling: here only the (-1,0) order propagates)."""
import sys
import numpy as np
import super_an as sa

d, f, S = sa.prep(sys.argv[1]); Lc = 2.0e-9; th = np.radians(d.get("theta", 0)); A = 30e-3
corners = [(Rs, C, L) for Rs in (5.2, 7.0) for C in (0.022e-12, 0.025e-12, 0.030e-12) for L in (0.2e-9, 0.3e-9, 0.4e-9)]
for i in range(len(f)):
    lam = 3e8 / f[i]; sr = np.sin(th) - lam / A; prop = abs(sr) < 1          # (-1,0) order propagating?
    tr = np.degrees(np.arcsin(sr)) if prop else np.nan; cr = np.cos(np.arcsin(sr)) / np.cos(th) if prop else 1.0
    rows = []
    for Rs, C, L in corners:
        R = {c: sa.close(S, i, sa.loads(f, i, c, Rs, C, L, Lc)) for c in ("AA", "BB", "AB")}
        p = {c: (np.sum(abs(R[c][:2]) ** 2), np.nansum(abs(R[c][2:]) ** 2)) for c in R}
        sw = np.linalg.norm((R["AA"][:2] - R["BB"][:2]) / 2) ** 2
        rows.append([1 - sum(p["AA"]), 1 - sum(p["BB"]), 1 - sum(p["AB"]), p["AB"][0], p["AB"][1], sw])
    r = np.array(rows); w = r[np.argmin(r[:, 4])] if prop else r[np.argmax(r[:, 2])]
    po = -10 * np.log10(w[5]) + 3.92
    if not prop:
        print(f"{f[i]/1e9:.2f} GHz  (-1,0) evanescent | worst corner: absorbed AA {w[0]:.3f} BB {w[1]:.3f} AB {w[2]:.3f} | max over corners: AB absorbed {r[:,2].max():.3f} vs uniform {r[:,:2].max():.3f}"); continue
    print(f"{f[i]/1e9:.2f} GHz  steered beam {tr:6.1f} deg  | worst corner: absorbed AA {w[0]:.3f} BB {w[1]:.3f} AB {w[2]:.3f}"
          f" | AB specular {w[3]:.3f} steered {w[4]:.3f} ({-10*np.log10(w[4]):.2f} dB)"
          f" | ideal 1-bit PO ref (2/pi)^2 = 3.92 dB; scaled by this cell's |R_sw|^2: {po:.2f} dB"
          f" | max over corners: AB absorbed {r[:,2].max():.3f} vs uniform {r[:,:2].max():.3f}")
