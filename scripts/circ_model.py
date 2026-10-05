"""Palace-calibrated circuit model of one diagonal chain of the X pol-rotation cell (design study, 2026-10-03).
Sheet admittance (normalised to 1/eta0):  y = eta0 * [ jw*C0 + 1 / (jw*La + 1/(jw*Cp) + Rl + n*Z_diode) ]  over the exact
RO4003C + air(G) shorted spacer; r_on = R(Z_on), r_off = R(Z_off); cross = |r_on - r_off|/2 (180 deg exact by symmetry).
Parameters C0, La, Cp, n come from eig_extract.py on a full-wave multiport.
  python circ_model.py closure  : model with params extracted from omp3 s5.49/G10 vs its full-wave losses
  python circ_model.py opt      : worst-case-optimal circuit targets per G (MADP-000907 tolerance box)"""
import sys
import numpy as np
from scipy.optimize import differential_evolution
from eig_extract import y_down

ETA0 = 376.730313
FB = np.array([2.4, 2.44, 2.4835, 5.15, 5.6, 6.1, 6.6, 7.125]) * 1e9
MADP = dict(Ron=5.2, C=0.025e-12, L=0.3e-9, Roff=60.0, Cb=(0.022e-12, 0.030e-12), Lb=(0.2e-9, 0.4e-9))


def refl(f, G, Zd, C0, La, Cp, n, Rl=0.1):
    w = 2 * np.pi * f
    y = ETA0 * (1j * w * C0 + 1 / (1j * w * La + 1 / (1j * w * Cp) + Rl + n * Zd)) + y_down(f, G)
    return (1 - y) / (1 + y)


def loss(f, G, p, d=MADP, worst=True):
    C0, La, Cp, n = p; w = 2 * np.pi * f; out = []
    Cs = d["Cb"] + (d["C"],) if worst else (d["C"],); Ls = d["Lb"] + (d["L"],) if worst else (d["L"],)
    for C in Cs:
        for L in Ls:
            for mm in ((0.95, 1.0, 1.05) if worst else (1.0,)):
                zon = d["Ron"] + 1j * w * L; zoff = d["Roff"] + 1 / (1j * w * C * mm) + 1j * w * L
                out.append(-20 * np.log10(abs(refl(f, G, zon, C0, La, Cp, n) - refl(f, G, zoff, C0, La, Cp, n)) / 2))
    return np.max(out, axis=0)


if __name__ == "__main__":
    if sys.argv[1] == "closure":
        f4 = np.array([2.44, 5.15, 6.1, 7.125]) * 1e9
        p = (0.080e-12, 3.52e-9, 0.458e-12, 1.58)                    # extracted from omp3 s5.49 G10 (C0 mean, n mean)
        print("model nominal", np.round(loss(f4, 10.0, p, worst=False), 2), " full-wave nominal [1.73 1.5 2.0 1.51]")
        print("model worst  ", np.round(loss(f4, 10.0, p), 2), " full-wave robust worst 2.33")
    else:
        for G in (8.0, 9.0, 10.0):
            for nfix in (1.6,):
                obj = lambda q: loss(FB, G, (q[0] * 1e-12, q[1] * 1e-9, q[2] * 1e-12, nfix)).max()
                r = differential_evolution(obj, [(0.01, 0.15), (0.5, 20.0), (0.05, 1.5)], seed=1, tol=1e-7, maxiter=300, polish=True)
                C0, La, Cp = r.x; fr = 1 / (2 * np.pi * np.sqrt(La * 1e-9 * Cp * 1e-12)) / 1e9
                print(f"G {G:4.1f} n {nfix}: worst {r.fun:.2f} dB  C0 {C0:.3f} pF  La {La:.2f} nH  Cp {Cp:.3f} pF  (series res {fr:.2f} GHz)"
                      f"  per-f {np.round(loss(FB, G, (C0 * 1e-12, La * 1e-9, Cp * 1e-12, nfix)), 2)}")
            # what the CURRENT family can do: La, Cp pinned near today's values, only C0 free
            obj2 = lambda q: loss(FB, G, (q[0] * 1e-12, 3.5e-9, 0.46e-12, 1.6)).max()
            r2 = differential_evolution(obj2, [(0.01, 0.15)], seed=1, tol=1e-8)
            print(f"        today's La/Cp (3.5 nH, 0.46 pF), best C0 {r2.x[0]:.3f} pF -> worst {r2.fun:.2f} dB")
