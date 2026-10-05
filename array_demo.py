"""Illustrative 1-bit beam steering with the simulated cell (array factor only, isotropic elements, plane-wave feed).
Each element reflects R_A or R_B (2x2 TE/TM). Coding c_n in {A,B}: steerable part +-R_sw, fixed specular part R_fix.
1-D aperture of N cells along x (15 mm pitch), incidence normal, TE in, observe cross-pol (TM) and co-pol (TE).
Writes figures/fig4_beam_steering.png and prints main-beam / quantisation-lobe / specular-residue levels."""
import os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "scripts")); os.chdir(HERE)
import via_an as v
from mp_pol_analyze import R0

N, P = 32, 15e-3; Lc = 2.5e-9; C0 = 299792458.0
d, f, S = v.prep("data/boss_d0.json")


def states(i):
    w = 2 * np.pi * f[i]; zl = 1j * w * Lc + w * Lc / 30; D = v.D
    on = D["Ron"] + 1j * w * D["L"] + zl; off = D["Roff"] + 1 / (1j * w * D["C"]) + 1j * w * D["L"] + zl
    A = v.te_in(S, i, np.array([on, on, off, off, 1e-3])); B = v.te_in(S, i, np.array([off, off, on, on, 1e-3]))
    return A, B


fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.3), sharey=True)
for ax, fg in zip(axes, (2.44, 5.5, 6.5)):
    i = int(np.argmin(abs(f - fg * 1e9))); A, B = states(i); lam = C0 / f[i]; k = 2 * np.pi / lam
    x = (np.arange(N) - (N - 1) / 2) * P; th = np.radians(np.linspace(-90, 90, 721))
    for ts, c in ((30, "#d95f02"), (-50, "#1b9e77")):
        code = np.cos(k * x * np.sin(np.radians(ts))) >= 0          # 1-bit quantised linear phase
        refl = np.where(code[:, None], A[None, :], B[None, :])         # (N, 2): TE, TM out per element
        AF = lambda comp: abs(np.exp(1j * k * np.outer(np.sin(th), x)) @ refl[:, comp]) / N
        cross, co = AF(1), AF(0)
        ax.plot(np.degrees(th), 20 * np.log10(cross + 1e-9), color=c, lw=1.4, label=f"cross-pol, steer {ts:+d} deg")
        ax.plot(np.degrees(th), 20 * np.log10(co + 1e-9), color=c, lw=0.8, ls=":", label=f"co-pol residue")
        pk = np.degrees(th[np.argmax(cross)])
        print(f"{f[i]/1e9:.2f} GHz steer {ts:+d}: main beam {pk:+.1f} deg at {20*np.log10(cross.max()):.1f} dB (re perfect mirror); "
              f"mirror lobe {20*np.log10(cross[np.argmin(abs(np.degrees(th)+pk))]):.1f} dB; specular co-pol {20*np.log10(co[np.argmin(abs(th))]):.1f} dB")
    ax.set_title(f"{f[i]/1e9:.2f} GHz", fontsize=9); ax.set_xlabel("angle (deg)"); ax.set_ylim(-40, 0); ax.grid(alpha=0.3)
axes[0].set_ylabel("level re a perfect mirror (dB)"); axes[0].legend(fontsize=6, loc="lower left")
fig.suptitle("1-bit steering, 32 cells x 15 mm, normal incidence, array factor only (illustrative)", fontsize=9)
fig.tight_layout(); fig.savefig("figures/fig4_beam_steering.png", dpi=200)
