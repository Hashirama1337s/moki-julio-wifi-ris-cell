"""Fig. 6: worst-case equivalent bit-phase error and loss vs mounting capacitance across each diode, 0 deg, both bands,
earlier v3 shape vs this (release) shape (both with boss + via). Full parts box; series inductor 2.3 nH (generic SRF 20 GHz model)."""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "scripts")); os.chdir(HERE)
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import real_box as rb
from real_all import chip_model
rb.chip = chip_model(20e9)
pads = np.array([0.0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03])
D = [("v3 shape (stub 0.3, gap 0.64 mm, spacer 12 mm)", "data/boss_d0.json", "#d95f02"),
     ("this shape (stub 2.3, gap 1.7 mm, spacer 12.75 mm)", "data/v4_d0.json", "#1b9e77")]
fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.4))
for lab, p, c in D:
    d, f, S = rb.prep(p); m = ((f >= 2.399e9) & (f <= 2.4836e9)) | ((f >= 5.149e9) & (f <= 7.126e9)); ph = []; ls = []
    for cp in pads:
        wv, pe = rb.worst(f, S, 2.3e-9, cp * 1e-12, phase=True); ph.append(pe[m].max()); ls.append(wv[m].max())
    print(lab, "phase", np.round(ph, 1), "loss", np.round(ls, 2))
    ax[0].plot(pads, ph, "o-", color=c, label=lab); ax[1].plot(pads, ls, "o-", color=c, label=lab)
ax[0].axhline(20, color="k", lw=0.8, ls="--"); ax[0].text(0.0, 21, "180 +- 20 deg target", fontsize=7)
ax[1].axhline(2, color="k", lw=0.8, ls="--"); ax[1].text(0.0, 2.05, "2 dB target", fontsize=7)
ax[0].set_xlabel("mounting capacitance across each diode (pF)"); ax[0].set_ylabel("worst equivalent bit-phase error (deg)")
ax[1].set_xlabel("mounting capacitance across each diode (pF)"); ax[1].set_ylabel("worst loss (dB)")
ax[0].legend(fontsize=6.5); ax[0].grid(alpha=0.3); ax[1].grid(alpha=0.3); ax[1].set_ylim(0, 2.5)
fig.suptitle("0 deg, both bands, full parts box, 2.3 nH inductor: phase is far more sensitive to mounting capacitance than loss", fontsize=8.5)
fig.tight_layout(); fig.savefig("figures/fig6_mounting_capacitance.png", dpi=200); print("saved")
