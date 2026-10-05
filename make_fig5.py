"""fig5: worst-case loss per incidence case for three series-inductor models (ideal / real high-SRF / common SRF-10 GHz part)."""
import os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "scripts")); os.chdir(HERE)
import real_box as rb
from real_all import CASES, chip_model
from phi_an import net
models = [("ideal 2.5 nH, Q30, Rs 5.2", lambda w, L, Q=14, srf=None: 1j * w * L + w * L / 30, 2.5, (5.2,)),
          ("real 2.0 nH, SRF 20 GHz, Rs 5.2+7", chip_model(20e9), 2.0, (5.2, 7.0)),
          ("common 2.0 nH, SRF 10 GHz, Rs 5.2+7", chip_model(10e9), 2.0, (5.2, 7.0))]
data = []
for name, te, tm in CASES:
    if tm: d, f, S = net(te, tm); full = True
    else: d, f, S = rb.prep(te); full = False
    data.append((name, f, S, full))
vals = np.zeros((len(models), len(data)))
for i, (lab, chip, L, rs) in enumerate(models):
    rb.chip = chip
    for j, (name, f, S, full) in enumerate(data):
        wv = rb.worst(f, S, L * 1e-9, 0.0, rs_list=rs, full=full)
        b = ((f >= 2.399e9) & (f <= 2.4836e9)) | ((f >= 5.149e9) & (f <= 7.126e9)); vals[i, j] = wv[b].max()
    print(lab, np.round(vals[i], 2), flush=True)
fig, ax = plt.subplots(figsize=(8.5, 3.4)); x = np.arange(len(data)); wd = 0.27
for i, (lab, *_r) in enumerate(models):
    ax.bar(x + (i - 1) * wd, vals[i], wd, label=lab, color=["#1b9e77", "#7570b3", "#d95f02"][i])
ax.axhline(2.0, color="#c0392b", ls="--", lw=1); ax.set_xticks(x); ax.set_xticklabels([d[0].replace("2.4@", "2.4 GHz @") for d in data], fontsize=7, rotation=20)
ax.set_ylabel("worst-case loss (dB)"); ax.set_ylim(0, 6); ax.legend(fontsize=7)
ax.set_title("Boss design: worst-case loss per case (e = edge plane, d = diagonal plane; 5-7 GHz unless marked 2.4 GHz)", fontsize=8)
fig.tight_layout(); fig.savefig("figures/fig5_inductor_requirement.png", dpi=200)
