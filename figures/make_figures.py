"""Figures for the Wi-Fi RIS cell repository (simulation results only), recomputed from the shipped network files.

usage:  python figures/make_figures.py [repo_root]      (repo_root defaults to the folder above this script; about 1-2 minutes)
writes: figures/loss_by_angle.svg, figures/bit_phase_error.svg  (next to this script)

The numbers come from the same functions and settings as reproduce.py (release design v4: 2.3 nH series inductor with a
20 GHz self-resonance, MADP-000907 diode box, mounting capacitance 0.013 / 0.015 / 0.017 pF, worst case over all of them).
The script prints the band maxima so they can be compared with the README table.
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else os.path.dirname(HERE)
OUT = HERE if os.path.basename(HERE) == "figures" else os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import real_box as rb                     # noqa: E402  (repository code)
from real_all import chip_model           # noqa: E402
from phi_an import net                    # noqa: E402

# ---- lab style (STYLE.md) -------------------------------------------------------------------------------------------
BG, BORDER, GRID, AX, TXT, ACC = "#1c2026", "#2f353d", "#22272e", "#9aa4ae", "#e8edf2", "#39ff14"
CAT = ["#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#0ca30c", "#9085e9", "#e66767"]
NEUTRAL = "#5b636c"
plt.rcParams.update({
    "svg.fonttype": "none", "svg.hashsalt": "moki-julio-wifi-ris-cell",
    "font.family": "sans-serif", "font.sans-serif": ["Segoe UI", "Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 14, "figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG,
    "axes.edgecolor": AX, "axes.labelcolor": TXT, "text.color": TXT, "xtick.color": AX, "ytick.color": AX,
    "axes.titlecolor": TXT, "grid.color": GRID, "grid.linewidth": 0.8, "lines.linewidth": 2,
})
W_IN = 12.5


def frame(fig):
    fig.patches.append(FancyBboxPatch((0.004, 0.004), 0.992, 0.992, boxstyle="round,pad=0,rounding_size=0.012",
                                      transform=fig.transFigure, fill=False, ec=BORDER, lw=1.5, zorder=-1))


def save(fig, name):
    fig.savefig(os.path.join(OUT, name), format="svg", metadata={"Date": None, "Creator": None})
    plt.close(fig)
    print("wrote figures/" + name)


def style_ax(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(True); ax.set_axisbelow(True)


# ---- the release settings, exactly as in reproduce.py -----------------------------------------------------------------
rb.chip = chip_model(20e9)
D = os.path.join(ROOT, "data")
PADS = (0.013e-12, 0.015e-12, 0.017e-12)
MADP = dict(ct=(0.022e-12, 0.030e-12, 0.025e-12), rs=(5.2, 7.0))
L_REL = 2.3e-9
BANDS = ((2.399e9, 2.4836e9), (5.149e9, 7.126e9))      # reproduce.py band masks
WIFI = ((2.4, 2.4835), (5.15, 7.125))                  # GHz, shading
XL1, XL2 = (2.3, 2.6), (5.0, 7.3)                       # GHz shown (simulated: 2.3-2.6 and 4.8-7.3)
# (label, 5-7 GHz files (TE[, TM]), 2.4 GHz files (TE[, TM]))
CASES = [("straight on (0 deg)", ("v4_d0", None), None),
         ("30 deg, along an edge", ("v4_d30", None), ("v4_24_t30", None)),
         ("45 deg, along an edge", ("v4_d45", None), ("v4_24_t45", None)),
         ("30 deg, along the diagonal", ("v4_p45t30_te", "v4_p45t30_tm"), ("v4_d24_p45t30_te", "v4_d24_p45t30_tm")),
         ("45 deg, along the diagonal", ("v4_p45_te", "v4_p45_tm"), ("v4_24_p45", "v4_24_p45_tm"))]


def worst_curve(te, tm):
    pte = os.path.join(D, te + ".json")
    d, f, S = net(pte, os.path.join(D, tm + ".json")) if tm else rb.prep(pte)
    lw = np.zeros(len(f)); pw = np.zeros(len(f))
    for cp in PADS:
        wv, pe = rb.worst(f, S, L_REL, cp, rs_list=MADP["rs"], full=bool(tm), phase=True, ct=MADP["ct"])
        lw = np.maximum(lw, wv); pw = np.maximum(pw, pe)
    return f, lw, pw


def typical_phase_0deg():
    """typical parts, straight on: the last block of reproduce.py (v4 design line)."""
    d, f, S = rb.prep(os.path.join(D, "v4_d0.json")); w = 2 * np.pi * f; zv = np.full(len(f), 1e-3 + 0j)
    cp = 0.015e-12; C0 = 0.025e-12; z = rb.chip(w, L_REL)
    dio = lambda x: 1 / (1 / x + 1j * w * cp)
    on = dio(5.2 + 1j * w * 0.3e-9) + z; off = dio(60 + 1 / (1j * w * C0) + 1j * w * 0.3e-9) + z; R = []
    for Z in (np.stack([on, on, off, off, zv], 1), np.stack([off, off, on, on, zv], 1)):
        g = (Z - rb.R0) / (Z + rb.R0)
        R.append(np.array([S[i, :2, 0] + S[i, :2, 1:] @ (np.diag(g[i]) @ np.linalg.solve(np.eye(5) - S[i, 2:, 1:] @ np.diag(g[i]), S[i, 2:, 0]))
                           for i in range(len(f))]))
    A, B = R
    return f, np.degrees(2 * np.arctan(np.linalg.norm((A + B) / 2, axis=1) / np.linalg.norm((A - B) / 2, axis=1)))


def band_max(f, v):
    return [float(v[(f >= a) & (f <= b)].max()) if ((f >= a) & (f <= b)).any() else None for a, b in BANDS]


def shown_max(curves, j):
    """largest plotted value inside the shown frequency windows (so nothing is clipped)."""
    m = 0.0
    for v in curves.values():
        g = v[0] / 1e9; y = v[j]
        w = ((g >= XL1[0]) & (g <= XL1[1])) | ((g >= XL2[0]) & (g <= XL2[1]))
        m = max(m, float(y[w].max()))
    return m


def broken_axes(fig, rect, ylab):
    """two panels on one shared y scale: 2.3-2.6 GHz and 4.8-7.3 GHz (the frequencies that were simulated)."""
    gs = fig.add_gridspec(1, 2, width_ratios=[0.5, 2.3], wspace=0.05, left=rect[0], right=rect[1], bottom=rect[2], top=rect[3])
    a1 = fig.add_subplot(gs[0]); a2 = fig.add_subplot(gs[1], sharey=a1)
    a1.set_xlim(*XL1); a2.set_xlim(*XL2)
    a1.set_xticks([2.4, 2.5]); a2.set_xticks([5, 5.5, 6, 6.5, 7])
    for a in (a1, a2):
        style_ax(a)
    a2.spines["left"].set_visible(False); a2.tick_params(axis="y", left=False, labelleft=False)
    a1.set_ylabel(ylab)
    for a, (lo, hi) in ((a1, WIFI[0]), (a2, WIFI[1])):
        a.axvspan(lo, hi, color=NEUTRAL, alpha=0.18, lw=0, zorder=0)
    # break marks
    for a, x, dx in ((a1, 1.0, 0.05), (a2, 0.0, 0.006)):        # "//" marks where the frequency axis is cut
        for off in (-1.2 * dx, 0.0):
            a.plot([x + off, x + off + dx], [-0.015, 0.015], transform=a.transAxes, color=AX, lw=1.2, clip_on=False)
    fig.text((rect[0] + rect[1]) / 2, rect[2] - 0.085, "frequency (GHz), axis cut between 2.6 and 5.0   ·   shaded: the Wi-Fi bands",
             ha="center", fontsize=13.5, color=TXT)
    return a1, a2


def fig_loss(curves):
    fig = plt.figure(figsize=(W_IN, 7.4))
    frame(fig)
    fig.text(0.03, 0.93, "Loss stays under the 2 dB budget at every simulated angle but one", fontsize=21, weight="bold")
    fig.text(0.03, 0.88, "Simulation only. Worst case over the whole parts box (diode spread, inductor tolerance, mounting "
                         "capacitance); all cells in the same state.", fontsize=13, color=AX)
    a1, a2 = broken_axes(fig, (0.08, 0.70, 0.17, 0.80), "worst-case loss (dB)")
    for k, (lab, hi, lo) in enumerate(CASES):
        c = CAT[k]
        runs = [(a1, curves[(k, "hi")]), (a2, curves[(k, "hi")])] + ([(a1, curves[(k, "lo")])] if lo else [])
        for a, (f, loss, _) in runs:
            sparse = len(f) <= 6                  # the 2.4 GHz runs at two frequencies only: dots, no joining line
            a.plot(f / 1e9, loss, color=c, lw=0 if sparse else 2, marker="o" if sparse else None, ms=9)
    for a in (a1, a2):
        a.axhline(2.0, color=TXT, lw=1.4, ls=(0, (5, 4)))
    a1.set_ylim(0, 1.08 * shown_max(curves, 1))
    # legend as direct labels on the right, in the order of the series
    hs = [Line2D([], [], color=CAT[k], lw=2.5, marker="o", ms=7) for k in range(len(CASES))]
    hs.append(Line2D([], [], color=TXT, lw=1.4, ls=(0, (5, 4))))
    fig.legend(hs, [c[0] for c in CASES] + ["2 dB loss budget"], loc="center left", bbox_to_anchor=(0.715, 0.50), frameon=False, fontsize=13,
               title="wave arrives", title_fontsize=13.5, alignment="left", labelspacing=0.9)
    fig.text(0.725, 0.24, "dots: runs at only two frequencies\n(2.4 and 2.4835 GHz). Lines join\nthe simulated frequencies.", fontsize=12, color=AX, va="top")
    save(fig, "loss_by_angle.svg")


def fig_phase(curves, ft, pt):
    fig = plt.figure(figsize=(W_IN, 7.0))
    frame(fig)
    fig.text(0.03, 0.93, "The 180-degree flip: strict target met only straight on, with typical parts, at 5-7 GHz",
             fontsize=19.5, weight="bold")
    fig.text(0.03, 0.88, "Simulation only. Error of the switched phase from a clean 180 degrees, wave arriving straight on (0 deg).",
             fontsize=13, color=AX)
    a1, a2 = broken_axes(fig, (0.08, 0.70, 0.17, 0.80), "bit-phase error (degrees)")
    f, _, pw = curves[(0, "hi")]
    for a in (a1, a2):
        a.plot(f / 1e9, pw, color=CAT[0], lw=2)
        a.plot(ft / 1e9, pt, color=CAT[1], lw=2)
        a.axhline(20.0, color=TXT, lw=1.4, ls=(0, (5, 4)))
    a1.set_ylim(0, 1.08 * shown_max({0: (f, pw), 1: (ft, pt)}, 1))
    hs = [Line2D([], [], color=CAT[0], lw=2.5), Line2D([], [], color=CAT[1], lw=2.5),
          Line2D([], [], color=TXT, lw=1.4, ls=(0, (5, 4)))]
    fig.legend(hs, ["worst case over the parts box", "typical parts", "target: 20 degrees or less"], loc="center left", bbox_to_anchor=(0.715, 0.55),
               frameon=False, fontsize=13, labelspacing=0.9)
    save(fig, "bit_phase_error.svg")


if __name__ == "__main__":
    curves = {}
    for k, (lab, hi, lo) in enumerate(CASES):
        curves[(k, "hi")] = worst_curve(*hi)
        if lo:
            curves[(k, "lo")] = worst_curve(*lo)
        f = np.concatenate([curves[(k, s)][0] for s in ("hi", "lo") if (k, s) in curves])
        L = np.concatenate([curves[(k, s)][1] for s in ("hi", "lo") if (k, s) in curves])
        P = np.concatenate([curves[(k, s)][2] for s in ("hi", "lo") if (k, s) in curves])
        lm, pm = band_max(f, L), band_max(f, P)
        print(f"  {lab:28s} loss 2.4 / 5-7: {lm[0]:.2f} / {lm[1]:.2f} dB   phase {pm[0]:.1f} / {pm[1]:.1f} deg")
    ft, pt = typical_phase_0deg()
    q = band_max(ft, pt)
    print(f"  0 deg typical parts: phase {q[0]:.1f} / {q[1]:.1f} deg")
    fig_loss(curves)
    fig_phase(curves, ft, pt)
