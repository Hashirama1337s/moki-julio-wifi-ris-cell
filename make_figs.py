"""Figures 1-3 (fig 1: release shape; figs 2-3: v3 shape). usage:  python make_figs.py"""
import os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "scripts")); os.chdir(HERE)
import via_an as v
import types
sys.modules.setdefault("gmsh", types.ModuleType("gmsh"))
sys.modules.setdefault("tile_lib", types.SimpleNamespace(make_config=0, run_palace=0, floquet_S=0, read_csv=0, deembed=0))
import palace_mp_pol as geo_mod   # geometry() only; gmsh/tile_lib stubbed (Windows side)

OUT = "figures"
BANDS = [(2.4, 2.4835), (5.15, 7.125)]


def shade(ax):
    for a, b in BANDS:
        ax.axvspan(a, b, color="#d8ecd8", zorder=0)
    ax.axhline(2.0, color="#c0392b", ls="--", lw=1); ax.text(7.32, 2.0, "spec 2 dB", va="center", fontsize=7, color="#c0392b")


def fig_geometry():
    pec, ports = geo_mod.geometry(2.0, 0.4, 1.275, 1.7, 2.3)        # release shape: stub 2.3 mm, diode gap 1.7 mm
    fig, ax = plt.subplots(figsize=(4.2, 4.2))
    for P in pec:
        xs, ys = zip(*(P + [P[0]])); ax.fill(xs, ys, color="#b87333")
    for i, P in enumerate(ports):
        xs, ys = zip(*(P + [P[0]])); ax.fill(xs, ys, color="#2c7fb8" if i < 2 else "#9ecae1")
    ax.add_patch(plt.Rectangle((-7.5, -7.5), 15, 15, fill=False, lw=1, ls=":"))
    ax.add_patch(plt.Rectangle((-2, -2), 4, 4, fill=False, lw=1, ls="--", color="k"))
    ax.set_aspect("equal"); ax.set_xlim(-7.8, 7.8); ax.set_ylim(-7.8, 7.8)
    ax.set_title("Cell (15 mm): copper; diode + inductor gaps ON (dark) / OFF (light),\nstate A shown; dashed = ground boss below (4 mm)", fontsize=8)
    ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)")
    fig.tight_layout(); fig.savefig(f"{OUT}/fig1_geometry.png", dpi=200); plt.close(fig)


def curve(path, Lc, via="short"):
    d, f, S = v.prep(path); return f / 1e9, v.worst(f, S, Lc * 1e-9, via)


def fig_loss():
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    runs = [("data/boss_d0.json", "0 deg", "#1b9e77"), ("data/boss_d45.json", "45 deg (plane along edge)", "#d95f02")]
    if os.path.exists("data/boss_d30.json"):
        runs.insert(1, ("data/boss_d30.json", "30 deg", "#7570b3"))
    for p, lab, c in runs:
        f, wv = curve(p, 2.5); ax.plot(f, wv, color=c, lw=1.6, label=lab)
    f, wv = curve("data/via_d45_g12.json", 2.5); ax.plot(f, wv, color="#999999", lw=1, ls=":", label="45 deg, full-height via (before)")
    shade(ax); ax.set_ylim(0, 6); ax.set_xlim(2.3, 7.3)
    ax.set_xlabel("frequency (GHz)"); ax.set_ylabel("worst-case loss (dB)")
    ax.set_title("v3 boss design, Lc 2.5 nH, G 12 mm: worst case over the diode box (C_T, L_s, mismatch)", fontsize=9)
    ax.legend(fontsize=7, loc="upper left"); fig.tight_layout(); fig.savefig(f"{OUT}/fig2_loss_vs_angle.png", dpi=200); plt.close(fig)


def fig_fix():
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    for p, via, lab, c in (("data/via_d45_g12.json", "open", "no via", "#e41a1c"), ("data/via_d45_g12.json", "short", "full-height via (12.5 mm post)", "#ff7f00"),
                           ("data/boss_d45.json", "short", "raised boss + 2.5 mm via neck", "#1b9e77")):
        f, wv = curve(p, 2.5, via); ax.plot(f, wv, color=c, lw=1.6, label=lab)
    shade(ax); ax.set_xlim(4.8, 7.3); ax.set_ylim(0, 11)
    ax.set_xlabel("frequency (GHz)"); ax.set_ylabel("worst-case loss (dB)")
    ax.set_title("45 deg incidence: the centre-charging mode and its fix", fontsize=9); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f"{OUT}/fig3_angle_fix.png", dpi=200); plt.close(fig)


if __name__ == "__main__":
    fig_geometry(); fig_loss(); fig_fix(); print("wrote fig1-fig3; figures/ now holds:", sorted(os.listdir(OUT)))
