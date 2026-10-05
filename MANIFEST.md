# Manifest (release 1.0.0)

## Scripts (`scripts/`)

| file | role |
|---|---|
| `palace_mp_pol.py` | geometry + gmsh mesh + Palace config; one run = full multiport (Floquet TE/TM + 4 diode ports [+ via port]). Env: `RIS_G` spacer, `RIS_A` period, `RIS_VIA` via-port width, `RIS_BOSS=b,neck`, `RIS_BIAS=wp,gp` corner bias posts, `PALACE_SWEEP` adaptive dense sweep, `PALACE_LIN` solver options (oblique: `{"PCMatShifted": true}`); `--direct R2,..` = passive-load validation run; `--theta/--phi/--pol` incidence |
| `tile_lib.py` | Palace config / run / CSV helpers |
| `mp_pol_analyze.py` | close the network with any loads; validation vs direct runs |
| `eig_extract.py`, `circ_model.py` | per-diagonal eigen-reflection, sheet admittance, La-Cp-C0 circuit read-out and circuit optimum |
| `sheet_scan.py`, `net6.py` | spacer-height transfer (validated to 1e-3 at normal incidence only) |
| `one_scan.py`, `lever_chip.py`, `robust_real.py` | spacer x inductor scans, tolerance box, mounting parasitics, inductor Q/SRF |
| `angle_an.py`, `te_only.py`, `via_an.py`, `via_joint.py`, `phi_an.py` | oblique incidence (reciprocity + C2 renormalisation), via port, both polarisations |
| `real_box.py`, `real_all.py` | realistic parts: chip inductor with self-resonance and Q(f), diode R_s 5.2/7.0 ohm, mounting capacitance; all angles |
| `super_mp.py`, `super_an.py`, `sup_energy.py` | 2-cell supercell (AA / BB / AB codes): build + run, network closure with realistic parts, energy accounting and steered-beam efficiency |
| `vendor_chip.py`, `coilcraft_to_csv.py` | measured Coilcraft inductor model; Touchstone .s2p -> series-impedance table |
| `bias_real.py` | bias posts with realistic parts (stacked 33 + 10 nH chip chokes, ideal choke, open) vs the cell without posts |
| `extra_checks.py` | prints the README text numbers not covered by reproduce.py (rotated-channel phase, validations, supercell checks, Foster check, loss-free switches, panel leak, inductor choice, option A band edge) |
| `macom_extract.py` | RF diode model numbers from MACOM's measured S-parameter files (download them yourself; see vendor/README.txt) |
| `shape_exact.py` | copper-shape scoring at normal incidence: exact 6-port network, analytic spacer transfer, per-diode mismatch, mounting capacitance and diode C_T box as options |
| `bias_an.py` | corner bias posts (ports 7-10): choke models (ideal, resistive, real chips, stacked chips) |

## Data (`data/`, Palace network JSON: frequencies, geometry, S matrix)

The design in this release is the `v4_*` shape (stub 2.3 mm, diode gap 1.7 mm, spacer 12.75 mm, boss + via). The `boss_*`, `sup_*`, `bias_*`
files are the earlier v3 shape (stub 0.3 mm, gap 0.644 mm, spacer 12 mm), kept for the comparison table and the sections run only on it.

| file | content |
|---|---|
| `v4_d0.json` | this shape (release design): 0 deg, 101 points 2.3-7.3 GHz |
| `v4_d30.json` | this shape (release design): 30 deg along the edge (phi 0), 51 points 4.8-7.3 GHz |
| `v4_d45.json` | this shape (release design): 45 deg along the edge, 51 points 4.8-7.3 GHz |
| `v4_p45t30_te.json` | this shape (release design): 30 deg along the diagonal (phi 45), TE multiport, 51 points |
| `v4_p45t30_tm.json` | this shape (release design): 30 deg along the diagonal, TM-in direct run, 51 points |
| `v4_p45_te.json` | this shape (release design): 45 deg along the diagonal, TE multiport, 51 points |
| `v4_p45_tm.json` | this shape (release design): 45 deg along the diagonal, TM-in direct run, 51 points |
| `v4_24_t30.json` | this shape (release design): 2.4 / 2.4835 GHz, 30 deg along the edge |
| `v4_24_t45.json` | this shape (release design): 2.4 / 2.4835 GHz, 45 deg along the edge |
| `v4_24_p45.json` | this shape (release design): 2.4 / 2.4835 GHz, 45 deg along the diagonal, TE multiport |
| `v4_24_p45_tm.json` | this shape (release design): 2.4 / 2.4835 GHz, 45 deg along the diagonal, TM-in direct |
| `v4_d24_p45t30_te.json` | this shape (release design): 2.3-2.6 GHz (31 points), 30 deg along the diagonal, TE multiport |
| `v4_d24_p45t30_tm.json` | this shape (release design): 2.3-2.6 GHz, 30 deg along the diagonal, TM-in direct |
| `boss_d0.json` | boss design, 0 deg, 101 points 2.3-7.3 GHz |
| `boss_d30.json`, `boss_d45.json` | boss design, 30 / 45 deg (phi 0), 51 points 4.8-7.3 GHz |
| `boss24_t30.json`, `boss24_t45.json` | boss design, 2.4 / 2.4835 GHz at 30 / 45 deg |
| `boss_n.json` | boss design, 0 deg, 4 points |
| `via_d45_g12*.json`, `via_d30_g12.json` | full-height via (before the boss), 45 / 30 deg |
| `dense_s2_g12_n.json` | no via / no boss, 0 deg, 101 points |
| `v12_s4.json`, `vdir_s4.json` | earlier design, 12 frequencies + direct-run validation |
| `ang_t45_p0_*.json` | 45 deg network + direct validation runs (TE and TM in) |
| `boss_p45_te/tm.json`, `boss24_p45*.json`, `boss24_p45t30_*.json` | boss design on the diagonal plane: 45 deg dense, 2.4 GHz at 30/45 deg (TE multiport + TM direct) |
| `bossval_t45_te.json` | direct validation run of the boss network at 45 deg |
| `bossfine_n.json` | mesh-convergence run (h0 0.18) |
| `bossn15_d45.json`, `boss6_d45.json`, `boss6_n.json` | boss variants (neck 1.5 mm, boss 6 mm) |
| `sup_n.json`, `sup_t45.json` | 2-cell supercell, 0 deg (4 points) and 45 deg phi 0 (6.0 / 6.5 / 7.1 GHz); 14-port rows TE00, TM00, TE(-1,0), TM(-1,0) + 10 lumped |
| `explore/shape_*.json` | copper-shape search runs at normal incidence (no boss, spacer 12 mm), scored with `shape_exact.py` |
| `explore/omp3_*.json` | earlier no-boss design runs used by the exploration scripts (spacer / inductor scans, circuit read-out) |
| `boss_p45t30_te.json`, `boss_p45t30_tm.json` | v3 shape, 30 deg along the diagonal (TE multiport + TM-in direct), 51 points |
| `d24_p45t30_te.json`, `d24_p45t30_tm.json` | v3 shape, 2.3-2.6 GHz dense (31 points), 30 deg along the diagonal (the `d24_` prefix predates the `boss_` naming) |
| `bias_n.json`, `bias_t45.json` | boss design + corner bias posts, 0 / 45 deg |

## Reproduce the headline table

`python reproduce.py` (top level) regenerates both README results tables (this shape and the v3 shape under the same rules), the
typical-parts rows, the lower-capacitance diode option, the network-vs-direct validation and the mesh check from `data/`.
Figures: `python make_figs.py` (figs 1-3), `python array_demo.py` (fig 4), `python make_fig5.py` (fig 5), `python make_fig6.py` (fig 6).
Individual scripts (run from `scripts/`; data paths are `../data/...`):

```
python via_an.py ../data/boss_d0.json 2.5      # 0 deg, via short vs open (identical at 0 deg); ideal Q30 inductor, nominal diodes
python via_an.py ../data/boss_d45.json 2.5     # 45 deg along the edge
python mp_pol_analyze.py ../data/v12_s4.json --direct ../data/vdir_s4.json   # network vs direct run (earlier design)
python angle_an.py 45 0 2.25 --val             # 45-deg validation (earlier no-boss design)
python real_all.py --srf 20 --L 2.0 --cp 0 0.02   # loss-only re-score (2-point 2.4 GHz files; the README table itself is reproduce.py)
python sup_energy.py ../data/sup_t45.json      # supercell energy accounting + steered beam
python bias_real.py ../data/bias_t45.json ../data/boss_d45.json 2.0 2.2   # bias posts with realistic parts, 45 deg (README bias section)
python bias_real.py ../data/bias_n.json ../data/boss_n.json 2.0 2.2       # same, 0 deg
python bias_an.py ../data/bias_n.json 2.5        # simple choke models (open, ideal 47 nH, 1 kohm feed, 22 nH chip), 0 deg
python bias_an.py ../data/bias_t45.json 2.5      # same, 45 deg
python sup_energy.py ../data/sup_n.json        # supercell absorption at 0 deg
python super_an.py ../data/sup_t45.json --single ../data/boss_d45.json   # supercell vs single cell (uniform codes)
python via_an.py ../data/bossn15_d45.json 2.5   # boss neck 1.5 mm (vs boss_d45: neck 2.5 mm)
python via_an.py ../data/boss6_d45.json 2.5     # 6 mm boss
python via_an.py ../data/boss_n.json 2.5        # mesh check, ideal inductor (compare bossfine_n.json)
python via_an.py ../data/bossfine_n.json 2.5
python via_an.py ../data/via_d45_g12.json 2.5   # full-height via (v3): via short vs open, the centre-charging notch
python extra_checks.py                         # the remaining README numbers (rotated-channel phase, loss-free absorption, R_s trend, ...)
python shape_exact.py ../data/explore/shape_ph_st23g170.json --pad 0.015 --padtol 0.002 --L 1.0:4.0:0.1 --G 9:13.5:0.25   # shape scoring (no boss, analytic spacer; suggests 2.6 nH, the full FEM with the boss prefers 2.3-2.4 nH: extra_checks 4g)
python shape_exact.py ../data/explore/shape_ph_st23g210.json --pad 0 --padtol 0 --L 1.0:4.0:0.1 --G 9:13.5:0.25   # without mounting capacitance (~18 deg)
python macom_extract.py                        # needs the MACOM .s2p files in vendor/macom/
python net6.py                                 # spacer transfer validated against a real run at another spacer (|dR| <= 1e-3)
python super_an.py ../data/sup_n.json --single ../data/boss_n.json     # supercell vs single cell at 0 deg
python vendor_chip.py                          # this shape scored with the measured Coilcraft inductor (needs vendor/coilcraft_02DS-2N3_Zseries.csv)
```

## Software (tested versions)

- Analysis (everything except the two solver drivers): Python 3.11.9, numpy 2.4.4, matplotlib 3.11.0 (`requirements.txt`). Windows or Linux.
- Field solver (`palace_mp_pol.py`, `super_mp.py`): Palace 0.18.1 (Spack, Ubuntu 24.04 under WSL2), gmsh 4.15.2, numpy 2.5.3, Python 3.12.3.
  `scripts/palace.sh` is the launcher (adapt it, or point `PALACE_SH` at your own). New runs write to `mp/` in the working directory.
- Memory: the 2-cell supercell at 45 deg peaked at 12.4 GB (16 GB machine); a 4-cell supercell would not fit.

## Release checks

- Two independent reviews of the README (final verdict: ready after minor fixes, all applied) and a clean-room audit: a fresh copy of this
  folder runs every command in this manifest, regenerates all figures byte-identical, and prints every computed number quoted in
  README.md (datasheet values, dimensions and literature figures are quoted inputs).
- Every reference opened at the primary source (two author/title errors in draft references were corrected before release).
- No absolute paths, no personal names, no vendor data files in the package (vendor files are rebuilt from the user's own downloads).
- The release shape was verified as full FEM runs at every angle (no analytic transfer behind any oblique or headline number).
