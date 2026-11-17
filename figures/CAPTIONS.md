# Captions - moki-julio-wifi-ris-cell

Both figures are **simulation results** (full-wave field solver + exact network algebra), recomputed from the shipped
network files with the repository's own functions and the release settings of `reproduce.py`. Nothing was built or
measured, and no figure shows hardware or a build.

Run: `python figures/make_figures.py .` from the repository root (Python 3.11, numpy, matplotlib; about 10 s). It prints
the band maxima; they equal the README table rows (0.91 / 0.81, 1.04 / 1.24, 2.21 / 2.28 dB; 27.2 / 30.1 deg worst case
and 22.6 / 19.3 deg typical, straight on).

Settings (from `reproduce.py`): release design v4; 2.3 nH series inductor, chip model with 20 GHz self-resonance
(`real_all.chip_model(20e9)`); MADP-000907 diode box (C_T 0.022 / 0.025 / 0.030 pF, R_s 5.2 / 7.0 ohm, L_s 0.2-0.4 nH,
+-5 % pair mismatch, inductor +-0.1 nH); mounting capacitance 0.013 / 0.015 / 0.017 pF; worst case = maximum over all of
these (`real_box.worst`, then max over the three mounting values).

## 1. `figures/loss_by_angle.svg` (hero)

**Point.** At every simulated angle but one, the worst-case loss stays under the 2 dB budget in both Wi-Fi bands; the
miss is a wave arriving at 45 degrees along the cell's diagonal.

**Caption.** Simulation only: worst-case loss of the designed cell against frequency, for five directions of the incoming
wave, over the whole parts box (diode spread, inductor tolerance, mounting capacitance), with all cells in the same
state. Shaded: the Wi-Fi bands 2.4-2.4835 GHz and 5.15-7.125 GHz; the frequency axis is cut between 2.6 and 5.0 GHz.
Every direction stays under the 2 dB budget in both bands except 45 degrees along the cell diagonal (up to 2.21 dB at
2.4 GHz and 2.28 dB at 5-7 GHz). Three of the 2.4 GHz runs have only two frequencies (dots).

**Alt text.** Line chart of worst-case loss in dB against frequency, in two panels (2.3-2.6 GHz and 5.0-7.3 GHz) with
the Wi-Fi bands shaded and a dashed 2 dB line. Four directions stay between about 0.6 and 1.4 dB inside the bands; the
45-degree diagonal direction rises above 2 dB from about 6.2 GHz and at 2.4835 GHz.

**Sources.**
| series | files (`data/`) | values |
|---|---|---|
| straight on | `v4_d0.json` (101 points, 2.3-7.3 GHz) | 0.91 / 0.81 dB band maxima |
| 30 deg along an edge | `v4_d30.json` (4.8-7.3 GHz), `v4_24_t30.json` (2.4, 2.4835 GHz) | 0.98 / 0.96 dB |
| 45 deg along an edge | `v4_d45.json`, `v4_24_t45.json` | 1.04 / 1.24 dB |
| 30 deg along the diagonal | `v4_p45t30_te/_tm.json`, `v4_d24_p45t30_te/_tm.json` (2.3-2.6 GHz, 31 points) | 1.37 / 1.24 dB |
| 45 deg along the diagonal | `v4_p45_te/_tm.json`, `v4_24_p45.json` + `v4_24_p45_tm.json` | 2.21 / 2.28 dB |
| 2 dB budget, band edges | `README.md`, `reproduce.py` (`band_max` masks) | |

Shown window: 2.3-2.6 GHz and 5.0-7.3 GHz (the 4.8-5.0 GHz part of the simulated range is left out; it lies outside the
Wi-Fi band). The y range is set from the data so no curve is clipped.

## 2. `figures/bit_phase_error.svg`

**Point.** The strict "clean 180-degree flip" target (error 20 degrees or less) is met only straight on, with typical
parts, in the 5-7 GHz band; the worst case over the parts box misses it.

**Caption.** Simulation only: how far the cell's switched phase is from a clean 180-degree flip, for a wave arriving
straight on. With typical parts the error stays under the usual 20-degree target across the 5.15-7.125 GHz band (at most
19.3 degrees) but not across 2.4-2.4835 GHz (up to 22.6 degrees); over the whole parts box the worst case is 27.2 and
30.1 degrees. At a 45-degree tilt the errors are larger (up to 46 and 70 degrees; see TECHNICAL.md).

**Alt text.** Line chart of bit-phase error in degrees against frequency at normal incidence, two panels with the Wi-Fi
bands shaded and a dashed 20-degree target line; the typical-parts curve stays just under 20 degrees in the 5-7 GHz band
while the worst-case curve rises to about 30 degrees.

**Sources.**
| series | file | how |
|---|---|---|
| worst case over the parts box | `data/v4_d0.json` | `real_box.worst(..., phase=True)`, max over mounting 0.013-0.017 pF: 27.2 / 30.1 deg |
| typical parts | `data/v4_d0.json` | the "typical parts, 0 deg" block of `reproduce.py` (C_T 0.025 pF, R_s 5.2 ohm, L_s 0.3 nH, 0.015 pF): 22.6 / 19.3 deg |
| 20-degree target; 46 / 70 deg at 45 deg tilt | `README.md` "What we do not claim" | |
