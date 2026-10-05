"""Every README number not printed by reproduce.py, from the shipped data. Run from scripts/: python extra_checks.py
1. rotated-channel A-vs-B phase in the edge plane (this shape and the v3 shape)
2. v3 boss design at 45 deg: dense reduced-order network vs a direct Palace run with passive loads
3. 2-cell supercell (v3 shape): loss-free absorption, AB absorption vs diode R_s, AB vs BA agreement
4. prose numbers: inductor self-resonance factor, Rayleigh frequency, Foster check, loss-free switches at 45 deg diagonal, panel leak,
   dense-vs-direct in dB, series-inductor choice for this shape, option A band-edge resonance"""
import os, json
import numpy as np
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data") + "/"
import real_box as rb
from real_all import chip_model
from mp_pol_analyze import R0
import super_an as sa

rb.chip = chip_model(20e9)

print("1. rotated-channel phase (edge plane; cross-pol A vs B, full box): max |arg(-A/B)| in band; v4 = this shape (2.3 nH, mounting 0.015 pF), boss = v3 shape (2.0 nH)")
for name, p in (("v4 0 deg", "v4_d0"), ("v4 30 deg edge", "v4_d30"), ("v4 45 deg edge", "v4_d45"), ("v4 2.4 GHz 30 deg edge", "v4_24_t30"),
                ("v4 2.4 GHz 45 deg edge", "v4_24_t45"), ("v3 0 deg", "boss_d0"), ("v3 30 deg edge", "boss_d30"), ("v3 45 deg edge", "boss_d45"),
                ("v3 2.4 GHz 30 deg edge", "boss24_t30"), ("v3 2.4 GHz 45 deg edge", "boss24_t45")):
    Lc, cpad = (2.3e-9, 0.015e-12) if p.startswith("v4") else (2.0e-9, 0.0)
    d, f, S = rb.prep(D + p + ".json"); w = 2 * np.pi * f; zv = np.full(len(f), 1e-3 + 0j); err = np.zeros(len(f))
    for Rs in (5.2, 7.0):
        for C in (0.022e-12, 0.025e-12, 0.030e-12):
            for Ls in (0.2e-9, 0.3e-9, 0.4e-9):
                for dl in (-0.1e-9, 0.0, 0.1e-9):
                    for mm in (0.95, 1.0, 1.05):
                        z = rb.chip(w, Lc + dl); dio = (lambda x: 1 / (1 / x + 1j * w * cpad)) if cpad else (lambda x: x)
                        on = dio(Rs + 1j * w * Ls) + z
                        o1 = dio(60 + 1 / (1j * w * C) + 1j * w * Ls) + z; o2 = dio(60 + 1 / (1j * w * C * mm) + 1j * w * Ls) + z
                        R = []
                        for Z in (np.stack([on, on, o2, o2, zv], 1), np.stack([o1, o1, on, on, zv], 1)):
                            g = (Z - R0) / (Z + R0)
                            R.append(np.array([S[i, 1, 0] + S[i, 1, 1:] @ (np.diag(g[i]) @ np.linalg.solve(np.eye(5) - S[i, 2:, 1:] @ np.diag(g[i]), S[i, 2:, 0]))
                                               for i in range(len(f))]))
                        err = np.maximum(err, abs(np.degrees(np.angle(-R[0] / R[1]))))
    m = ((f >= 2.399e9) & (f <= 2.4836e9)) | ((f >= 5.149e9) & (f <= 7.126e9))
    print(f"   {name:26s} {err[m].max():.2f} deg")

print("2. boss design, 45 deg: network (dense adaptive sweep) closed with the direct run's passive loads vs the direct run")
v = json.load(open(D + "bossval_t45_te.json")); d, f, S = rb.prep(D + "boss_d45.json")
R = np.array([v["direct"][k] for k in ("2", "3", "4", "5", "6")], float)
for fi, te, tm in zip(v["f"], v["TE"], v["TM"]):
    i = int(np.argmin(abs(f - fi * (1e9 if fi < 100 else 1)))); g = (R - R0) / (R + R0)
    out = S[i, :2, 0] + S[i, :2, 1:] @ (np.diag(g) @ np.linalg.solve(np.eye(5) - S[i, 2:, 1:] @ np.diag(g), S[i, 2:, 0]))
    print(f"   {f[i]/1e9:.2f} GHz: |dTE| {abs(out[0] - complex(*te)):.4f}  |dTM| {abs(out[1] - complex(*tm)):.4f}  (|TM| {abs(complex(*tm)):.2f})")

print("3. 2-cell supercell")
real = sa.chip
for p in ("sup_n", "sup_t45"):
    d, f, S = sa.prep(D + p + ".json")
    for i in range(len(f)):
        w = 2 * np.pi * f[i]; sa.chip = lambda w_, L: 1j * w_ * L
        on = 1j * w * 0.3e-9 + 1j * w * 2e-9; off = 1 / (1j * w * 0.025e-12) + 1j * w * 0.3e-9 + 1j * w * 2e-9
        zAB = np.array([on, on, off, off, off, off, on, on, 1e-3, 1e-3], complex); zAA = np.array([on, on, off, off] * 2 + [1e-3, 1e-3], complex)
        rAB = sa.close(S, i, zAB); rAA = sa.close(S, i, zAA); sa.chip = real
        ab = sa.close(S, i, sa.loads(f, i, "AB", 5.2, 0.025e-12, 0.3e-9, 2e-9)); ba = sa.close(S, i, sa.loads(f, i, "BA", 5.2, 0.025e-12, 0.3e-9, 2e-9))
        print(f"   {p:8s} {f[i]/1e9:.2f} GHz: loss-free parts absorb AA {1 - np.sum(abs(rAA[:2]) ** 2):.3f}, AB {1 - np.nansum(abs(rAB) ** 2):.3f} | "
              f"AB vs BA max |dR| {np.nanmax(abs(abs(ab) - abs(ba))):.4f}")
    if p == "sup_t45":
        i = int(np.argmin(abs(f - 6.5e9)))
        for Rs in (5.2, 2.0, 0.5):
            sa.chip = lambda w_, L: 1j * w_ * L
            r = sa.close(S, i, sa.loads(f, i, "AB", Rs, 0.025e-12, 0.3e-9, 2e-9)); sa.chip = real
            print(f"   AB at 6.5 GHz, 45 deg, ideal inductor, typical diodes, R_s {Rs}: absorbed {1 - np.nansum(abs(r) ** 2):.3f}")

# ---- 4. remaining README prose numbers (added for release 1.0.0) ----
from phi_an import net
from mp_pol_analyze import load as _load


def _close5(S, Z):
    nf = len(Z); g = (Z - R0) / (Z + R0); G = np.zeros((nf, 5, 5), complex); G[:, range(5), range(5)] = g
    if S.shape[2] == 7:
        X = np.linalg.solve(np.eye(5) - S[:, 2:, 2:] @ G, S[:, 2:, :2]); return S[:, :2, :2] + S[:, :2, 2:] @ G @ X
    b = np.linalg.solve(np.eye(5) - S[:, 2:, 1:] @ G, S[:, 2:, 0][:, :, None]); return S[:, :2, 0][:, :, None] + S[:, :2, 1:] @ G @ b


def _eq(A, B):
    return np.degrees(2 * np.arctan(np.linalg.norm((A + B) / 2, axis=1) / np.linalg.norm((A - B) / 2, axis=1)))


print("4a. single-pole chip inductor at its 10 GHz minimum self-resonance: L_eff / L at 7 GHz =", round(1 / (1 - (7 / 10) ** 2), 2))
print("4b. Rayleigh (grazing) frequency of the (-1,0) order, 30 mm supercell, 45 deg: %.2f GHz" % (2.998e8 / (0.03 * (1 + np.sin(np.radians(45)))) / 1e9))

print("4c. Foster check (v3 shape, 0 deg, 8-corner box): best equivalent phase error with an IDEAL series reactance X chosen per frequency")
d, f, S = rb.prep(D + "boss_d0.json"); m = (f >= 5.149e9) & (f <= 7.126e9); f, S = f[m], S[m]; w = 2 * np.pi * f; zv = np.full(len(f), 1e-3 + 0j)
X = np.arange(-1500, 1501, 2.0); floor = np.full(len(f), 1e9); xbest = np.zeros(len(f))
for x in X:
    pe = np.zeros(len(f)); ls = np.zeros(len(f))
    for Rs in (5.2, 7.0):
        for C in (0.022e-12, 0.030e-12):
            for Ls in (0.2e-9, 0.4e-9):
                on = Rs + 1j * w * Ls + 1j * x; off = 60 + 1 / (1j * w * C) + 1j * w * Ls + 1j * x
                A = _close5(S, np.stack([on, on, off, off, zv], 1)); B = _close5(S, np.stack([off, off, on, on, zv], 1))
                pe = np.maximum(pe, _eq(A, B).max(1)); ls = np.maximum(ls, (-20 * np.log10(np.linalg.norm((A - B) / 2, axis=1))).max(1))
    pe = np.where(ls > 2.0, 1e9, pe); upd = pe < floor; floor[upd] = pe[upd]; xbest[upd] = x
print(f"    floor {floor.min():.1f}-{floor.max():.1f} deg over 5.15-7.125 GHz; needed X: " +
      " ".join(f"{f[i]/1e9:.2f}:{xbest[i]:.0f}" for i in range(0, len(f), 4)) + " ohm (a passive reactance must rise with frequency)")

print("4d. ideal switches (ON = jwL, OFF = open), this shape, 45 deg along the diagonal, best L: worst equivalent phase error")
d, f, S = net(D + "v4_p45_te.json", D + "v4_p45_tm.json"); m = (f >= 5.149e9) & (f <= 7.126e9); f, S = f[m], S[m]; w = 2 * np.pi * f; zv = np.full(len(f), 1e-3 + 0j)
best = min((_eq(_close5(S, np.stack([1j * w * L * 1e-9 + 1e-6, 1j * w * L * 1e-9 + 1e-6, np.full(len(f), 1e9 + 0j), np.full(len(f), 1e9 + 0j), zv], 1)),
                _close5(S, np.stack([np.full(len(f), 1e9 + 0j), np.full(len(f), 1e9 + 0j), 1j * w * L * 1e-9 + 1e-6, 1j * w * L * 1e-9 + 1e-6, zv], 1))).max(), L)
           for L in np.arange(0.5, 6.01, 0.1))
print(f"    {best[0]:.1f} deg at L {best[1]:.1f} nH.  Sheet scaling at 45 deg: TE 1/cos = {1/np.cos(np.radians(45)):.2f}, TM cos = {np.cos(np.radians(45)):.2f}")

print("4e. panel leak (this shape, 0 deg, typical parts, 2.3 nH, 0.015 pF): 32-cell array factor, 1-bit code steering to 40 deg, E along a diagonal")
d, f, S = rb.prep(D + "v4_d0.json"); x = (np.arange(32) - 15.5) * 15e-3; th = np.radians(np.linspace(-89, 89, 1781)); worst_change = 0
for fg in (2.44e9, 5.5e9, 6.5e9, 7.0e9):
    i = int(np.argmin(abs(f - fg))); ww = 2 * np.pi * f[i]; z = rb.chip(np.array([ww]), 2.3e-9)[0]; cp = 0.015e-12
    dio = lambda q: 1 / (1 / q + 1j * ww * cp); on = dio(5.2 + 1j * ww * 0.3e-9) + z; off = dio(60 + 1 / (1j * ww * 0.025e-12) + 1j * ww * 0.3e-9) + z
    A = _close5(S[i:i + 1], np.array([[on, on, off, off, 1e-3]]))[0, :, 0]; B = _close5(S[i:i + 1], np.array([[off, off, on, on, 1e-3]]))[0, :, 0]
    ru, rv = A[0] + A[1], A[0] - A[1]; fix, sw = (ru + rv) / 2, (ru - rv) / 2; k = 2 * np.pi * f[i] / 2.998e8
    code = np.where(np.cos(k * x * np.sin(np.radians(40))) >= 0, 1.0, -1.0); E = np.exp(1j * k * np.outer(np.sin(th), x))
    m40 = np.argmin(abs(np.degrees(th) - 40)); with_leak = abs(E @ (fix + code * sw))[m40] / 32; no_leak = abs(E @ (code * sw))[m40] / 32
    spec = abs(E @ np.full(32, fix))[np.argmin(abs(th))] / 32; worst_change = max(worst_change, abs(20 * np.log10(with_leak / no_leak)))
    e = np.degrees(2 * np.arctan(abs(fix) / abs(sw)))
    print(f"    {f[i]/1e9:.2f} GHz: eq. error {e:4.1f} deg, specular lobe {20*np.log10(spec/no_leak):+.1f} dB re steered beam")
print(f"    steered beam with vs without the leak: max change {worst_change:.2f} dB")
for e in (20, 28, 30):
    print(f"    rule of thumb, equal-amplitude bit: eq. error {e} deg -> specular lobe {20*np.log10(np.tan(np.radians(e/2))*np.pi/2):.1f} dB re steered beam")

print("4f. boss design at 45 deg, dense network vs direct run: |dTM| 0.0101 on |TM| 0.44 = %.2f dB" % (20 * np.log10(1 + 0.0101 / 0.44)))

print("4g. series inductor for this shape (full FEM with the boss, 0 deg, mounting 0.015 +- 0.002 pF, full box): worst equivalent phase error")
d, f, S = rb.prep(D + "v4_d0.json"); mb = ((f >= 2.399e9) & (f <= 2.4836e9)) | ((f >= 5.149e9) & (f <= 7.126e9))
for L in (2.0, 2.2, 2.3, 2.4, 2.6):
    pw = np.zeros(len(f))
    for cp in (0.013e-12, 0.015e-12, 0.017e-12):
        pw = np.maximum(pw, rb.worst(f, S, L * 1e-9, cp, phase=True)[1])
    print(f"    {L:.1f} nH: {pw[mb].max():.1f} deg")

print("4h. option A (MA4AGFCP910, 1.9 nH), 45 deg along the diagonal, loss at the lower band edge vs mid-band")
d, f, S = net(D + "v4_p45_te.json", D + "v4_p45_tm.json"); lw = np.zeros(len(f))
for cp in (0.013e-12, 0.015e-12, 0.017e-12):
    lw = np.maximum(lw, rb.worst(f, S, 1.9e-9, cp, rs_list=(5.2, 6.0), full=True, ct=(0.016e-12, 0.021e-12, 0.018e-12)))
print("    " + " ".join(f"{f[k]/1e9:.2f}:{lw[k]:.2f}" for k in range(len(f)) if 5.1e9 <= f[k] <= 5.6e9) + " dB")

print("4i. generic inductor Q model 14 sqrt(f / 0.5 GHz) at 1.7 GHz = %.1f vs 0201DS datasheet Q 64 at 1.7 GHz (quoted input): %.0f %%" % (14 * np.sqrt(1.7 / 0.5), 100 * 14 * np.sqrt(1.7 / 0.5) / 64))
