"""Close the X-cell multiport (from palace_mp_pol.py) with any port loads; validate vs a direct Palace run; evaluate
real diodes in state A (NE+SW ON) and state B (mirror: NW+SE ON), incl. tolerance corners and D-to-D mismatch.
Palace convention exp(+jwt): Z_L = jwL, Z_C = 1/(jwC). Diode: ON = R_on + jwL_s; OFF = R_off + 1/(jwC_off) + jwL_s.
usage: python mp_pol_analyze.py <multiport.json> [--direct direct.json] [--diode MPP4203|MADP008120]"""
import json, sys
import numpy as np

R0 = 50.0
DIODES = {  # parts lane (experts/parts/PARTS.md): nominal + tolerance box; R_off series-equivalent
    "MPP4203": dict(R_on=2.0, C=0.09e-12, L=0.2e-9, R_off=20.0, Cbox=(0.07e-12, 0.10e-12), Lbox=(0.15e-9, 0.30e-9)),
    "MADP008120": dict(R_on=2.2, C=0.12e-12, L=0.6e-9, R_off=67.0, Cbox=(0.105e-12, 0.15e-12), Lbox=(0.5e-9, 0.7e-9)),
}


def load(path):
    d = json.load(open(path)); f = np.array(d["f"]) * 1e9
    S = np.array([[[complex(*c) for c in row] for row in Sf] for Sf in d["S"]]) if "S" in d else None
    return d, f, S


def close(S, ZL):
    """S (nf,6,5); ZL (nf,4) -> (TE, TM) reflected for unit TE incidence."""
    gl = (ZL - R0) / (ZL + R0)
    out_te, out_tm = [], []
    for i in range(S.shape[0]):
        Spp = S[i, 2:, 1:]; SpF = S[i, 2:, 0]; Gm = np.diag(gl[i])
        ap = Gm @ np.linalg.solve(np.eye(4) - Spp @ Gm, SpF)          # a_p = Gamma b_p, b_p = (I - Spp Gamma)^-1 SpF
        out_te.append(S[i, 0, 0] + S[i, 0, 1:] @ ap); out_tm.append(S[i, 1, 0] + S[i, 1, 1:] @ ap)
    return np.array(out_te), np.array(out_tm)


def zd(f, state, d, C=None, L=None):
    w = 2 * np.pi * f; C = d["C"] if C is None else C; L = d["L"] if L is None else L
    return d["R_on"] + 1j * w * L if state == "on" else d["R_off"] + 1 / (1j * w * C) + 1j * w * L


def states(S, f, d, C1=None, C2=None, L=None):
    """state A: ports 2,3 (NE,SW) ON with diodes 'D1'; ports 4,5 OFF with diodes 'D2' (capacitance C1/C2 allow mismatch)."""
    on1, off2 = zd(f, "on", d, C1, L), zd(f, "off", d, C2, L)      # A: D1 on, D2 off
    on2, off1 = zd(f, "on", d, C2, L), zd(f, "off", d, C1, L)      # B: D2 on, D1 off
    _, A = close(S, np.stack([on1, on1, off2, off2], 1))
    _, B = close(S, np.stack([off1, off1, on2, on2], 1))
    return A, B


if __name__ == "__main__":
    d, f, S = load(sys.argv[1])
    if "--direct" in sys.argv:
        dd, f2, _ = load(sys.argv[sys.argv.index("--direct") + 1])
        R = np.array([dd["direct"][str(k)] for k in (2, 3, 4, 5)], float)
        te, tm = close(S, np.tile(R, (len(f), 1)).astype(complex))
        TEd = np.array([complex(*z) for z in dd["TE"]]); TMd = np.array([complex(*z) for z in dd["TM"]])
        print("VALIDATION model vs direct (R loads", R, "):")
        print("  TE dphase", np.round(np.degrees(np.angle(te / TEd)), 2), " d|TE|", np.round(abs(te) - abs(TEd), 4))
        print("  TM dphase", np.round(np.degrees(np.angle(tm / TMd)), 2), " d|TM|", np.round(abs(tm) - abs(TMd), 4))
    name = sys.argv[sys.argv.index("--diode") + 1] if "--diode" in sys.argv else "MPP4203"
    dio = DIODES[name]
    A, B = states(S, f, dio)
    err = 180 - np.abs(np.degrees(np.angle(A / B)))
    print(f"\n{name} nominal: f GHz {np.round(f/1e9,3)}")
    print("  |TM| A", np.round(abs(A), 3), " B", np.round(abs(B), 3), " loss dB", np.round(-20 * np.log10(np.minimum(abs(A), abs(B))), 2))
    print("  180-|dphi|", np.round(err, 2))
    worst_loss, worst_err = 0, 0
    for C in dio["Cbox"] + (dio["C"],):
        for L in dio["Lbox"]:
            for mm in (0.95, 1.0, 1.05):                         # D1 vs D2 capacitance mismatch +-5 % (within a lot)
                A, B = states(S, f, dio, C1=C, C2=C * mm, L=L)
                worst_loss = max(worst_loss, float(np.max(-20 * np.log10(np.minimum(abs(A), abs(B))))))
                worst_err = max(worst_err, float(np.max(180 - np.abs(np.degrees(np.angle(A / B))))))
    print(f"  tolerance box (C {dio['Cbox']}, L {dio['Lbox']}, D1/D2 mismatch +-5%): worst loss {worst_loss:.2f} dB, worst 180-err {worst_err:.2f} deg")
