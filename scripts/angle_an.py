"""Oblique-incidence analysis of the X cell (2026-10-03). Inputs per angle: ang_t<th>_p<ph>_mp.json (TE multiport) and
ang_..._tm.json (TM in, ports on 50 ohm, raw port V). Port<-Floquet coupling is normalised by reciprocity + C2
(S[port k, pol in] = c * S[pol out, port rot180(k)]; c printed with its spread over ports = self-check).
For any loads: 2x2 reflection R (TE/TM). States A/B: R_A = R_fix + R_sw, R_B = R_fix - R_sw; the switched part flips
sign EXACTLY for any angle/polarisation; its size is the steerable efficiency, R_fix e is the specular residue.
usage: python angle_an.py <theta> <phi> [Lc_nH] [--val]"""
import json, os, sys
import numpy as np


def wmean(k, a):
    """reciprocity factor averaged with weights |S|^2: ports not excited by symmetry (e.g. phi 45) carry ~0 and garbage ratios"""
    wt = abs(a) ** 2; return (k * wt).sum(1) / wt.sum(1)
from mp_pol_analyze import load, R0
DATA = os.environ.get("RIS_DATA", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"))   # network JSON files

ROT = [1, 0, 3, 2]                                     # rot180: NE<->SW, NW<->SE
D = dict(Ron=5.2, C=0.025e-12, L=0.3e-9, Roff=60.0, Cb=(0.022e-12, 0.030e-12), Lb=(0.2e-9, 0.4e-9)); Q = 30.0


def net(th, ph):
    d, f, S = load(f"{DATA}/ang_t{th}_p{ph}_mp.json"); t = json.load(open(f"{DATA}/ang_t{th}_p{ph}_tm.json")); nf = len(f)
    c = lambda a: np.array([complex(*z) for z in a])
    Vtm = np.stack([c(v) for v in t["V"]], 1) / np.sqrt(2 * R0)             # (nf,4)
    S6 = np.zeros((nf, 6, 6), complex); S6[:, :2, 0] = S[:, :2, 0]; S6[:, :, 2:] = S[:, :, 1:]
    S6[:, 0, 1] = c(t["TE"]); S6[:, 1, 1] = c(t["TM"])
    kte = S[:, 0, 1:][:, ROT] / S[:, 2:, 0]; ktm = S[:, 1, 1:][:, ROT] / Vtm
    # the C2 partner relation carries a -1 (in-plane vectors flip under rot180): physical factor = -k (TE checked: -k = +1.000)
    sg = float(os.environ.get("TM_SIGN", "1"))
    S6[:, 2:, 0] = S[:, 2:, 0] * (-wmean(kte, S[:, 2:, 0]))[:, None]; S6[:, 2:, 1] = Vtm * (-sg * wmean(ktm, Vtm))[:, None]
    info = (np.round(kte.mean(1), 3), np.max(abs(kte - kte.mean(1)[:, None])), np.round(ktm.mean(1), 3), np.max(abs(ktm - ktm.mean(1)[:, None])))
    return f, S6, info


def y_stub(f, G, th, pol):
    """normalised (to the mode impedance) admittance looking down from the sheet at oblique incidence: RO4003C over air G."""
    k = 2 * np.pi * f / 299792458.0 * 1e-3; s2 = np.sin(np.radians(th)) ** 2; kz0 = k * np.sqrt(1 - s2); kzs = k * np.sqrt(3.55 - s2)
    zs = kz0 / kzs if pol == "TE" else (kzs / 3.55) / kz0
    za = 1j * np.tan(kz0 * G); zd = zs * (za + 1j * zs * np.tan(kzs * 0.508)) / (zs + 1j * za * np.tan(kzs * 0.508))
    return 1 / zd


def regap(S6, f, th, G0, G1):
    k = 2 * np.pi * f / 299792458.0 * 1e-3 * np.cos(np.radians(th)); out = np.zeros_like(S6); I = np.eye(6)
    for i in range(len(f)):
        Dm = np.diag([np.exp(1j * k[i] * 20.0)] * 2 + [1] * 4); s = Dm @ S6[i] @ Dm
        y = (I - s) @ np.linalg.inv(I + s)
        for j, pol in enumerate(("TE", "TM")):
            y[j, j] += y_stub(f[i:i + 1], G1, th, pol)[0] - y_stub(f[i:i + 1], G0, th, pol)[0]
        out[i] = np.linalg.inv(Dm) @ ((I - y) @ np.linalg.inv(I + y)) @ np.linalg.inv(Dm)   # back to the port plane
    return out


def close6(S6, ZL):
    gl = (ZL - R0) / (ZL + R0); R = np.zeros((S6.shape[0], 2, 2), complex)
    for i in range(S6.shape[0]):
        Gm = np.diag(gl[i]); R[i] = S6[i, :2, :2] + S6[i, :2, 2:] @ Gm @ np.linalg.solve(np.eye(4) - S6[i, 2:, 2:] @ Gm, S6[i, 2:, :2])
    return R


def metrics(f, S6, Lc):
    w = 2 * np.pi * f; acc = {k: [] for k in ("TE", "TM", "resTE", "resTM", "xTE", "xTM", "perr")}
    for C in D["Cb"] + (D["C"],):
        for L in D["Lb"] + (D["L"],):
            for dl in ((-0.1e-9, 0.0, 0.1e-9) if Lc > 0 else (0.0,)):
                for mm in (0.95, 1.0, 1.05):
                    lc = Lc + dl; zl = 1j * w * lc + w * lc / Q
                    on = D["Ron"] + 1j * w * L + zl; o1 = D["Roff"] + 1 / (1j * w * C) + 1j * w * L + zl; o2 = D["Roff"] + 1 / (1j * w * C * mm) + 1j * w * L + zl
                    A = close6(S6, np.stack([on, on, o2, o2], 1)); B = close6(S6, np.stack([o1, o1, on, on], 1))
                    sw, fx = (A - B) / 2, (A + B) / 2
                    for j, p in enumerate(("TE", "TM")):
                        acc[p].append(-20 * np.log10(np.linalg.norm(sw[:, :, j], axis=1)))        # steerable loss, incident pol p
                        acc["res" + p].append(20 * np.log10(np.linalg.norm(fx[:, :, j], axis=1)))  # fixed specular residue
                        acc["x" + p].append(-20 * np.log10(np.minimum(abs(A[:, 1 - j, j]), abs(B[:, 1 - j, j]))))  # plain cross-pol
                    acc["perr"].append(180 - abs(np.degrees(np.angle(A[:, 1, 0] / B[:, 1, 0]))))
    return {k: np.max(v, axis=0) for k, v in acc.items()}


if __name__ == "__main__":
    th, ph = sys.argv[1], sys.argv[2]; Lc = float(sys.argv[3]) * 1e-9 if len(sys.argv) > 3 and not sys.argv[3].startswith("-") else 2.25e-9
    f, S6, info = net(th, ph)
    print(f"theta {th} phi {ph}  f {np.round(f / 1e9, 3)}  Lc {Lc * 1e9:.2f} nH")
    print(f"  reciprocity factor TE {info[0]} spread {info[1]:.1e} | TM {info[2]} spread {info[3]:.1e}")
    if "--val" in sys.argv:
        for pol, j in (("te", 0), ("tm", 1)):
            v = json.load(open(f"{DATA}/ang_t{th}_p{ph}_v{pol}.json")); c = lambda a: np.array([complex(*z) for z in a])
            Rl = np.array([v["direct"][str(k)] for k in (2, 3, 4, 5)], float); R = close6(S6, np.tile(Rl, (len(f), 1)).astype(complex))
            print(f"  VALIDATION {pol.upper()} in, R loads {Rl}: max |model - direct| TE out {np.max(abs(R[:, 0, j] - c(v['TE']))):.4f}  TM out {np.max(abs(R[:, 1, j] - c(v['TM']))):.4f}")
    m = metrics(f, S6, Lc)
    print("  worst-case steerable loss dB: TE in", np.round(m["TE"], 2), " TM in", np.round(m["TM"], 2))
    print("  plain cross-pol loss dB     : TE in", np.round(m["xTE"], 2), " TM in", np.round(m["xTM"], 2))
    print("  specular residue dB (max)   : TE in", np.round(m["resTE"], 1), " TM in", np.round(m["resTM"], 1))
    print("  cross-pol A/B phase error deg (max):", np.round(m["perr"], 2))
    best = min(((max(metrics(f, S6, l * 1e-9)[k].max() for k in ("TE", "TM")), l) for l in np.arange(0.5, 4.01, 0.25)))
    print(f"  best Lc at this angle: {best[1]:.2f} nH -> worst steerable loss {best[0]:.2f} dB")
