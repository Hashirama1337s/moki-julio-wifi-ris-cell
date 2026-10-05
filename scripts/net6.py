"""Full 6-port network of the X cell at normal incidence from one multiport run (2026-10-03): ports [TE, TM, NE, SW, NW, SE].
The TM-in column comes from reciprocity + the diagonal mirror (checked on the TE column). Spacer change G0 -> G1 is an
isotropic shunt stub change at the sheet plane (exact for a thin sheet): y' = y + (y_down(G1) - y_down(G0)) on the two
polarisation ports. Works for ANY loads (asymmetric too, e.g. one diode per diagonal)."""
import numpy as np
from mp_pol_analyze import load, R0
from eig_extract import y_down, C0, HUP
import os
DATA = os.environ.get("RIS_DATA", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"))   # network JSON files


def s6(path):
    d, f, S = load(path); nf = len(f); S6 = np.zeros((nf, 6, 6), complex)
    S6[:, :, 0] = S[:, :, 0]; S6[:, :, 2:] = S[:, :, 1:]
    S6[:, 0, 1] = S[:, 1, 0]; S6[:, 1, 1] = S[:, 0, 0]          # mirror x<->y: TM-in Floquet block
    S6[:, 2:, 1] = S[:, 1, 1:]                                   # reciprocity: S[port, TM in] = S[TM out, port]
    rec = np.max(abs(S[:, 2:, 0] - S[:, 0, 1:]))                  # check on TE: S[port, TE in] vs S[TE out, port]
    return d, f, S6, rec


def regap(S6, f, G0, G1):
    k = 2 * np.pi * f / C0 * 1e-3; out = np.zeros_like(S6); I = np.eye(6)
    for i in range(len(f)):
        Dm = np.diag([np.exp(1j * k[i] * HUP)] * 2 + [1] * 4); s = Dm @ S6[i] @ Dm
        y = (I - s) @ np.linalg.inv(I + s)
        dy = y_down(f[i:i + 1], G1)[0] - y_down(f[i:i + 1], G0)[0]; y[0, 0] += dy; y[1, 1] += dy
        out[i] = np.linalg.inv(Dm) @ ((I - y) @ np.linalg.inv(I + y)) @ np.linalg.inv(Dm)   # back to the port plane
    return out


def close6(S6, ZL):
    """ZL (nf,4) on ports NE,SW,NW,SE -> 2x2 reflection matrix (nf,2,2) in the [TE,TM] basis."""
    gl = (ZL - R0) / (ZL + R0); R = np.zeros((S6.shape[0], 2, 2), complex)
    for i in range(S6.shape[0]):
        Gm = np.diag(gl[i]); Spp = S6[i, 2:, 2:]
        R[i] = S6[i, :2, :2] + S6[i, :2, 2:] @ Gm @ np.linalg.solve(np.eye(4) - Spp @ Gm, S6[i, 2:, :2])
    return R


if __name__ == "__main__":
    a = f"{DATA}/explore/omp3_5p891_0p400_1p275_0p644_0p300_10p000.json"; b = f"{DATA}/explore/omp3_5p891_0p400_1p275_0p644_0p300_9p000.json"
    d, f, Sa, rec = s6(a); _, _, Sb, rb = s6(b); print("reciprocity residual (TE column):", rec, rb)
    w = 2 * np.pi * f; on = 5.2 + 1j * w * 0.3e-9; off = 60 + 1 / (1j * w * 0.025e-12) + 1j * w * 0.3e-9; sh = 0 * on
    for name, Z in (("two/diag", np.stack([on, on, off, off], 1)), ("one/diag", np.stack([on, sh, sh, off], 1))):
        Rp = close6(regap(Sa, f, 10.0, 9.0), Z); Rr = close6(regap(Sb, f, 9.0, 9.0), Z)
        print(name, "G10->9 predicted |cross|", np.round(abs(Rp[:, 1, 0]), 4), " real G9", np.round(abs(Rr[:, 1, 0]), 4),
              " max |dR|", np.round(np.max(abs(Rp - Rr)), 4))
