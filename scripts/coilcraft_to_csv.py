"""Build vendor/coilcraft_02DS-2N3_Zseries.csv (used by vendor_chip.py) from the manufacturer's measured S-parameter file.
The vendor file is not redistributed here: download the 0201DS S-parameter archive from the Coilcraft 0201DS product page,
extract 02DS-2N3.s2p, then:   python coilcraft_to_csv.py path/to/02DS-2N3.s2p
Series impedance of a 2-port series element: Z = 2 * R0 * (1 - S21) / S21."""
import os, sys
import numpy as np

UNIT = {"HZ": 1e-9, "KHZ": 1e-6, "MHZ": 1e-3, "GHZ": 1.0}


def read_s2p(path):
    unit, fmt, r0, rows = 1e-9, "MA", 50.0, []
    for line in open(path, encoding="latin-1"):
        line = line.split("!")[0].strip()
        if not line:
            continue
        if line.startswith("#"):
            tok = line[1:].upper().split()
            unit = next((UNIT[t] for t in tok if t in UNIT), unit)
            fmt = next((t for t in tok if t in ("MA", "DB", "RI")), fmt)
            if "R" in tok:
                r0 = float(tok[tok.index("R") + 1])
            continue
        rows.append([float(x) for x in line.split()])
    a = np.array(rows); f = a[:, 0] * unit; p, q = a[:, 3], a[:, 4]      # S21 = columns 3, 4 (order S11 S21 S12 S22)
    s21 = {"MA": p * np.exp(1j * np.radians(q)), "DB": 10 ** (p / 20) * np.exp(1j * np.radians(q)), "RI": p + 1j * q}[fmt]
    return f, s21, r0


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit("usage: python coilcraft_to_csv.py path/to/02DS-2N3.s2p")
    f, s21, r0 = read_s2p(sys.argv[1]); z = 2 * r0 * (1 - s21) / s21
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "vendor", "coilcraft_02DS-2N3_Zseries.csv")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    hdr = ("Coilcraft 0201DS 2.3 nH (02DS-2N3.s2p, coilcraft.com 0201DS product page)\n"
           "series impedance Z = 2*R0*(1-S21)/S21 from the measured 2-port S-parameters; columns: f_GHz, ReZ_ohm, ImZ_ohm")
    np.savetxt(out, np.c_[f, z.real, z.imag], delimiter=",", fmt="%.3f", header=hdr)
    print("wrote", out, len(f), "points")
