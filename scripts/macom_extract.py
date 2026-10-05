"""MACOM measured diode S-parameters -> the RF diode numbers quoted in README "Diode model". Vendor files are not redistributed:
download MA4AGP907_SPAR.zip and MA4AGFCP910_SPAR.zip from https://cdn.macom.com/s-parameters/ and unzip the .s2p files into vendor/macom/.
Each file is a series-mounted diode on a probe fixture. The series branch of a pi network is B of the ABCD matrix, so the fixture's shunt
pad capacitance drops out. OFF: effective capacitance -1 / (w Im B) over 2-8 GHz. ON: R = Re B, L = slope of Im B vs w over 2-8 GHz."""
import glob, os
import numpy as np

R0 = 50.0
HERE = os.path.dirname(os.path.abspath(__file__))


def load(p):
    rows = [list(map(float, l.split())) for l in open(p) if l.strip() and l.lstrip()[0] not in "!#"]
    a = np.array(rows); f = a[:, 0] * 1e9                       # files are '# GHZ S RI R 50'
    S11, S21, S12, S22 = (a[:, k] + 1j * a[:, k + 1] for k in (1, 3, 5, 7))
    return f, R0 * ((1 + S11) * (1 + S22) - S12 * S21) / (2 * S21)


files = sorted(glob.glob(os.path.join(HERE, "..", "vendor", "macom", "*.s2p")))
if not files:
    raise SystemExit("put the MACOM .s2p files (MA4AGP907_SPAR.zip, MA4AGFCP910_SPAR.zip from cdn.macom.com/s-parameters/, unzipped) in vendor/macom/")
for p in files:
    f, B = load(p); m = (f >= 2e9) & (f <= 8e9); w = 2 * np.pi * f[m]; name = os.path.basename(p)
    if "Reverse" in name:
        C = -1 / (w * B[m].imag)
        print(f"{name:40s} OFF: C {C.min() * 1e15:.1f}-{C.max() * 1e15:.1f} fF over 2-8 GHz")
    else:
        L = np.polyfit(w, B[m].imag, 1)[0]
        print(f"{name:40s} ON : R {np.median(B[m].real):.2f} ohm (2-8 GHz {B[m].real.min():.2f}-{B[m].real.max():.2f}), L {L * 1e12:.0f} pH")
