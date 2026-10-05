"""Palace config / run / CSV helpers (Floquet port + lumped diode ports).  Run on Linux / WSL (numpy only).
mm / GHz.  Palace is exp(+j w t); Gamma = reflected/incident tangential E.  Codes using exp(-i w t)
(e.g. openEMS) -> compare with conj(Gamma_palace).  S11 is referenced at the Floquet port plane z = zt + hup and de-embedded to
the pattern plane z = zt:  Gamma(zt) = S11 * exp(+2 j kz hup),  kz = k0 cos(theta)."""
import json, math, os, re, subprocess, time
import numpy as np

C0 = 299792458.0
FREQS = [2.4, 2.44, 2.4835, 5.15, 5.3969, 5.6438, 5.8906, 6.1375, 6.3844, 6.6312, 6.8781, 7.125]
HERE = os.path.dirname(os.path.abspath(__file__))
PALACE_SH = os.environ.get("PALACE_SH", os.path.join(HERE, "palace.sh"))   # launcher: loads Palace, then exec palace "$@"
NP = int(os.environ.get("PALACE_NP", "6"))
KAPPA = 2 * math.pi * 5e9 * 8.854e-12 * 3.55 * 0.0027      # S/m, RO4003C loss as a conductivity equal to tan d 0.0027 at 5 GHz (used for the 0.508 mm board)


def k0_mm(f_ghz):
    return 2 * math.pi * f_ghz * 1e9 / C0 * 1e-3


def make_config(mesh, out, freqs, theta, phi, pol, loads, hup=20.0, order=2, tol=1e-8, excite=None,
                lossmodel="kappa", maxits=400, extra_solver=None):
    """loads: dict {2: {"R":..,"L":..,"C":..}, 3: {...}} lumped ports on the diode rectangles (attributes 8, 9).
    excite: None -> only the Floquet port excited; or dict {port_index: excitation_index} for multi-excitation."""
    th = math.radians(theta); ph = math.radians(phi)
    kref = k0_mm(1.0)
    kx, ky = kref * math.sin(th) * math.cos(ph), kref * math.sin(th) * math.sin(ph)
    mat_sub = {"Attributes": [2], "Permeability": 1.0, "Permittivity": 3.55}
    if lossmodel == "kappa":
        mat_sub["Conductivity"] = KAPPA
    elif lossmodel == "tand":
        mat_sub["LossTan"] = 0.0027
    bc = {
        "Periodic": {"FloquetWaveVector": [kx, ky, 0.0], "FloquetReferenceFrequency": 1.0,
                     "BoundaryPairs": [{"DonorAttributes": [1], "ReceiverAttributes": [2]},
                                       {"DonorAttributes": [3], "ReceiverAttributes": [4]}]},
        "PEC": {"Attributes": [5, 7]},
        "FloquetPort": [{"Index": 1, "Attributes": [6], "Excitation": (True if excite is None else int(excite.get(1, 0))),
                         "IncidentPolarization": pol, "MaxOrder": 0}],
        "LumpedPort": [],
    }
    for idx, att in ((2, 8), (3, 9)):
        if idx in loads:
            d = dict(Index=idx, Attributes=[att], Direction="+Y")
            d.update(loads[idx])
            if excite is not None:
                d["Excitation"] = int(excite.get(idx, 0))
            bc["LumpedPort"].append(d)
    if not bc["LumpedPort"]:
        del bc["LumpedPort"]
    cfg = {
        "Problem": {"Type": "Driven", "Verbose": 2, "Output": out},
        "Model": {"Mesh": mesh, "L0": 1e-3},
        "Domains": {"Materials": [{"Attributes": [1], "Permeability": 1.0, "Permittivity": 1.0}, mat_sub]},
        "Boundaries": bc,
        "Solver": {"Order": order,
                   "Driven": {"Samples": [{"Type": "Point", "Freq": list(freqs)}]},
                   "Linear": {"Type": "Default", "KSPType": "GMRES", "Tol": tol, "MaxIts": maxits}},
    }
    if extra_solver:
        cfg["Solver"].update(extra_solver)
    return cfg


def run_palace(cfg, workdir, name, np_=None):
    os.makedirs(workdir, exist_ok=True)
    path = os.path.join(workdir, name + ".json")
    cfg = json.loads(json.dumps(cfg))
    cfg["Problem"]["Output"] = os.path.join(workdir, "out_" + name)
    cfg["Model"]["Mesh"] = os.path.abspath(cfg["Model"]["Mesh"])
    with open(path, "w") as fh:
        json.dump(cfg, fh, indent=1)
    log = os.path.join(workdir, name + ".log")
    t0 = time.time()
    cmd = f'export HWLOC_COMPONENTS=-gl; source ~/spack/share/spack/setup-env.sh && spack load palace && palace -np {np_ or NP} "{path}" > "{log}" 2>&1'
    r = subprocess.run(["bash", "-c", cmd])
    dt = time.time() - t0
    if r.returncode != 0:
        raise RuntimeError(f"palace failed, see {log}")
    return dt, cfg["Problem"]["Output"], log


def read_csv(fn):
    with open(fn) as fh:
        header = [h.strip() for h in fh.readline().split(",")]
    data = np.atleast_2d(np.genfromtxt(fn, delimiter=",", skip_header=1))
    return header, data


def floquet_S(outdir, pol_out, exc=1):
    """Complex S of specular order, outgoing polarisation pol_out ('TE'/'TM'), for excitation index exc."""
    h, d = read_csv(os.path.join(outdir, "port-floquet-S.csv"))
    norm = lambda s: s.replace(" ", "")
    pat_m = f"|S[P1(0;0){pol_out}][{exc}]|"; pat_p = f"arg(S[P1(0;0){pol_out}][{exc}])"
    im = [i for i, x in enumerate(h) if norm(x).startswith(pat_m)]
    ip = [i for i, x in enumerate(h) if norm(x).startswith(pat_p)]
    if not im:
        raise KeyError(f"{pat_m} not in {h}")
    return d[:, 0], 10 ** (d[:, im[0]] / 20) * np.exp(1j * np.radians(d[:, ip[0]]))


def deembed(f, S, theta, hup):
    return np.array([s * np.exp(2j * k0_mm(fi) * math.cos(math.radians(theta)) * hup) for fi, s in zip(f, S)])


def mobius(Gs, Zs):
    """Gamma(Z) = (a Z + b)/(c Z + 1) per frequency from three loads. Gs: list of 3 arrays (nf,)"""
    coef = []
    for k in range(len(Gs[0])):
        A = np.array([[Z, 1.0, -Z * G[k]] for Z, G in zip(Zs, Gs)], dtype=complex)
        coef.append(np.linalg.solve(A, np.array([G[k] for G in Gs], dtype=complex)))
    coef = np.array(coef)
    return coef[:, 0], coef[:, 1], coef[:, 2]


def gamma_of(coef, Z):
    a, b, c = coef
    return (a * Z + b) / (c * Z + 1.0)


def diode_Z_palace(freqs, state, d):
    """exp(+j w t): Z_L = +j w L, Z_C = 1/(j w C) = -j/(wC)."""
    w = 2 * np.pi * np.asarray(freqs) * 1e9
    if state == "on":
        return d["R_on"] + 1j * w * d["L_s"]
    return -1j / (w * d["C_off"]) + 1j * w * d["L_s"]


def wrap180(x):
    return (np.asarray(x) + 180.0) % 360.0 - 180.0
