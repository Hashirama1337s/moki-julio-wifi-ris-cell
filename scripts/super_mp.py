"""Two-cell supercell (2 x 1 along x) of the boss design, multiport (2026-10-04: opposite-code neighbour test).
Each cell = the validated single-cell geometry (palace_mp_pol.geometry) shifted to x = -A/2 and +A/2, each with its own boss and
centre-via port. Ports: 2-5 cell 0 diodes (NE, SW, NW, SE), 6-9 cell 1 diodes, 10 / 11 cell 0 / 1 via. Floquet port with MaxOrder 1
so the (-1, 0) order (the beam an alternating A/B code steers to) is recorded.
S (nf, 4 + 10, 1 + 10): rows [TE00, TM00, TE(-1,0), TM(-1,0), b2..b11]; cols [TE in, a2..a11]. Normalisation as palace_mp_pol.
usage (WSL): RIS_G=12 python super_mp.py <tag> <f1,..> [--theta 45 --phi 0] [--direct R2,..,R11]"""
import argparse, json, math, os, re, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import gmsh
from tile_lib import make_config, run_palace, read_csv
from palace_mp_pol import geometry, cplx, ORDER, DIRS

A = 15.0; NX = 2; AX, AY = NX * A, A; G = float(os.environ.get("RIS_G", "12.0")); T = 0.508; R0 = 50.0
BOSS = (4.0, 2.5); VIA = 1.0; XC = [-A / 2, A / 2]


def make_mesh(fn, pec_all, port_all, h0=0.25, gr=0.6, hmax=3.0, hup=20.0):
    zt = round(G + T, 6); zb = round(G - BOSS[1], 6)
    gmsh.initialize(); gmsh.option.setNumber("General.Verbosity", 0); gmsh.model.add("super"); occ = gmsh.model.occ
    vols = [occ.addBox(-AX / 2, -AY / 2, 0, AX, AY, G), occ.addBox(-AX / 2, -AY / 2, G, AX, AY, T), occ.addBox(-AX / 2, -AY / 2, zt, AX, AY, hup)]
    b = BOSS[0]
    cut, _ = occ.cut([(3, vols[0])], [(3, occ.addBox(xc - b / 2, -b / 2, 0, b, b, zb)) for xc in XC]); vols[0] = cut[0][1]

    def surf(P, z):
        pts = [occ.addPoint(x, y, z) for x, y in P]
        return occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(pts[i], pts[(i + 1) % len(pts)]) for i in range(len(pts))])])
    ps = [surf(P, zt) for P in pec_all]
    out, _ = occ.fuse([(2, ps[0])], [(2, x) for x in ps[1:]]); pec_s = [t for d, t in out]
    port_s = [surf(P, zt) for P in port_all]
    via_s = []
    for xc in XC:
        for (dx, dy) in ((VIA / 2, 0), (-VIA / 2, 0), (0, VIA / 2), (0, -VIA / 2)):
            p = [occ.addPoint(xc, 0, zb), occ.addPoint(xc + dx, dy, zb), occ.addPoint(xc + dx, dy, zt), occ.addPoint(xc, 0, zt)]
            via_s.append(occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(p[i], p[(i + 1) % 4]) for i in range(4)])]))
    tools = [(2, t) for t in pec_s] + [(2, t) for t in port_s] + [(2, t) for t in via_s]
    _, outmap = occ.fragment([(3, v) for v in vols], tools); occ.synchronize()
    pieces = outmap[len(vols):]; n0 = len(pec_s); n1 = n0 + len(port_s)
    pec_tags = sorted({t for m in pieces[:n0] for d, t in m if d == 2})
    port_tags = [sorted({t for d, t in m if d == 2}) for m in pieces[n0:n1]]
    via_tags = [sorted({t for d, t in m if d == 2}) for m in pieces[n1:]]
    pec_tags = [t for t in pec_tags if all(t not in pt for pt in port_tags + via_tags)]
    top = zt + hup; EPS = 1e-6
    groups = {k: [] for k in range(1, 8)}; xs, Xs, ys, Ys = [], [], [], []
    for dim, tag in gmsh.model.getEntities(2):
        x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(2, tag); cx, cy, cz = occ.getCenterOfMass(2, tag)
        if x1 < -AX / 2 + EPS: groups[1].append(tag); xs.append((tag, cy, cz))
        elif x0 > AX / 2 - EPS: groups[2].append(tag); Xs.append((tag, cy, cz))
        elif y1 < -AY / 2 + EPS: groups[3].append(tag); ys.append((tag, cx, cz))
        elif y0 > AY / 2 - EPS: groups[4].append(tag); Ys.append((tag, cx, cz))
        elif z1 < EPS: groups[5].append(tag)
        elif abs(z0 - top) < EPS and abs(z1 - top) < EPS: groups[6].append(tag)
    groups[7] = pec_tags
    vt_all = sum(via_tags, [])
    for dim, tag in gmsh.model.getEntities(2):                     # boss walls + tops -> PEC
        x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(2, tag)
        for xc in XC:
            if z1 > EPS and z1 < zb + EPS and max(abs(x0 - xc), abs(x1 - xc), abs(y0), abs(y1)) < b / 2 + EPS and tag not in vt_all and tag not in groups[7]:
                groups[7].append(tag)
    for i, pt in enumerate(port_tags):
        groups[20 + i] = pt
    for i, vt in enumerate(via_tags):
        groups[40 + i] = vt
    for k, v in groups.items():
        if v: gmsh.model.addPhysicalGroup(2, v, k)
    air, sub = [], []
    for dim, tag in gmsh.model.getEntities(3):
        cx, cy, cz = occ.getCenterOfMass(3, tag); (sub if (G < cz < zt) else air).append(tag)
    gmsh.model.addPhysicalGroup(3, air, 1); gmsh.model.addPhysicalGroup(3, sub, 2)
    key = lambda s: (round(s[1], 4), round(s[2], 4))
    gmsh.model.mesh.setPeriodic(2, [s[0] for s in sorted(Xs, key=key)], [s[0] for s in sorted(xs, key=key)], [1, 0, 0, AX, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1])
    gmsh.model.mesh.setPeriodic(2, [s[0] for s in sorted(Ys, key=key)], [s[0] for s in sorted(ys, key=key)], [1, 0, 0, 0, 0, 1, 0, AY, 0, 0, 1, 0, 0, 0, 0, 1])
    f = gmsh.model.mesh.field; f.add("MathEval", 1)
    f.setString(1, "F", f"Min({hmax}, {h0} + {gr}*(((z-{zt})+Abs(z-{zt}))/2 + (({G}-z)+Abs({G}-z))/2))"); f.setAsBackgroundMesh(1)
    for o in ("Mesh.MeshSizeFromPoints", "Mesh.MeshSizeExtendFromBoundary", "Mesh.MeshSizeFromCurvature"):
        gmsh.option.setNumber(o, 0)
    gmsh.model.mesh.generate(3); gmsh.option.setNumber("Mesh.MshFileVersion", 2.2); gmsh.option.setNumber("Mesh.Binary", 1)
    n = len(gmsh.model.mesh.getElementsByType(4)[0]); gmsh.write(fn); gmsh.finalize()
    assert all(len(p) > 0 for p in port_tags) and len(port_tags) == 8 and len(via_tags) == 8
    return n


def config(mesh, freqs, theta, phi, loads=None):
    cfg = make_config(mesh, "x", freqs, theta, phi, "TE", {}, excite=None if loads else {1: 1})
    cfg["Boundaries"]["FloquetPort"][0]["MaxOrder"] = 1
    lp = []
    for c in range(2):
        for i, k in enumerate(ORDER):
            sx, sy = DIRS[k]; idx = 2 + 4 * c + i
            d = {"Index": idx, "Attributes": [20 + 4 * c + i], "Direction": [sx / math.sqrt(2), sy / math.sqrt(2), 0.0]}
            if loads: d["R"] = float(loads[idx])
            else: d["R"] = R0; d["Excitation"] = idx
            lp.append(d)
    for c in range(2):
        idx = 10 + c; d = {"Index": idx, "Elements": [{"Attributes": [40 + 4 * c + j], "Direction": "+Z"} for j in range(4)]}
        if loads: d["R"] = float(loads[idx])
        else: d["R"] = R0; d["Excitation"] = idx
        lp.append(d)
    cfg["Boundaries"]["LumpedPort"] = lp
    if os.environ.get("PALACE_LIN"):
        cfg["Solver"]["Linear"].update(json.loads(os.environ["PALACE_LIN"]))
    return cfg


def order_col(h, d, order, pol, e):
    """complex S for Floquet order string like '0;0' or '-1;0'; NaN if Palace did not write it."""
    norm = lambda s: s.replace(" ", "")
    name = f"S[P1({order}){pol}][{e}]"
    mi = [i for i, x in enumerate(h) if norm(x).startswith("|" + name + "|")]
    if not mi:
        return np.full(d.shape[0], np.nan + 0j)
    ai = [i for i, x in enumerate(h) if norm(x).startswith("arg(" + name + ")")][0]
    return 10 ** (d[:, mi[0]] / 20) * np.exp(1j * np.radians(d[:, ai]))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("tag"); ap.add_argument("freqs")
    ap.add_argument("--theta", type=float, default=0.0); ap.add_argument("--phi", type=float, default=0.0); ap.add_argument("--direct", default="")
    ap.add_argument("--h0", type=float, default=0.25)
    a = ap.parse_args(); fq = [float(x) for x in a.freqs.split(",")]
    os.makedirs("mp", exist_ok=True)
    pec1, ports1 = geometry(2.0, 0.4, 1.275, 0.644, 0.3)
    pec_all = [[(x + xc, y) for x, y in P] for xc in XC for P in pec1]
    port_all = [[(x + xc, y) for x, y in P] for xc in XC for P in ports1]
    os.makedirs("mp", exist_ok=True); mesh = f"mp/{a.tag}.msh"; nt = make_mesh(mesh, pec_all, port_all, h0=a.h0)
    loads = {2 + i: float(r) for i, r in enumerate(a.direct.split(","))} if a.direct else None
    cfg = config(mesh, fq, a.theta, a.phi, loads)
    dt, od, log = run_palace(cfg, "mp/" + a.tag, a.tag)
    fh, fd = read_csv(od + "/port-floquet-S.csv"); f = fd[:, 0]
    orders = sorted({m.group(1) for x in fh for m in [re.search(r"P1\(([-0-9;]+)\)", x.replace(" ", ""))] if m})
    res = dict(tag=a.tag, f=list(map(float, f)), tets=nt, wall=dt, theta=a.theta, phi=a.phi, direct=loads, orders=orders, NX=NX, A=A, G=G, boss=BOSS)
    NP = 10
    if loads:
        for o in ("0;0", "-1;0"):
            for pol in ("TE", "TM"):
                res[f"{pol}_{o}"] = [[float(z.real), float(z.imag)] for z in order_col(fh, fd, o, pol, 1)]
    else:
        vh, vd = read_csv(od + "/port-V.csv")
        S = np.zeros((len(f), 4 + NP, 1 + NP), complex)
        for e in range(1, 2 + NP):
            for r, (o, pol) in enumerate((("0;0", "TE"), ("0;0", "TM"), ("-1;0", "TE"), ("-1;0", "TM"))):
                S[:, r, e - 1] = order_col(fh, fd, o, pol, e)
            vinc = None if e == 1 else vd[:, [i for i, x in enumerate(vh) if x.replace(" ", "").startswith(f"V_inc[{e}][{e}]")][0]]
            for k in range(2, 2 + NP):
                V = cplx(vh, vd, f"V[{k}][{e}]")
                S[:, 2 + k, e - 1] = V / math.sqrt(2 * R0) if e == 1 else (V - (vinc if k == e else 0)) / vinc
        res["S"] = [[[[float(S[i, r, c].real), float(S[i, r, c].imag)] for c in range(1 + NP)] for r in range(4 + NP)] for i in range(len(f))]
    json.dump(res, open(f"mp/{a.tag}.json", "w"))
    print(f"{a.tag}: tets {nt}, {dt:.0f}s, orders {orders}, {'direct' if loads else 'multiport'}", flush=True)
