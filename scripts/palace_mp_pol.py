"""Exact-diode engine for the switchable X pol-rotation cell (Palace, normal incidence, TE in).
Geometry 'xsw2': axis-aligned centre pad (cp) + 4 diagonal stubs (width w, length st) + gap g + arms to corner pads (s).
Each gap holds a lumped port (attrs 8-11, R = 50 ohm, own excitation 2-5, direction along its diagonal). One Palace run
(excitations: 1 = Floquet TE, 2-5 = each port) gives the network S: rows [TE out, TM out, b2..b5], cols [TE in, a2..a5],
in Palace's normalisation (lumped: S_ke = V_k / V_inc,e, S_ee = (V_e - V_inc)/V_inc; Floquet->lumped: V_k / sqrt(R)).
Any diode (series R-L-C, impossible as a Palace lumped port, which is parallel RLC) is then exact algebra.
usage (WSL): python palace_mp_pol.py <tag> <f1,..> --cap s --w w --cpad cp --gap g [--st 0.3]
             [--direct R2,R3,R4,R5]  (validation: same geometry, ports as passive R loads, Floquet excitation only)"""
import argparse, json, math, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import gmsh
from tile_lib import make_config, run_palace, floquet_S, read_csv, deembed

BOSS = [float(v) for v in os.environ.get("RIS_BOSS", "0,0").split(",")]   # raised ground boss under the centre: side b (mm), neck height hn (mm, air gap boss top -> substrate)
VIA = float(os.environ.get("RIS_VIA", "0"))      # centre-via port width (mm); 0 = no via port
A = float(os.environ.get("RIS_A", "15.0")); G = float(os.environ.get("RIS_G", "8.0")); T = 0.508; R0 = 50.0   # RIS_G: air-spacer height probe (mm)
DIRS = {"NE": (1, 1), "SW": (-1, -1), "NW": (-1, 1), "SE": (1, -1)}
ORDER = ["NE", "SW", "NW", "SE"]           # ports 2,3,4,5 (state A: NE+SW ON, NW+SE OFF; state B = mirror)


def rot(pts, ang):
    c, s = math.cos(ang), math.sin(ang); return [(c * x - s * y, s * x + c * y) for x, y in pts]


def rect(cx, cy, lx, ly):
    return [(cx - lx / 2, cy - ly / 2), (cx + lx / 2, cy - ly / 2), (cx + lx / 2, cy + ly / 2), (cx - lx / 2, cy + ly / 2)]


def geometry(s, w, cp, g, st, ring=0.0, taper=0.0):
    """ring > 0: corner pads become square frames of that width (concave); taper > 0: arm width grows linearly from w at the
    diode to `taper` at the pad. Both keep full D4 symmetry (each shape is symmetric about its own diagonal)."""
    c = A / 2 - 0.1 - s / 2; inner = c - s / 2; wt = taper if taper > 0 else w   # RIS_A: cell period (mm)
    assert ring == 0 or (ring >= 0.2 and s - 2 * ring >= 0.2), "ring width out of range"
    assert inner >= cp / 2 + 0.2 and g >= 0.2 and w >= 0.2, "geometry violates minimum gaps"
    assert cp / math.sqrt(2) + st + g + 0.2 <= inner * math.sqrt(2), "diode gap reaches the corner pad"
    pec = [rect(0, 0, cp, cp)]; ports = []
    d_pad = cp / 2                                     # stub starts inside the pad (overlap along the diagonal)
    d_st = cp / math.sqrt(2) + st                      # stub end
    for k in ORDER:
        sx, sy = DIRS[k]; ang = math.atan2(sy, sx)
        pec.append(rot(rect((d_pad + d_st) / 2, 0, d_st - d_pad, w), ang))          # stub
        ports.append(rot(rect(d_st + g / 2, 0, g, w), ang))                          # gap = port
        d1 = (c if ring == 0 else inner + ring / 2) * math.sqrt(2); d0 = d_st + g
        pec.append(rot([(d0, -w / 2), (d1, -wt / 2), (d1, wt / 2), (d0, w / 2)], ang))     # arm (tapered if wt != w)
        if ring == 0:
            pec.append(rect(sx * c, sy * c, s, s))                                  # corner pad
        else:                                                                       # concave frame pad
            X, Y, o = sx * c, sy * c, s / 2 - ring / 2
            pec += [rect(X, Y - o, s, ring), rect(X, Y + o, s, ring), rect(X - o, Y, ring, s), rect(X + o, Y, ring, s)]
    ext = max(max(abs(x), abs(y)) for P in pec + ports for x, y in P)
    assert ext <= A / 2 - 0.1 + 1e-9
    return pec, ports


def make_mesh(fn, pec, ports, h0=0.25, gr=0.6, hmax=3.0, hup=20.0):
    zt = round(G + T, 6)
    gmsh.initialize(); gmsh.option.setNumber("General.Verbosity", 0); gmsh.model.add("xsw2"); occ = gmsh.model.occ
    vols = [occ.addBox(-A / 2, -A / 2, 0, A, A, G), occ.addBox(-A / 2, -A / 2, G, A, A, T), occ.addBox(-A / 2, -A / 2, zt, A, A, hup)]
    zb = 0.0
    if BOSS[0] > 0:                                      # metal boss = hole in the air spacer, its faces become PEC (group 7)
        zb = round(G - BOSS[1], 6); b = BOSS[0]
        cut, _ = occ.cut([(3, vols[0])], [(3, occ.addBox(-b / 2, -b / 2, 0, b, b, zb))]); vols[0] = cut[0][1]

    def surf(P):
        pts = [occ.addPoint(x, y, zt) for x, y in P]
        return occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(pts[i], pts[(i + 1) % len(pts)]) for i in range(len(pts))])])
    ps = [surf(P) for P in pec]
    out, _ = occ.fuse([(2, ps[0])], [(2, x) for x in ps[1:]]); pec_s = [t for d, t in out]
    port_s = [surf(P) for P in ports]
    via_s = []
    if VIA > 0:                                          # centre via as a port: 4 vertical half-strips (C4-symmetric cross), ground -> centre pad
        for (x0, y0, x1, y1) in ((0, 0, VIA / 2, 0), (0, 0, -VIA / 2, 0), (0, 0, 0, VIA / 2), (0, 0, 0, -VIA / 2)):
            p = [occ.addPoint(x0, y0, zb), occ.addPoint(x1, y1, zb), occ.addPoint(x1, y1, zt), occ.addPoint(x0, y0, zt)]
            via_s.append(occ.addPlaneSurface([occ.addCurveLoop([occ.addLine(p[i], p[(i + 1) % 4]) for i in range(4)])]))
    tools = [(2, t) for t in pec_s] + [(2, t) for t in port_s] + [(2, t) for t in via_s]
    _, outmap = occ.fragment([(3, v) for v in vols], tools); occ.synchronize()
    pieces = outmap[len(vols):]
    pec_tags = sorted({t for m in pieces[:len(pec_s)] for d, t in m if d == 2})
    port_tags = [sorted({t for d, t in m if d == 2}) for m in pieces[len(pec_s):len(pec_s) + len(port_s)]]
    via_tags = [sorted({t for d, t in m if d == 2}) for m in pieces[len(pec_s) + len(port_s):]]
    pec_tags = [t for t in pec_tags if all(t not in pt for pt in port_tags + via_tags)]
    top = zt + hup; EPS = 1e-6
    groups = {k: [] for k in range(1, 16)}; xs, Xs, ys, Ys = [], [], [], []
    for dim, tag in gmsh.model.getEntities(2):
        x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(2, tag); cx, cy, cz = occ.getCenterOfMass(2, tag)
        if x1 < -A / 2 + EPS: groups[1].append(tag); xs.append((tag, cy, cz))
        elif x0 > A / 2 - EPS: groups[2].append(tag); Xs.append((tag, cy, cz))
        elif y1 < -A / 2 + EPS: groups[3].append(tag); ys.append((tag, cx, cz))
        elif y0 > A / 2 - EPS: groups[4].append(tag); Ys.append((tag, cx, cz))
        elif z1 < EPS: groups[5].append(tag)
        elif abs(z0 - top) < EPS and abs(z1 - top) < EPS: groups[6].append(tag)
    groups[7] = pec_tags
    if BOSS[0] > 0:
        for dim, tag in gmsh.model.getEntities(2):
            x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(2, tag)
            if z1 > EPS and z1 < zb + EPS and max(abs(x0), abs(x1), abs(y0), abs(y1)) < BOSS[0] / 2 + EPS and tag not in sum(via_tags, []):
                groups[7].append(tag)               # boss top + side walls
    for i, pt in enumerate(port_tags):
        groups[8 + i] = pt
    for i, vt in enumerate(via_tags):
        groups[12 + i] = vt
    for k, v in groups.items():
        if v: gmsh.model.addPhysicalGroup(2, v, k)
    air, sub = [], []
    for dim, tag in gmsh.model.getEntities(3):
        cx, cy, cz = occ.getCenterOfMass(3, tag); (sub if (G < cz < zt) else air).append(tag)
    gmsh.model.addPhysicalGroup(3, air, 1); gmsh.model.addPhysicalGroup(3, sub, 2)
    key = lambda s: (round(s[1], 4), round(s[2], 4))
    gmsh.model.mesh.setPeriodic(2, [s[0] for s in sorted(Xs, key=key)], [s[0] for s in sorted(xs, key=key)], [1, 0, 0, A, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1])
    gmsh.model.mesh.setPeriodic(2, [s[0] for s in sorted(Ys, key=key)], [s[0] for s in sorted(ys, key=key)], [1, 0, 0, 0, 0, 1, 0, A, 0, 0, 1, 0, 0, 0, 0, 1])
    f = gmsh.model.mesh.field; f.add("MathEval", 1)
    f.setString(1, "F", f"Min({hmax}, {h0} + {gr}*(((z-{zt})+Abs(z-{zt}))/2 + (({G}-z)+Abs({G}-z))/2))"); f.setAsBackgroundMesh(1)
    for o in ("Mesh.MeshSizeFromPoints", "Mesh.MeshSizeExtendFromBoundary", "Mesh.MeshSizeFromCurvature"):
        gmsh.option.setNumber(o, 0)
    gmsh.model.mesh.generate(3); gmsh.option.setNumber("Mesh.MshFileVersion", 2.2); gmsh.option.setNumber("Mesh.Binary", 1)
    n = len(gmsh.model.mesh.getElementsByType(4)[0]); gmsh.write(fn); gmsh.finalize()
    assert all(len(p) > 0 for p in port_tags), "a diode port surface is missing"
    return n


def config(mesh, freqs, loads=None, theta=0.0, phi=0.0, pol="TE"):
    """loads None -> multiport (each lumped port its own excitation); else dict {port_idx: R} passive, Floquet-only."""
    cfg = make_config(mesh, "x", freqs, theta, phi, pol, {}, excite=None if loads else {1: 1})
    lp = []
    for i, k in enumerate(ORDER):
        sx, sy = DIRS[k]; d = {"Index": 2 + i, "Attributes": [8 + i], "Direction": [sx / math.sqrt(2), sy / math.sqrt(2), 0.0]}
        if loads:
            d["R"] = float(loads[2 + i])
        else:
            d["R"] = R0; d["Excitation"] = 2 + i
        lp.append(d)
    if VIA > 0:                                          # port 6 = centre via (4 parallel elements, +Z)
        d = {"Index": 6, "Elements": [{"Attributes": [12 + j], "Direction": "+Z"} for j in range(4)]}
        if loads:
            d["R"] = float(loads[6])
        else:
            d["R"] = R0; d["Excitation"] = 6
        lp.append(d)
    cfg["Boundaries"]["LumpedPort"] = lp
    if os.environ.get("PALACE_SWEEP"):                        # "fmin,fmax,step,tol": adaptive (reduced-order) dense sweep
        a0, a1, st, tol = [float(v) for v in os.environ["PALACE_SWEEP"].split(",")]
        cfg["Solver"]["Driven"] = {"Samples": [{"Type": "Linear", "MinFreq": a0, "MaxFreq": a1, "FreqStep": st}], "AdaptiveTol": tol}
    if os.environ.get("PALACE_LIN"):                          # JSON dict merged into Solver.Linear (oblique incidence needs help)
        cfg["Solver"]["Linear"].update(json.loads(os.environ["PALACE_LIN"]))
    return cfg


def cplx(h, d, name):
    norm = lambda s: s.replace(" ", "")
    mi = [i for i, x in enumerate(h) if norm(x).startswith("|" + name + "|")]
    if mi:
        ai = [i for i, x in enumerate(h) if norm(x).startswith("arg(" + name + ")")][0]
        return 10 ** (d[:, mi[0]] / 20) * np.exp(1j * np.radians(d[:, ai]))
    re = [i for i, x in enumerate(h) if norm(x).startswith("Re{" + name + "}")][0]
    im = [i for i, x in enumerate(h) if norm(x).startswith("Im{" + name + "}")][0]
    return d[:, re] + 1j * d[:, im]


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("tag"); ap.add_argument("freqs")
    for k, v in (("cap", 5.9), ("w", 0.4), ("cpad", 1.27), ("gap", 0.64), ("st", 0.3), ("h0", 0.25), ("ring", 0.0), ("taper", 0.0)):
        ap.add_argument("--" + k, type=float, default=v)
    ap.add_argument("--direct", default=""); ap.add_argument("--theta", type=float, default=0.0)
    ap.add_argument("--phi", type=float, default=0.0); ap.add_argument("--pol", default="TE")
    a = ap.parse_args(); fq = [float(x) for x in a.freqs.split(",")]
    os.makedirs("mp", exist_ok=True)
    pec, ports = geometry(a.cap, a.w, a.cpad, a.gap, a.st, a.ring, a.taper)
    os.makedirs("mp", exist_ok=True); mesh = f"mp/{a.tag}.msh"; nt = make_mesh(mesh, pec, ports, h0=a.h0)
    loads = None
    if a.direct:
        loads = {2 + i: float(r) for i, r in enumerate(a.direct.split(","))}
    cfg = config(mesh, fq, loads, a.theta, a.phi, a.pol)
    t0 = time.time(); dt, od, log = run_palace(cfg, "mp/" + a.tag, a.tag)
    fh, fd = read_csv(od + "/port-floquet-S.csv"); f = fd[:, 0]
    res = dict(tag=a.tag, f=list(map(float, f)), geom=dict(A=A, boss=BOSS, via=VIA, s=a.cap, w=a.w, cp=a.cpad, g=a.gap, st=a.st, ring=a.ring, taper=a.taper, G=G), tets=nt, wall=dt, direct=loads, theta=a.theta, phi=a.phi, pol=a.pol)
    if loads:
        res["TE"] = [[float(z.real), float(z.imag)] for z in cplx(fh, fd, "S[P1(0;0)TE][1]")]
        res["TM"] = [[float(z.real), float(z.imag)] for z in cplx(fh, fd, "S[P1(0;0)TM][1]")]
        vh, vd = read_csv(od + "/port-V.csv")                # raw port voltages (for the TM-in column at oblique incidence)
        res["V"] = [[[float(z.real), float(z.imag)] for z in cplx(vh, vd, f"V[{k}]")] for k in range(2, 6 + (1 if VIA > 0 else 0))]
    else:
        vh, vd = read_csv(od + "/port-V.csv")
        NP = 4 + (1 if VIA > 0 else 0)
        S = np.zeros((len(f), 2 + NP, 1 + NP), complex)     # rows TE,TM,b2..b(1+NP) ; cols TE_in,a2..a(1+NP)
        for e in range(1, 2 + NP):
            S[:, 0, e - 1] = cplx(fh, fd, f"S[P1(0;0)TE][{e}]"); S[:, 1, e - 1] = cplx(fh, fd, f"S[P1(0;0)TM][{e}]")
            vinc = None if e == 1 else vd[:, [i for i, x in enumerate(vh) if x.replace(" ", "").startswith(f"V_inc[{e}][{e}]")][0]]
            for k in range(2, 2 + NP):
                V = cplx(vh, vd, f"V[{k}][{e}]")
                if e == 1:
                    S[:, k, 0] = V / math.sqrt(2 * R0)          # Floquet-driven port waves: 1/sqrt(2R) (validated: kappa = 0.7071 at all f)
                else:
                    S[:, k, e - 1] = (V - (vinc if k == e else 0)) / vinc
        res["S"] = [[[[float(S[i, r, c].real), float(S[i, r, c].imag)] for c in range(1 + NP)] for r in range(2 + NP)] for i in range(len(f))]
    json.dump(res, open(f"mp/{a.tag}.json", "w"))
    print(f"{a.tag}: tets {nt}, {dt:.0f}s, {'direct' if loads else 'multiport'}", flush=True)
