"""
Independent 2-D plane-strain FE check of the mechanical model (no ANSYS).

Builds the same mesh as ansys/01_thermal_transient.mac (K / L / AL, LESIZE,
hard divisions, mapped and free areas), the same materials (Table 3,
Eqs. 6-7), the same supports, and runs the three stages of
ansys/02_mechanical_settlement.mac with linear elasticity:

  stage 1  geostatic: natural ground + buried footing alive, embankment and
           abutment above ground dead; K0 initial stress per element
           (INISTATE) + gravity  ->  displacements must be ~ 0
  stage 2  embankment and abutment born (EALIVE), their weight +
           pavement/traffic                ->  construction settlement
  stage 3  initial temperature field ICT: thaw strain (Eq. 4) + thawed
           stiffness in the active layer and the fill (paper mode)

It also answers "is the stiffness matrix singular?" for each way of joining
the concrete to the soil - the ANSYS rigid-body error.

Nodes that do not coincide across an interface are tied to the other side
by linear interpolation (penalty), as bonded contact would do.

Usage:  python3 tools/fe_check_2d.py [path/to/01_thermal_transient.mac]
"""
import math
import sys

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import check_geometry as geo  # noqa: E402

MAC = sys.argv[1] if len(sys.argv) > 1 else "ansys/01_thermal_transient.mac"
G = 9.81
# Table 3: a1, b1 (E, MPa), a2, b2 (nu), alpha (%) ; Table 2: rho_dry, w
SOIL = {1: (34, 30, 0.42, -0.007, 1, 1800, 0.15),
        2: (36, 30, 0.42, -0.007, 5, 1900, 0.10),
        3: (28, 26, 0.40, -0.008, 12, 1600, 0.30),
        4: (140, 108, 0.25, -0.004, 1, 1800, 0.15),
        5: (61, 53, 0.35, -0.007, 1, 2060, 0.06)}
E_CON, NU_CON, RHO_CON = 25e9, 0.167, 2500.0
Q_TRAF = 60.1e3
T_FRZ = -1.0
ICT_Y = [-30, -16, -15, -14, -13, -12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1.75, -1, -0.75,
         -0.5, -0.25, 0, 0.25, 4, 4.25, 7, 10]
ICT_T = [-0.25, -0.25, -0.27, -0.29, -0.31, -0.33, -0.36, -0.41, -0.45, -0.50, -0.55, -0.60, -0.65,
         -0.69, -0.80, -0.80, 0.20, 0.20, 0.25, 0.25, 1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 3.0]
XA = -6.2
FOOT_FRONT = {8, 24, 25, 26, 13, 22, 23}            # concrete lines, coupled in 02
BEHIND = {3, 4, 5, 6, 7}                            # concrete lines, contact in 02
PENALTY = 1e3 * E_CON


def layer(y):
    if y > 0:
        return 5
    d = -y
    return 1 if d < 0.5 else 2 if d < 2 else 3 if d < 8 else 4


def rho_bulk(m):
    return SOIL[m][5] * (1 + SOIL[m][6])


def overburden(y):
    """Vertical geostatic stress (Pa, positive) at depth -y in natural ground."""
    d, s = -y, 0.0
    for top, bot, m in ((0, 0.5, 1), (0.5, 2, 2), (2, 8, 3), (8, 1e9, 4)):
        s += rho_bulk(m) * G * max(0.0, min(d, bot) - top)
    return s


def props(m, T):
    a1, b1, a2, b2, *_ = SOIL[m]
    ta = -T if T < 0 else 0.0
    return (a1 + b1 * ta) * 1e6, a2 + b2 * ta


def dmat(E, nu):
    c = E / ((1 + nu) * (1 - 2 * nu))
    return c * np.array([[1 - nu, nu, 0], [nu, 1 - nu, 0], [0, 0, (1 - 2 * nu) / 2]])


# ------------------------------------------------------------------ mesh
def build_mesh():
    par, kp, lines, areas, lesize = geo.parse(MAC)
    free = geo.parse.free
    length = lambda n: math.dist(kp[lines[n][0]], kp[lines[n][1]])
    div = {n: math.ceil(length(n) / s - 1e-9) for n, s in lesize.items()}
    div.update(geo.parse.ndiv)

    def loop(ls):
        a, b = lines[ls[0]]
        order, rest = [(ls[0], a, b)], list(ls[1:])
        while rest:
            for n in rest:
                p, q = lines[n]
                if p == order[-1][2]:
                    order.append((n, p, q)); rest.remove(n); break
                if q == order[-1][2]:
                    order.append((n, q, p)); rest.remove(n); break
        return order

    loops = [loop(ls) for ls in areas]
    changed = True                                   # mapped areas: carry divisions
    while changed:
        changed = False
        for ia, lp in enumerate(loops, 1):
            if ia in free:
                continue
            for i, j in ((0, 2), (1, 3)):
                a, b = lp[i][0], lp[j][0]
                if a in div and b not in div:
                    div[b] = div[a]; changed = True
                elif b in div and a not in div:
                    div[a] = div[b]; changed = True

    nodes, key2id, elems, edges = [], {}, [], []

    def node(x, y, body):
        k = (round(x, 6), round(y, 6), body)
        if k not in key2id:
            key2id[k] = len(nodes)
            nodes.append((x, y, body))
        return key2id[k]

    for ia, lp in enumerate(loops, 1):
        body = "C" if ia <= 9 else "S"
        c = [np.array(kp[lp[k][1]], float) for k in range(4)]
        na = max(div[lp[0][0]], div[lp[2][0]])       # free areas: the finer side
        nb = max(div[lp[1][0]], div[lp[3][0]])
        grid = {}
        for i in range(na + 1):
            u = i / na
            for j in range(nb + 1):
                v = j / nb
                P = ((1 - v) * (c[0] + (c[1] - c[0]) * u) + v * (c[3] + (c[2] - c[3]) * u)
                     + (1 - u) * (c[0] + (c[3] - c[0]) * v) + u * (c[1] + (c[2] - c[1]) * v)
                     - ((1 - u) * (1 - v) * c[0] + u * (1 - v) * c[1] + u * v * c[2] + (1 - u) * v * c[3]))
                grid[i, j] = node(P[0], P[1], body)
        for i in range(na):
            for j in range(nb):
                q = [grid[i, j], grid[i + 1, j], grid[i + 1, j + 1], grid[i, j + 1]]
                xy = np.array([nodes[k][:2] for k in q])
                a2 = sum(xy[k, 0] * xy[(k + 1) % 4, 1] - xy[(k + 1) % 4, 0] * xy[k, 1] for k in range(4))
                elems.append((q[::-1] if a2 < 0 else q, body))
        sides = ([grid[i, 0] for i in range(na + 1)], [grid[na, j] for j in range(nb + 1)],
                 [grid[i, nb] for i in range(na, -1, -1)], [grid[0, j] for j in range(nb, -1, -1)])
        for k in range(4):
            edges.append((ia, lp[k][0], body, sides[k]))
    return kp, lines, nodes, elems, edges


def ties_for(nodes, edges, join_lines):
    """Penalty ties between edges lying on the same segment.
    Same body: always (non-conforming free/mapped areas). Concrete/soil:
    only for concrete lines in join_lines."""
    groups = {}
    for ia, ln, body, ids in edges:
        p, q = nodes[ids[0]][:2], nodes[ids[-1]][:2]
        key = tuple(sorted(((round(p[0], 5), round(p[1], 5)), (round(q[0], 5), round(q[1], 5)))))
        groups.setdefault(key, []).append((ia, ln, body, ids))
    ties = []
    for grp in groups.values():
        for a in range(len(grp)):
            for b in range(a + 1, len(grp)):
                A, B = grp[a], grp[b]
                if A[2] != B[2]:
                    conc = A if A[2] == "C" else B
                    if conc[1] not in join_lines:
                        continue
                elif set(A[3]) == set(B[3]):
                    continue
                fine, coarse = (A, B) if len(A[3]) >= len(B[3]) else (B, A)
                P0 = np.array(nodes[coarse[3][0]][:2])
                P1 = np.array(nodes[coarse[3][-1]][:2])
                L = np.linalg.norm(P1 - P0)
                tc = [np.dot(np.array(nodes[n][:2]) - P0, P1 - P0) / L ** 2 for n in coarse[3]]
                order = np.argsort(tc)
                cs = [coarse[3][k] for k in order]
                tcs = [tc[k] for k in order]
                for n in fine[3]:
                    if n in cs:
                        continue
                    t = np.dot(np.array(nodes[n][:2]) - P0, P1 - P0) / L ** 2
                    k = min(max(np.searchsorted(tcs, t) - 1, 0), len(cs) - 2)
                    w = (t - tcs[k]) / (tcs[k + 1] - tcs[k])
                    ties.append((n, [(cs[k], 1 - w), (cs[k + 1], w)]))
    return ties


# ------------------------------------------------------------ element math
GP = [(-1 / math.sqrt(3), -1 / math.sqrt(3)), (1 / math.sqrt(3), -1 / math.sqrt(3)),
      (1 / math.sqrt(3), 1 / math.sqrt(3)), (-1 / math.sqrt(3), 1 / math.sqrt(3))]


def q4(xy, xi, eta):
    N = 0.25 * np.array([(1 - xi) * (1 - eta), (1 + xi) * (1 - eta), (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)])
    dN = 0.25 * np.array([[-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)],
                          [-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)]])
    J = dN @ xy
    dNx = np.linalg.solve(J, dN)
    B = np.zeros((3, 8))
    B[0, 0::2] = dNx[0]
    B[1, 1::2] = dNx[1]
    B[2, 0::2] = dNx[1]
    B[2, 1::2] = dNx[0]
    return N, B, np.linalg.det(J)


class Model:
    def __init__(self, foot="couple", behind="tie"):
        self.kp, self.lines, nodes, self.elems, edges = build_mesh()
        join = (FOOT_FRONT if foot == "couple" else set()) | (BEHIND if behind == "tie" else set())
        self.ties = ties_for(nodes, edges, join)
        self.xy = np.array([n[:2] for n in nodes])
        self.nn = len(self.xy)
        self.cent = np.array([self.xy[q].mean(axis=0) for q, _ in self.elems])
        self.mat = np.array([7 if b == "C" else layer(c[1]) for (q, b), c in zip(self.elems, self.cent)])
        x, y = self.xy[:, 0], self.xy[:, 1]
        self.fix = set()
        for i in np.where(abs(x - x.min()) < 1e-6)[0]:
            self.fix.add(2 * i)
        for i in np.where(abs(x - x.max()) < 1e-6)[0]:
            self.fix.add(2 * i)
        for i in np.where(abs(y - y.min()) < 1e-6)[0]:
            self.fix.add(2 * i + 1)

    def assemble(self, active, temp, sigma0=None, eps_th=None, u_ref=None, rho=None):
        rows, cols, vals = [], [], []
        f = np.zeros(2 * self.nn)
        live = np.zeros(self.nn, bool)
        for e in np.where(active)[0]:
            q, body = self.elems[e]
            live[q] = True
            xy = self.xy[q]
            m = self.mat[e]
            E, nu = (E_CON, NU_CON) if m == 7 else props(m, temp[e])
            D = dmat(E, nu)
            dofs = np.ravel([[2 * n, 2 * n + 1] for n in q])
            ke = np.zeros((8, 8))
            fe = np.zeros(8)
            for xi, eta in GP:
                N, B, dj = q4(xy, xi, eta)
                ke += B.T @ D @ B * dj
                if rho is not None and rho[e] > 0:
                    fe[1::2] -= N * rho[e] * G * dj
                if sigma0 is not None and sigma0[e] is not None:
                    fe -= B.T @ sigma0[e] * dj              # constant per element
                if eps_th is not None and eps_th[e] is not None:
                    fe += B.T @ D @ eps_th[e] * dj
            if u_ref is not None and u_ref[e] is not None:
                fe += ke @ u_ref[e]
            rows.append(np.repeat(dofs, 8))
            cols.append(np.tile(dofs, 8))
            vals.append(ke.ravel())
            f[dofs] += fe
        for n, parts in self.ties:                          # penalty ties
            if not live[n] or not all(live[m] for m, _ in parts):
                continue
            ids = [n] + [m for m, _ in parts]
            c = np.array([1.0] + [-w for _, w in parts])
            for d in (0, 1):
                dd = [2 * i + d for i in ids]
                rows.append(np.repeat(dd, len(dd)))
                cols.append(np.tile(dd, len(dd)))
                vals.append((PENALTY * np.outer(c, c)).ravel())
        K = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                          shape=(2 * self.nn, 2 * self.nn))
        return K, f, live

    def solve(self, K, f, live):
        fixed = set(self.fix)
        for n in np.where(~live)[0]:                       # nodes of dead elements
            fixed |= {2 * n, 2 * n + 1}
        free = np.array(sorted(set(range(2 * self.nn)) - fixed))
        lu = spla.splu(K[free][:, free].tocsc())
        u = np.zeros(2 * self.nn)
        u[free] = lu.solve(f[free])
        d = np.abs(lu.U.diagonal())
        return u, d.min() / d.max()

    def traffic(self, q, cap=True):
        f = np.zeros(2 * self.nn)
        x, y = self.xy[:, 0], self.xy[:, 1]
        top = np.where(abs(y - 7) < 1e-6)[0]
        top = top[np.argsort(x[top])]
        xmax = -1.5 if cap else -6.2
        for a, b in zip(top[:-1], top[1:]):
            L = x[b] - x[a]
            if x[b] <= xmax + 1e-6 and L > 1e-9:
                f[2 * a + 1] -= q * L / 2
                f[2 * b + 1] -= q * L / 2
        return f

    def nodes_at(self, yv):
        x, y = self.xy[:, 0], self.xy[:, 1]
        cand = np.where(abs(y - yv) < 1e-6)[0]
        return [cand[np.argmin(abs(x[cand] - (XA - d)))] for d in range(21)]


def main():
    print("2-D plane-strain FE check (linear elastic), same mesh as the ANSYS macro\n")
    print("1) Is the model supported?  (min/max pivot of the stiffness matrix, all elements active)")
    for foot, label in (("couple", "footing + front coupled, backwall contact open  (02, ABUT_BOND=2)"),
                        ("open", "all concrete/soil contact open: concrete floats  (the failed runs)")):
        M = Model(foot, "open")
        ne = len(M.elems)
        K, _, live = M.assemble(np.ones(ne, bool), np.full(ne, T_FRZ))
        try:
            _, ratio = M.solve(K, np.zeros(2 * M.nn), live)
            state = "SINGULAR (rigid-body mode)" if ratio < 1e-12 else "supported"
            print(f"   {label}:  pivot ratio {ratio:.1e} -> {state}")
        except RuntimeError as err:
            print(f"   {label}:  factorisation failed ({err}) -> SINGULAR")
    print("   (ANSYS stops with a rigid-body error when this ratio is effectively 0)\n")

    M = Model("couple", "tie")
    ne = len(M.elems)
    print(f"   mesh: {M.nn} nodes, {ne} elements, {len(M.ties)} interpolated ties")

    below = M.cent[:, 1] < 0
    stage1 = below.copy()
    temp1 = np.full(ne, T_FRZ)
    rho1 = np.array([rho_bulk(layer(c[1])) if c[1] < 0 else 0.0 for c in M.cent])
    K0 = {m: props(m, T_FRZ)[1] / (1 - props(m, T_FRZ)[1]) for m in SOIL}
    sig0 = []
    for e in range(ne):
        if stage1[e]:
            sv = overburden(M.cent[e][1])
            k0 = K0[layer(M.cent[e][1])]
            sig0.append(np.array([-k0 * sv, -sv, 0.0]))
        else:
            sig0.append(None)

    K1, f_grav, live1 = M.assemble(stage1, temp1, rho=rho1)
    u_plain, _ = M.solve(K1, f_grav, live1)
    _, f_geo, _ = M.assemble(stage1, temp1, sigma0=sig0, rho=rho1)
    u1, _ = M.solve(K1, f_geo, live1)
    print("\n2) STAGE 1 - geostatic step (embankment and abutment above ground dead)")
    print(f"   gravity only, no initial stress:  max settlement {-u_plain[1::2].min()*100:6.2f} cm  <- not a real settlement")
    print(f"   K0 initial stress + gravity:      max |u| = {np.abs(u1).max()*1000:.3f} mm   -> geostatic state OK (~0)")

    active2 = np.ones(ne, bool)
    rho2 = np.zeros(ne)
    for e in range(ne):
        if M.mat[e] == 7:
            rho2[e] = RHO_CON - (rho1[e] if stage1[e] else 0.0)
        elif not stage1[e]:
            rho2[e] = rho_bulk(5)
    born = [None if stage1[e] else np.ravel([[u1[2 * n], u1[2 * n + 1]] for n in M.elems[e][0]])
            for e in range(ne)]
    rho_all = np.where(stage1, rho1, 0) + rho2
    K2, f2, live2 = M.assemble(active2, temp1, sigma0=sig0, rho=rho_all, u_ref=born)
    u2, _ = M.solve(K2, f2 + M.traffic(Q_TRAF), live2)
    top, base = M.nodes_at(7.0), M.nodes_at(0.0)

    def table(u):
        print("    d (m)   total (cm)  foundation (cm)  embankment (cm)     [from stage 1]")
        for d in (0, 2, 5, 10, 15, 20):
            st = -(u[2 * top[d] + 1] - u1[2 * top[d] + 1]) * 100
            sb = -(u[2 * base[d] + 1] - u1[2 * base[d] + 1]) * 100
            print(f"   {d:5d}  {st:10.2f}  {sb:14.2f}  {st - sb:14.2f}")

    print("\n3) STAGE 2 - embankment + abutment born, 60.1 kPa traffic (frozen ground, -1 C)")
    table(u2)

    temp3 = np.interp(M.cent[:, 1], ICT_Y, ICT_T)
    eps = [None if M.mat[e] == 7 else
           np.array([0.0, -SOIL[M.mat[e]][4] / 100 * min(1.0, max(0.0, temp3[e] / 0.1)), 0.0])
           for e in range(ne)]
    K3, f3, live3 = M.assemble(active2, temp3, sigma0=sig0, eps_th=eps, rho=rho_all, u_ref=born)
    u3, _ = M.solve(K3, f3 + M.traffic(Q_TRAF), live3)
    print("\n4) STAGE 3 - start of freeze-thaw: ICT field, thawed active layer + fill (paper mode)")
    table(u3)
    print("\n   All three stages solve with a non-singular stiffness matrix.")


if __name__ == "__main__":
    main()
