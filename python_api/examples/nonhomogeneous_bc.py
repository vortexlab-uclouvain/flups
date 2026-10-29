#!/usr/bin/env python3
"""
Non-homogeneous Dirichlet / Neumann boundary conditions with FLUPS: RHS correction

FLUPS only knows homogeneous symmetric conditions: ODD (u = 0 on the wall) and EVEN (du/dn = 0).
The lattice kernels (LGF_p, MEHR_*L) are the exact inverse of a stencil (Solver.get_stencil), so any
ghost closure of the form "reflection + known data" turns into an exact correction of the rhs:

    L_h R(u) = M_h f - L_h D    on the unknown nodes,   u = R(u) + D

where R(u) is the symmetric (ODD/EVEN) extension of u, the one FLUPS solves for, and D = u - R(u)
is zero on the unknown nodes and known data on the wall/ghost nodes (pyflups.bc, and
doc/symmetric_conditions.md for the derivation).

The solve only uses what a user has: f sampled on the nodes of the domain, and g (Dirichlet, u = g)
or q (Neumann, du/dn = q with n the inward normal) sampled on the nodes of each wall. A manufactured
solution with values and slopes that vary along every wall provides these samples and the error
measurement, nothing else. Two modes:
    - "data" (default): the ghosts are PDE Taylor expansions about the wall,
          u(n) = sum_k a_k n^k / k!,  a_0 = g (Dirichlet) or a_1 = q (Neumann),  a_{k+2} = d_n^k f - Lap_t a_k,
      keeping the terms k < q (--ghost-order, default q = p), with every a_k computed by finite
      differences of the samples: one-sided in the normal direction, along the wall for Lap_t.
      Order p needs a Dirichlet ghost O(h^p) and a Neumann ghost O(h^{p+1}), so a_k only needs
      O(h^(P-k)) with P = p (Dirichlet) or p + 1 (Neumann): each finite difference is taken with the
      accuracy the next coefficient requires. The rhs is M_h f. MEHR_*L also needs f beyond the walls
      (in M_h f) and the wall data beyond the face edges (the edge/corner ghosts read by the diagonal
      terms of the stencil): both are extrapolated from the samples. The cross-shaped LGF_p stencils
      need neither. This checks the convergence order;
    - "exact": an algebra check, not a solve from data. The ghosts and the rhs = L_h u come from the
      manufactured solution itself, so the error must be at round-off level for every kernel and BC
      set (stencil returned by FLUPS == kernel symbol, node layout, reflections).

The data are node-centred: N + 1 nodes per direction, both wall nodes included.

The BC string lists the 6 faces x-, x+, y-, y+, z-, z+ with D (Dirichlet), N (Neumann),
P (periodic, both sides of a direction). Without any Dirichlet face the solution is known up to a
constant only: FLUPS drops the k = 0 mode, so the weighted mean is matched before measuring the
error and the compatibility defect of the corrected rhs is reported.

Examples (the FLUPS shared library is found through FLUPS_LIB_PATH):
    export FLUPS_LIB_PATH=/path/to/flups/lib
    mpirun -n 1 python nonhomogeneous_bc.py --kernel=LGF_2 --bc=DDNNPP --res=16,32,64
    mpirun -n 2 python nonhomogeneous_bc.py --kernel=MEHR_4L --bc=DNDDND --mode=exact --nproc=1,1,2
    mpirun -n 1 python nonhomogeneous_bc.py --suite
"""

import argparse
import math
import os
import sys

import numpy as np
from mpi4py import MPI

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import pyflups as pf  # noqa: E402
from pyflups import bc as pbc  # noqa: E402

KERNELS = {
    "LGF_2": (pf.GreenType.LGF_2, 2),
    "LGF_4": (pf.GreenType.LGF_4, 4),
    "LGF_6": (pf.GreenType.LGF_6, 6),
    "LGF_8": (pf.GreenType.LGF_8, 8),
    "MEHR_4L": (pf.GreenType.MEHR_4L, 4),
    "MEHR_6L": (pf.GreenType.MEHR_6L, 6),
}
BCTYPE = {"D": pf.BoundaryType.ODD, "N": pf.BoundaryType.EVEN, "P": pf.BoundaryType.PER}


def unit(d, m=1):
    e = [0, 0, 0]
    e[d] = m
    return tuple(e)


# ---------------------------------------------------------------------------------------------
# manufactured solution u = X(x) Y(y) Z(z), each factor Im(exp(lam s + i beta)): only used to
# sample the data (f, g, q) and to measure the error, plus the "exact" algebra check
# ---------------------------------------------------------------------------------------------
class Manufactured:
    # non-symmetric wall factors: neither the value nor the slope vanishes on any wall
    WALL = [(0.8 + 2.6j, 0.3), (-0.6 + 3.1j, 1.1), (0.5 + 2.2j, -0.7)]

    def __init__(self, bc, L):
        self.lam, self.beta = [], []
        for d in range(3):
            if bc[d][0] == "P":
                self.lam.append(2j * math.pi * (d + 1) / L)
                self.beta.append(0.4 * (d + 1))
            else:
                self.lam.append(self.WALL[d][0])
                self.beta.append(self.WALL[d][1])

    def phi(self, d, k, s):
        return np.imag(self.lam[d] ** k * np.exp(self.lam[d] * s + 1j * self.beta[d]))

    def du(self, alpha, xs):
        """d^alpha u on the tensor grid xs = [x, y, z] (1D arrays)"""
        return np.multiply.outer(np.multiply.outer(self.phi(0, alpha[0], xs[0]), self.phi(1, alpha[1], xs[1])),
                                 self.phi(2, alpha[2], xs[2]))

    def f(self, xs):
        """f = lap u"""
        return sum(self.du(unit(d, 2), xs) for d in range(3))


def sample_data(sol, bc, x, L):
    """what the user knows: f on the nodes, g or q = du/dn (n inward) on the nodes of each wall"""
    f = sol.f(x)
    walls = {}
    for d in range(3):
        if bc[d][0] == "P":
            continue
        for side in range(2):
            xs = list(x)
            xs[d] = np.array([0.0 if side == 0 else L])
            if bc[d][side] == "D":
                walls[(d, side)] = sol.du((0, 0, 0), xs).squeeze(axis=d)
            else:
                walls[(d, side)] = (1.0 if side == 0 else -1.0) * sol.du(unit(d), xs).squeeze(axis=d)
    return f, walls


# ---------------------------------------------------------------------------------------------
# finite differences on the samples
# ---------------------------------------------------------------------------------------------
def fd_weights(offsets, k):
    """weights of the k-th derivative at 0 from samples at the offsets (in units of h), exact for the
    polynomials of degree len(offsets) - 1"""
    m = len(offsets)
    b = np.zeros(m)
    b[k] = math.factorial(k)
    return np.linalg.solve(np.vander(np.asarray(offsets, dtype=float), m, increasing=True).T, b)


def fd_at(a, axis, x, t, k, npts, h):
    """k-th derivative along axis at the position t (k = 0: interpolation/extrapolation), from the npts
    samples a(x) closest to t; the axis is removed"""
    j0 = int(np.clip(np.searchsorted(x, t) - npts // 2, 0, len(x) - npts))
    idx = np.arange(j0, j0 + npts)
    w = fd_weights((x[idx] - t) / h, k) / h**k
    return np.tensordot(w, np.take(a, idx, axis=axis), axes=(0, axis))


def d2(a, axis, x, h, acc, period):
    """second derivative along axis on every sample, O(h^acc): centred, one-sided near the ends of a
    bounded direction, wrapped with a period of `period` points in a periodic one"""
    if period:
        m = (acc + 2) // 2
        offs = np.arange(-m, m + 1)
        w = fd_weights(offs, 2) / h**2
        idx = np.arange(a.shape[axis])
        return sum(wj * np.take(a, (idx + o) % period, axis=axis) for wj, o in zip(w, offs))
    return np.stack([fd_at(a, axis, x, t, 2, acc + 2, h) for t in x], axis=axis)


def extend(a, axis, x, r, h, npts, period):
    """a extended by r points on each side of axis: wrapped if periodic, else extrapolated through the
    npts closest samples"""
    if period:
        return np.take(a, np.arange(-r, a.shape[axis] + r) % period, axis=axis)
    lo = [fd_at(a, axis, x, x[0] - m * h, 0, npts, h) for m in range(r, 0, -1)]
    hi = [fd_at(a, axis, x, x[-1] + m * h, 0, npts, h) for m in range(1, r + 1)]
    return np.concatenate([np.stack(lo, axis=axis), a, np.stack(hi, axis=axis)], axis=axis)


# ---------------------------------------------------------------------------------------------
# wall/ghost data for pbc.ghost_remainder: ghost(d, side, m) on the extended face
# ---------------------------------------------------------------------------------------------
class DataGhosts:
    """Taylor ghosts from the samples only: a_0 = g or a_1 = q, a_{k+2} = d_n^k f - Lap_t a_k, each
    finite difference taken with the accuracy O(h^(P-k-2)) that a_{k+2} needs"""

    def __init__(self, f, walls, bc, x, L, h, r, order, ghost_order, period):
        self.f, self.walls, self.bc, self.x, self.L, self.h, self.r = f, walls, bc, x, L, h, r
        self.order, self.K, self.period = order, ghost_order or order, period
        self.coefs = {}

    def _coefs(self, d, side):
        neumann = self.bc[d][side] == "N"
        P = self.order + neumann
        sg = 1.0 if side == 0 else -1.0
        xw = 0.0 if side == 0 else self.L
        tdirs = [t for t in range(3) if t != d]
        k = 1 if neumann else 0
        a = {k: self.walls[(d, side)]}
        while k + 2 < self.K:
            acc = max(P - k - 2, 1)
            dnf = sg**k * fd_at(self.f, d, self.x[d], xw, k, k + acc, self.h)
            lapt = sum(d2(a[k], i, self.x[t], self.h, acc, self.period[t]) for i, t in enumerate(tdirs))
            a[k + 2] = dnf - lapt
            k += 2
        # beyond the face edges, only read by the diagonal terms of MEHR_*
        for j in a:
            for i, t in enumerate(tdirs):
                a[j] = extend(a[j], i, self.x[t], self.r, self.h, max(P - j, 2), self.period[t])
        return a

    def __call__(self, d, side, m):
        if (d, side) not in self.coefs:
            self.coefs[(d, side)] = self._coefs(d, side)
        return pbc.taylor_ghost(self.coefs[(d, side)], m, self.h, BCTYPE[self.bc[d][side]])


class ExactData:
    """algebra check only: the ghost data of the manufactured solution itself"""

    def __init__(self, sol, bc, xs, L, h):
        self.sol, self.bc, self.xs, self.L, self.h = sol, bc, xs, L, h

    def _u(self, d, s):
        xs = list(self.xs)
        xs[d] = np.array([s])
        return self.sol.du((0, 0, 0), xs).squeeze(axis=d)

    def __call__(self, d, side, m):
        sw, sg = (0.0, 1.0) if side == 0 else (self.L, -1.0)
        if m == 0:
            return self._u(d, sw)
        uin, uout = self._u(d, sw + sg * m * self.h), self._u(d, sw - sg * m * self.h)
        return uout + uin if self.bc[d][side] == "D" else uout - uin


# ---------------------------------------------------------------------------------------------
# one solve
# ---------------------------------------------------------------------------------------------
def bcs(bc):
    return "".join(a + b for a, b in bc)


def weights_1d(bc_d, n):
    """quadrature of the k = 0 mode: trapezoidal for EVEN-EVEN (DCT-I), duplicate excluded if PER"""
    w = np.ones(n)
    if bc_d[0] == "P":
        w[-1] = 0.0
    else:
        w[0] = 0.5 if bc_d[0] == "N" else 0.0
        w[-1] = 0.5 if bc_d[1] == "N" else 0.0
    return w


def solve_case(kernel, bc, N, mode, ghost_order, nproc, comm, L=1.0):
    green, order = KERNELS[kernel]
    h = L / N
    n = N + 1
    fbc = [[BCTYPE[bc[d][s]] for s in range(2)] for d in range(3)]

    topo = pf.Topology(axis=0, lda=1, nglob=[n] * 3, nproc=nproc, is_complex=False, comm=comm)
    solver = pf.Solver(topo, fbc, [h] * 3, [L] * 3, center_type=[pf.CenterType.NODE_CENTER] * 3)
    solver.set_green_type(green)
    solver.setup()
    lhs = solver.get_stencil(pf.StencilType.LHS)
    rhs_op = solver.get_stencil(pf.StencilType.RHS)
    r = max(pbc.stencil_reach(lhs[0]), pbc.stencil_reach(rhs_op[0]))

    # corrected rhs on the global grid
    x = [np.arange(n) * h for _ in range(3)]
    sol = Manufactured(bc, L)
    if mode == "exact":
        xext = [np.arange(-r, n + r) * h for _ in range(3)]
        D = pbc.ghost_remainder(fbc, n, r, ExactData(sol, bc, xext, L, h))
        rhs = pbc.apply_stencil(*lhs, sol.du((0, 0, 0), xext), r) / h**2
    else:
        f, walls = sample_data(sol, bc, x, L)
        period = [None if bc[d][0] != "P" else N for d in range(3)]
        ghosts = DataGhosts(f, walls, bc, x, L, h, r, order, ghost_order, period)
        D = pbc.ghost_remainder(fbc, n, r, ghosts)
        fext = f
        for d in range(3):
            fext = extend(fext, d, x[d], r, h, order + 1, period[d])
        rhs = pbc.apply_stencil(*rhs_op, fext, r)
    rhs = pbc.corrected_rhs(*lhs, rhs, D, h, fbc)

    # FLUPS solve on the local pencil, gathered back on every rank
    nloc = [topo.get_nloc(d) for d in range(3)]
    nmem = [topo.get_nmem(d) for d in range(3)]
    i0 = topo.get_istart()
    loc = tuple(slice(i0[d], i0[d] + nloc[d]) for d in range(3))
    frhs = np.zeros(topo.memsize)
    fsol = np.zeros(topo.memsize)
    frhs.reshape(nmem[2], nmem[1], nmem[0])[:nloc[2], :nloc[1], :nloc[0]] = rhs[loc].T
    solver.solve(fsol, frhs, pf.SolverType.STD)
    v = np.zeros((n, n, n))
    v[loc] = fsol.reshape(nmem[2], nmem[1], nmem[0])[:nloc[2], :nloc[1], :nloc[0]].T
    comm.Allreduce(MPI.IN_PLACE, v, op=MPI.SUM)
    del solver, topo

    # error against the manufactured solution
    uh = pbc.recover(v, D, fbc)
    uex = sol.du((0, 0, 0), x)
    w = np.multiply.outer(np.multiply.outer(weights_1d(bc[0], n), weights_1d(bc[1], n)), weights_1d(bc[2], n))
    has_dirichlet = "D" in bcs(bc)
    shift = 0.0 if has_dirichlet else np.sum(w * (uh - uex)) / np.sum(w)
    compat = 0.0 if has_dirichlet else np.sum(w * rhs) / np.sum(w)
    err = uh - shift - uex
    return np.max(np.abs(err)) / np.max(np.abs(uex)), math.sqrt(np.mean(err**2)) / np.max(np.abs(uex)), compat


def run(kernel, bc, res, mode, ghost_order, nproc, comm, out=print):
    order = KERNELS[kernel][1]
    out(f"# {kernel:8s} bc={bcs(bc)} mode={mode}"
        + (f" ghost_order={ghost_order or order}" if mode == "data" else "")
        + f"  (expected order {order})")
    out(f"#   {'N':>4s} {'err_inf':>11s} {'err_2':>11s} {'ord_inf':>7s} {'ord_2':>7s} {'compat':>10s}")
    prev, rows = None, []
    for N in res:
        einf, e2, compat = solve_case(kernel, bc, N, mode, ghost_order, nproc, comm)
        oi = o2 = float("nan")
        if prev is not None:
            q = math.log(N / prev[0])
            oi, o2 = math.log(prev[1] / einf) / q, math.log(prev[2] / e2) / q
        out(f"    {N:4d} {einf:11.3e} {e2:11.3e} {oi:7.2f} {o2:7.2f} {compat:10.2e}")
        prev = (N, einf, e2)
        rows.append((N, einf, e2, oi, o2))
    return rows


def parse_bc(s):
    s = s.upper()
    if len(s) != 6 or any(c not in "DNP" for c in s):
        raise argparse.ArgumentTypeError("--bc needs 6 letters among D, N, P (x-, x+, y-, y+, z-, z+)")
    bc = [(s[2 * d], s[2 * d + 1]) for d in range(3)]
    if any((a == "P") != (b == "P") for a, b in bc):
        raise argparse.ArgumentTypeError("P must be set on both sides of a direction")
    return bc


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kernel", default="LGF_2", choices=list(KERNELS))
    ap.add_argument("--bc", default="DNDDNN", type=parse_bc, help="6 letters D/N/P for x-, x+, y-, y+, z-, z+")
    ap.add_argument("--res", default="16,32,64", help="comma-separated N (cells per direction)")
    ap.add_argument("--mode", default="data", choices=["data", "exact"],
                    help="data: solve from the samples of f, g, q; exact: algebra check with the exact ghosts")
    ap.add_argument("--ghost-order", type=int, default=0, help="keep the Taylor ghost terms a_k with k < q, default q = order")
    ap.add_argument("--nproc", default="", help="ranks per direction, e.g. 1,1,2 (default 1,1,size)")
    ap.add_argument("--suite", action="store_true", help="every kernel on 3 BC sets, exact + data")
    args = ap.parse_args()

    comm = MPI.COMM_WORLD
    rank, size = comm.Get_rank(), comm.Get_size()
    nproc = [int(v) for v in args.nproc.split(",")] if args.nproc else [1, 1, size]
    if int(np.prod(nproc)) != size:
        raise SystemExit(f"--nproc={nproc} does not match {size} ranks")
    out = print if rank == 0 else (lambda *a, **k: None)
    res = [int(v) for v in args.res.split(",")]

    if not args.suite:
        run(args.kernel, args.bc, res, args.mode, args.ghost_order, nproc, comm, out)
        return

    bcsets = ["DDDDDD", "NNNNNN", "DNPPND"]
    worst_exact, fails = 0.0, []
    for kernel in KERNELS:
        for b in bcsets:
            bc = parse_bc(b)
            (_, einf, _, _, _), = run(kernel, bc, [32], "exact", 0, nproc, comm, out)
            worst_exact = max(worst_exact, einf)
            rows = run(kernel, bc, res, "data", 0, nproc, comm, out)
            ordr = rows[-1][3]
            if not ordr > KERNELS[kernel][1] - 0.5:
                fails.append(f"{kernel} {b}: order {ordr:.2f}")
    out(f"# suite: worst exact-mode error {worst_exact:.2e}; data-mode order failures: {fails or 'none'}")


if __name__ == "__main__":
    main()
