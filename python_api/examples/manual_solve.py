#!/usr/bin/env python3
"""
Advanced mode: the steps of Solver.solve done one by one, and what can be done in between

Solver.solve(field, rhs) is the sequence

    do_copy(topo, rhs, FLUPS_FORWARD)      rhs -> inner buffer of FLUPS (host -> device on a GPU build)
    do_fft(sign=FLUPS_FORWARD)             physical -> spectral, in place in the inner buffer
    do_mult(solver_type=STD)               multiplication by the Green's function and the FFT normalisation
    do_fft(sign=FLUPS_BACKWARD)            spectral -> physical, in place
    do_copy(topo, field, FLUPS_BACKWARD)   inner buffer -> field

This example runs it in three ways, on a right-hand side made of a few Laplacian eigenmodes, so that
every result is known exactly:

    1. the five steps above: the same result as solve(), bit for bit;
    2. a Gaussian filter exp(-sigma^2 |k|^2 / 2) applied to the spectral solution, between do_mult and
       the backward FFT: a solve and a filter for the price of a solve;
    3. a user Green's function instead of do_mult: the screened Poisson equation
       lap(u) - alpha^2 u = f, i.e. u_hat = -f_hat / (|k|^2 + alpha^2).

Spectral data: between do_fft(FLUPS_FORWARD) and do_fft(FLUPS_BACKWARD) the inner buffer holds the
spectrum in the layout of solver.get_inner_topo_spectral(). It is real when every direction is
symmetric (DST/DCT) and complex as soon as one is periodic. It is read and written back with do_copy and
that topology (works on both FFTW and GPU builds), viewed with pf.field_view(topo_spec, array) as
[component, i0, i1, i2], and solver.get_wavenumbers() gives the wave numbers of its entries. On an FFTW
build, pf.field_view(topo_spec, solver.get_inner_buffer()) avoids the two copies.

FFT normalisation: the transforms of FLUPS are unnormalised, do_mult applies the factor. A user
multiplier must apply it too: part 3 measures it once with a forward + backward FFT of a constant field.

The checks use the spectral kernel CHAT_2, for which an eigenmode is solved exactly (the lattice
kernels LGF_* / MEHR_* invert a finite-difference Laplacian instead: --kernel shows that part 1 still
reproduces solve() exactly, the analytic errors then measure the discretisation error).

The BC string lists the 6 faces x-, x+, y-, y+, z-, z+ with D (homogeneous Dirichlet, ODD),
N (homogeneous Neumann, EVEN), P (periodic, both sides of a direction). Non-homogeneous walls: see
nonhomogeneous_bc.py.

Examples (the FLUPS shared library is found through FLUPS_LIB_PATH):
    export FLUPS_LIB_PATH=/path/to/flups/lib
    mpirun -n 1 python manual_solve.py                      # Dirichlet in x and y, Neumann in z
    mpirun -n 2 python manual_solve.py --bc=PPDNNN --center=cell --nproc=1,1,2
"""

import argparse
import os
import sys

import numpy as np
from mpi4py import MPI

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import pyflups as pf  # noqa: E402

BCTYPE = {"D": pf.BoundaryType.ODD, "N": pf.BoundaryType.EVEN, "P": pf.BoundaryType.PER}

# right-hand side: sum of amplitude * phi_m0(x) phi_m1(y) phi_m2(z)
MODES = [(1.0, (1, 1, 1)), (0.5, (2, 1, 3)), (0.25, (3, 2, 1))]


def eigenmode(bc_d, m, L):
    """1D Laplacian eigenmode number m >= 1 for the walls bc_d of a direction: (phi, k)"""
    if bc_d == ("P", "P"):
        k = 2.0 * np.pi * m / L
        return (lambda x: np.sin(k * x + 0.3)), k
    k = (m - 0.5 * (bc_d[0] != bc_d[1])) * np.pi / L
    if bc_d[0] == "D":
        return (lambda x: np.sin(k * x)), k
    return (lambda x: np.cos(k * x)), k


def analytic(bc, L, x, multiplier):
    """sum over the modes of amplitude * multiplier(|k|^2) * phi(x, y, z) on the local grid x"""
    u = 0.0
    for a, m in MODES:
        phi, k = zip(*[eigenmode(bc[d], m[d], L) for d in range(3)])
        u = u + a * multiplier(sum(kd**2 for kd in k)) * np.multiply.outer(np.multiply.outer(phi[0](x[0]), phi[1](x[1])), phi[2](x[2]))
    return u


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
    ap.add_argument("--bc", default="DDDDNN", type=parse_bc, help="6 letters D/N/P for x-, x+, y-, y+, z-, z+")
    ap.add_argument("--N", type=int, default=32, help="cells per direction")
    ap.add_argument("--center", default="node", choices=["node", "cell"], help="node- or cell-centred data")
    ap.add_argument("--kernel", default="CHAT_2", choices=[g.name for g in pf.GreenType])
    ap.add_argument("--nproc", default="", help="ranks per direction, e.g. 1,1,2 (default 1,1,size)")
    args = ap.parse_args()

    comm = MPI.COMM_WORLD
    rank, size = comm.Get_rank(), comm.Get_size()
    nproc = [int(v) for v in args.nproc.split(",")] if args.nproc else [1, 1, size]
    if int(np.prod(nproc)) != size:
        raise SystemExit(f"--nproc={nproc} does not match {size} ranks")
    out = print if rank == 0 else (lambda *a, **k: None)

    bc, N, L = args.bc, args.N, 1.0
    cell = args.center == "cell"
    center = pf.CenterType.CELL_CENTER if cell else pf.CenterType.NODE_CENTER
    h = L / N
    n = N if cell else N + 1  # node-centred data includes both wall nodes

    topo = pf.Topology(axis=0, lda=1, nglob=[n] * 3, nproc=nproc, is_complex=False, comm=comm)
    solver = pf.Solver(topo, [[BCTYPE[bc[d][s]] for s in range(2)] for d in range(3)], [h] * 3, [L] * 3,
                       center_type=[center] * 3)
    solver.set_green_type(pf.GreenType[args.kernel])
    solver.setup()

    # local grid and rhs, written through a [component, i0, i1, i2] view of the FLUPS layout
    istart = topo.get_istart()
    x = [(istart[d] + np.arange(topo.get_nloc(d)) + 0.5 * cell) * h for d in range(3)]
    rhs = np.zeros(topo.memsize)
    field = np.zeros(topo.memsize)
    pf.field_view(topo, rhs)[0] = analytic(bc, L, x, lambda k2: 1.0)

    def error(uex):
        """max error of the solution in field, relative to max |uex|"""
        e = comm.allreduce(np.max(np.abs(pf.field_view(topo, field)[0] - uex), initial=0.0), op=MPI.MAX)
        return e / comm.allreduce(np.max(np.abs(uex), initial=0.0), op=MPI.MAX)

    def spectrum_apply(op):
        """apply op(spectrum view, kx, ky, kz) in place to the spectral data held by the inner buffer"""
        topo_spec = solver.get_inner_topo_spectral()
        spec = np.zeros(topo_spec.memsize)
        solver.do_copy(topo_spec, spec, pf.FLUPS_BACKWARD)  # inner buffer -> spec
        op(pf.field_view(topo_spec, spec), *solver.get_wavenumbers())
        solver.do_copy(topo_spec, spec, pf.FLUPS_FORWARD)  # spec -> inner buffer
        return topo_spec

    out(f"# {args.kernel}, BC {''.join(a + b for a, b in bc)}, {args.center}-centred, {n}^3 points, {size} rank(s)")

    # reference: solve
    solver.solve(field, rhs)
    ref = field.copy()
    poisson = analytic(bc, L, x, lambda k2: -1.0 / k2)
    out(f"solve():                     error {error(poisson):.2e}")

    # 1. the steps of solve
    field[:] = 0.0
    solver.do_copy(topo, rhs, pf.FLUPS_FORWARD)
    solver.do_fft(sign=pf.FLUPS_FORWARD)
    solver.do_mult(solver_type=pf.SolverType.STD)
    solver.do_fft(sign=pf.FLUPS_BACKWARD)
    solver.do_copy(topo, field, pf.FLUPS_BACKWARD)
    diff = comm.allreduce(np.max(np.abs(field - ref)), op=MPI.MAX)
    out(f"1. step by step:             error {error(poisson):.2e}, max |field - solve()| = {diff:.1e}")

    # 2. Gaussian filter of the solution, in spectral space
    sigma = 2.0 * h

    def gaussian(s, kx, ky, kz):
        s *= np.exp(-0.5 * sigma**2 * (kx**2 + ky**2 + kz**2))

    solver.do_copy(topo, rhs, pf.FLUPS_FORWARD)
    solver.do_fft(sign=pf.FLUPS_FORWARD)
    solver.do_mult(solver_type=pf.SolverType.STD)
    topo_spec = spectrum_apply(gaussian)
    solver.do_fft(sign=pf.FLUPS_BACKWARD)
    solver.do_copy(topo, field, pf.FLUPS_BACKWARD)
    kind = "complex" if topo_spec.is_complex else "real"
    out(f"2. filtered solution:        error {error(analytic(bc, L, x, lambda k2: -np.exp(-0.5 * sigma**2 * k2) / k2)):.2e}"
        f" ({kind} spectrum, axis {topo_spec.axis})")

    # 3. user Green's function for lap(u) - alpha^2 u = f, the FFT normalisation measured once
    pf.field_view(topo, field)[0] = 1.0
    solver.do_copy(topo, field, pf.FLUPS_FORWARD)
    solver.do_fft(sign=pf.FLUPS_FORWARD)
    solver.do_fft(sign=pf.FLUPS_BACKWARD)
    solver.do_copy(topo, field, pf.FLUPS_BACKWARD)  # 1 / normfact on every point, 0 on the ODD wall nodes
    normfact = 1.0 / comm.allreduce(np.max(field), op=MPI.MAX)
    alpha = 10.0

    def screened(s, kx, ky, kz):
        s *= -normfact / (kx**2 + ky**2 + kz**2 + alpha**2)

    solver.do_copy(topo, rhs, pf.FLUPS_FORWARD)
    solver.do_fft(sign=pf.FLUPS_FORWARD)
    spectrum_apply(screened)
    solver.do_fft(sign=pf.FLUPS_BACKWARD)
    solver.do_copy(topo, field, pf.FLUPS_BACKWARD)
    out(f"3. screened Poisson:         error {error(analytic(bc, L, x, lambda k2: -1.0 / (k2 + alpha**2))):.2e}"
        f" (alpha = {alpha:g}, normfact = {normfact:.6e})")


if __name__ == "__main__":
    main()
