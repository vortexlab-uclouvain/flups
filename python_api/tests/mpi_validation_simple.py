import os
import sys
import numpy as np
from mpi4py import MPI

# Ensure Python can find pyflups
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
os.environ.setdefault('PYTHONPATH', ROOT)
import pyflups as pf
from pyflups.enums import BoundaryType, GreenType, SolverType, DiffType, CenterType


def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size_comm = comm.Get_size()

    # Domain setup similar to samples/validation
    nglob = [64, 64, 64]
    L = [2.0 * np.pi, 2.0 * np.pi, 2.0 * np.pi]
    h = [L[0] / nglob[0], L[1] / nglob[1], L[2] / nglob[2]]

    # Split processes evenly (user can adapt); default 1,1,1 if single rank
    def factorize(n):
        # naive factorization into three dims as close as possible
        a = int(round(n ** (1/3)))
        npx = max(1, a)
        npy = max(1, a)
        npz = max(1, n // (npx * npy))
        while npx * npy * npz != n:
            if npx * npy * npz < n:
                npx += 1
            else:
                npz = max(1, npz-1)
        return [npx, npy, npz]

    nproc = factorize(size_comm)

    topo = pf.Topology(axis=0, lda=1, nglob=nglob, nproc=nproc, is_complex=False, comm=comm)

    bc = [
        [BoundaryType.PER, BoundaryType.PER],
        [BoundaryType.PER, BoundaryType.PER],
        [BoundaryType.PER, BoundaryType.PER],
    ]

    solver = pf.Solver(topo, bc, h, L, order_diff=DiffType.NOD, center_type=[CenterType.CELL_CENTER]*3)
    solver.set_green_type(GreenType.LGF_2)
    solver.setup()

    # Build RHS: f(x)=sin(x)+sin(y)+sin(z)
    nloc = [topo.get_nloc(0), topo.get_nloc(1), topo.get_nloc(2)]
    nmem = [topo.get_nmem(0), topo.get_nmem(1), topo.get_nmem(2)]
    # Allocate with padded memory sizes (nmem) to respect internal padding
    rhs = np.zeros((nmem[0]*nmem[1]*nmem[2],), dtype=np.float64)
    field = np.zeros_like(rhs)

    istart = topo.get_istart()
    ax = topo.axis
    size = [topo.get_nmem(0), topo.get_nmem(1), topo.get_nmem(2)]
    nf = 1

    for i0 in range(nloc[0]):
        x = (istart[0]+i0) * h[0]
        for i1 in range(nloc[1]):
            y = (istart[1]+i1) * h[1]
            for i2 in range(nloc[2]):
                z = (istart[2]+i2) * h[2]
                val = np.sin(x) + np.sin(y) + np.sin(z)
                lid = pf.locID(ax, i0, i1, i2, 0, ax, size, nf)
                rhs[lid] = val

    solver.solve(field, rhs, solver_type=SolverType.STD)

    # Analytical solution phi = -rhs
    phi_analytical = -rhs
    # Compute local error
    num = np.linalg.norm(field - phi_analytical)
    den = max(1e-12, np.linalg.norm(phi_analytical))
    rel_err_local = num / den

    # Reduce to global max error
    rel_err_global = comm.allreduce(rel_err_local, op=MPI.MAX)

    if rank == 0:
        print(f"MPI ranks: {size_comm}, nproc grid: {nproc}, relative error: {rel_err_global:.3e}")
        if rel_err_global < 1e-2:
            print("OK: validation passed")
            sys.exit(0)
        else:
            print("FAIL: validation error too high")
            sys.exit(1)


if __name__ == "__main__":
    main()
