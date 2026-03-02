import os
import numpy as np
from mpi4py import MPI

# Ensure Python can find pyflups
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
os.environ.setdefault('PYTHONPATH', ROOT)
import pyflups as pf
from pyflups.enums import BoundaryType, GreenType, SolverType, DiffType, CenterType


def test_simple_poisson_periodic():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nproc = [1, 1, 1]

    # Domain setup similar to samples/validation
    nglob = [64, 64, 64]
    L = [2.0 * np.pi, 2.0 * np.pi, 2.0 * np.pi]
    h = [L[0] / nglob[0], L[1] / nglob[1], L[2] / nglob[2]]

    topo = pf.Topology(axis=0, lda=1, nglob=nglob, nproc=nproc, is_complex=False, comm=comm)

    bc = [
        [BoundaryType.PER, BoundaryType.PER],
        [BoundaryType.PER, BoundaryType.PER],
        [BoundaryType.PER, BoundaryType.PER],
    ]

    solver = pf.Solver(topo, bc, h, L, order_diff=DiffType.NOD, center_type=[CenterType.CELL_CENTER]*3)
    solver.set_green_type(GreenType.LGF_2)
    solver.setup()

    # Build RHS: f(x)=sin(x)+sin(y)+sin(z) (scaled)
    nloc = [topo.get_nloc(0), topo.get_nloc(1), topo.get_nloc(2)]
    rhs = np.zeros((nloc[0]*nloc[1]*nloc[2],), dtype=np.float64)
    field = np.zeros_like(rhs)

    # Simple initialization: use indices mapping like in validation
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

    # Solve
    solver.solve(field, rhs, solver_type=SolverType.STD)

    # Basic sanity: solution finite and non-trivial
    assert np.all(np.isfinite(field)), "Non-finite values in solution"
    assert np.max(np.abs(field)) > 1e-6, "Solution seems zero; check solver path"

    # Analytical solution for periodic Poisson with rhs = sin(x)+sin(y)+sin(z)
    # Since ∇^2(sin) = -sin, the solution φ satisfies φ = -rhs (up to constant)
    phi_analytical = -rhs

    err = np.linalg.norm(field - phi_analytical) / max(1e-12, np.linalg.norm(phi_analytical))
    # Allow modest tolerance due to discretization/Green kernel choice
    assert err < 1e-2, f"Relative error too high: {err}"
