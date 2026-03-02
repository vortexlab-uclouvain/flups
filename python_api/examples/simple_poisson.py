#!/usr/bin/env python3
"""
Simple Poisson Solver Example using PyFLUPS

This example demonstrates basic usage of PyFLUPS to solve
a simple 3D Poisson equation with periodic boundary conditions.

Equation: ∇²φ = f
where f is a point source at the center of the domain.

Usage:
    mpirun -np 4 python simple_poisson.py
"""

import numpy as np
from mpi4py import MPI
import pyflups as pf


def main():
    # MPI setup
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    if rank == 0:
        print("="*60)
        print("PyFLUPS - Simple Poisson Solver Example")
        print("="*60)
        print(f"Running on {size} MPI processes")
        print()

    # =========================================================================
    # DOMAIN PARAMETERS
    # =========================================================================
    nglob = [64, 64, 64]  # Global grid size
    L = [2.0 * np.pi, 2.0 * np.pi, 2.0 * np.pi]  # Domain size
    h = [L[i] / nglob[i] for i in range(3)]  # Grid spacing

    if rank == 0:
        print(f"Domain size: {L}")
        print(f"Grid points: {nglob}")
        print(f"Grid spacing: {h}")
        print()

    # =========================================================================
    # MPI DECOMPOSITION
    # =========================================================================
    # Determine process decomposition
    if size == 1:
        nproc = [1, 1, 1]
    elif size == 2:
        nproc = [2, 1, 1]
    elif size == 4:
        nproc = [2, 2, 1]
    elif size == 8:
        nproc = [2, 2, 2]
    else:
        # Simple decomposition for other sizes
        nproc = [size, 1, 1]

    if rank == 0:
        print(f"Process decomposition: {nproc}")
        print()

    # =========================================================================
    # BOUNDARY CONDITIONS
    # =========================================================================
    # Periodic in all directions
    bc = [
        [pf.BoundaryType.PER, pf.BoundaryType.PER],  # x
        [pf.BoundaryType.PER, pf.BoundaryType.PER],  # y
        [pf.BoundaryType.PER, pf.BoundaryType.PER],  # z
    ]

    if rank == 0:
        print("Boundary conditions:")
        print("  x: PERIODIC")
        print("  y: PERIODIC")
        print("  z: PERIODIC")
        print()

    # =========================================================================
    # CREATE TOPOLOGY
    # =========================================================================
    if rank == 0:
        print("Creating topology...")

    topo = pf.Topology(
        axis=0,           # Fastest rotating index
        lda=1,            # Scalar field (1 component)
        nglob=nglob,
        nproc=nproc,
        is_complex=False,
        comm=comm
    )

    # Get local sizes
    nloc = [topo.get_nloc(i) for i in range(3)]
    nmem = [topo.get_nmem(i) for i in range(3)]
    istart = topo.get_istart()

    if rank == 0:
        print("  Topology created successfully")
        print(f"  Local size on rank 0: {nloc}")
        print(f"  Memory size on rank 0: {nmem}")
        print(f"  Start indices on rank 0: {istart}")
        print()

    # =========================================================================
    # CREATE SOLVER
    # =========================================================================
    if rank == 0:
        print("Creating solver...")

    solver = pf.Solver(
        topology=topo,
        boundary_conditions=bc,
        h=h,
        L=L,
        order_diff=pf.DiffType.NOD,
        center_type=[pf.CenterType.CELL_CENTER] * 3
    )

    # Set Green's function type
    solver.set_green_type(pf.GreenType.LGF_2)

    # Setup solver
    if rank == 0:
        print("  Setting up solver (allocating memory)...")

    solver.setup(change_comm=False)

    alloc_size = solver.get_alloc_size()
    if rank == 0:
        print("  Solver setup complete")
        print(f"  Allocated memory: {alloc_size * 8 / 1024**2:.2f} MB")
        print()

    # =========================================================================
    # PREPARE DATA
    # =========================================================================
    if rank == 0:
        print("Preparing data...")

    # Calculate local padded memory size (respect internal padding)
    local_size_mem = nmem[0] * nmem[1] * nmem[2]

    # Allocate arrays with padded sizes
    rhs = np.zeros(local_size_mem, dtype=np.float64)
    field = np.zeros(local_size_mem, dtype=np.float64)

    # Create a point source at the center of the domain
    # Find which rank contains the center
    center_global = [n // 2 for n in nglob]

    # Check if this rank contains the center
    contains_center = True
    for i in range(3):
        if not (istart[i] <= center_global[i] < istart[i] + nloc[i]):
            contains_center = False
            break

    if contains_center:
        # Calculate local indices of center
        center_local = [center_global[i] - istart[i] for i in range(3)]
        center_idx = pf.locID(
            axsrc=0,
            i0=center_local[0],
            i1=center_local[1],
            i2=center_local[2],
            lia=0,
            axtrg=topo.axis,
            size=nmem,
            nf=1
        )
        rhs[center_idx] = 1.0
        print(f"Rank {rank}: Point source placed at local index {center_idx}")

    if rank == 0:
        print()

    # =========================================================================
    # SOLVE
    # =========================================================================
    if rank == 0:
        print("Solving Poisson equation...")

    comm.Barrier()
    t_start = MPI.Wtime()

    solver.solve(field, rhs, pf.SolverType.STD)

    comm.Barrier()
    t_end = MPI.Wtime()

    if rank == 0:
        print(f"  Solution computed in {t_end - t_start:.4f} seconds")
        print()

    # =========================================================================
    # ANALYZE RESULTS
    # =========================================================================
    if rank == 0:
        print("Analyzing results...")

    # Compute statistics
    field_min = np.min(field)
    field_max = np.max(field)
    field_mean = np.mean(field)

    # Gather global statistics
    global_min = comm.allreduce(field_min, op=MPI.MIN)
    global_max = comm.allreduce(field_max, op=MPI.MAX)
    # Use padded local size for reduction; divide by physical points for readability
    global_mean = comm.allreduce(field_mean * local_size_mem, op=MPI.SUM) / np.prod(nglob)

    if rank == 0:
        print(f"  Field min:  {global_min:.6e}")
        print(f"  Field max:  {global_max:.6e}")
        print(f"  Field mean: {global_mean:.6e}")
        print()

    # =========================================================================
    # EXPORT (Optional - requires HDF5 support in FLUPS)
    # =========================================================================
    # Uncomment to export results
    # if rank == 0:
    #     print("Exporting results...")
    # try:
    #     pf.hdf5_dump(topo, "poisson_solution", field)
    #     if rank == 0:
    #         print("  Results exported to data/poisson_solution.h5")
    # except Exception as e:
    #     if rank == 0:
    #         print(f"  Export failed: {e}")

    # =========================================================================
    # CLEANUP
    # =========================================================================
    if rank == 0:
        print("="*60)
        print("Simulation completed successfully!")
        print("="*60)


if __name__ == "__main__":
    main()
