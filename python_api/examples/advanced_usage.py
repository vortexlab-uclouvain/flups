#!/usr/bin/env python3
"""
Advanced Usage Example - PyFLUPS

Demonstrates advanced features:
- Using internal FLUPS buffers
- Manual FFT operations
- Spectral space manipulation
- Profiling

Usage:
    mpirun -np 4 python advanced_usage.py
"""

import numpy as np
from mpi4py import MPI
import pyflups as pf


def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    if rank == 0:
        print("="*60)
        print("PyFLUPS - Advanced Usage Example")
        print("="*60)
        print()

    # Setup
    nglob = [32, 32, 32]
    L = [1.0, 1.0, 1.0]
    h = [L[i] / nglob[i] for i in range(3)]
    nproc = [2, 2, 1] if comm.Get_size() == 4 else [comm.Get_size(), 1, 1]

    bc = [
        [pf.BoundaryType.PER, pf.BoundaryType.PER],
        [pf.BoundaryType.PER, pf.BoundaryType.PER],
        [pf.BoundaryType.PER, pf.BoundaryType.PER],
    ]

    # =========================================================================
    # EXAMPLE 1: Using a Profiler
    # =========================================================================
    if rank == 0:
        print("EXAMPLE 1: Using Profiler")
        print("-" * 40)

    # Create profiler
    prof = pf.Profiler("advanced_example")

    # Create topology
    topo = pf.Topology(
        axis=0, lda=1, nglob=nglob, nproc=nproc,
        is_complex=False, comm=comm
    )

    # Create solver with profiler
    solver = pf.Solver(
        topology=topo,
        boundary_conditions=bc,
        h=h, L=L,
        profiler=prof
    )

    solver.set_green_type(pf.GreenType.LGF_4)
    solver.setup()

    # Prepare data
    nloc = [topo.get_nloc(i) for i in range(3)]
    size = nloc[0] * nloc[1] * nloc[2]

    rhs = np.random.random(size)
    field = np.zeros(size, dtype=np.float64)

    # Solve
    solver.solve(field, rhs)

    # Display profiler - all ranks must call this (it's a collective MPI operation)
    if rank == 0:
        print("\nProfiling results:")
    prof.display()  # ALL ranks must call this

    # =========================================================================
    # EXAMPLE 2: Getting Spectral Information
    # =========================================================================
    if rank == 0:
        print("EXAMPLE 2: Spectral Information")
        print("-" * 40)

    kfact, koffset, symstart = solver.get_spectral_info()

    if rank == 0:
        print(f"k-factor:  {kfact}")
        print(f"k-offset:  {koffset}")
        print(f"symstart:  {symstart}")
        print()

    # =========================================================================
    # EXAMPLE 3: Using Internal Buffer
    # =========================================================================
    if rank == 0:
        print("EXAMPLE 3: Internal Buffer Access")
        print("-" * 40)

    # Get internal buffer size
    alloc_size = solver.get_alloc_size()

    if rank == 0:
        print(f"Internal buffer size: {alloc_size} doubles")
        print(f"Memory: {alloc_size * 8 / 1024**2:.2f} MB")

    # solver.get_inner_buffer() wraps this buffer in a numpy array without copy, on FFTW builds only:
    # on a GPU build it lives on the device. solver.do_copy reads/writes it on every build, see
    # manual_solve.py.
    if rank == 0:
        print()

    # =========================================================================
    # EXAMPLE 4: Manual FFT Operations
    # =========================================================================
    if rank == 0:
        print("EXAMPLE 4: Manual FFT Operations")
        print("-" * 40)

    # solve() is exactly these five steps, on the inner buffer of FLUPS (see manual_solve.py for
    # spectral-space operations between them)
    rhs = np.random.random(size)
    result = np.zeros(size, dtype=np.float64)
    reference = np.zeros(size, dtype=np.float64)
    solver.solve(reference, rhs)

    solver.do_copy(topo, rhs, pf.FLUPS_FORWARD)          # 1. rhs -> inner buffer
    solver.do_fft(sign=pf.FLUPS_FORWARD)                 # 2. forward FFT, in place
    solver.do_mult(solver_type=pf.SolverType.STD)        # 3. Green's function and normalisation
    solver.do_fft(sign=pf.FLUPS_BACKWARD)                # 4. backward FFT, in place
    solver.do_copy(topo, result, pf.FLUPS_BACKWARD)      # 5. inner buffer -> result

    diff = comm.allreduce(np.max(np.abs(result - reference)), op=MPI.MAX)
    if rank == 0:
        print(f"max |steps - solve()| = {diff:.1e}")
        print()

    # =========================================================================
    # EXAMPLE 5: Different Green Functions Comparison
    # =========================================================================
    if rank == 0:
        print("EXAMPLE 5: Comparing Green Functions")
        print("-" * 40)

    green_types = [
        pf.GreenType.LGF_2,
        pf.GreenType.LGF_4,
        pf.GreenType.HEJ_2,
        pf.GreenType.CHAT_2,
    ]

    test_rhs = np.random.random(size)

    results = {}

    for green_type in green_types:
        solver_test = pf.Solver(
            topology=topo,
            boundary_conditions=bc,
            h=h, L=L
        )
        solver_test.set_green_type(green_type)
        solver_test.setup()

        field_test = np.zeros(size, dtype=np.float64)

        comm.Barrier()
        t_start = MPI.Wtime()
        solver_test.solve(field_test, test_rhs)
        comm.Barrier()
        t_end = MPI.Wtime()

        results[green_type] = {
            'time': t_end - t_start,
            'max': np.max(np.abs(field_test)),
            'mean': np.mean(field_test)
        }

        # Cleanup to free memory
        del solver_test

    if rank == 0:
        print("Green function comparison:")
        print(f"{'Type':<12} {'Time (s)':<12} {'Max':<12} {'Mean':<12}")
        print("-" * 50)
        for gt in green_types:
            name = gt.name
            r = results[gt]
            print(f"{name:<12} {r['time']:<12.4f} {r['max']:<12.4e} {r['mean']:<12.4e}")

    print()

    # =========================================================================
    # EXAMPLE 6: Memory Layout and Index Computation
    # =========================================================================
    if rank == 0:
        print("EXAMPLE 6: Memory Layout")
        print("-" * 40)

    nmem = [topo.get_nmem(i) for i in range(3)]
    axis = topo.axis
    lda = topo.lda

    if rank == 0:
        print(f"Axis (FRI): {axis}")
        print(f"Local size: {nloc}")
        print(f"Memory size: {nmem}")
        print()
        print("Computing linear index for point (1,2,3) with component 0:")

    # Compute linear index
    linear_idx = pf.locID(
        axsrc=0,
        i0=1, i1=2, i2=3,
        lia=0,
        axtrg=axis,
        size=nmem,
        nf=lda
    )

    if rank == 0 and linear_idx < size:
        print(f"  Linear index: {linear_idx}")

    print()

    # =========================================================================
    # Cleanup
    # =========================================================================
    if rank == 0:
        print("="*60)
        print("Advanced examples completed!")
        print("="*60)


if __name__ == "__main__":
    main()
