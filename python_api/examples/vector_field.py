#!/usr/bin/env python3
"""
Vector Field Example using PyFLUPS

This example demonstrates solving for a vector field using the
rotational (Biot-Savart) formulation.

Equation: ∇²φ = ∇ × ω
where ω is a vorticity field.

Usage:
    mpirun -np 4 python vector_field.py
"""

import numpy as np
from mpi4py import MPI
import pyflups as pf


def create_vortex_ring(x, y, z, R=0.5, r0=0.1):
    """
    Create a vortex ring vorticity distribution
    
    Parameters
    ----------
    x, y, z : arrays
        Grid coordinates
    R : float
        Ring radius
    r0 : float
        Core radius
    """
    # Distance from z-axis
    rho = np.sqrt(x**2 + y**2)
    
    # Distance from ring center
    dist = np.sqrt((rho - R)**2 + z**2)
    
    # Gaussian vorticity distribution
    omega = np.exp(-(dist / r0)**2)
    
    # Convert to toroidal vorticity (azimuthal direction)
    omega_x = -y / (rho + 1e-10) * omega
    omega_y = x / (rho + 1e-10) * omega
    omega_z = np.zeros_like(omega)
    
    return omega_x, omega_y, omega_z


def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    if rank == 0:
        print("="*60)
        print("PyFLUPS - Vector Field Example")
        print("="*60)
        print(f"Running on {size} MPI processes")
        print()
    
    # Domain parameters
    nglob = [64, 64, 64]
    L = [2.0, 2.0, 2.0]
    h = [L[i] / nglob[i] for i in range(3)]
    
    # Process decomposition
    if size == 4:
        nproc = [2, 2, 1]
    elif size == 8:
        nproc = [2, 2, 2]
    else:
        nproc = [size, 1, 1]
    
    # Boundary conditions - periodic
    bc = [
        [pf.BoundaryType.PER, pf.BoundaryType.PER],
        [pf.BoundaryType.PER, pf.BoundaryType.PER],
        [pf.BoundaryType.PER, pf.BoundaryType.PER],
    ]
    
    if rank == 0:
        print(f"Grid: {nglob}")
        print(f"Domain: {L}")
        print()
    
    # Create topology for vector field (lda=3)
    topo = pf.Topology(
        axis=0,
        lda=3,  # 3-component vector field
        nglob=nglob,
        nproc=nproc,
        is_complex=False,
        comm=comm
    )
    
    # Get local information
    nloc = [topo.get_nloc(i) for i in range(3)]
    nmem = [topo.get_nmem(i) for i in range(3)]
    istart = topo.get_istart()
    
    if rank == 0:
        print(f"Topology created (3-component vector)")
        print(f"Local size: {nloc}")
        print()
    
    # Create solver with ROT formulation
    solver = pf.Solver(
        topology=topo,
        boundary_conditions=bc,
        h=h,
        L=L,
        order_diff=pf.DiffType.SPE,  # Spectral derivatives
        center_type=[pf.CenterType.CELL_CENTER] * 3
    )
    
    solver.set_green_type(pf.GreenType.HEJ_2)
    solver.set_alpha(2.0)
    
    if rank == 0:
        print("Setting up solver with ROT formulation...")
    
    solver.setup()
    
    if rank == 0:
        print("Setup complete")
        print()
    
    # Prepare data arrays
    local_size = nloc[0] * nloc[1] * nloc[2] * 3  # 3 components
    
    omega = np.zeros(local_size, dtype=np.float64)  # Vorticity (RHS)
    velocity = np.zeros(local_size, dtype=np.float64)  # Velocity (solution)
    
    # Create coordinate arrays for this rank
    x_local = np.zeros((nloc[0], nloc[1], nloc[2]))
    y_local = np.zeros((nloc[0], nloc[1], nloc[2]))
    z_local = np.zeros((nloc[0], nloc[1], nloc[2]))
    
    for i2 in range(nloc[2]):
        for i1 in range(nloc[1]):
            for i0 in range(nloc[0]):
                x_local[i0, i1, i2] = (istart[0] + i0) * h[0] - L[0] / 2
                y_local[i0, i1, i2] = (istart[1] + i1) * h[1] - L[1] / 2
                z_local[i0, i1, i2] = (istart[2] + i2) * h[2] - L[2] / 2
    
    # Create vortex ring
    omega_x, omega_y, omega_z = create_vortex_ring(x_local, y_local, z_local)
    
    # Fill omega array with proper memory layout
    for i2 in range(nloc[2]):
        for i1 in range(nloc[1]):
            for i0 in range(nloc[0]):
                for comp in range(3):
                    idx = pf.locID(0, i0, i1, i2, comp, topo.axis, nmem, 3)
                    
                    if comp == 0:
                        omega[idx] = omega_x[i0, i1, i2]
                    elif comp == 1:
                        omega[idx] = omega_y[i0, i1, i2]
                    else:
                        omega[idx] = omega_z[i0, i1, i2]
    
    if rank == 0:
        print("Vorticity field initialized")
        print(f"  Max vorticity: {np.max(np.abs(omega)):.4f}")
        print()
    
    # Solve for velocity field from vorticity
    if rank == 0:
        print("Computing velocity field from vorticity...")
    
    comm.Barrier()
    t_start = MPI.Wtime()
    
    solver.solve(velocity, omega, pf.SolverType.ROT)
    
    comm.Barrier()
    t_end = MPI.Wtime()
    
    if rank == 0:
        print(f"  Solved in {t_end - t_start:.4f} seconds")
        print()
    
    # Analyze results
    vel_magnitude = np.zeros(nloc[0] * nloc[1] * nloc[2])
    
    for i2 in range(nloc[2]):
        for i1 in range(nloc[1]):
            for i0 in range(nloc[0]):
                vx_idx = pf.locID(0, i0, i1, i2, 0, topo.axis, nmem, 3)
                vy_idx = pf.locID(0, i0, i1, i2, 1, topo.axis, nmem, 3)
                vz_idx = pf.locID(0, i0, i1, i2, 2, topo.axis, nmem, 3)
                
                flat_idx = i0 + nloc[0] * (i1 + nloc[1] * i2)
                vel_magnitude[flat_idx] = np.sqrt(
                    velocity[vx_idx]**2 + 
                    velocity[vy_idx]**2 + 
                    velocity[vz_idx]**2
                )
    
    # Statistics
    vel_max_local = np.max(vel_magnitude)
    vel_max_global = comm.allreduce(vel_max_local, op=MPI.MAX)
    
    if rank == 0:
        print("Results:")
        print(f"  Max velocity magnitude: {vel_max_global:.6e}")
        print()
        print("="*60)
        print("Simulation completed successfully!")
        print("="*60)


if __name__ == "__main__":
    main()
