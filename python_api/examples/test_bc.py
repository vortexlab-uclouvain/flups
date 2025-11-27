#!/usr/bin/env python3
"""
Test boundary conditions passing
"""

import numpy as np
from mpi4py import MPI
from ctypes import c_int, POINTER, cast
import pyflups as pf

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if rank == 0:
    print("Testing BC creation...")

# Create BC values
bc_vals = [
    [int(pf.BoundaryType.PER), int(pf.BoundaryType.PER)],
    [int(pf.BoundaryType.PER), int(pf.BoundaryType.PER)],
    [int(pf.BoundaryType.PER), int(pf.BoundaryType.PER)],
]

# Method 1: Create 2D array
bc_rows = []
for i in range(3):
    row = (c_int * 2)()
    row[0] = bc_vals[i][0]
    row[1] = bc_vals[i][1]
    bc_rows.append(row)
    if rank == 0:
        print(f"Row {i}: {row[0]}, {row[1]} at {id(row)}")

# Create array of pointers
bc_ptrs = (POINTER(c_int) * 3)(
    cast(bc_rows[0], POINTER(c_int)),
    cast(bc_rows[1], POINTER(c_int)),
    cast(bc_rows[2], POINTER(c_int))
)

if rank == 0:
    print("BC array created successfully")
    for i in range(3):
        print(f"bc_ptrs[{i}][0] = {bc_ptrs[i][0]}, bc_ptrs[{i}][1] = {bc_ptrs[i][1]}")
