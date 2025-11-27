#!/usr/bin/env python3
"""
Test simple topology only
"""

import numpy as np
from mpi4py import MPI
import pyflups as pf

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if rank == 0:
    print("Testing topology creation...")

# Parameters
nglob = [64, 64, 64]
nproc = [2, 2, 1]

# Create topology
topo = pf.Topology(
    axis=0,
    lda=1,
    nglob=nglob,
    nproc=nproc,
    is_complex=False,
    comm=comm
)

if rank == 0:
    print(f"Topology created successfully")
    print(f"  axis: {topo.axis}")
    print(f"  lda: {topo.lda}")
    print(f"  nglob: [{topo.get_nglob(0)}, {topo.get_nglob(1)}, {topo.get_nglob(2)}]")
    print(f"  nloc: [{topo.get_nloc(0)}, {topo.get_nloc(1)}, {topo.get_nloc(2)}]")
    print(f"  isComplex: {topo.is_complex}")

if rank == 0:
    print("Test passed!")
