"""
FLUPS Python API - Low-level ctypes wrapper
Based on src/flups.h

Copyright © UCLouvain 2020
Licensed under the Apache License, Version 2.0
"""

import ctypes
import os
import sys
import numpy as np
from ctypes import (
    c_void_p, c_int, c_bool, c_double, c_size_t, c_char_p,
    POINTER, Structure, byref
)
from mpi4py import MPI

from .enums import BoundaryType, GreenType, SolverType, DiffType, CenterType


# Find and load the FLUPS library
def find_flups_library():
    """Find the FLUPS shared library"""
    possible_names = [
        'libflups.so',       # Linux
        'libflups.dylib',    # macOS
        'libflups.dll',      # Windows
        'flups.so',
        'flups.dylib',
    ]
    
    # Check in current directory and common paths
    search_paths = [
        os.getcwd(),
        os.path.join(os.getcwd(), 'lib'),
        os.path.join(os.getcwd(), '..', 'lib'),
        '/usr/local/lib',
        '/usr/lib',
    ]
    
    # Check environment variable
    if 'FLUPS_LIB_PATH' in os.environ:
        search_paths.insert(0, os.environ['FLUPS_LIB_PATH'])
    
    for path in search_paths:
        for name in possible_names:
            full_path = os.path.join(path, name)
            if os.path.exists(full_path):
                return full_path
    
    raise FileNotFoundError(
        "Could not find FLUPS library. Please set FLUPS_LIB_PATH environment variable "
        "or place the library in a standard location."
    )


# Load the library
try:
    _lib_path = find_flups_library()
    _lib = ctypes.CDLL(_lib_path) # Load the shared library
except FileNotFoundError as e:
    print(f"Warning: {e}")
    _lib = None


# MPI communicator conversion
def mpi_comm_to_int(comm):
    """Convert MPI4Py communicator to integer handle"""
    if comm is None:
        return MPI.COMM_WORLD.py2f()
    return comm.py2f()


# Opaque structure pointers
class FLUPS_Solver(Structure):
    pass

class FLUPS_Topology(Structure):
    pass

class FLUPS_Profiler(Structure):
    pass


# Define function signatures if library is loaded
if _lib is not None:
    # ========================================================================
    # MEMORY MANAGEMENT
    # ========================================================================
    
    _lib.flups_malloc.argtypes = [c_size_t]
    _lib.flups_malloc.restype = c_void_p
    
    _lib.flups_free.argtypes = [c_void_p]
    _lib.flups_free.restype = None
    
    _lib.flups_info.argtypes = [c_int, POINTER(c_char_p)]
    _lib.flups_info.restype = None
    
    # ========================================================================
    # TOPOLOGY
    # ========================================================================
    
    _lib.flups_topo_new.argtypes = [
        c_int,                    # axis
        c_int,                    # lda
        POINTER(c_int),           # nglob[3]
        POINTER(c_int),           # nproc[3]
        c_bool,                   # isComplex
        POINTER(c_int),           # axproc[3]
        c_int,                    # alignment
        c_int                     # MPI_Comm
    ]
    _lib.flups_topo_new.restype = POINTER(FLUPS_Topology)
    
    _lib.flups_topo_free.argtypes = [POINTER(FLUPS_Topology)]
    _lib.flups_topo_free.restype = None
    
    _lib.flups_topo_get_isComplex.argtypes = [POINTER(FLUPS_Topology)]
    _lib.flups_topo_get_isComplex.restype = c_bool
    
    _lib.flups_topo_get_lda.argtypes = [POINTER(FLUPS_Topology)]
    _lib.flups_topo_get_lda.restype = c_int
    
    _lib.flups_topo_get_axis.argtypes = [POINTER(FLUPS_Topology)]
    _lib.flups_topo_get_axis.restype = c_int
    
    _lib.flups_topo_get_nglob.argtypes = [POINTER(FLUPS_Topology), c_int]
    _lib.flups_topo_get_nglob.restype = c_int
    
    _lib.flups_topo_get_nloc.argtypes = [POINTER(FLUPS_Topology), c_int]
    _lib.flups_topo_get_nloc.restype = c_int
    
    _lib.flups_topo_get_nmem.argtypes = [POINTER(FLUPS_Topology), c_int]
    _lib.flups_topo_get_nmem.restype = c_int
    
    _lib.flups_topo_get_nproc.argtypes = [POINTER(FLUPS_Topology), c_int]
    _lib.flups_topo_get_nproc.restype = c_int
    
    _lib.flups_topo_get_istartGlob.argtypes = [POINTER(FLUPS_Topology), POINTER(c_int)]
    _lib.flups_topo_get_istartGlob.restype = None
    
    _lib.flups_topo_get_locsize.argtypes = [POINTER(FLUPS_Topology)]
    _lib.flups_topo_get_locsize.restype = c_size_t
    
    _lib.flups_topo_get_memsize.argtypes = [POINTER(FLUPS_Topology)]
    _lib.flups_topo_get_memsize.restype = c_size_t
    
    _lib.flups_topo_cmpt_rank_fromid.argtypes = [POINTER(FLUPS_Topology), c_int, c_int]
    _lib.flups_topo_cmpt_rank_fromid.restype = c_int
    
    _lib.flups_topo_cmpt_start_id_from_rank.argtypes = [POINTER(FLUPS_Topology), c_int, c_int]
    _lib.flups_topo_cmpt_start_id_from_rank.restype = c_int
    
    _lib.flups_topo_get_comm.argtypes = [POINTER(FLUPS_Topology)]
    _lib.flups_topo_get_comm.restype = c_int
    
    _lib.flups_topo_ranksplit.argtypes = [POINTER(FLUPS_Topology), c_int, POINTER(c_int)]
    _lib.flups_topo_ranksplit.restype = None
    
    _lib.flups_topo_rankindex.argtypes = [POINTER(FLUPS_Topology), POINTER(c_int)]
    _lib.flups_topo_rankindex.restype = c_int
    
    # ========================================================================
    # SOLVER
    # ========================================================================
    
    _lib.flups_init.argtypes = [
        POINTER(FLUPS_Topology),      # t
        POINTER(POINTER(c_int)),      # bc[3][2]
        POINTER(c_double),            # h[3]
        POINTER(c_double),            # L[3]
        c_int,                        # orderDiff (DiffType)
        POINTER(c_int)                # center_type[3]
    ]
    _lib.flups_init.restype = POINTER(FLUPS_Solver)
    
    _lib.flups_init_timed.argtypes = [
        POINTER(FLUPS_Topology),
        POINTER(POINTER(c_int)),
        POINTER(c_double),
        POINTER(c_double),
        c_int,
        POINTER(c_int),
        POINTER(FLUPS_Profiler)
    ]
    _lib.flups_init_timed.restype = POINTER(FLUPS_Solver)
    
    _lib.flups_cleanup.argtypes = [POINTER(FLUPS_Solver)]
    _lib.flups_cleanup.restype = None
    
    _lib.flups_cleanup_solver.argtypes = [POINTER(FLUPS_Solver)]
    _lib.flups_cleanup_solver.restype = None
    
    _lib.flups_cleanup_all.argtypes = [POINTER(FLUPS_Solver)]
    _lib.flups_cleanup_all.restype = None
    
    _lib.flups_cleanup_fftw.argtypes = []
    _lib.flups_cleanup_fftw.restype = None
    
    _lib.flups_set_greenType.argtypes = [POINTER(FLUPS_Solver), c_int]
    _lib.flups_set_greenType.restype = None
    
    _lib.flups_setup.argtypes = [POINTER(FLUPS_Solver), c_bool]
    _lib.flups_setup.restype = None
    
    _lib.flups_solve.argtypes = [
        POINTER(FLUPS_Solver),
        POINTER(c_double),
        POINTER(c_double),
        c_int
    ]
    _lib.flups_solve.restype = None
    
    # ========================================================================
    # SOLVER (Advanced)
    # ========================================================================
    
    _lib.flups_get_allocSize.argtypes = [POINTER(FLUPS_Solver)]
    _lib.flups_get_allocSize.restype = c_size_t
    
    _lib.flups_get_spectralInfo.argtypes = [
        POINTER(FLUPS_Solver),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double)
    ]
    _lib.flups_get_spectralInfo.restype = None
    
    _lib.flups_set_alpha.argtypes = [POINTER(FLUPS_Solver), c_double]
    _lib.flups_set_alpha.restype = None
    
    _lib.flups_get_innerBuffer.argtypes = [POINTER(FLUPS_Solver)]
    _lib.flups_get_innerBuffer.restype = POINTER(c_double)
    
    _lib.flups_get_innerTopo_physical.argtypes = [POINTER(FLUPS_Solver)]
    _lib.flups_get_innerTopo_physical.restype = POINTER(FLUPS_Topology)
    
    _lib.flups_get_innerTopo_spectral.argtypes = [POINTER(FLUPS_Solver)]
    _lib.flups_get_innerTopo_spectral.restype = POINTER(FLUPS_Topology)
    
    _lib.flups_skip_firstSwitchtopo.argtypes = [POINTER(FLUPS_Solver)]
    _lib.flups_skip_firstSwitchtopo.restype = None
    
    _lib.flups_do_copy.argtypes = [
        POINTER(FLUPS_Solver),
        POINTER(FLUPS_Topology),
        POINTER(c_double),
        c_int
    ]
    _lib.flups_do_copy.restype = None
    
    _lib.flups_do_FFT.argtypes = [POINTER(FLUPS_Solver), POINTER(c_double), c_int]
    _lib.flups_do_FFT.restype = None
    
    _lib.flups_do_mult.argtypes = [POINTER(FLUPS_Solver), POINTER(c_double), c_int]
    _lib.flups_do_mult.restype = None
    
    _lib.flups_pencilDirs.argtypes = [POINTER(POINTER(c_int)), POINTER(c_int)]
    _lib.flups_pencilDirs.restype = None
    
    _lib.flups_switchtopo_info.argtypes = [POINTER(FLUPS_Solver)]
    _lib.flups_switchtopo_info.restype = None
    
    # ========================================================================
    # PROFILER
    # ========================================================================
    
    _lib.flups_profiler_new.argtypes = []
    _lib.flups_profiler_new.restype = POINTER(FLUPS_Profiler)
    
    _lib.flups_profiler_new_n.argtypes = [c_char_p]
    _lib.flups_profiler_new_n.restype = POINTER(FLUPS_Profiler)
    
    _lib.flups_profiler_free.argtypes = [POINTER(FLUPS_Profiler)]
    _lib.flups_profiler_free.restype = None
    
    _lib.flups_profiler_disp.argtypes = [POINTER(FLUPS_Profiler)]
    _lib.flups_profiler_disp.restype = None
    
    # ========================================================================
    # HDF5
    # ========================================================================
    
    _lib.flups_hdf5_dump.argtypes = [
        POINTER(FLUPS_Topology),
        c_char_p,
        POINTER(c_double)
    ]
    _lib.flups_hdf5_dump.restype = None
    
    _lib.flups_print_data.argtypes = [POINTER(FLUPS_Topology), POINTER(c_double)]
    _lib.flups_print_data.restype = None


# Export the library for direct access if needed
def get_lib():
    """Get the underlying ctypes library object"""
    return _lib
