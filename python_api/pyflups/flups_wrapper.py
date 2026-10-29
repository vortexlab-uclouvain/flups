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
    POINTER, Structure, byref, cast, CFUNCTYPE
)
from mpi4py import MPI

from .enums import BoundaryType, GreenType, SolverType, DiffType, CenterType

# SwitchTopo communication backends, each built as its own library by `make python`
BACKENDS = ('a2a', 'nb', 'isr')


# Find and load the FLUPS library
def find_flups_library():
    """
    Find the FLUPS shared library

    The environment variable FLUPS_LIB gives either the full path to the library or the
    communication backend to use (a2a, nb or isr). The libraries are searched for in
    FLUPS_LIB_PATH first, then in the lib folder of this FLUPS checkout, the current directory,
    ./lib, ../lib, /usr/local/lib and /usr/lib. Without FLUPS_LIB, libflups is preferred, then
    libflups_nb, libflups_a2a and libflups_isr.
    """
    choice = os.environ.get('FLUPS_LIB', '')
    if choice and choice not in BACKENDS:
        if os.path.isfile(choice):
            return choice
        raise ImportError(f"FLUPS_LIB={choice} is neither an existing file nor one of the backends {BACKENDS}")

    stems = [f'libflups_{choice}'] if choice else ['libflups', 'libflups_nb', 'libflups_a2a', 'libflups_isr']
    exts = {'darwin': ['.dylib', '.so'], 'win32': ['.dll']}.get(sys.platform, ['.so'])
    search_paths = [
        os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'lib'),
        os.getcwd(),
        os.path.join(os.getcwd(), 'lib'),
        os.path.join(os.getcwd(), '..', 'lib'),
        '/usr/local/lib',
        '/usr/lib',
    ]
    if 'FLUPS_LIB_PATH' in os.environ:
        search_paths.insert(0, os.environ['FLUPS_LIB_PATH'])

    for stem in stems:
        for path in search_paths:
            for ext in exts:
                full_path = os.path.join(path, stem + ext)
                if os.path.isfile(full_path):
                    return os.path.normpath(full_path)

    raise ImportError(
        "Could not find the FLUPS library " + " or ".join(s + exts[0] for s in stems) + " in:\n  "
        + "\n  ".join(os.path.normpath(p) for p in search_paths)
        + "\nBuild it with `make python` and set FLUPS_LIB_PATH to its folder, or FLUPS_LIB to its path."
    )


# Load the library
_lib_path = find_flups_library()
_lib = ctypes.CDLL(_lib_path)


# MPI communicator conversion: MPI_Comm is a pointer in Open MPI and an int in MPICH
MPI_Comm_t = c_void_p if MPI._sizeof(MPI.Comm) == ctypes.sizeof(c_void_p) else c_int


def mpi_comm_handle(comm):
    """C handle (MPI_Comm) of an mpi4py communicator, to be passed to FLUPS (default: MPI.COMM_WORLD)"""
    if comm is None:
        comm = MPI.COMM_WORLD
    return MPI_Comm_t(MPI._handleof(comm))


def mpi_comm_from_handle(handle):
    """mpi4py communicator wrapping a C handle (MPI_Comm) returned by FLUPS"""
    if hasattr(MPI.Comm, 'fromhandle'):
        return MPI.Comm.fromhandle(handle)
    comm = MPI.Comm()
    MPI_Comm_t.from_address(MPI._addressof(comm)).value = handle
    return comm


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
        MPI_Comm_t                # MPI_Comm
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
    _lib.flups_topo_get_comm.restype = MPI_Comm_t
    
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
    
    _lib.flups_cleanup_backend.argtypes = []
    _lib.flups_cleanup_backend.restype = None
    
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

    _lib.flups_get_stencil.argtypes = [
        POINTER(FLUPS_Solver),
        c_int,                        # StencilType
        POINTER(c_int),               # offsets[3*n] or NULL
        POINTER(c_double)             # weights[n] or NULL
    ]
    _lib.flups_get_stencil.restype = c_int

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
