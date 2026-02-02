"""
FLUPS Python API - High-level Object-Oriented Interface

Copyright © UCLouvain 2020
Licensed under the Apache License, Version 2.0
"""

import numpy as np
from ctypes import c_int, c_double, c_bool, POINTER, byref, cast, c_void_p
from typing import Optional, Tuple, List, Union
from mpi4py import MPI

from . import flups_wrapper as fw
from .enums import (
    BoundaryType, GreenType, SolverType, DiffType, CenterType,
    FLUPS_FORWARD, FLUPS_BACKWARD, FLUPS_ALIGNMENT
)

# Global counter for active solvers
_active_solvers_count = 0


class Topology:
    """
    FLUPS Topology class
    
    Manages the domain decomposition and memory layout for parallel FFT operations.
    """
    
    def __init__(self, 
                 axis: int,
                 lda: int,
                 nglob: Union[List[int], np.ndarray],
                 nproc: Union[List[int], np.ndarray],
                 is_complex: bool = False,
                 axproc: Optional[Union[List[int], np.ndarray]] = None,
                 alignment: int = FLUPS_ALIGNMENT,
                 comm: Optional[MPI.Comm] = None,
                 _owned: bool = True):
        """
        Create a new topology
        
        Parameters
        ----------
        axis : int
            The direction aligned with the fastest rotating index (0, 1, or 2)
        lda : int
            Leading dimension of the array (number of components for vector field)
        nglob : array-like
            Global number of points in each direction [nx, ny, nz]
        nproc : array-like
            Number of processors per direction [npx, npy, npz]
        is_complex : bool, optional
            Whether the topology handles complex numbers (default: False)
        axproc : array-like, optional
            Correspondence between physical dimensions and rank decomposition
        alignment : int, optional
            Memory alignment constant (default: FLUPS_ALIGNMENT)
        comm : MPI.Comm, optional
            MPI communicator (default: MPI.COMM_WORLD)
        _owned : bool, optional
            Internal flag: if True, this object owns the topology and will free it.
            If False, the topology is owned by another object (e.g., Solver)
        """
        self._topo = None
        self._owned = _owned
        
        # Convert inputs to ctypes arrays
        nglob_arr = (c_int * 3)(*nglob)
        nproc_arr = (c_int * 3)(*nproc)
        
        if axproc is not None:
            axproc_arr = (c_int * 3)(*axproc)
            axproc_ptr = byref(axproc_arr)
        else:
            # Pass NULL pointer when axproc is not provided
            axproc_ptr = cast(None, POINTER(c_int))
        
        # Convert MPI communicator
        if comm is None:
            comm = MPI.COMM_WORLD
        mpi_comm = fw.mpi_comm_to_int(comm)
        
        # Create topology - with error handling for MPI incompatibilities
        try:
            self._topo = fw._lib.flups_topo_new(
                axis, lda, nglob_arr, nproc_arr, 
                is_complex, axproc_ptr, alignment, mpi_comm
            )
        except Exception as e:
            raise RuntimeError(
                f"Failed to create FLUPS topology. This may be due to MPI version mismatch.\n"
                f"FLUPS was likely compiled with a different MPI than the current environment.\n"
                f"Original error: {e}"
            ) from e
        
        if not self._topo:
            raise RuntimeError(
                "Failed to create FLUPS topology. Null pointer returned. "
                "This usually indicates an MPI communicator issue."
            )
    
    def __del__(self):
        """Destructor to free topology only if owned by this object"""
        if self._topo is not None and self._owned:
            fw._lib.flups_topo_free(self._topo)
            self._topo = None
    
    @property
    def is_complex(self) -> bool:
        """Check if topology works on complex numbers"""
        return fw._lib.flups_topo_get_isComplex(self._topo)
    
    @property
    def lda(self) -> int:
        """Get the leading dimension (number of vector components)"""
        return fw._lib.flups_topo_get_lda(self._topo)
    
    @property
    def axis(self) -> int:
        """Get the physical direction aligned in memory"""
        return fw._lib.flups_topo_get_axis(self._topo)
    
    def get_nglob(self, dim: int) -> int:
        """Get global number of points in a given direction"""
        return fw._lib.flups_topo_get_nglob(self._topo, dim)
    
    def get_nloc(self, dim: int) -> int:
        """Get local number of points in a given direction"""
        return fw._lib.flups_topo_get_nloc(self._topo, dim)
    
    def get_nmem(self, dim: int) -> int:
        """Get local memory size per direction"""
        return fw._lib.flups_topo_get_nmem(self._topo, dim)
    
    def get_nproc(self, dim: int) -> int:
        """Get number of processes in a given direction"""
        return fw._lib.flups_topo_get_nproc(self._topo, dim)
    
    def get_istart(self) -> np.ndarray:
        """Get start index of this process in all 3 directions"""
        istart = (c_int * 3)()
        fw._lib.flups_topo_get_istartGlob(self._topo, istart)
        return np.array([istart[0], istart[1], istart[2]])
    
    @property
    def locsize(self) -> int:
        """Get scalar local size on this rank"""
        return fw._lib.flups_topo_get_locsize(self._topo)
    
    @property
    def memsize(self) -> int:
        """Get memory size on this rank (in bytes)"""
        return fw._lib.flups_topo_get_memsize(self._topo)
    
    def cmpt_rank_from_id(self, global_id: int, direction: int) -> int:
        """Compute rank associated to a scalar global id"""
        return fw._lib.flups_topo_cmpt_rank_fromid(self._topo, global_id, direction)
    
    def cmpt_start_id_from_rank(self, rank_id: int, direction: int) -> int:
        """Compute starting index of a rank"""
        return fw._lib.flups_topo_cmpt_start_id_from_rank(self._topo, rank_id, direction)
    
    def get_comm(self) -> MPI.Comm:
        """Get the MPI communicator"""
        comm_int = fw._lib.flups_topo_get_comm(self._topo)
        return MPI.Comm.f2py(comm_int)
    
    def ranksplit(self, rank: int) -> np.ndarray:
        """Split rank into rank per dimensions"""
        rankd = (c_int * 3)()
        fw._lib.flups_topo_ranksplit(self._topo, rank, rankd)
        return np.array([rankd[0], rankd[1], rankd[2]])
    
    def rankindex(self, rankd: Union[List[int], np.ndarray]) -> int:
        """Get rank from rank per dimension"""
        rankd_arr = (c_int * 3)(*rankd)
        return fw._lib.flups_topo_rankindex(self._topo, rankd_arr)
    
    @property
    def _ptr(self):
        """Get the underlying C pointer (for internal use)"""
        return self._topo


class Profiler:
    """FLUPS Profiler for timing measurements"""
    
    def __init__(self, name: str = "default"):
        """
        Create a new profiler
        
        Parameters
        ----------
        name : str, optional
            Name of the profiler (default: "default")
        """
        if name == "default":
            self._prof = fw._lib.flups_profiler_new()
        else:
            self._prof = fw._lib.flups_profiler_new_n(name.encode('utf-8'))
        
        if not self._prof:
            raise RuntimeError("Failed to create FLUPS profiler")
    
    def __del__(self):
        """Destructor to free profiler"""
        if self._prof is not None:
            fw._lib.flups_profiler_free(self._prof)
            self._prof = None
    
    def display(self):
        """Display profiler information"""
        fw._lib.flups_profiler_disp(self._prof)
    
    @property
    def _ptr(self):
        """Get the underlying C pointer (for internal use)"""
        return self._prof


class Solver:
    """
    FLUPS Solver class
    
    Main class for solving Poisson equations using FFT-based methods.
    """
    
    def __init__(self,
                 topology: Topology,
                 boundary_conditions: List[List[BoundaryType]],
                 h: Union[List[float], np.ndarray],
                 L: Union[List[float], np.ndarray],
                 order_diff: DiffType = DiffType.NOD,
                 center_type: Union[List[CenterType], np.ndarray] = None,
                 profiler: Optional[Profiler] = None):
        """
        Create a new FLUPS solver
        
        Parameters
        ----------
        topology : Topology
            User-determined topology of data in physical space
        boundary_conditions : list of lists
            Boundary conditions [[xmin, xmax], [ymin, ymax], [zmin, zmax]]
            where each element is a BoundaryType
        h : array-like
            Physical space increment in each direction [hx, hy, hz]
        L : array-like
            Physical length of the domain in each direction [Lx, Ly, Lz]
        order_diff : DiffType, optional
            Order of derivatives for ROT solver (default: NOD)
        center_type : array-like, optional
            Center type for each direction (default: [CELL_CENTER]*3)
        profiler : Profiler, optional
            Profiler for timing measurements
        """
        self._solver = None
        self._topology = topology
        
        # Prepare boundary conditions array
        # bc is BoundaryType* [3][2], which is a 3x2 array of pointers
        # Each pointer points to an array of lda BoundaryType values
        lda = topology.lda
        
        # Create 6 separate arrays (3x2), each of size lda
        bc_arrays = []
        for i in range(3):
            for j in range(2):
                # Create an array of lda elements
                bc_array = (c_int * lda)()
                # Fill with the boundary condition value
                for comp in range(lda):
                    bc_array[comp] = int(boundary_conditions[i][j])
                bc_arrays.append(bc_array)
        
        # Create array of pointers to these arrays
        bc_ptrs_flat = []
        for bc_array in bc_arrays:
            ptr = cast(bc_array, POINTER(c_int))
            bc_ptrs_flat.append(ptr)
        
        # Cast to array of 6 pointers
        bc_c = (POINTER(c_int) * 6)(*bc_ptrs_flat)
        
        # Prepare other arrays
        h_arr = (c_double * 3)(*h)
        L_arr = (c_double * 3)(*L)
        
        if center_type is None:
            center_type = [CenterType.CELL_CENTER] * 3
        center_arr = (c_int * 3)(*[int(ct) for ct in center_type])
        
        # Create solver
        if profiler is None:
            self._solver = fw._lib.flups_init(
                topology._ptr, bc_c, h_arr, L_arr,
                int(order_diff), center_arr
            )
        else:
            self._solver = fw._lib.flups_init_timed(
                topology._ptr, bc_c, h_arr, L_arr,
                int(order_diff), center_arr, profiler._ptr
            )
        
        if not self._solver:
            raise RuntimeError("Failed to create FLUPS solver")
        
        self._bc_arrays = bc_arrays  # Keep reference to prevent garbage collection
        self._bc_c = bc_c
        
        # Increment active solver count
        global _active_solvers_count
        _active_solvers_count += 1
    
    def __del__(self):
        """Destructor to free solver
        
        When the last solver is destroyed, also clean up backend resources
        to allow multiple solvers to coexist safely.
        """
        global _active_solvers_count
        
        if self._solver is not None:
            # Delete the solver object
            fw._lib.flups_cleanup_solver(self._solver)
            self._solver = None
            
            # Decrement counter
            _active_solvers_count -= 1
            
            # Clean up backend when last solver is destroyed
            if _active_solvers_count <= 0:
                _active_solvers_count = 0  # Prevent negative count
                try:
                    fw._lib.flups_cleanup_backend()
                except:
                    pass  # Ignore errors during backend cleanup
    
    def set_green_type(self, green_type: GreenType):
        """
        Set the type of Green's function
        
        Must be called before setup()
        
        Parameters
        ----------
        green_type : GreenType
            Type of Green's function to use
        """
        fw._lib.flups_set_greenType(self._solver, int(green_type))
    
    def set_alpha(self, alpha: float):
        """
        Set alpha factor for regularized Hejlesen kernels
        
        Must be called before setup()
        
        Parameters
        ----------
        alpha : float
            Number of grid points in the smoothing Gaussian (default: 2.0)
        """
        fw._lib.flups_set_alpha(self._solver, c_double(alpha))
    
    def setup(self, change_comm: bool = False):
        """
        Setup the solver and allocate memory
        
        After this call, the solver cannot be changed anymore.
        
        Parameters
        ----------
        change_comm : bool, optional
            Allow FLUPS to change the communicator (default: False)
        """
        fw._lib.flups_setup(self._solver, change_comm)
    
    def solve(self, 
              field: np.ndarray,
              rhs: np.ndarray,
              solver_type: SolverType = SolverType.STD) -> np.ndarray:
        """
        Solve the Poisson equation
        
        Parameters
        ----------
        field : np.ndarray
            Output array for the solution (will be modified in-place if same as rhs)
        rhs : np.ndarray
            Right-hand side of the equation
        solver_type : SolverType, optional
            Type of solver to use (default: STD)
        
        Returns
        -------
        np.ndarray
            Solution field (same as input field array)
        """
        # Ensure arrays are contiguous and double precision
        if not field.flags['C_CONTIGUOUS']:
            field = np.ascontiguousarray(field, dtype=np.float64)
        if not rhs.flags['C_CONTIGUOUS']:
            rhs = np.ascontiguousarray(rhs, dtype=np.float64)
        
        field_ptr = field.ctypes.data_as(POINTER(c_double))
        rhs_ptr = rhs.ctypes.data_as(POINTER(c_double))
        
        fw._lib.flups_solve(self._solver, field_ptr, rhs_ptr, int(solver_type))
        
        return field
    
    def get_alloc_size(self) -> int:
        """Get maximum amount of memory required by FLUPS"""
        return fw._lib.flups_get_allocSize(self._solver)
    
    def get_spectral_info(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Get information for computing spectral modes
        
        Returns
        -------
        kfact : np.ndarray
            Multiplication factor for spectral modes
        koffset : np.ndarray
            Spectral offset based on boundary conditions
        symstart : np.ndarray
            First point which is symmetrized
        """
        kfact = (c_double * 3)()
        koffset = (c_double * 3)()
        symstart = (c_double * 3)()
        
        fw._lib.flups_get_spectralInfo(self._solver, kfact, koffset, symstart)
        
        return (np.array([kfact[0], kfact[1], kfact[2]]),
                np.array([koffset[0], koffset[1], koffset[2]]),
                np.array([symstart[0], symstart[1], symstart[2]]))
    
    def get_inner_buffer(self) -> np.ndarray:
        """Get the inner buffer used by FLUPS"""
        ptr = fw._lib.flups_get_innerBuffer(self._solver)
        size = self.get_alloc_size()
        # Create numpy array from pointer
        return np.ctypeslib.as_array(ptr, shape=(size,))
    
    def get_inner_topo_physical(self) -> Topology:
        """
        Get the first pencil topology in physical space
        
        Warning: This topology is owned by the solver and will be freed when the solver is destroyed.
        Do not call __del__ on this topology.
        """
        topo_ptr = fw._lib.flups_get_innerTopo_physical(self._solver)
        # Create a Topology wrapper with _owned=False (don't free this as it's owned by solver)
        topo = Topology.__new__(Topology)
        topo._topo = topo_ptr
        topo._owned = False  # Mark as not owned so __del__ won't free it
        return topo
    
    def get_inner_topo_spectral(self) -> Topology:
        """
        Get the spectral topology
        
        Warning: This topology is owned by the solver and will be freed when the solver is destroyed.
        Do not call __del__ on this topology.
        """
        topo_ptr = fw._lib.flups_get_innerTopo_spectral(self._solver)
        topo = Topology.__new__(Topology)
        topo._topo = topo_ptr
        topo._owned = False  # Mark as not owned so __del__ won't free it
        return topo
    
    def skip_first_switchtopo(self):
        """Instruct FLUPS to skip the first SwitchTopology operation"""
        fw._lib.flups_skip_firstSwitchtopo(self._solver)
    
    def do_copy(self, topology: Topology, data: np.ndarray, sign: int):
        """Copy data from user to FLUPS arrays"""
        data_ptr = data.ctypes.data_as(POINTER(c_double))
        fw._lib.flups_do_copy(self._solver, topology._ptr, data_ptr, sign)
    
    def do_fft(self, data: np.ndarray, sign: int):
        """Compute FFT (physical to spectral space)"""
        data_ptr = data.ctypes.data_as(POINTER(c_double))
        fw._lib.flups_do_FFT(self._solver, data_ptr, sign)
    
    def do_mult(self, data: np.ndarray, solver_type: SolverType):
        """Multiply with Green's function"""
        data_ptr = data.ctypes.data_as(POINTER(c_double))
        fw._lib.flups_do_mult(self._solver, data_ptr, int(solver_type))
    
    def switchtopo_info(self):
        """Print information about SwitchTopo"""
        fw._lib.flups_switchtopo_info(self._solver)
    
    @staticmethod
    def cleanup_backend():
        """Free data employed by backend (call after all solvers are freed)"""
        fw._lib.flups_cleanup_backend()


def hdf5_dump(topology: Topology, filename: str, data: np.ndarray):
    """
    Dump data in HDF5 format with XDMF description
    
    Parameters
    ----------
    topology : Topology
        Topology describing the data layout
    filename : str
        Output filename (without extension)
    data : np.ndarray
        Data to dump
    """
    data_ptr = data.ctypes.data_as(POINTER(c_double))
    fw._lib.flups_hdf5_dump(topology._ptr, filename.encode('utf-8'), data_ptr)


def print_data(topology: Topology, data: np.ndarray):
    """
    Print data to console
    
    Parameters
    ----------
    topology : Topology
        Topology describing the data layout
    data : np.ndarray
        Data to print
    """
    data_ptr = data.ctypes.data_as(POINTER(c_double))
    fw._lib.flups_print_data(topology._ptr, data_ptr)


def locID(axsrc: int, i0: int, i1: int, i2: int, lia: int, 
          axtrg: int, size: List[int], nf: int) -> int:
    """
    Compute memory local index for a point
    
    Parameters
    ----------
    axsrc : int
        FRI, reference axis aligned with index i0
    i0, i1, i2 : int
        Indices in respective directions
    lia : int
        Index of vector component
    axtrg : int
        Topology FRI (memory alignment)
    size : list
        Size of memory (in 012-order)
    nf : int
        Number of unknowns in one element
    
    Returns
    -------
    int
        Linear memory index
    """
    i = [i0, i1, i2]
    dax0 = (3 + axtrg - axsrc) % 3
    dax1 = (dax0 + 1) % 3
    dax2 = (dax0 + 2) % 3
    ax0 = axtrg
    ax1 = (ax0 + 1) % 3
    ax2 = (ax0 + 2) % 3
    
    return i[dax0] * nf + size[ax0] * nf * (i[dax1] + size[ax1] * (i[dax2] + lia * size[ax2]))
