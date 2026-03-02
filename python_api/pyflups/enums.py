"""
FLUPS Python API - Enumerations
Based on flups_interface.h

Copyright © UCLouvain 2020
Licensed under the Apache License, Version 2.0
"""

from enum import IntEnum


# Constants
FLUPS_FORWARD = -1  # equivalent to FFTW_FORWARD
FLUPS_BACKWARD = 1  # equivalent to FFTW_BACKWARD
FLUPS_BACKWARD_DIFF = 2
FLUPS_ALIGNMENT = 16


class BoundaryType(IntEnum):
    """List of supported boundary conditions"""
    EVEN = 0  # EVEN boundary condition = zero flux
    ODD = 1   # ODD boundary condition = zero value
    PER = 3   # PERiodic boundary conditions
    UNB = 4   # UNBounded boundary condition
    NONE = 9  # No boundary condition = dimension not used


class GreenType(IntEnum):
    """The type of Green's function used for the Poisson solver"""
    CHAT_2 = 0   # quadrature in zero, order 2, Chatelain et al. (2010)
    LGF_2 = 1    # Lattice Green's function, order 2, Gillis et al. (2018)
    HEJ_2 = 2    # regularized, order 2, Hejlesen et al. (2015)
    HEJ_4 = 3    # regularized, order 4, Hejlesen et al. (2015)
    HEJ_6 = 4    # regularized, order 6, Hejlesen et al. (2015)
    HEJ_8 = 5    # regularized, order 8, Hejlesen et al. (2015)
    HEJ_10 = 6   # regularized, order 10, Hejlesen et al. (2015)
    HEJ_0 = 7    # Fourier cutoff, spectral-like, Hejlesen et al. (2019)
    LGF_4 = 8    # Lattice Green's function, order 4
    LGF_6 = 9    # Lattice Green's function, order 6
    LGF_8 = 10   # Lattice Green's function, order 8
    MEHR_4L = 11 # Spectral equivalent of the left-hand side of the Mehrstellen HOC stencil, order 4
    MEHR_6L = 12 # Spectral equivalent of the left-hand side of the Mehrstellen HOC stencil, order 6
    MEHR_4F = 13 # Spectral equivalent of the full Mehrstellen HOC stencil, order 4
    MEHR_6F = 14 # Spectral equivalent of the full Mehrstellen HOC stencil, order 6


class SolverType(IntEnum):
    """The type of possible solvers"""
    STD = 0  # the standard poisson solver: ∇²(φ) = (rhs)
    ROT = 1  # the Bio-Savart poisson solver: ∇²(φ) = ∇ × (rhs)


class DiffType(IntEnum):
    """The type of derivative to be used with SolverType ROT"""
    NOD = 0  # Default parameter to be used with the STD type solve
    SPE = 1  # Spectral derivation
    FD2 = 2  # Spectral equivalent of 2nd order finite difference
    FD4 = 4  # Spectral equivalent of 4th order finite difference
    FD6 = 6  # Spectral equivalent of 6th order finite difference


class CenterType(IntEnum):
    """List of supported data center"""
    NODE_CENTER = 0  # NODE centered data
    CELL_CENTER = 1  # CELL centered data
