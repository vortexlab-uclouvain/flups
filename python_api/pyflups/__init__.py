"""
PyFLUPS - Python API for FLUPS Library

FLUPS is a Fourier-based Library of Unbounded Poisson Solvers.

Copyright © UCLouvain 2020
Licensed under the Apache License, Version 2.0
"""

__version__ = "1.0.0"
__author__ = "Pierre Balty, Thomas Gillis, Denis-Gabriel Caprace, and contributors"

from .enums import (
    BoundaryType,
    GreenType,
    SolverType,
    DiffType,
    CenterType,
    FLUPS_FORWARD,
    FLUPS_BACKWARD,
    FLUPS_BACKWARD_DIFF,
    FLUPS_ALIGNMENT,
)

from .flups import (
    Topology,
    Solver,
    Profiler,
    hdf5_dump,
    print_data,
    locID,
)

__all__ = [
    # Classes
    'Topology',
    'Solver',
    'Profiler',

    # Enums
    'BoundaryType',
    'GreenType',
    'SolverType',
    'DiffType',
    'CenterType',

    # Constants
    'FLUPS_FORWARD',
    'FLUPS_BACKWARD',
    'FLUPS_BACKWARD_DIFF',
    'FLUPS_ALIGNMENT',

    # Functions
    'hdf5_dump',
    'print_data',
    'locID',
]
