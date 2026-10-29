"""
FLUPS Python API - non-homogeneous Dirichlet / Neumann walls through a right-hand-side correction

FLUPS solves homogeneous symmetric problems: ODD (u = 0 on the wall) and EVEN (du/dn = 0). The
lattice kernels (LGF_*, MEHR_*) invert a known stencil exactly (Solver.get_stencil), so a wall with
data, u = g or du/dn = q, becomes a correction of the right-hand side. Write the field extended to
the ghost nodes as u = R(u) + D, with R(u) the symmetric extension FLUPS solves for and D zero on
every unknown node. Then

    L_h R(u) = M_h f - L_h D    on the unknown nodes,   u = R(u) + D,

and D only involves the boundary data: with n the inward distance to the wall,
    ODD  (Dirichlet): D_0 = g,  D_{-m} = u(-mh) + u(mh)
    EVEN (Neumann):             D_{-m} = u(-mh) - u(mh)
which taylor_ghost() evaluates from the wall Taylor coefficients a_k = d_n^k u. They follow from the
data and the PDE: a_0 = g or a_1 = q, a_{k+2} = d_n^k f - Lap_t a_k (Lap_t: Laplacian along the wall).
Order p needs the terms k < p (Dirichlet ghost O(h^p), Neumann ghost O(h^(p+1))).
See doc/symmetric_conditions.md and examples/nonhomogeneous_bc.py.

The helpers work on node-centred global arrays indexed [x, y, z], with the n = N + 1 nodes per
direction (both wall nodes included, h = L/N), every rank holding the whole domain. Extended arrays
carry r ghost layers on each side (r = stencil_reach), i.e. have n + 2r points per direction.

Copyright © UCLouvain 2020
Licensed under the Apache License, Version 2.0
"""

import math
from typing import Callable, Dict, List

import numpy as np

from .enums import BoundaryType


def stencil_reach(offsets: np.ndarray) -> int:
    """Number of ghost layers a stencil reads beyond a node (0 for an empty stencil)"""
    return int(np.max(np.abs(offsets))) if len(offsets) else 0


def apply_stencil(offsets: np.ndarray, weights: np.ndarray, ext: np.ndarray, r: int) -> np.ndarray:
    """
    Apply a stencil to an extended array

    Parameters
    ----------
    offsets, weights : np.ndarray
        The stencil, as returned by Solver.get_stencil (the weights are used as is: divide by h^2
        for the left-hand side stencil)
    ext : np.ndarray
        Array with r ghost layers on each side
    r : int
        Number of ghost layers of ext, at least the stencil reach

    Returns
    -------
    np.ndarray
        The stencil applied on the core of ext, of shape ext.shape - 2r
    """
    n = [s - 2 * r for s in ext.shape]
    out = np.zeros(n)
    for o, w in zip(offsets, weights):
        out += w * ext[r + o[0]:r + o[0] + n[0], r + o[1]:r + o[1] + n[1], r + o[2]:r + o[2] + n[2]]
    return out


def _check_bc(bc: List[List[BoundaryType]]):
    for d in range(3):
        for side in range(2):
            if bc[d][side] not in (BoundaryType.ODD, BoundaryType.EVEN, BoundaryType.PER):
                raise ValueError(f"direction {d}: only ODD, EVEN and PER walls take a correction")
        if (bc[d][0] == BoundaryType.PER) != (bc[d][1] == BoundaryType.PER):
            raise ValueError(f"direction {d}: PER must be set on both sides")


def _plane(d: int, i: int) -> tuple:
    s = [slice(None)] * 3
    s[d] = i
    return tuple(s)


def fixed_mask(bc: List[List[BoundaryType]], n: int) -> np.ndarray:
    """Nodes that are not solved for: the ODD wall nodes and the duplicated last node of PER directions"""
    _check_bc(bc)
    fixed = np.zeros((n, n, n), dtype=bool)
    for d in range(3):
        if bc[d][0] == BoundaryType.PER:
            fixed[_plane(d, n - 1)] = True
        for side in range(2):
            if bc[d][side] == BoundaryType.ODD:
                fixed[_plane(d, 0 if side == 0 else n - 1)] = True
    return fixed


def ghost_remainder(bc: List[List[BoundaryType]], n: int, r: int,
                    ghost: Callable[[int, int, int], np.ndarray]) -> np.ndarray:
    """
    Build D = u - R(u) on the extended grid from the boundary data

    D is built one direction at a time (x, then y, then z), each pass reading the values written by
    the previous ones: D_{-m} = -D_{+m} + G_m on ODD walls and D_{-m} = D_{+m} + N_m on EVEN walls.
    This gives the edge and corner ghosts as double/triple reflections.

    Parameters
    ----------
    bc : list of lists
        [[xmin, xmax], [ymin, ymax], [zmin, zmax]] with ODD (Dirichlet), EVEN (Neumann) or PER
    n : int
        Number of nodes per direction (N + 1)
    r : int
        Number of ghost layers (at least the reach of the left-hand side stencil)
    ghost : callable
        ghost(d, side, m) returns the data on the wall (side 0: min, 1: max) normal to d, as an
        array broadcastable to the extended face, of shape (n + 2r, n + 2r) over the two other
        directions in increasing order:
          - m = 0, ODD walls only: g, the wall value;
          - m >= 1, ODD walls: G_m = u(-mh) + u(mh); EVEN walls: N_m = u(-mh) - u(mh),
            with n the inward distance to the wall (see taylor_ghost).
        The face covers the ghost positions of the neighbouring walls. Only stencils with diagonal
        terms (MEHR_*) read them, and then use the data slightly beyond the face; the cross-shaped
        LGF_* stencils never do.

    Returns
    -------
    np.ndarray
        D on the extended grid, of shape (n + 2r,)*3, zero on the unknown nodes
    """
    _check_bc(bc)
    N = n - 1
    face = (n + 2 * r, n + 2 * r)
    D = np.zeros((n + 2 * r,) * 3)
    for d in range(3):
        if bc[d][0] == BoundaryType.PER:
            for m in range(1, r + 1):
                D[_plane(d, r - m)] = D[_plane(d, r + N - m)]
                D[_plane(d, r + N + m)] = D[_plane(d, r + m)]
            continue
        for side in range(2):
            odd = bc[d][side] == BoundaryType.ODD
            iw, step = (r, -1) if side == 0 else (r + N, 1)
            if odd:
                D[_plane(d, iw)] = np.broadcast_to(ghost(d, side, 0), face)
            for m in range(1, r + 1):
                refl = -1.0 if odd else 1.0
                D[_plane(d, iw + step * m)] = refl * D[_plane(d, iw - step * m)] + np.broadcast_to(ghost(d, side, m), face)
    return D


def taylor_ghost(coefs: Dict[int, np.ndarray], m: int, h: float, wall: BoundaryType) -> np.ndarray:
    """
    Ghost data from the Taylor coefficients of u on the wall

    Parameters
    ----------
    coefs : dict
        {k: a_k}, a_k = d_n^k u on the wall (n inward): the even k for an ODD wall (a_0 = g), the odd
        k for an EVEN wall (a_1 = q). Order p needs the terms k < p.
    m : int
        Ghost layer (m = 0 returns g for an ODD wall)
    h : float
        Grid spacing
    wall : BoundaryType
        ODD (Dirichlet) or EVEN (Neumann)

    Returns
    -------
    np.ndarray
        G_m = 2 sum_{k even} (mh)^k/k! a_k for ODD, N_m = -2 sum_{k odd} (mh)^k/k! a_k for EVEN
    """
    if wall == BoundaryType.ODD and m == 0:
        return coefs[0]
    parity, sign = (0, 2.0) if wall == BoundaryType.ODD else (1, -2.0)
    out = 0.0
    for k, a in coefs.items():
        if k % 2 == parity:
            out = out + sign * (m * h) ** k / math.factorial(k) * a
    return out


def corrected_rhs(offsets: np.ndarray, weights: np.ndarray, rhs: np.ndarray, remainder: np.ndarray, h: float,
                  bc: List[List[BoundaryType]]) -> np.ndarray:
    """
    Right-hand side to give to FLUPS: rhs - L_h D on the unknown nodes, zero on the fixed ones

    Parameters
    ----------
    offsets, weights : np.ndarray
        The left-hand side stencil (Solver.get_stencil(StencilType.LHS))
    rhs : np.ndarray
        M_h f on the n^3 nodes (f itself for the LGF kernels)
    remainder : np.ndarray
        D from ghost_remainder
    h : float
        Grid spacing
    bc : list of lists
        The boundary conditions given to the solver
    """
    n = rhs.shape[0]
    r = (remainder.shape[0] - n) // 2
    out = rhs - apply_stencil(offsets, weights, remainder, r) / h**2
    out[fixed_mask(bc, n)] = 0.0
    return out


def recover(v: np.ndarray, remainder: np.ndarray, bc: List[List[BoundaryType]]) -> np.ndarray:
    """
    Solution u = R(u) + D from the FLUPS output v = R(u), on the n^3 nodes

    The ODD wall nodes are set to the wall data and the duplicated PER node to the first one.
    """
    n = v.shape[0]
    r = (remainder.shape[0] - n) // 2
    u = np.where(fixed_mask(bc, n), 0.0, v)
    for d in range(3):
        if bc[d][0] == BoundaryType.PER:
            u[_plane(d, n - 1)] = u[_plane(d, 0)]
    return u + remainder[r:r + n, r:r + n, r:r + n]
