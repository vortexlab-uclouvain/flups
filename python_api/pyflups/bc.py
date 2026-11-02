"""
FLUPS Python API - non-homogeneous Dirichlet / Neumann walls through a right-hand-side correction

FLUPS solves homogeneous symmetric problems: ODD (u = 0 on the wall) and EVEN (du/dn = 0). The
lattice kernels (LGF_*, MEHR_*) invert a known stencil exactly (Solver.get_lhs_stencil), so a wall with
data, u = g or du/dn = q, becomes a correction of the right-hand side. Write the field extended to
the ghost nodes as u = R(u) + D, with R(u) the symmetric extension FLUPS solves for and D zero on
every unknown node. Then

    L_h R(u) = M_h f - L_h D    on the unknown nodes,   u = R(u) + D,

and D only involves the boundary data: with n the inward distance to the wall and s_m the distance
of the ghost layer m (s_m = m h for node-centred data, (m - 1/2) h for cell-centred data),
    ODD  (Dirichlet): D_0 = g (node-centred only),  D_{-m} = u(-s_m) + u(s_m)
    EVEN (Neumann):                                 D_{-m} = u(-s_m) - u(s_m)
which taylor_ghost() evaluates from the wall Taylor coefficients a_k = d_n^k u. They follow from the
data and the PDE: a_0 = g or a_1 = q, a_{k+2} = d_n^k f - Lap_t a_k (Lap_t: Laplacian along the wall).
Order p needs the terms k < p (Dirichlet ghost O(h^p), Neumann ghost O(h^(p+1))).
See doc/nonhomogeneous_conditions.md and examples/nonhomogeneous_bc.py.

The helpers work on global arrays indexed [x, y, z], every rank holding the whole domain, with
h = L/N and n points per direction:
    - node-centred data (default): n = N + 1 nodes at x_i = i h, both wall nodes included. The ODD
      wall nodes are not solved for (FLUPS returns 0 there) and hold g in the solution;
    - cell-centred data (center=CenterType.CELL_CENTER): n = N nodes at x_i = (i + 1/2) h, every one
      unknown. The walls lie half-way between the first node and the first ghost, so the ghost
      layer m mirrors the node m - 1, and a Dirichlet wall always needs Taylor ghosts.
Extended arrays carry r ghost layers on each side (r = stencil_reach), i.e. have n + 2r points per
direction.

Copyright © UCLouvain 2020
Licensed under the Apache License, Version 2.0
"""

import math
from typing import Callable, Dict, List

import numpy as np

from .enums import BoundaryType, CenterType


def stencil_reach(offsets: np.ndarray) -> int:
    """Number of ghost layers a stencil reads beyond a node (0 for an empty stencil)"""
    return int(np.max(np.abs(offsets))) if len(offsets) else 0


def apply_stencil(offsets: np.ndarray, weights: np.ndarray, ext: np.ndarray, r: int) -> np.ndarray:
    """
    Apply a stencil to an extended array

    Parameters
    ----------
    offsets, weights : np.ndarray
        The stencil, as returned by Solver.get_lhs_stencil or get_rhs_stencil (the weights are used as is: divide by h^2
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


def _is_cell(center: CenterType) -> bool:
    if center not in (CenterType.NODE_CENTER, CenterType.CELL_CENTER):
        raise ValueError(f"unknown center type {center}")
    return center == CenterType.CELL_CENTER


def fixed_mask(bc: List[List[BoundaryType]], n: int, center: CenterType = CenterType.NODE_CENTER) -> np.ndarray:
    """Nodes that are not solved for: for node-centred data the ODD wall nodes and the duplicated last
    node of PER directions, none for cell-centred data"""
    _check_bc(bc)
    fixed = np.zeros((n, n, n), dtype=bool)
    if _is_cell(center):
        return fixed
    for d in range(3):
        if bc[d][0] == BoundaryType.PER:
            fixed[_plane(d, n - 1)] = True
        for side in range(2):
            if bc[d][side] == BoundaryType.ODD:
                fixed[_plane(d, 0 if side == 0 else n - 1)] = True
    return fixed


def ghost_remainder(bc: List[List[BoundaryType]], n: int, r: int, ghost: Callable[[int, int, int], np.ndarray],
                    center: CenterType = CenterType.NODE_CENTER) -> np.ndarray:
    """
    Build D = u - R(u) on the extended grid from the boundary data

    D is built one direction at a time (x, then y, then z), each pass reading the values written by
    the previous ones: D_{-m} = -D_{+m} + G_m on ODD walls and D_{-m} = D_{+m} + N_m on EVEN walls,
    with +m the mirror image of the ghost -m (the node m for node-centred data, m - 1 for cell-centred
    data). This gives the edge and corner ghosts as double/triple reflections.

    Parameters
    ----------
    bc : list of lists
        [[xmin, xmax], [ymin, ymax], [zmin, zmax]] with ODD (Dirichlet), EVEN (Neumann) or PER
    n : int
        Number of nodes per direction (N + 1 node-centred, N cell-centred)
    r : int
        Number of ghost layers (at least the reach of the left-hand side stencil)
    ghost : callable
        ghost(d, side, m) returns the data on the wall (side 0: min, 1: max) normal to d, as an
        array broadcastable to the extended face, of shape (n + 2r, n + 2r) over the two other
        directions in increasing order:
          - m = 0, ODD walls of node-centred data only: g, the wall value;
          - m >= 1, ODD walls: G_m = u(-s_m) + u(s_m); EVEN walls: N_m = u(-s_m) - u(s_m), with n
            the inward distance to the wall and s_m = m h (node) or (m - 1/2) h (cell), see
            taylor_ghost.
        The face covers the ghost positions of the neighbouring walls. Only stencils with diagonal
        terms (MEHR_*) read them, and then use the data slightly beyond the face; the cross-shaped
        LGF_* stencils never do.
    center : CenterType, optional
        NODE_CENTER (default) or CELL_CENTER

    Returns
    -------
    np.ndarray
        D on the extended grid, of shape (n + 2r,)*3, zero on the unknown nodes
    """
    _check_bc(bc)
    cell = _is_cell(center)
    period = n if cell else n - 1
    first, last = r, r + n - 1  # extended indices of the first and last nodes
    face = (n + 2 * r, n + 2 * r)
    D = np.zeros((n + 2 * r,) * 3)
    for d in range(3):
        if bc[d][0] == BoundaryType.PER:
            for m in range(1, r + 1):
                D[_plane(d, first - m)] = D[_plane(d, first - m + period)]
                D[_plane(d, last + m)] = D[_plane(d, last + m - period)]
            continue
        for side in range(2):
            odd = bc[d][side] == BoundaryType.ODD
            iw, step = (first, -1) if side == 0 else (last, 1)
            if odd and not cell:
                D[_plane(d, iw)] = np.broadcast_to(ghost(d, side, 0), face)
            for m in range(1, r + 1):
                refl = -1.0 if odd else 1.0
                mirror = iw - step * (m - 1 if cell else m)
                D[_plane(d, iw + step * m)] = refl * D[_plane(d, mirror)] + np.broadcast_to(ghost(d, side, m), face)
    return D


def taylor_ghost(coefs: Dict[int, np.ndarray], m: int, h: float, wall: BoundaryType,
                 center: CenterType = CenterType.NODE_CENTER) -> np.ndarray:
    """
    Ghost data from the Taylor coefficients of u on the wall

    Parameters
    ----------
    coefs : dict
        {k: a_k}, a_k = d_n^k u on the wall (n inward): the even k for an ODD wall (a_0 = g), the odd
        k for an EVEN wall (a_1 = q). Order p needs the terms k < p.
    m : int
        Ghost layer (m = 0 returns g for an ODD wall of node-centred data)
    h : float
        Grid spacing
    wall : BoundaryType
        ODD (Dirichlet) or EVEN (Neumann)
    center : CenterType, optional
        NODE_CENTER (default): the layer m lies at s = m h; CELL_CENTER: at s = (m - 1/2) h

    Returns
    -------
    np.ndarray
        G_m = 2 sum_{k even} s^k/k! a_k for ODD, N_m = -2 sum_{k odd} s^k/k! a_k for EVEN
    """
    cell = _is_cell(center)
    if m == 0:
        if cell or wall != BoundaryType.ODD:
            raise ValueError("the layer m = 0 only exists on the ODD walls of node-centred data")
        return coefs[0]
    s = (m - 0.5 * cell) * h
    parity, sign = (0, 2.0) if wall == BoundaryType.ODD else (1, -2.0)
    out = 0.0
    for k, a in coefs.items():
        if k % 2 == parity:
            out = out + sign * s**k / math.factorial(k) * a
    return out


def corrected_rhs(offsets: np.ndarray, weights: np.ndarray, rhs: np.ndarray, remainder: np.ndarray, h: float,
                  bc: List[List[BoundaryType]], center: CenterType = CenterType.NODE_CENTER) -> np.ndarray:
    """
    Right-hand side to give to FLUPS: rhs - L_h D on the unknown nodes, zero on the fixed ones

    Parameters
    ----------
    offsets, weights : np.ndarray
        The left-hand side stencil (Solver.get_lhs_stencil())
    rhs : np.ndarray
        M_h f on the n^3 nodes (f itself for the LGF kernels)
    remainder : np.ndarray
        D from ghost_remainder
    h : float
        Grid spacing
    bc : list of lists
        The boundary conditions given to the solver
    center : CenterType, optional
        NODE_CENTER (default) or CELL_CENTER
    """
    n = rhs.shape[0]
    r = (remainder.shape[0] - n) // 2
    out = rhs - apply_stencil(offsets, weights, remainder, r) / h**2
    out[fixed_mask(bc, n, center)] = 0.0
    return out


def recover(v: np.ndarray, remainder: np.ndarray, bc: List[List[BoundaryType]],
            center: CenterType = CenterType.NODE_CENTER) -> np.ndarray:
    """
    Solution u = R(u) + D from the FLUPS output v = R(u), on the n^3 nodes

    For node-centred data the ODD wall nodes are set to the wall data and the duplicated PER node to
    the first one. For cell-centred data D vanishes on every node and u = v.
    """
    n = v.shape[0]
    r = (remainder.shape[0] - n) // 2
    cell = _is_cell(center)
    u = np.where(fixed_mask(bc, n, center), 0.0, v)
    for d in range(3):
        if bc[d][0] == BoundaryType.PER and not cell:
            u[_plane(d, n - 1)] = u[_plane(d, 0)]
    return u + remainder[r:r + n, r:r + n, r:r + n]
