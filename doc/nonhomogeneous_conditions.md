# Non-homogeneous boundary conditions: how to modify the right-hand side

FLUPS only solves homogeneous symmetric problems: an ODD wall imposes $u = 0$ and an EVEN wall imposes
$\partial_n u = 0$. This page explains how to impose a wall with data, $u = g$ (Dirichlet) or $\partial_n u = q$
(Neumann, $n$ the inward distance to the wall), with the same solver: the boundary data is moved into the right-hand
side, and the solve itself is unchanged. The method is exact for the lattice kernels (`LGF_*`, `MEHR_*`). It does not
apply to `CHAT_2` and `HEJ_*`.

FLUPS solves $\nabla^2 u = f$. Section 1 explains the method, section 2 how to get the stencil coefficients from FLUPS,
section 3 gives the recipe and section 4 the Python helpers.

## 1. The method

### 1.1 A one-dimensional example

Take $u'' = f$ on $[0, L]$ with node-centred data, $x_i = i h$, and the `LGF_2` kernel. `LGF_2` is the exact inverse
of the 3-point stencil: the field FLUPS returns satisfies, at every unknown node $i$,

$$ \frac{u_{i-1} - 2 u_i + u_{i+1}}{h^2} = f_i . $$

At the first node next to a wall, this stencil reads one value on the other side: $u_{i-1}$ for $i = 0$, or the wall
node $u_0$ for $i = 1$. **FLUPS does not store ghost points**: these values only exist inside the stencil, and since
FLUPS solves on the symmetric extension of the domain, they are always the mirror image of the field. With ODD the wall
node is $u_0 = 0$ and $u_{-m} = -u_m$. With EVEN, $u_{-m} = u_m$.

We call **ghost layer** $m$ the plane of points, parallel to the wall, that lies $m$ planes outside it: $x = -m h$ for
node-centred data, $x = -(m - 1/2) h$ for cell-centred data. A stencil of reach $r$ reads the layers $1$ to $r$. They
are not part of the domain and FLUPS has no memory for them: they only exist in the user's own arrays, extended by $r$
layers beyond each wall, where the user builds the correction described below.

- **Dirichlet, $u(0) = g$ (ODD wall).** The true equation at node 1 has $u_0 = g$, while FLUPS uses $u_0 = 0$. Move the
  difference to the right-hand side: give FLUPS $f_1 - g/h^2$ instead of $f_1$, keep $f_i$ for the other nodes, and set
  $u_0 = g$ after the solve. The discrete problem is the same, so the solution is exactly the one of the stencil.
- **Neumann, $u'(0) = q$ (EVEN wall).** The wall node is an unknown. Its equation reads the ghost $u_{-1}$, which is
  $u_1$ for FLUPS and $u_1 - 2hq$ to second order for the true field. Give FLUPS $f_0 + 2q/h$ instead of $f_0$.

With cell-centred data, $x_i = (i + 1/2) h$, there is no wall node, and the classical ghosts are $u_{-1} = 2g - u_0$
(Dirichlet) and $u_{-1} = u_0 - hq$ (Neumann), against $-u_0$ and $u_0$ for FLUPS. The corrections are $f_0 - 2g/h^2$
and $f_0 + q/h$.

### 1.2 General form

The lattice kernels are the exact inverse of a stencil on the symmetric extension of the domain: the field FLUPS
returns solves $L_h u = M_h f$, where $L_h$ is a discrete Laplacian and $M_h$ is the identity for `LGF_*`:

| kernel              | $h^2 L_h$ ($\delta^2_{d,m} u = u_{+m} - 2u + u_{-m}$ along $d$, $\delta^2_d = \delta^2_{d,1}$) | $M_h$ | reach |
| ------              | ------                                                                                    | ------ | :---: |
| `LGF_p`             | $\sum_d \sum_m c_m \delta^2_{d,m}$, $c = (1)$, $(4/3, -1/12)$, $(3/2, -3/20, 1/90)$, $(8/5, -1/5, 8/315, -1/560)$ | $1$ | $p/2$ |
| `MEHR_4L`/`MEHR_4F` | $\sum_d \delta^2_d + \frac{1}{6} \sum_{i<j} \delta^2_i \delta^2_j$                          | $1 + \frac{1}{12} \sum_d \delta^2_d$ | 1 |
| `MEHR_6L`/`MEHR_6F` | `MEHR_4L` $+ \frac{1}{30} \delta^2_0 \delta^2_1 \delta^2_2$                                    | $1 + \frac{1}{12} \sum_d \delta^2_d + \frac{1}{90} \sum_{i<j} \delta^2_i \delta^2_j - \frac{1}{240} \sum_d \delta^4_d$ | 1 |

The reach $r$ is the number of points the stencil reads on each side, so the nodes within $r$ of a wall read up to $r$
ghost layers. Let $R(u)$ be the mirror extension FLUPS works with, and $u$ the true field extended to the ghosts. Their
difference $D = u - R(u)$ vanishes on every unknown node. It is non-zero only on the ghost layers and, for node-centred
data, on the ODD wall nodes ($D_0 = g$). Since $L_h u = M_h f$ on the unknown nodes,

$$ L_h R(u) = M_h f - L_h D \quad \text{on the unknown nodes}, \qquad u = R(u) + D . $$

FLUPS solves the first equation exactly (ODD for Dirichlet walls, EVEN for Neumann walls), and adding $D$ back amounts
to setting the Dirichlet wall nodes to $g$. The correction $L_h D$ is non-zero only on the $r$ layers next to a wall.
Section 1.1 is the case $r = 1$, with $D_0 = g$ (Dirichlet) or $D_{-1} = -2hq$ (Neumann).

`CHAT_2` and `HEJ_*` are not the inverse of a stencil, so this argument does not apply to them.

### 1.3 Ghost values from the boundary data

The ghost layer $m \geq 1$ lies at the distance $s_m$ outside the wall and mirrors the node at $s_m$ inside:

- node-centred data: $N+1$ nodes at $x_i = i h$, $s_m = m h$, the ghost $-m$ mirrors the node $m$. The wall nodes are
  data on ODD walls (FLUPS ignores them on input and returns 0) and unknowns on EVEN walls;
- cell-centred data: $N$ nodes at $x_i = (i + 1/2) h$, $s_m = (m - 1/2) h$, the ghost $-m$ mirrors the node $m - 1$.
  There is no wall node and every node is an unknown.

$D$ follows from a Taylor expansion of $u$ about the wall, with the coefficients $a_k = \partial_n^k u$:

| wall           | on the wall | ghost layer $m$ |
| ------         | ------      | ------          |
| ODD, $u = g$   | $D_0 = g$ (node-centred only) | $D_{-m} = u(-s_m) + u(s_m) = 2 \sum_{k\ \mathrm{even}} \frac{s_m^k}{k!} a_k$ |
| EVEN, $\partial_n u = q$ | (unknown) | $D_{-m} = u(-s_m) - u(s_m) = -2 \sum_{k\ \mathrm{odd}} \frac{s_m^k}{k!} a_k$ |

The coefficients follow from the data and the PDE ($\Delta_t$ is the Laplacian along the wall):

$$ a_0 = g \ \text{(Dirichlet)}, \qquad a_1 = q \ \text{(Neumann)}, \qquad a_{k+2} = \partial_n^k f - \Delta_t a_k . $$

They only involve $g$ or $q$ and their tangential derivatives, and $f$ and its normal derivatives on the wall: $D$ is
known before the solve.

- **Accuracy.** Order $p$ needs the terms $k < p$: a Dirichlet ghost $O(h^p)$, a Neumann ghost $O(h^{p+1})$. With one
  term fewer, the order drops by two. With node-centred data and Dirichlet walls, `LGF_2` and `MEHR_*L` only read the
  wall node, so $D_0 = g$ is enough and the correction is exact. With cell-centred data, a Dirichlet wall always needs
  Taylor ghosts: $a_0$ for order 2, up to $a_2$ for order 4 and $a_4$ for order 6.
- **Edges and corners.** Build $D$ one direction at a time, each pass reading the values written by the previous one:
  $D_{-m} = -D_{+m} + (\text{ghost data})$ on ODD walls and $D_{-m} = D_{+m} + (\text{ghost data})$ on EVEN walls. Near an
  edge, the diagonal terms of the Mehrstellen stencils read ghosts outside two walls at once, which use the boundary
  data slightly beyond the face. The cross-shaped `LGF_*` stencils never read them. The rule is the same for both
  centrings, the mirror of the ghost $-m$ being the node $m$ (node-centred) or $m-1$ (cell-centred).

## 2. Getting the stencil coefficients from FLUPS

The weights of $L_h$ and $M_h$ are returned by `flups_get_stencil`:

```c
int flups_get_stencil(FLUPS_Solver* s, const FLUPS_StencilType type, int* offsets, double* weights);
```

- `type` is `STENCIL_LHS` for $L_h$ or `STENCIL_RHS` for $M_h$;
- the stencil is $L_h u(x) = \frac{1}{h^2} \sum_i w_i \, u(x + o_i h)$ and $M_h f(x) = \sum_i w_i \, f(x + o_i h)$, with
  the integer offsets $o_i$ in `offsets[3*i + d]` and the weights $w_i$ in `weights[i]`;
- the return value is the number of points $n$ (at most 27). Call the function with `NULL` arrays to get $n$, then
  with arrays of $3n$ and $n$ entries;
- it returns 0 for `CHAT_2` and `HEJ_*`, which are not the inverse of a stencil;
- it can be called once the Green's function type is set, before or after `flups_setup`.

In Python, `Solver.get_stencil(pf.StencilType.LHS)` (or `RHS`) returns the arrays `(offsets, weights)`.

## 3. Imposing the conditions with FLUPS, step by step

1. **Boundary conditions.** Declare the Dirichlet walls ODD and the Neumann walls EVEN. The other directions keep
   their usual condition.
2. **Kernel.** Use `LGF_*` or `MEHR_*L`. With `MEHR_*F`, FLUPS applies $M_h$ to the whole right-hand side, including
   the correction $-L_h D$, and undoing it would need the non-local $M_h^{-1}$.
3. **Stencils.** Get $L_h$ and, for `MEHR_*L`, $M_h$ with `flups_get_stencil` (section 2). Their reach $r$ is the
   number of ghost layers to build.
4. **Ghost values.** On each wall, compute the coefficients $a_k$ needed for the target order and build $D$ on the $r$
   ghost layers, direction by direction (section 1.3). For node-centred data, also set $D_0 = g$ on the ODD walls.
5. **Right-hand side.** On every unknown node, compute $M_h f - L_h D$. For `LGF_*`, $M_h f = f$. For `MEHR_*L`,
   $M_h f$ reads $f$ one (order 4) or two (order 6) layers beyond the wall.
6. **Solve** with the modified right-hand side.
7. **Recover.** For node-centred data, set the ODD wall nodes to $g$ (FLUPS returns 0 there). Cell-centred data has
   no wall node, and the solution is used as returned.

**No Dirichlet wall.** With only EVEN and PER directions, $u$ is known up to a constant and FLUPS drops the $k=0$ mode,
i.e. it removes the mean of the right-hand side (trapezoidal for node-centred data, plain for cell-centred data). The
modified right-hand side satisfies the compatibility condition up to the discretisation error.

## 4. Python helpers and validation

`pyflups.bc` implements these steps for node- and cell-centred data (`ghost_remainder`, `taylor_ghost`,
`corrected_rhs`, `recover`, `apply_stencil`, all with a `center` argument; usage in `python_api/README.md`), and
`python_api/examples/nonhomogeneous_bc.py` validates them on a manufactured solution with values and slopes that vary
along every wall (`--suite`, with `--center=node` or `--center=cell`). For both centrings, with ghost data taken from
the exact solution the error is at round-off level ($\leq 10^{-13}$) for every lattice kernel and for Dirichlet,
Neumann and mixed walls, which also checks `flups_get_stencil` against the kernels. Solving from the samples only
($f$ on the nodes, $g$ or $q$ on the wall nodes), with every $a_k$ computed by finite differences of accuracy
$O(h^{P-k})$ ($P = p$ for Dirichlet, $p+1$ for Neumann walls) and the values `MEHR_*L` needs beyond the walls and face
edges extrapolated, the observed order is the design order of the kernel (2, 4, 6, 8 for `LGF_*`, 4 and 6 for
`MEHR_*L`).
