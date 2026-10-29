Considering symetric boundary conditions, here are all the possible cases: 

![test-2](/uploads/3cdcfee9c2e981ba77e1585a6fa2f622/test-2.png)

This will impact the parameter in the FFTW_Plan_Dim class. 

## Cell centered RHS 
|parameter \ BCs       | ODD - ODD              | ODD - EVEN        | EVEN - ODD        |  EVEN - EVEN          |
| ------               | :------:               | :------:          | :------:          | :------:              |
| Transform type       | RODFT10 & RODFT01      | RODFT11           | REDFT11           | REDFT10 & REDFT01     |
| Normalisation factor | $\frac{1}{2N}$         | $\frac{1}{2N}$    | $\frac{1}{2N}$    | $\frac{1}{2N}$        |
| Offset of the modes  | 0.0                    | 0.5               | 0.5               | 0.0                   |
| # of points out      | N + 1                  | N                 | N                 | N + 1                 |
| Output correction    | The 0 mode = 0         | /                 | /                 | The fli-flop mode = 0 |

## Vertex centered 
|parameter \ BCs       | ODD - ODD                   | ODD - EVEN        | EVEN - ODD        | EVEN - EVEN           |
| ------               | :------:                    | :------:          | :------:          | :------:              |
|Transform type        | RODFT00                     | RODFT01 & RODFT10 | REDFT01 & REDFT10 | REDFT00               |
| Normalisation factor | $\frac{1}{2(N+1)}$          | $\frac{1}{2N}$    | $\frac{1}{2N}$    | $\frac{1}{2(N-1)}$    |
| Offset of the modes  | 0.0                         | 0.0               | 0.0               | 0.0                   |
| # of points out      | N - 2                       | N - 1             | N - 1             | N                     |
| Input correction     | First and last data ignored | First data ignored| Last data ignored | /                     |
The ignored nodes of an ODD side are the wall nodes; FLUPS returns 0 there.

## Non-homogeneous Dirichlet and Neumann conditions

FLUPS only solves homogeneous symmetric problems: ODD imposes $u = 0$ on the wall and EVEN imposes $\partial_n u = 0$.
A wall with data, $u = g$ (Dirichlet) or $\partial_n u = q$ (Neumann, $n$ the inward distance to the wall), is imposed
by the user through a correction of the right-hand side. No change to the solve is needed.

### Why it works

The lattice kernels (`LGF_*`, `MEHR_*`) are the exact inverse of a finite-difference stencil on the symmetric
(mirrored) extension of the domain: the returned field solves $L_h u = M_h f$, with $M_h$ the identity for `LGF_*`.
The stencils are returned by `flups_get_stencil` (`Solver.get_stencil` in pyflups):

| kernel              | $h^2 L_h$ ($\delta^2_{d,m} u = u_{+m} - 2u + u_{-m}$ along $d$, $\delta^2_d = \delta^2_{d,1}$) | $M_h$ | reach |
| ------              | ------                                                                                    | ------ | :---: |
| `LGF_p`             | $\sum_d \sum_m c_m \delta^2_{d,m}$, $c = (1)$, $(4/3, -1/12)$, $(3/2, -3/20, 1/90)$, $(8/5, -1/5, 8/315, -1/560)$ | $1$ | $p/2$ |
| `MEHR_4L`/`MEHR_4F` | $\sum_d \delta^2_d + \frac{1}{6} \sum_{i<j} \delta^2_i \delta^2_j$                          | $1 + \frac{1}{12} \sum_d \delta^2_d$ | 1 |
| `MEHR_6L`/`MEHR_6F` | `MEHR_4L` $+ \frac{1}{30} \delta^2_0 \delta^2_1 \delta^2_2$                                    | $1 + \frac{1}{12} \sum_d \delta^2_d + \frac{1}{90} \sum_{i<j} \delta^2_i \delta^2_j - \frac{1}{240} \sum_d \delta^4_d$ | 1 |

`CHAT_2` and `HEJ_*` are not the inverse of a stencil (`flups_get_stencil` returns 0 points): the correction below does
not apply to them.

Near a wall the stencil reads values on the wall and beyond it (the ghosts). FLUPS fills them with the mirror image of
the field, $R(u)$, which corresponds to $g = 0$ or $q = 0$. Write the true field, extended to the ghosts, as
$u = R(u) + D$, where $D$ is zero on every unknown node. Then

$$ L_h R(u) = M_h f - L_h D \quad \text{on the unknown nodes}, \qquad u = R(u) + D. $$

FLUPS solves the first equation with the corrected right-hand side $M_h f - L_h D$ (ODD for Dirichlet, EVEN for
Neumann walls), and the solution is recovered by adding $D$ back, i.e. by setting the Dirichlet wall nodes to $g$.

### Building D from the data (node-centred data)

On a wall, with $m \geq 1$ the ghost layer at distance $mh$ outside:

| wall           | on the wall | ghost layer $m$ |
| ------         | ------      | ------          |
| ODD, $u = g$   | $D_0 = g$   | $D_{-m} = u(-mh) + u(mh) = 2 \sum_{k\ \mathrm{even}} \frac{(mh)^k}{k!} a_k$ |
| EVEN, $\partial_n u = q$ | (unknown) | $D_{-m} = u(-mh) - u(mh) = -2 \sum_{k\ \mathrm{odd}} \frac{(mh)^k}{k!} a_k$ |

The wall Taylor coefficients $a_k = \partial_n^k u$ follow from the data and the PDE ($\Delta_t$ is the Laplacian
along the wall):

$$ a_0 = g \ \text{(Dirichlet)}, \qquad a_1 = q \ \text{(Neumann)}, \qquad a_{k+2} = \partial_n^k f - \Delta_t a_k. $$

They only involve $g$ or $q$ and their tangential derivatives, and $f$ and its derivatives on the wall.

- **Accuracy.** Order $p$ needs the terms $k < p$: a Dirichlet ghost $O(h^p)$, a Neumann ghost $O(h^{p+1})$. With one
  term fewer the order drops by two. For Dirichlet walls, `LGF_2` and `MEHR_*L` only read the wall node, so
  $D_0 = g$ is enough and the correction is exact.
- **Edges and corners.** Build $D$ one direction at a time, each pass reading the values written by the previous one:
  $D_{-m} = -D_{+m} + (\text{ghost data})$ on ODD walls and $D_{-m} = D_{+m} + (\text{ghost data})$ on EVEN walls. Near an
  edge, the diagonal terms of the Mehrstellen stencils read ghosts outside two walls at once, which use the boundary
  data slightly beyond the face; the cross-shaped `LGF_*` stencils never read them.
- **Mehrstellen.** Use `MEHR_*L` and apply $M_h$ to $f$ yourself: $M_h f$ reads $f$ one (order 4) or two (order 6)
  layers beyond the wall. With `MEHR_*F`, FLUPS applies $M_h$ to the whole right-hand side, including $-L_h D$, and
  undoing it would need the non-local $M_h^{-1}$.
- **No Dirichlet wall.** With only EVEN and PER directions, $u$ is known up to a constant and FLUPS drops the $k=0$
  mode, i.e. it projects out the (trapezoidal) mean of the right-hand side. The corrected right-hand side satisfies
  the compatibility condition up to the discretisation error.
- **Cell-centred data.** The same principle holds, but the wall lies half-way between nodes: the ghosts sit at
  $(m - 1/2) h$ and a Dirichlet wall always needs Taylor ghosts. This case has not been validated.

### Python helpers and validation

`pyflups.bc` implements these steps for node-centred data (`ghost_remainder`, `taylor_ghost`, `corrected_rhs`,
`recover`, `apply_stencil`), and `python_api/examples/nonhomogeneous_bc.py` validates them on a manufactured solution
with values and slopes that vary along every wall. With ghost data taken from the exact solution the error is at
round-off level ($\leq 10^{-13}$) for every lattice kernel and for Dirichlet, Neumann and mixed walls, which also
checks `flups_get_stencil` against the kernels. Solving from the samples only ($f$ on the nodes, $g$ or $q$ on the wall
nodes), with every $a_k$ computed by finite differences of accuracy $O(h^{P-k})$ ($P = p$ for Dirichlet, $p+1$ for
Neumann walls) and the values `MEHR_*L` needs beyond the walls and face edges extrapolated, the observed order is the
design order of the kernel (2, 4, 6, 8 for `LGF_*`, 4 and 6 for `MEHR_*L`).
