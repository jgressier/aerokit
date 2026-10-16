# Chebyshev spectral discretisation

`aerokit.common.numspectral` supplies dense collocation differentiation
matrices for smooth one-dimensional fields.  It is intended for the stability
operators in `aerokit.stability` and for small spectral boundary-value or
eigenvalue problems; it is not a finite-volume discretisation.

## One domain

`ChebCollocation(npts, xmin, xmax)` uses `npts` Chebyshev--Lobatto nodes.  Its
default affine mapping returns the physical points in increasing order:

```python
import numpy as np
from aerokit.common.numspectral import ChebCollocation

grid = ChebCollocation(41, xmin=-1.0, xmax=1.0)
f = np.exp(grid.x)
df = grid.matder(1) @ f
assert np.allclose(df, f)
```

`matder(m)` is the physical-coordinate matrix for the `m`-th derivative.  The
reference matrices are computed once and cached.  When a non-affine mapping is
given, the mapping applies the chain rule before the matrix is returned.

`extrapol(fk, x)` evaluates the polynomial interpolant through values `fk` at
the collocation nodes.  It also permits extrapolation, but that use is usually
ill-conditioned.  `fit_to_gauss(x, f)` is a least-squares Chebyshev fit to
arbitrary data; despite its historical name, it does not perform Gaussian
quadrature.

## Two domains and their interface

`MultiDomainChebCollocation(npts, xmin, interface, xmax)` creates a left grid
on `[xmin, interface]` and a right grid on `[interface, xmax]`.  The interface
node is stored only once:

```text
global index:  0 ... npts - 2  npts - 1  npts ... 2*npts - 2
                  left interior  interface       right interior
```

Thus `grid.npts == 2*npts - 1` and `grid.interface_index == npts - 1`.  A
single unknown at that index makes the *field value* continuous.

### One-sided rows versus the assembled row

Both subdomains retain a legitimate one-sided derivative at the shared point.
For derivative order `m`, the internal `_embedded_interface_row` operation
places each local row in global-column ordering:

```text
D^-_Gamma = [ D_left[-1, :] | 0 ... 0 ]
D^+_Gamma = [ 0 ... 0 | D_right[0, :] ]
```

The common interface value belongs to both local rows but appears once in the
global vector.  `matder(m)` installs their arithmetic mean at the interface,
`(D^-_Gamma + D^+_Gamma) / 2`.  This gives a useful derivative for a smooth
field.  It does **not** enforce derivative continuity and should not be used
as an implicit transmission condition.

Use `interface_derivative_jumps(maxorder)` when assembling an equation that
requires smooth matching.  Its row `jumps[m - 1]` represents

```text
(d^m f / dx^m)_left - (d^m f / dx^m)_right.
```

Replacing suitable system rows by these vectors imposes zero derivative jump.
For example, the two-domain Orr--Sommerfeld assembly uses the first three
jump rows.  This construction is suitable for a smooth solution across the
interface.  A material or flux discontinuity requires a problem-specific
interface flux or jump law instead.
