# Ideal blade rows

`aerokit.aero.blade.BladeRow` is a deliberately small, two-dimensional model
for preliminary compressor or turbine velocity triangles.  It uses `(x, y)`
components; `y` can also be read as the local tangential `theta` component.

```python
from aerokit.aero.blade import BladeRow
from aerokit.aero.model2D import State2d

inlet = State2d(rho=1.2, u=120., v=20., p=1.e5)
rotor = BladeRow(omega=800., radius=0.25)  # Uy = 200 m/s

# The specified 25 degrees is the outlet angle in the relative frame.
outlet = rotor.downstream(inlet, beta=25., conservation="massflow")
```

## Reference frames

The absolute velocity is `C` and the relative velocity is `W`.  With blade
velocity `U = (0, Uy)`, the convention is

\[
 C = W + U.
\]

`row.relative(absolute_state)` and `row.absolute(relative_state)` make this
conversion explicit.  `downstream` takes and returns absolute states by
default; use `frame="relative"` to work entirely in the blade frame.

## Flow closure

The model turns the relative flow to the requested `beta` while conserving
relative total pressure and relative total temperature.  Two choices are
available:

- `conservation="vx"` (default) keeps `Vx` constant.  It is a useful first
  velocity-triangle approximation.
- `conservation="massflow"` conserves `rho * Vx * area`.  `area_ratio` is the
  outlet-to-inlet area ratio.  The outlet is solved isentropically on the
  inlet Mach branch, unless `branch="subsonic"` or `branch="supersonic"` is
  selected.

The mass-flow closure detects a choked requested outlet.  It is a single
streamtube model: radial equilibrium, losses, and the planned stack of radial
layers are intentionally outside this first implementation.

For a rotor, `specific_work(inlet, outlet)` returns Euler work
`Uy * (Ctheta_out - Ctheta_in)`.

## Compressor-map example

[`examples/engine/compressor_map.py`](../../examples/engine/compressor_map.py)
solves a family of ideal operating points by varying blade speed `U` and the
imposed ratio `Pt_in / p_out`.  Each point is adjusted until the predicted
outlet static pressure equals `p_out`; the result is shown as conventional
speed lines in `Vx`--`Pt_in/p_out` coordinates and in the equivalent
`psi`--`phi` diagram:

```console
python examples/engine/compressor_map.py --save compressor-map.png
```
