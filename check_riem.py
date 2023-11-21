import aerokit.instance.riemann as riem
import aerokit.aero.unsteady1D as uq
import numpy as np

qL = uq.unsteady_state(rho=2., u=0., p=10.)
qR = uq.unsteady_state(rho=1., u=0., p=1.)
pb = riem.riemann_pb(qL, qR)