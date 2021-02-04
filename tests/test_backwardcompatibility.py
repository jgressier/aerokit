import aerokit.aero.Isentropic as Is
import aerokit.aero.riemann     as oldRiem
import aerokit.instance.riemann as Riem
import aerokit.aero.nozzle      as oldNoz
import aerokit.instance.nozzle  as Noz
import aerokit.aero.unsteady1D  as uq
import numpy as np
import pytest

def stagnation_i2t():
    assert Is.TtTs_Mach(2.) == Is.TiTs_Mach(2.)
    assert Is.Ptps_Mach(2.) == Is.PiPs_Mach(2.)
    assert Mach_PiPs(3.)    == Mach_PtPs(3.)
    assert Mach_TiTs(3.)    == Mach_TtTs(3.)
    assert Velocity_MachTi(.8, 300.) == Velocity_MachTt(.8, 300.)

def riemann():
    pbriem = oldRiem.riemann_pb(uq.unsteady_state(1., -10., 100.),
                                uq.unsteady_state(1., -20., 1000.) )
    assert True

def nozzle():
    x = np.linspace(0., 6., 100)
    noz = oldNoz.nozzle(x, 2.-np.sin(x))
    assert True # if no error