import aerokit.aero.Isentropic as Is
import aerokit.aero.ShockWave as sw
import aerokit.aero.riemann     as oldRiem
import aerokit.instance.riemann as Riem
import aerokit.aero.nozzle      as oldNoz
import aerokit.instance.nozzle  as Noz
import aerokit.aero.unsteady1D  as uq
import numpy as np
import pytest

def test_stagnation_i2t():
    assert Is.TtTs_Mach(2.) == Is.TiTs_Mach(2.)
    assert Is.PtPs_Mach(2.) == Is.PiPs_Mach(2.)
    assert Is.Mach_PiPs(3.) == Is.Mach_PtPs(3.)
    assert Is.Mach_TiTs(3.) == Is.Mach_TtTs(3.)
    assert Is.Velocity_MachTi(.8, 300.) == Is.Velocity_MachTt(.8, 300.)

def test_Ptshock()
    assert sw.Pt_ratio(2.) == sw.Pi_ratio(2.)
    
def test_riemann():
    pbriem = oldRiem.riemann_pb(uq.unsteady_state(1., -10., 100.),
                                uq.unsteady_state(1., -20., 1000.) )
    assert True

def test_nozzle():
    x = np.linspace(0., 6., 100)
    noz = oldNoz.nozzle(x, 2.-np.sin(x))
    assert True # if no error