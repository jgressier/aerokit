import importlib
import warnings

import aerokit.aero.Isentropic as Is
import aerokit.aero.ShockWave as sw
import aerokit.aero.riemann as oldRiem
import aerokit.aero.nozzle as oldNoz
import aerokit.aero.unsteady1D as uq
import numpy as np


def test_deprecated_aero_instance_modules_warn():
    for module in (oldRiem, oldNoz):
        with warnings.catch_warnings(record=True) as recorded:
            warnings.simplefilter("always", DeprecationWarning)
            importlib.reload(module)
        assert any(item.category is DeprecationWarning for item in recorded)


def test_stagnation_i2t():
    assert Is.TtTs_Mach(2.0) == Is.TiTs_Mach(2.0)
    assert Is.PtPs_Mach(2.0) == Is.PiPs_Mach(2.0)
    assert Is.Mach_PiPs(3.0) == Is.Mach_PtPs(3.0)
    assert Is.Mach_TiTs(3.0) == Is.Mach_TtTs(3.0)
    assert Is.Velocity_MachTi(0.8, 300.0) == Is.Velocity_MachTt(0.8, 300.0)


def test_Ptshock():
    assert sw.Pt_ratio(2.0) == sw.Pi_ratio(2.0)


def test_riemann():
    pbriem = oldRiem.riemann_pb(
        uq.unsteady_state(1.0, -10.0, 100.0), uq.unsteady_state(1.0, -20.0, 1000.0)
    )
    assert pbriem.pstar() > 0.0


def test_nozzle():
    x = np.linspace(0.0, 6.0, 100)
    oldNoz.nozzle(x, 2.0 - np.sin(x))
    assert True  # if no error
