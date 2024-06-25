# python module: viscmodel


class ViscModel():

    _needed_args = { }

    def __init__(self, viscosity: float):
        self._viscosity = viscosity

    def compute(self, **kwargs):
        assert set(kwargs.keys()) == self._needed_args
        return self._viscosity


class Cst(ViscModel):
    pass

    
class PowerLaw(ViscModel):

    _needed_args = { 'T' }

    def __init__(self, visc_ref: float, Tref: float, power: float):
        self._visc_ref = visc_ref
        self._Tref = Tref
        self._power = power

    def compute(self, **kwargs):
        assert set(kwargs.keys()) == self._needed_args
        return self._visc_ref * (kwargs['T']/self._Tref)**self._power
