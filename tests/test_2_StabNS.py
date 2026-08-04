import numpy as np

from aerokit.common.mapping import SemiInfiniteAlgebraicMapping
from aerokit.stability.NS import NSaxi


def test_nsaxi_semi_infinite_mapping():
    operator = NSaxi(20, basestate=None, mapping=SemiInfiniteAlgebraicMapping(10.0))
    assert operator.r[0] == 0.0
    assert np.isinf(operator.r[-1])
    assert np.all(np.isfinite(operator._radial_matder(1)))
