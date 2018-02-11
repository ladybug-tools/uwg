import pytest
from UWG import psychrometrics
from math import log10

import pprint
import decimal

dd = lambda x: decimal.Decimal.from_float(x)
pp = lambda x: pprint.pprint(x)



def test_psychrometric():

    # Input values
    Tdb_in = 297.5311337413935
    w_in = 0.018576773131376
    P = 10090

    Tdb, w, phi, h, Tdp, v = psychrometrics.psychrometrics (Tdb_in, w_in, P)

    assert Tdb == pytest.approx(24.381133741393512, abs=1e-15)  # dry bulb temperature [C]
    assert w == pytest.approx(0.018576773131376, abs=1e-15)     # Humidity ratio [kgv/kgd]
    assert phi == pytest.approx(9.581555922752541, abs=1e-15)   # Relative Humidity Pw/Pws*100
    assert h == pytest.approx(7.183036653518451e+04, abs=1e-15) # enthalpy [J/kgd]
    assert Tdp == pytest.approx(-10.012150181172135, abs=1e-15) # Wet bulb temp [C]
    assert v == pytest.approx(8.717031296493113, abs=1e-15)     # specific volume [m3/kga]s

    """
    print dd(Tdb), '\n24.381133741393512'
    print dd(w), '\n0.018576773131376'
    print dd(phi), '\n9.581555922752541'
    print dd(h), '\n71830.36653518451'
    print dd(Tdp), '\n-10.012150181172135'
    print dd(v), '\n8.717031296493113'
    """
if __name__ == "__main__":
    test_psychrometric()
