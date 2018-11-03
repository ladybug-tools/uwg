import pytest
from uwg import psychrometrics
from math import log10

import pprint
import decimal

dd = decimal.Decimal.from_float
pp = pprint.pprint


def test_psychrometric_float_point():

    # Input values
    Tdb_in = 297.5311337413935
    w_in = 0.018576773131376
    P = 10090

    Tdb, w, phi, h, Tdp, v = psychrometrics.psychrometrics(Tdb_in, w_in, P)

    assert Tdb == pytest.approx(24.381133741393512, abs=1e-15)  # dry bulb temperature [C]
    assert w == pytest.approx(0.018576773131376, abs=1e-15)     # Humidity ratio [kgv/kgd]
    assert phi == pytest.approx(9.581555922752541, abs=1e-15)   # Relative Humidity Pw/Pws*100
    assert h == pytest.approx(7.183036653518451e+04, abs=1e-15) # enthalpy [J/kgd]
    assert Tdp == pytest.approx(-10.012150181172135, abs=1e-15) # Wet bulb temp [C]
    assert v == pytest.approx(8.717031296493113, abs=1e-15)     # specific volume [m3/kga]s


def test_psychrometric_simple_1():
    # Really simple, coarse tests
    # Input values
    Tdb_in = 20.0 + 273.15
    w_in = 0.002
    P = 101325.0

    Tdb, w, phi, h, Tdp, v = psychrometrics.psychrometrics(Tdb_in, w_in, P)

    #Tests for 20, 0.002, atmosphere
    assert Tdb+273.15 == pytest.approx(Tdb_in,True,1e-14)       # Tdb [C]
    assert w == pytest.approx(w_in,1e-14)                       # W [kgv/kga]
    assert phi == pytest.approx(13.,1)                          # RH [%}
    assert h/1000. == pytest.approx(25.,1)                      # Enthalpy [J/kga]
    assert v == pytest.approx(0.83,1e-2)                        # Spec. vol [m^3 kg-1]

def test_psychrometric_simple_2():
    # Really simple, coarse tests
    # Input values
    Tdb_in = 40.0 + 273.15
    w_in = 0.009
    P = 101325.0

    Tdb, w, phi, h, Tdp, v = psychrometrics.psychrometrics(Tdb_in, w_in, P)

    #Tests for 40, 0.009, atmosphere
    assert Tdb+273.15 == pytest.approx(Tdb_in,True,1e-14)       # Tdb [C]
    assert w == pytest.approx(w_in,1e-14)                       # W [kgv/kga]
    assert phi == pytest.approx(19.5,1e-1)                      # RH [%}
    assert h/1000. == pytest.approx(63.,1)                      # Enthalpy [J/kga]
    assert v == pytest.approx(0.9,1e-1)                         # Spec. vol [m^3 kg-1]

if __name__ == "__main__":
    pass
    #test_psychrometric_simple_1()
    #test_psychrometric_simple_2()
    #test_psychrometric_float_point()
