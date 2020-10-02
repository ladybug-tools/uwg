"""Functions for sychrometric calculations."""
from __future__ import division

from math import log, pow, exp


def psychrometrics(Tdb_in, w_in, P):
    """Modified version of Psychometrics by Tea Zakula MIT Building Technology Lab.

    Args:
        Tdb_in: Number for dry bulb temperature (K).
        w_in: Number for Humidity Ratio (kgv/kgda).
        P: Atmospheric Station Pressure (P).

    Returns:
        A tuple with six values

        - Tdb: Number for dry bulb temperature in (C).

        - w: Number for Humidity Ratio (kgv/kgda).

        - phi: Number for relative humidity (Pw/Pws*100).

        - h: Number for enthalpy (J/kga).

        - Tdp: Number for dew point temperature (C).

        - v: Number for specific volume (m3/kga).
    """
    # Change units
    c_air = 1006.  # [J/kg] air heat capacity, value from ASHRAE Fundamentals
    hlg = 2501000.  # [J/kg] latent heat, value from ASHRAE Fundamentals
    cw = 1860.  # [J/kg] value from ASHRAE Fundamentals
    P = P / 1000.  # convert from Pa to kPa

    Tdb = Tdb_in - 273.15
    w = w_in

    # phi (RH) calculation from Tdb and w
    Pw = (w * P) / (0.621945 + w)  # partial pressure of water vapor
    Pws = saturation_pressure(Tdb)  # Get saturation pressure for given Tdb
    phi = Pw / Pws * 100.0

    # enthalpy calculation from Tdb and w
    h = c_air * Tdb + w * (hlg + cw * Tdb)  # [J kga-1]

    # specific volume calculation from Tdb and w
    v = 0.287042 * (Tdb + 273.15) * (1 + 1.607858 * w) / P

    # dew point calculation from w
    _pw = (w * P) / (0.621945 + w)  # water vapor partial pressure in kPa
    alpha = log(_pw)

    Tdp = 6.54 + 14.526 * alpha + pow(alpha, 2) * 0.7389 + pow(alpha, 3) * 0.09486 + \
        pow(_pw, 0.1984) * 0.4569  # valid for Tdp between 0 C and 93 C

    return Tdb, w, phi, h, Tdp, v


def saturation_pressure(Tdb_):
    """Saturation pressure from dry bulb pressure.

    Args:
        Tdb_: Dry bulb temperature (K).

    Returns:
        Saturation pressure.
    """
    T = Tdb_ + 273.15

    # N.B In Matlab, negative values are converted to complex values.
    # log(-x) = log(x) + log(-1) = log(x) + i*pi
    # Python will throw an exception. Negative value occurs here if simulation
    # timestep (dtSim) is large, i.e 3600s.
    _Pws = exp(-1 * 5.8002206e3 / T + 1.3914993 + 4.8640239e-2 * T * -1.0 +
               4.1764768e-5 * pow(T, 2) - 1.4452093e-8 * pow(T, 3) +
               6.5459673 * log(T))  # Pa

    return _Pws / 1000.  # kPa


def moist_air_density(P, Tdb, H):
    """Moist air density (kgv m-3) from dry bulb temp, HR, and pressure[1][2][3].

    Note:
        [1] ASHRAE Fundamentals (2005) ch. 6 eqn. 28
        [2] ASHRAE Fundamentals (2009) ch. 1 eqn. 28
        [3] https://github.com/psychrometrics/Libraries/blob/master/Psychrometrics_SI.cpp

    Args:
        P: Number for Pressure (P).
        Tdb: Number for dry bulb temperature in (C).
        H: Number for Humidity Ratio (kgv/kgda).

    Returns:
        Moist air density (kgv m-3).
    """

    return P / (1000 * 0.287042 * Tdb * (1. + 1.607858 * H))


def hum_from_rhum_temp(RH, T, P):
    """Specific Humidity (kgh20/kgn202) from RH, T and Pa.

    Args:
        RH: Number for Relative Humidity.
        T: Number for dry bulb temperature (C).
        P: Number for air pressure (P).

    Returns:
        Specific Humidity (kgh20/kgn202).
    """
    # Saturation vapour pressure from ASHRAE
    C8 = -5.8002206e3
    C9 = 1.3914993
    C10 = -4.8640239e-2
    C11 = 4.1764768e-5
    C12 = -1.4452093e-8
    C13 = 6.5459673

    T += 273.15

    PWS = exp(C8 / T + C9 + C10 * T + C11 * pow(T, 2) + C12 * pow(T, 3) + C13 * log(T))
    PW = RH * PWS / 100.0  # Vapour pressure

    return 0.62198 * PW / (P - PW)  # 4. Specific humidity
