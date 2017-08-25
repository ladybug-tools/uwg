from math import log, pow, exp

def Psychrometrics (Tdb_in, w_in, P):
    """
    Modified version of Psychometrics by Tea Zakula
    MIT Building Technology Lab
    Input: Tdb_in, w_in, P
    Output: Tdb, w, phi, h, Tdp, v

    where:

    Tdb (dry bulb temperature) and Tdp(dew point temperature) in C
    w (humidity ratio) in kg/kg of dry air
    phi (relative humidity) in %
    h (enthalpy) in J/kg of dry air
    v (specific volume) in m3/kg of dry air
    P (Atmospheric Station Pressure) in Pa
    """

    #TODO: Test this!

    c_air = 1006.   # [J/kg] air heat capacity, value from ASHRAE Fundamentals
    hlg = 2501000.  # [J/kg] latent heat, value from ASHRAE Fundamentals
    cw  = 1860.     # [J/kg] value from ASHRAE Fundamentals
    P = P/1000.     # convert from Pa to kPa

    Tdb = Tdb_in - 273.15
    w = w_in

    # phi (RH) calculation from Tdb and w

    Pw = w*P/(0.621945 + w)                             # partial pressure of water vapor [1]
    Pws = saturation_pressure(Tdb)                      # Get saturation pressure for given Tdb
    phi = Pw/Pws*100.0

    # enthalpy calculation from Tdb and w
    h = c_air*Tdb + w*(hlg+cw*Tdb)                      # [J kga-1] [2]

    # specific volume calculation from Tdb and w
    v = 0.287042 * (Tdb+273.15)*(1+1.607858*w)/P        # ?

    # dew point calculation from w
    pw = (P*w)/(0.621945+w) # water vapor partial pressure in kPa
    alpha = log(pw)
    Tdp = 6.54 + 14.526*alpha + 0.7389*pow(alpha,2) +
          0.09486*pow(alpha,3) + 0.4569*pow(pw,0.1984)  # valid for Tdp between 0 C and 93 C

    # [1] Derivation of Pw
    # w = 0.621945 * Pw / (P - Pw)
    # w / 0.622 = Pw / (P - Pw)
    # 0.622 / w = (P - Pw) / Pw
    # 0.622/w + 1/1 = P/Pw
    # (0.622 + w)/w = P/Pw
    # Pw = w*P / (0.622 + w)

    # [2] Derivation for h
    # (J*K-1*kga-1)*K + kgv*kga-1 * (J*K-1*kgv-1 + J*kga-1 * K)

    return  Tdb, w, phi, h, Tdp, v

def saturation_pressure(Tdb):
    T = Tdb+273.15
    Pws = exp(-(5.8002206e3)/T+1.3914993+-(4.8640239e-2)*T +
          (4.1764768e-5)*math(T,2) - (1.4452093e-8)*pow(T,3) + 6.5459673*log(T))  #in Pa
    Pws = Pws/1000.                                                               # in kPa
    return Pws
