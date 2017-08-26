from test import Test
from math import log, pow, exp

def psychrometrics (Tdb_in, w_in, P):
    """
    Modified version of Psychometrics by Tea Zakula
    MIT Building Technology Lab
    Input: Tdb_in, w_in, P
    Output: Tdb, w, phi, h, Tdp, v

    where:
    Tdb_in = [K] dry bulb temperature
    w_in = [kgv/kgda] Humidity Ratio
    P = [P] Atmospheric Station Pressure

    Tdb: [C] dry bulb temperature
    w: [kgv/kgda] Humidity Ratio
    phi: [Pw/Pws*100] relative humidity
    Tdp: [C] dew point temperature
    h: [J/kga] enthalpy
    v: [m3/kga] specific volume
    """
    # Change units
    c_air = 1006.   # [J/kg] air heat capacity, value from ASHRAE Fundamentals
    hlg = 2501000.  # [J/kg] latent heat, value from ASHRAE Fundamentals
    cw  = 1860.     # [J/kg] value from ASHRAE Fundamentals
    P = P/1000.     # convert from Pa to kPa

    Tdb = Tdb_in - 273.15
    w = w_in

    # phi (RH) calculation from Tdb and w
    Pw = (w*P)/(0.621945 + w)                             # partial pressure of water vapor [1]
    Pws = saturation_pressure(Tdb)                      # Get saturation pressure for given Tdb
    phi = Pw/Pws*100.0

    # enthalpy calculation from Tdb and w
    h = c_air*Tdb + w*(hlg+cw*Tdb)                      # [J kga-1] [2]

    # specific volume calculation from Tdb and w
    v = 0.287042 * (Tdb+273.15)*(1+1.607858*w)/P        # ?

    # dew point calculation from w
    pw = (w*P)/(0.621945 + w) # water vapor partial pressure in kPa
    alpha = log(pw)
    Tdp = 6.54 + 14.526*alpha + pow(alpha,2)*0.7389 + pow(alpha,3)*0.09486 + pow(pw,0.1984)*0.4569  # valid for Tdp between 0 C and 93 C

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

def saturation_pressure(Tdb_):
    T = Tdb_ + 273.15
    Pws = exp(-1*(5.8002206e3) / T+1.3914993 + (4.8640239e-2)*T*(-1.) + (4.1764768e-5)*pow(T,2) - (1.4452093e-8)*pow(T,3) + 6.5459673*log(T))  #in Pa
    Pws = Pws/1000.                                                               # in kPa
    return Pws

if __name__ == "__main__":

    psychro_test = Test("test_psychrometrics",run_test=False)
    #Test 1
    Tdb_in = 20.0 + 273.15
    w_in = 0.002
    P = 101325.0
    Tdb, w, phi, h, Tdp, v = psychrometrics(Tdb_in, w_in, P)

    #Tests for 20, 0.002, atmosphere
    psychro_test.test_equality_tol(Tdb+273.15,Tdb_in,True,1e-10) #Tdb [C]
    psychro_test.test_equality_tol(w,w_in,True,1e-10) #W [kgv/kga]
    psychro_test.test_equality_tol(phi,12.,True,2.) #RH [%}
    psychro_test.test_equality_tol(h/1000.,25,True,2.) #enthalpy [J/kga]
    psychro_test.test_equality_tol(v,0.84,True,0.03) #spec. vol [m^3 kg-1]

    # Test 2
    Tdb_in = 40.0 + 273.15
    w_in = 0.009
    P = 101325.0
    Tdb, w, phi, h, Tdp, v = psychrometrics(Tdb_in, w_in, P)

    #Tests for 40, 0.009, atmosphere
    psychro_test.test_equality_tol(Tdb+273.15,Tdb_in,True,1e-10) #Tdb [C]
    psychro_test.test_equality_tol(w,w_in,True,1e-10) #W [kgv/kga]
    psychro_test.test_equality_tol(phi,20.,True,2.) #RH [%}
    psychro_test.test_equality_tol(h/1000.,65,True,2.) #enthalpy [J/kga]
    psychro_test.test_equality_tol(v,0.9,True,0.03) #spec. vol [m^3 kg-1]

    print psychro_test.test_results()
