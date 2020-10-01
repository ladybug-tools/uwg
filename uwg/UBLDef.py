"""Class definition for the Urban Boundary Layer (UBL)."""
from __future__ import division

try:
    range = xrange
except NameError:
    pass

import logging
from .utilities import is_near_zero


class UBLDef(object):
    """Urban Boundary Layer (UBL) calculations.

    Args:
        location: Text string for relative location within a city. Choose from "N", "NE",
            "E", "SE", "S", "SW", "W", "NW", or "C".
        charLength: Value for characteristic length of the urban area [m].
        initialTemp: Value for initial temperature [K].
        maxdx: Value for maximum discretization length for the UBL model [m].
        dayBLHeight: Value for daytime mixing height [m].
        nightBLHeight: Value for nighttime boundary-layer height [m].

    Properties:
        * location -- relative location within a city (N, NE, E, SE, S, SW, W, NW, C).
        * charLength -- characteristic length of the urban area (m)
        * perimeter -- horizontal urban area (m2)
        * urbArea -- length of the side of the urban area orthogonal
        * orthLength -- length to side of urban area orthogonal to wind direction (m)
        * paralLength -- length of side of urban area parallel to wind direction (m)
        * ublTemp -- urban boundary layer temperature (K)
        * ublTempdx -- urban boundary layer temperature discretization (K)
        * dayBLHeight -- daytime mixing height, orig = 700
        * nightBLHeight -- nighttime boundary-layer ht (m), Sing: 80, Bub-Cap: 50, orig 80
    """

    def __init__(self, location, charLength, initialTemp, maxdx, dayBLHeight,
                 nightBLHeight):
        self.location = location
        self.charLength = charLength
        self.perimeter = 4. * charLength
        self.urbArea = charLength ** 2
        self.orthLength = charLength

        numdx = round(charLength / min(charLength, maxdx))
        self.paralLength = charLength / numdx
        self.ublTemp = initialTemp
        self.ublTempdx = [initialTemp for x in range(int(numdx))]
        self.dayBLHeight = dayBLHeight
        self.nightBLHeight = nightBLHeight

        # Logger will be disabled by default unless explicitly called in tests
        self.logger = logging.getLogger(__name__)

    def ublmodel(self, UCM, RSM, rural, forc, parameter, simTime):
        # Note that only one urban canyon area is considered
        self.sensHeat = UCM.sensHeat
        heatDif = max(self.sensHeat - rural.sens, 0)
        Cp = parameter.cp  # Heat capacity of air (J/kg.K)
        k_w = parameter.circCoeff  # k_w per Bueno 'the uwg', eq 8
        g = parameter.g  # Gravity
        v_wind = max(forc.wind, parameter.windMin)  # wind velocity

        # Air density
        refDens = 0.
        for iz in range(RSM.nzref):
            refDens = \
                refDens + RSM.densityProfC[iz] * RSM.dz[iz] / \
                (RSM.z[RSM.nzref - 1] + RSM.dz[RSM.nzref - 1] / 2.)

        forDens = 0.
        for iz in range(RSM.nzfor):
            forDens = \
                forDens + RSM.densityProfC[iz] * RSM.dz[iz] / \
                (RSM.z[RSM.nzfor - 1] + RSM.dz[RSM.nzfor - 1] / 2.)

        # ---------------------------------------------------------------------
        # Day
        # ---------------------------------------------------------------------
        time = simTime.secDay / 3600.
        noon = 12.
        daylimit = parameter.dayThreshold      # sunlight threshold for day (~150W/m^2)
        nightlimit = parameter.dayThreshold    # sunlight threshold for night (~50W/m^2)
        sunlight = forc.dir + forc.dif

        # If dir & dif light is greater than threshold, use day
        is_day = (sunlight > daylimit) and (time < noon or is_near_zero(time - noon)) \
            or (sunlight > nightlimit) and (time > noon) or (self.sensHeat > 150.0)
        if is_day:
            # Circulation velocity per Bueno 'the uwg', eq 8
            self.logger.debug("{} Day ubl calcs".format(__name__))
            h_UBL = self.dayBLHeight            # Day boundary layer height
            eqTemp = RSM.tempProf[RSM.nzref - 1]
            eqWind = RSM.windProf[RSM.nzref - 1]

            Csurf = UCM.Q_ubl * simTime.dt / (h_UBL * refDens * Cp)
            u_circ = k_w * (g * heatDif / Cp / refDens / eqTemp * h_UBL) ** (1. / 3.)

            if v_wind > u_circ:  # Forced problem (usually this)
                advCoef = self.orthLength*eqWind * simTime.dt / self.urbArea * 1.4
                self.ublTemp = (Csurf + advCoef * eqTemp + self.ublTemp) / (1. + advCoef)
                self.ublTempdx = [self.ublTemp for x in range(len(self.ublTempdx))]
            else:  # Convective problem
                advCoef = self.perimeter * u_circ * simTime.dt / self.urbArea * 1.4
                self.ublTemp = (Csurf + advCoef * eqTemp + self.ublTemp) / (1 + advCoef)
                self.ublTempdx = [self.ublTemp for x in range(len(self.ublTempdx))]

        # ---------------------------------------------------------------------
        # Night
        # ---------------------------------------------------------------------
        else:
            self.logger.debug("{} Night ubl calcs".format(__name__))
            h_UBL = self.nightBLHeight  # Night boundary layer height
            Csurf = UCM.Q_ubl * simTime.dt / (h_UBL * refDens * Cp)
            self.ublTemp, self.ublTempdx = \
                UBLDef.nightforc(self.ublTempdx, simTime.dt, h_UBL, self.paralLength,
                                 self.charLength, RSM, Csurf)

        self.logger.debug("ublTemp = {}".format(self.ublTemp))

    @staticmethod
    def nightforc(ublTempdx, dt, h_UBL, paralLength, charLength, RSM, Csurf):
        # Night forcing (RSM.nzfor = number of layers of forcing)
        # Average potential temperature & wind speed of the profile
        intAdv1 = 0.
        for iz in range(RSM.nzfor):
            intAdv1 = intAdv1 + RSM.windProf[iz] * RSM.tempProf[iz] * RSM.dz[iz]

        advCoef1 = 1.4 * dt / paralLength / h_UBL * intAdv1

        intAdv2 = 0
        for iz in range(RSM.nzfor):
            intAdv2 = intAdv2 + RSM.windProf[iz] * RSM.dz[iz]

        advCoef2 = 1.4 * dt / paralLength / h_UBL * intAdv2

        ublTempdx[0] = (Csurf + advCoef1 + ublTempdx[0]) / (1 + advCoef2)
        ublTemp = ublTempdx[0]

        for i in range(1, int(charLength) // int(paralLength)):
            eqTemp = ublTempdx[i - 1]
            ublTempdx[i] = (Csurf + advCoef2*eqTemp + ublTempdx[i]) / (1 + advCoef2)
            ublTemp = ublTemp + ublTempdx[i]

        ublTemp = ublTemp / float(charLength) * float(paralLength)

        return ublTemp, ublTempdx

    def __repr__(self):
        return 'UBL,\n location: {}\n charLength: {}\n perimeter: {}\n urbArea: {}\n ' \
            'orthLength: {}\n paralLength: {}\n ublTemp: {}\n dayBLHeight: {}\n ' \
            'nightBLHeight: {}'.format(
                self.location, self.charLength, self.perimeter, self.urbArea,
                self.orthLength, self.paralLength, self.ublTemp, self.dayBLHeight,
                self.nightBLHeight)
