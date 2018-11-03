from __future__ import division

try:
    range = xrange
except NameError:
    pass

import logging


class UBLDef(object):
    """
    %   Class definition for the Urban Boundary Layer (UBL)

    properties
        location;       % relative location within a city (N,NE,E,SE,S,SW,W,NW,C)
        charLength;     % characteristic length of the urban area (m)
        perimeter;      % urban area perimeter (m)
        urbArea;        % horizontal urban area (m2)
        orthLength;     % length of the side of the urban area orthogonal
                        % to the wind direction (m)
        paralLength;    % length of the side of the urban area parallel
                        % to the wind direction (m)
        ublTemp;        % urban boundary layer temperature (K)
        ublTempdx;      % urban boundary layer temperature discretization (K)
        advHeat;        % advection heat flux (W m-2)
        sensHeat;       % sensible heat flux (W m-2)
        dayBLHeight;    % daytime mixing height, orig = 700
        nightBLHeight;  % Sing: 80, Bub-Cap: 50, nighttime boundary-layer height (m); orig 80
    end
    """

    def __init__(self,location,charLength,initialTemp,maxdx,dayBLHeight,nightBLHeight):
        self.location = location                                    # relative location within a city (N,NE,E,SE,S,SW,W,NW,C)
        self.charLength = charLength                                # characteristic length of the urban area (m)
        self.perimeter = 4. * charLength
        self.urbArea = charLength**2                                # horizontal urban area (m2)
        self.orthLength = charLength                                # length of the side of the urban area orthogonal
                                                                    # to the wind direction (m)
        numdx = round(charLength/min(charLength,maxdx))
        self.paralLength = charLength/numdx                         # length of the side of the urban area parallel
                                                                    # to the wind direction (m)
        self.ublTemp = initialTemp                                  # urban boundary layer temperature (K)
        self.ublTempdx = [initialTemp for x in range(int(numdx))]   # urban boundary layer temperature discretization (K)
        self.dayBLHeight = dayBLHeight                              # daytime mixing height, orig = 700
        self.nightBLHeight = nightBLHeight                          # Sing: 80, Bub-Cap: 50, nighttime boundary-layer height (m); orig 80

        # Logger will be disabled by default unless explicitly called in tests
        self.logger = logging.getLogger(__name__)

    def __repr__(self):
        return "UBL: urbArea {}m2, charLength {}m".format(
            self.urbArea,
            self.charLength
            )

    def is_near_zero(self, num, eps=1e-14):
        return abs(float(num)) < eps

    def UBLModel(self,UCM,RSM,rural,forc,parameter,simTime):
        # Note that only one urban canyon area is considered
        self.sensHeat = UCM.sensHeat
        heatDif = max(self.sensHeat - rural.sens,0)
        Cp = parameter.cp                           # Heat capacity of air (J/kg.K)
        k_w = parameter.circCoeff                   # k_w per Bueno 'the uwg', eq 8
        g = parameter.g                             # Gravity
        v_wind = max(forc.wind,parameter.windMin)   # wind velocity

        # Air density
        refDens = 0.
        for iz in range(RSM.nzref):
            refDens = refDens + RSM.densityProfC[iz] * RSM.dz[iz] / (RSM.z[RSM.nzref-1] + RSM.dz[RSM.nzref-1]/2.)

        forDens = 0.
        for iz in range(RSM.nzfor):
            forDens = forDens + RSM.densityProfC[iz] * RSM.dz[iz] / (RSM.z[RSM.nzfor-1] + RSM.dz[RSM.nzfor-1]/2.)

        # ---------------------------------------------------------------------
        # Day
        # ---------------------------------------------------------------------
        time = simTime.secDay/3600.
        noon = 12.
        daylimit = parameter.dayThreshold      # sunlight threshold for day (~150W/m^2)
        nightlimit = parameter.dayThreshold    # sunlight threshold for night (~50W/m^2)
        sunlight = forc.dir + forc.dif

        # If dir & dif light is greater than threshold, use day
        is_day = (sunlight > daylimit) and (time < noon or self.is_near_zero(time-noon)) \
                or (sunlight > nightlimit) and (time > noon) or (self.sensHeat > 150.0)
        if is_day:
            # Circulation velocity per Bueno 'the uwg', eq 8
            self.logger.debug("{} Day ubl calcs".format(__name__))
            h_UBL = self.dayBLHeight            # Day boundary layer height
            eqTemp = RSM.tempProf[RSM.nzref-1]
            eqWind = RSM.windProf[RSM.nzref-1]

            Csurf = UCM.Q_ubl*simTime.dt/(h_UBL*refDens*Cp)
            u_circ = k_w*(g*heatDif/Cp/refDens/eqTemp*h_UBL)**(1./3.)

            if v_wind > u_circ:   # Forced problem (usually this)
                advCoef  = self.orthLength*eqWind*simTime.dt/self.urbArea*1.4
                self.ublTemp = (Csurf + advCoef * eqTemp + self.ublTemp)/(1. + advCoef)
                self.ublTempdx = [self.ublTemp for x in range(len(self.ublTempdx))]
            else:                   # Convective problem
                advCoef  = self.perimeter*u_circ*simTime.dt/self.urbArea*1.4
                self.ublTemp = (Csurf+advCoef*eqTemp + self.ublTemp)/(1 + advCoef)
                self.ublTempdx = [self.ublTemp for x in range(len(self.ublTempdx))]



        # ---------------------------------------------------------------------
        # Night
        # ---------------------------------------------------------------------
        else:
            self.logger.debug("{} Night ubl calcs".format(__name__))
            h_UBL = self.nightBLHeight      # Night boundary layer height
            Csurf = UCM.Q_ubl*simTime.dt/(h_UBL*refDens*Cp)
            self.ublTemp, self.ublTempdx = self.NightForc(self.ublTempdx,simTime.dt, \
                h_UBL,self.paralLength,self.charLength,RSM,Csurf)

        self.logger.debug("ublTemp = {}".format(self.ublTemp))

    def NightForc(self,ublTempdx,dt,h_UBL,paralLength,charLength,RSM,Csurf):
        # Night forcing (RSM.nzfor = number of layers of forcing)
        # Average potential temperature & wind speed of the profile
        intAdv1 = 0.
        for iz in range(RSM.nzfor):
            intAdv1 = intAdv1 + RSM.windProf[iz] * RSM.tempProf[iz] * RSM.dz[iz]

        advCoef1 = 1.4*dt/paralLength/h_UBL*intAdv1

        intAdv2 = 0
        for iz in range(RSM.nzfor):
            intAdv2 = intAdv2 + RSM.windProf[iz]*RSM.dz[iz]

        advCoef2 = 1.4*dt/paralLength/h_UBL*intAdv2

        ublTempdx[0] = (Csurf + advCoef1 + ublTempdx[0])/(1 + advCoef2)
        ublTemp = ublTempdx[0]

        for i in range(1,int(charLength)//int(paralLength)):
            eqTemp = ublTempdx[i-1]
            ublTempdx[i] = (Csurf + advCoef2*eqTemp + ublTempdx[i])/(1 + advCoef2)
            ublTemp = ublTemp + ublTempdx[i]

        # ublTemp/charLength*paralLength;
        ublTemp = ublTemp/float(charLength)*float(paralLength)

        return ublTemp, ublTempdx
