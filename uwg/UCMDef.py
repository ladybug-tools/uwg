from __future__ import division, print_function

try:
    range = xrange
except NameError:
    pass

from math import sqrt, pow
import copy


class UCMDef(object):
    """
    Definition of Urban Canopy - Building Energy Model Class

    properties
        road;          % Road element class (moved from BEM)

        % Urban Canyon Parameters
        bldHeight;     % average building height (m)
        bldDensity;    % horizontal building density (footprint)
        verToHor;      % vertical-to-horizontal urban area ratio (facade area/urban area)
        treeCoverage;  % horizontal tree density (footprint)
        sensAnthrop;   % sensible anthropogenic heat (other than from buildings) (W m-2)
        latAnthrop;    % latent anthropogenic heat (other than from buildings) (W m-2)
        z0u;           % urban roughness length (m)
        l_disp;          % urban displacement length (m)
        roadShad;      % shadowing of roads
        canWidth;      % canyon width (m)
        bldWidth;      % bld width (m)
        canAspect;     % canyon aspect ratio
        roadConf;      % road-sky configuration factors (sky view factor)
        alb_wall;      % average wall albedo
        wallConf;      % wall-sky configuration factors (sky view factor)
        VFwallroad;    % wall-road view factor
        VFroadwall;    % road-wall view factor
        facArea;       % facade area (m2)
        roadArea;      % road area (m2)
        roofArea;      % roof area (m2) (also building area)
        facAbsor;      % average facade absortivity
        roadAbsor;     % average road absortivity
        h_mix;         % waste heat mix into canyon ratio

        % Urban Canyon Variables
        canTemp;       % canyon air temperature (db) (K)
        Tdp;           % dew point temperature
        Twb;           % wetbulb temperature
        canHum;        % canyon specific humidity (kg kg-1)
        canRHum;       % canyon relative humidity (%)
        canWind;       % urban canyon wind velocity (m s-1)
        turbU;         % canyon turbulent velocities (m s-1)
        turbV;         % canyon turbulent velocities (m s-1)
        turbW;         % canyon turbulent velocities (m s-1)
        ublTemp;       % urban boundary layer temperature (K)
        ublTempdx;     % urban boundary layer temperature discretization (K)
        ublWind;       % urban boundary layer wind velocity (m s-1)
        ustar;         % friction velocity (m s-1)
        ustarMod;      % modified friction velocity (m s-1)
        uExch;         % exchange velocity (m s-1)
        treeLatHeat;   % latent heat from trees (W m-2)
        treeSensHeat;  % sensible heat from trees (W m-2)
        sensHeat;      % urban sensible heat (W m-2)
        latHeat;       % urban latent heat (W m-2)
        windProf;      % urban wind profile
        Q_roof;        % sensible heat flux from building roof (convective)
        Q_wall;        % sensible heat flux from building wall (convective)
        Q_window;      % sensible heat flux from building window (via U-factor)
        Q_road;        % sensible heat flux from road (convective)
        Q_hvac;        % sensible heat flux from HVAC waste
        Q_traffic;     % sensible heat flux from traffic (net)
        Q_ubl;         % Convective heat exchange with UBL layer
        Q_vent;        % Convective heat exchange from ventilation/infiltration
        SolRecWall;    % Solar received by wall
        SolRecRoof;    % Solar received by roof
        SolRecRoad;    % Solar received by road
        roadTemp;      % average road temperature (K)
        roofTemp;      % average roof temperature (K)
        wallTemp;      % average wall temperature (K)
        ElecTotal;     % Total Electricity consumption of urban area
        GasTotal;      % Total Gas consumption of the urban area
    """
    CANYON_TEMP_BOUND_ERROR = "Something obviously went wrong (UCMDef.py)... "

    def __init__(self,bldHeight,bldDensity,verToHor,treeCoverage,sensAnthrop,latAnthrop,
        initialTemp,initialHum,initialWind,parameter,r_glaze,SHGC,alb_wall,road):
        self.road = road                                            # Road element class (moved from BEM)
        self.bldHeight = bldHeight                                  # average building height (m)
        self.verToHor = verToHor                                    # vertical-to-horizontal urban area ratio (facade area/urban area)
        self.bldDensity = bldDensity                                # horizontal building density (footprint/total_area)
        self.treeCoverage = treeCoverage                            # horizontal tree density (footprint)
        self.sensAnthrop = sensAnthrop                              # sensible anthropogenic heat (other than from buildings) (W m-2)
        self.latAnthrop = latAnthrop                                # latent anthropogenic heat (other than from buildings) (W m-2)
        self.roadShad = min(treeCoverage/(1-bldDensity),1)          # fraction of road not building shadowed
        # Key to understanding next few formulas is that UWG_Matlab assumes bld_area = square, so sqrt(bld_area) = side length
        self.bldWidth = 4*bldHeight*bldDensity/verToHor             # bld width (side length) derived from bldDensity and verToHor (m)
        d = self.bldWidth/(sqrt(bldDensity))                        # urban area width == sqrt(bldDensity) == ratio of bld footprint_width / urban area footprint
        self.canWidth = d - self.bldWidth                           # canyon width (m) = urban area width - building width
        self.canAspect = bldHeight/self.canWidth                    # canyon aspect ratio
        self.roadConf = pow(pow(self.canAspect,2)+1,0.5) - self.canAspect   # road-sky configuration factor (sky view factor SVF)
        self.wallConf = 0.5 * (self.canAspect + 1 -
            pow(pow(self.canAspect,2)+1,0.5)) / (self.canAspect)    # wall-sky configuration factor (sky view factor SVF)
        self.facArea = 4*self.bldWidth*bldHeight                    # bld width [m]
        self.roadArea = d*d - pow(self.bldWidth,2)                  # road area [m2] = urban_area_ - bld_area
        self.roofArea = pow(self.bldWidth,2)                        # roof area [m2]
        self.canTemp = initialTemp                                  # canyon air temperature (db) (K)
        self.roadTemp = initialTemp                                 # average road temperature (K)
        self.canHum = initialHum                                    # canyon specific humidity (kg kg-1)
        self.ublWind = max(initialWind,parameter.windMin)           # urban boundary layer wind volocity (m s-1)
        self.canWind = initialWind                                  # urban canyon wind velocity (m s-1)
        self.ustar = 0.1*initialWind                                # friction velocity (m s-1)
        self.ustarMod = 0.1*initialWind                             # modified friction velocity (m s-1)

        # Calculate z0u = urban roughness length (m)
        frontDens = verToHor/4.                                     # density of just street facing facade
        if frontDens < 0.15:
          self.z0u = frontDens * self.bldHeight
        else:
          self.z0u = 0.15 * self.bldHeight

        # Calculate l_dsp = urban displacement length (m)
        if frontDens < 0.05:
          self.l_disp = 3 * frontDens * self.bldHeight
        elif frontDens < 0.15:
          self.l_disp = (0.15+5.5*(frontDens-0.05))*self.bldHeight
        elif frontDens < 1:
          self.l_disp = (0.7+0.35*(frontDens-0.15))*self.bldHeight
        else:
          self.l_disp = 0.5*self.bldHeight

        self.alb_wall = alb_wall                                    # average wall albedo (-)
        self.facAbsor = (1-r_glaze)*(1-alb_wall) + r_glaze*(1-0.75*SHGC) # avg facade absorptivity (-) == wall_mat_fraction * absorption + window_frac * non_solar_heat_gain
        self.roadAbsor = (1-road.vegCoverage)*(1-road.albedo)       # average road absorptivity
        self.sensHeat = 0.0                                         # urban sensible heat [W m-2]
        # Variables set in urbflux()
        self.latHeat = None                                         # urban latent heat [W m-2]
        self.windProf = []                                          # wind profile
        self.canRHum = None
        self.Tdp = None

    def __repr__(self):
        return "UCMDef: ver2Hor = {}, bldDens = {}, canyon H/W = {}/{}, canyon aspect ratio = {}, bld width = {}, roof Area = {}".format(
            self.verToHor,
            self.bldDensity,
            int(self.bldHeight),
            int(self.canWidth),
            round(self.canAspect,1),
            self.facArea,
            self.roofArea
            )

    def UCModel(self,BEM,T_ubl,forc,parameter):
        # Calculate the urban canyon temperature per The uwg (2012) Eq. 10
        dens = forc.pres/(1000*0.287042*self.canTemp*(1.+1.607858*self.canHum))     # air density
        dens_ubl = forc.pres/(1000*0.287042*T_ubl*(1.+1.607858*forc.hum))           # air density
        Cp_air = parameter.cp

        self.Q_wall = 0.
        self.Q_window = 0.
        self.Q_road = 0.
        self.Q_hvac = 0.
        self.Q_traffic = 0.
        self.Q_vent = 0.
        self.Q_ubl = 0.
        self.ElecTotal = 0.
        self.GasTotal = 0.
        self.roofTemp = 0.
        self.wallTemp = 0.

        # Road to Canyon
        T_road = self.road.layerTemp[0]
        h_conv = self.road.aeroCond
        H1 = T_road*h_conv*self.roadArea       # Heat (Sens) from road surface
        H2 = h_conv*self.roadArea

        H1 = H1 + T_ubl*self.roadArea*self.uExch*Cp_air*dens_ubl # Heat from UBL
        H2 = H2 + self.roadArea*self.uExch*Cp_air*dens_ubl
        Q = (self.roofArea+self.roadArea)*(self.sensAnthrop + self.treeSensHeat*self.treeCoverage)

        # Building energy output to canyon, in terms of absolute (total) values
        for j in range(len(BEM)):
            # Re-naming variable for readability
            building = BEM[j].building
            wall = BEM[j].wall
            T_indoor = building.indoorTemp
            T_wall = wall.layerTemp[0]
            R_glazing = building.glazingRatio
            A_wall = (1.-R_glazing)*self.facArea
            A_window = R_glazing*self.facArea
            U_window = building.uValue

            H1 = H1 + BEM[j].frac*(                 # fraction of the urban floor space of this typology
                T_indoor*A_window*U_window +        # window U
                T_wall*A_wall*h_conv +              # Wall conv
                T_indoor*self.roofArea*BEM[j].building.vent*BEM[j].building.nFloor*Cp_air*dens +            # Vent
                T_indoor*self.roofArea*BEM[j].building.infil*self.bldHeight/3600.0*Cp_air*dens)             # Infil

            H2 = H2 + BEM[j].frac*(
                A_window*U_window +
                A_wall*h_conv +
                self.roofArea*BEM[j].building.vent*BEM[j].building.nFloor*Cp_air*dens +                     # Vent
                self.roofArea*BEM[j].building.infil*self.bldHeight/3600.0*Cp_air*dens)                      # Infil

            Q = Q + BEM[j].frac*(
                self.roofArea*building.sensWaste*self.h_mix +           # HVAC waste heat
                A_window*BEM[j].wall.solRec*(1.0-BEM[j].building.shgc))   # heat that didn't make it to inside

            self.wallTemp = self.wallTemp + BEM[j].frac*T_wall
            self.roofTemp = self.roofTemp + BEM[j].frac*BEM[j].roof.layerTemp[0]
            self.Q_ubl = self.Q_ubl + BEM[j].frac*self.bldDensity*(BEM[j].roof.sens + BEM[j].building.sensWaste*(1.-self.h_mix)) # Changed by Jiachen Mao in March 2017

        # Solve for canyon temperature
        self.canTemp = (H1 + Q)/H2

        # Heat flux based per m^2 of urban area
        self.Q_road = h_conv*(T_road-self.canTemp)*(1.-self.bldDensity)  # Sensible heat from road (W/m^2 of urban area)
        self.Q_ubl = self.Q_ubl + self.uExch*Cp_air*dens*(self.canTemp-T_ubl)*(1.-self.bldDensity)
        self.Q_wall = h_conv*(self.wallTemp-self.canTemp)*(self.verToHor)
        self.Q_traffic = self.sensAnthrop

        # Building energy output to canyon, per m^2 of urban area
        T_can = copy.copy(self.canTemp)
        for j in range(len(BEM)):
            V_vent = BEM[j].building.vent*BEM[j].building.nFloor  # ventilation volume per m^2 of building
            V_infil = BEM[j].building.infil*self.bldHeight/3600.0
            T_indoor = BEM[j].building.indoorTemp
            U_window = BEM[j].building.uValue                     # Added by Jiachen Mao in March 2017
            R_glazing = BEM[j].building.glazingRatio              # Changed by Jiachen Mao in March 2017

            self.Q_window = self.Q_window + BEM[j].frac*self.verToHor*R_glazing*U_window*(T_indoor-T_can)
            self.Q_window = self.Q_window + BEM[j].frac*self.verToHor*R_glazing*BEM[j].wall.solRec*(1.-BEM[j].building.shgc)
            self.Q_vent = self.Q_vent + BEM[j].frac*self.bldDensity*Cp_air*dens*(V_vent + V_infil)*(T_indoor-T_can)
            self.Q_hvac = self.Q_hvac + BEM[j].frac*self.bldDensity*BEM[j].building.sensWaste*self.h_mix

            self.Q_roof = self.Q_roof + BEM[j].frac*self.bldDensity*BEM[j].roof.sens

            # Total Electrical & Gas power in MW
            self.ElecTotal = self.ElecTotal + BEM[j].fl_area*BEM[j].building.ElecTotal/1.e6
            self.GasTotal = self.GasTotal + BEM[j].fl_area*BEM[j].building.GasTotal/1.e6

        # Sensible Heat
        # N.B In the current uwg code, latent heat from evapotranspiration, stagnant water,
        # or anthropogenic sources is not modelled due to the difficulty of validation, and
        # lack of reliability of precipitation data from EPW files.
        self.sensHeat = self.Q_wall + self.Q_road + self.Q_vent + self.Q_window + self.Q_hvac + self.Q_traffic + self.treeSensHeat + self.Q_roof

        # Error checking
        if self.canTemp > 350. or self.canTemp < 250:
            raise Exception(self.CANYON_TEMP_BOUND_ERROR)
