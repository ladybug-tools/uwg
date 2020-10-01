"""Class for Urban Canopy Model."""
from __future__ import division

try:
    range = xrange
except NameError:
    pass

from math import sqrt, pow
import copy


class UCMDef(object):
    """Definition of Urban Canopy - Building Energy Model Class

    Args:
        bldHeight: Average building height (m).
        bldDensity: Horizontal building density (footprint/total_area).
        verToHor: Vertical-to-horizontal urban area ratio (facade area/urban area).
        treeCoverage: Horizontal tree density (footprint).
        sensAnthrop: Sensible anthropogenic heat (other than from buildings) (W m-2).
        latAnthrop: Latent anthropogenic heat (other than from buildings) (W m-2).
        initialTemp: Canyon air temperature (db) (K).
        initialHum: Canyon specific humidity (kg kg-1).
        initialWind: Urban canyon wind velocity (m s-1).
        parameter: Param object.
        r_glaze: Area-weighted average of glazing ratio from urban building stock.
        SHGC: Area-weighted average of SHGC from urban building stock.
        alb_wall: Area-weighted average of wall albedo from urban building stock.
        road: Road element class (moved from BEM).

    Properties
        * road -- Road element class (moved from BEM)
        * # Urban Canyon Parameters
        * bldHeight -- average building height (m)
        * bldDensity -- horizontal building density (footprint)
        * verToHor -- vertical-to-horizontal urban area ratio (facade area/urban area)
        * treeCoverage -- horizontal tree density (footprint)
        * sensAnthrop -- sensible anthropogenic heat (other than from buildings) (W m-2)
        * latAnthrop -- latent anthropogenic heat (other than from buildings) (W m-2)
        * z0u -- urban roughness length (m)
        * l_disp -- urban displacement length (m)
        * roadShad -- shadowing of roads
        * canWidth -- canyon width (m)
        * bldWidth -- bld width (m)
        * canAspect -- canyon aspect ratio
        * roadConf -- road-sky configuration factors (sky view factor)
        * alb_wall -- average wall albedo
        * wallConf -- wall-sky configuration factors (sky view factor)
        * VFwallroad -- wall-road view factor
        * VFroadwall -- road-wall view factor
        * facArea -- facade area (m2)
        * roadArea -- road area (m2)
        * roofArea -- roof area (m2) (also building area)
        * facAbsor -- average facade absortivity
        * roadAbsor -- average road absortivity
        * h_mix -- waste heat mix into canyon ratio
        * # Urban Canyon Variables
        * canTemp -- canyon air temperature (db) (K)
        * Tdp -- dew point temperature
        * Twb -- wetbulb temperature
        * canHum -- canyon specific humidity (kg kg-1)
        * canRHum -- canyon relative humidity ()
        * canWind -- urban canyon wind velocity (m s-1)
        * turbU -- canyon turbulent velocities (m s-1)
        * turbV --  canyon turbulent velocities (m s-1)
        * turbW -- canyon turbulent velocities (m s-1)
        * ublTemp -- urban boundary layer temperature (K)
        * ublTempdx -- urban boundary layer temperature discretization (K)
        * ublWind -- urban boundary layer wind velocity (m s-1)
        * ustar -- friction velocity (m s-1)
        * ustarMod -- modified friction velocity (m s-1)
        * uExch -- exchange velocity (m s-1)
        * treeLatHeat -- latent heat from trees (W m-2)
        * treeSensHeat -- sensible heat from trees (W m-2)
        * sensHeat -- urban sensible heat (W m-2)
        * latHeat -- urban latent heat (W m-2)
        * windProf -- urban wind profile
        * Q_roof -- sensible heat flux from building roof (convective)
        * Q_wall -- sensible heat flux from building wall (convective)
        * Q_window -- sensible heat flux from building window (via U-factor)
        * Q_road -- sensible heat flux from road (convective)
        * Q_hvac -- sensible heat flux from HVAC waste
        * Q_traffic -- sensible heat flux from traffic (net)
        * Q_ubl -- Convective heat exchange with UBL layer
        * Q_vent -- Convective heat exchange from ventilation/infiltration
        * SolRecWall -- Solar received by wall
        * SolRecRoof -- Solar received by roof
        * SolRecRoad -- Solar received by road
        * roadTemp -- average road temperature (K)
        * roofTemp -- average roof temperature (K)
        * wallTemp -- average wall temperature (K)
        * ElecTotal -- Total Electricity consumption of urban area
        * GasTotal -- Total Gas consumption of the urban area
    """

    def __init__(self, bldHeight, bldDensity, verToHor, treeCoverage, sensAnthrop,
                 latAnthrop, initialTemp, initialHum, initialWind, parameter,
                 r_glaze, SHGC, alb_wall, road):

        self.road = road
        self.bldHeight = bldHeight
        self.verToHor = verToHor
        self.bldDensity = bldDensity
        self.treeCoverage = treeCoverage
        self.vegcover = (1 - self.bldDensity) * self.road.vegcoverage
        self.sensAnthrop = sensAnthrop
        self.latAnthrop = latAnthrop
        self.roadShad = min(treeCoverage / (1 - bldDensity), 1)
        # UWG_Matlab assumes bld_area = square, so sqrt(bld_area) = side length
        # bld width: (side length) derived from bldDensity and verToHor (m)
        self.bldWidth = 4 * bldHeight * bldDensity / verToHor
        # urban area width == sqrt(bldDensity) == \
        # ratio of bld footprint_width / urban area footprint
        d = self.bldWidth / (sqrt(bldDensity))
        # canyon width (m) = urban area width - building width
        self.canWidth = d - self.bldWidth
        self.canAspect = bldHeight / self.canWidth
        self.roadConf = pow(pow(self.canAspect, 2) + 1, 0.5) - self.canAspect
        self.wallConf = (
            0.5 * (self.canAspect + 1 - pow(pow(self.canAspect, 2) + 1, 0.5)) /
            self.canAspect)
        self.facArea = 4 * self.bldWidth * bldHeight
        self.roadArea = d * d - pow(self.bldWidth, 2)
        self.roofArea = pow(self.bldWidth, 2)
        self.canTemp = initialTemp
        self.roadTemp = initialTemp
        self.canHum = initialHum
        self.ublWind = max(initialWind, parameter.windMin)
        self.canWind = initialWind
        self.ustar = 0.1 * initialWind
        self.ustarMod = 0.1 * initialWind

        # Calculate z0u = urban roughness length (m)
        frontDens = verToHor / 4.  # density of just street facing facade
        if frontDens < 0.15:
            self.z0u = frontDens * self.bldHeight
        else:
            self.z0u = 0.15 * self.bldHeight

        # Calculate l_dsp = urban displacement length (m)
        if frontDens < 0.05:
            self.l_disp = 3 * frontDens * self.bldHeight
        elif frontDens < 0.15:
            self.l_disp = (0.15 + 5.5 * (frontDens - 0.05)) * self.bldHeight
        elif frontDens < 1:
            self.l_disp = (0.7 + 0.35 * (frontDens - 0.15)) * self.bldHeight
        else:
            self.l_disp = 0.5 * self.bldHeight

        self.alb_wall = alb_wall
        self.facAbsor = (1 - r_glaze) * (1 - alb_wall) + r_glaze * (1 - 0.75 * SHGC)
        self.roadAbsor = (1 - road.vegcoverage) * (1 - road.albedo)
        self.sensHeat = 0.0

        # Variables set in urbflux function
        self.latHeat = None
        self.windProf = []
        self.canRHum = None
        self.Tdp = None

    def UCModel(self, BEM, T_ubl, forc, parameter):
        """Calculate the urban canyon temperature per The uwg (2012) Eq. 10."""

        # air density
        dens = \
            forc.pres / (1000 * 0.287042 * self.canTemp * (1. + 1.607858 * self.canHum))
        dens_ubl = forc.pres / (1000 * 0.287042 * T_ubl * (1. + 1.607858 * forc.hum))
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
        H1 = T_road * h_conv * self.roadArea  # Heat (Sens) from road surface
        H2 = h_conv * self.roadArea

        H1 = H1 + T_ubl * self.roadArea * self.uExch * Cp_air * dens_ubl  # Heat from UBL
        H2 = H2 + self.roadArea * self.uExch * Cp_air * dens_ubl
        # W = m2 * W/m2
        Q = (self.roofArea + self.roadArea) * (self.sensAnthrop + self.treeSensHeat)

        # Building energy output to canyon, in terms of absolute (total) values
        for j in range(len(BEM)):
            # Re-naming variable for readability
            building = BEM[j].building
            wall = BEM[j].wall
            T_indoor = building.indoor_temp
            T_wall = wall.layerTemp[0]
            R_glazing = building.glazing_ratio
            A_wall = (1. - R_glazing) * self.facArea
            A_window = R_glazing*self.facArea
            U_window = building.u_value

            H1 = H1 + BEM[j].frac*(
                T_indoor * A_window * U_window +  # window U
                T_wall * A_wall * h_conv +  # wall conv
                (T_indoor * self.roofArea * BEM[j].building.vent *
                 BEM[j].building.nFloor * Cp_air * dens) +  # Vent
                (T_indoor * self.roofArea * BEM[j].building.infil * self.bldHeight /
                 3600.0 * Cp_air * dens))  # Infil

            H2 = H2 + BEM[j].frac * (
                A_window * U_window +  # window U
                A_wall * h_conv +  # wall conv
                (self.roofArea * BEM[j].building.vent * BEM[j].building.nFloor *
                 Cp_air * dens) +  # Vent
                (self.roofArea * BEM[j].building.infil * self.bldHeight / 3600.0 *
                 Cp_air * dens))  # Infil

            Q = Q + BEM[j].frac * (
                self.roofArea * building.sensWaste * self.h_mix +  # HVAC waste heat
                # heat that didn't make it to inside
                A_window * BEM[j].wall.solRec * (1.0 - BEM[j].building.shgc))

            self.wallTemp = self.wallTemp + BEM[j].frac * T_wall
            self.roofTemp = self.roofTemp + BEM[j].frac * BEM[j].roof.layerTemp[0]
            self.Q_ubl = (self.Q_ubl + BEM[j].frac * self.bldDensity *
                          (BEM[j].roof.sens + BEM[j].building.sensWaste *
                          (1. - self.h_mix)))  # Changed by Jiachen Mao in March 2017

        # Solve for canyon temperature
        self.canTemp = (H1 + Q) / H2

        # Heat flux based per m^2 of urban area
        # Sensible heat from road (W/m^2 of urban area)
        self.Q_road = h_conv * (T_road - self.canTemp) * (1. - self.bldDensity)
        self.Q_ubl = (
            self.Q_ubl + self.uExch * Cp_air * dens * (self.canTemp - T_ubl) *
            (1. - self.bldDensity))
        self.Q_wall = h_conv * (self.wallTemp - self.canTemp) * self.verToHor
        self.Q_traffic = self.sensAnthrop

        # Building energy output to canyon, per m^2 of urban area
        T_can = copy.copy(self.canTemp)
        for j in range(len(BEM)):
            # ventilation volume per m^2 of building
            V_vent = BEM[j].building.vent * BEM[j].building.nFloor
            V_infil = BEM[j].building.infil * self.bldHeight / 3600.0
            T_indoor = BEM[j].building.indoor_temp
            U_window = BEM[j].building.u_value  # Added by Jiachen Mao in 03/17
            R_glazing = BEM[j].building.glazing_ratio  # Changed by Jiachen Mao in 03/17

            self.Q_window = (self.Q_window + BEM[j].frac * self.verToHor * R_glazing *
                             U_window * (T_indoor - T_can))
            self.Q_window = (self.Q_window + BEM[j].frac * self.verToHor * R_glazing *
                             BEM[j].wall.solRec * (1. - BEM[j].building.shgc))
            self.Q_vent = (self.Q_vent + BEM[j].frac * self.bldDensity * Cp_air * dens *
                           (V_vent + V_infil) * (T_indoor - T_can))
            self.Q_hvac = (self.Q_hvac + BEM[j].frac * self.bldDensity *
                           BEM[j].building.sensWaste * self.h_mix)

            self.Q_roof = self.Q_roof + BEM[j].frac * self.bldDensity * BEM[j].roof.sens

            # Total Electrical & Gas power in MW
            self.ElecTotal = \
                self.ElecTotal + BEM[j].fl_area * BEM[j].building.ElecTotal / 1.e6
            self.GasTotal = \
                self.GasTotal + BEM[j].fl_area * BEM[j].building.GasTotal / 1.e6

        # ------------------------------------------------------------------------------
        # Sensible Heat
        # Note 1: In the current uwg code, latent heat from evapotranspiration, stagnant
        # water, or anthropogenic sources is not modeled due to the difficulty of
        # validation, and lack of reliability of precipitation data from EPW files.
        # Note 2: Changed sensHeat to multiply treeSensHeat by treeCoverage fraction.
        # Since treeSensHeat calculated in solarcalcs does not multiply sensheat by tree
        # fraction. For example, if there were 0 trees, this value should be 0.
        # ------------------------------------------------------------------------------
        self.sensHeat = (
            self.Q_wall + self.Q_road + self.Q_vent + self.Q_window + self.Q_hvac +
            self.Q_traffic + self.treeSensHeat + self.Q_roof)

        # Error checking
        if self.canTemp > 350. or self.canTemp < 250:
            raise Exception('Got canyon temperature at {} C. Something obviously went '
                            'wrong (UCMDef.py).'.format(self.canTemp - 273.15))

    def __repr__(self):
        return 'UCMDef,\n verToHor: {}\n bldDensity: {}\n bldHeight: {}\n canWidth: ' \
            '{}\n canAspect: {}\n facArea: {}\n roofArea: {}'.format(
                self.verToHor, self.bldDensity, self.bldHeight, self.canWidth,
                self.canAspect, self.facArea, self.roofArea)
