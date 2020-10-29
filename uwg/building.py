"""Class of specified building characteristics."""
from __future__ import division

from .psychrometrics import psychrometrics, moist_air_density
from .utilities import is_near_zero, float_in_range, float_positive

try:
    str = basestring
except NameError:
    pass


class Building(object):
    """Specified building characteristics.

    Args:
        floor_height: Floor height in meters.
        int_heat_night: Nighttime internal sensible heat gain [W/m2].
        int_heat_day: Daytime internal sensible heat gain [W/m2].
        int_heat_frad: Value between 0 and 1 for radiant fraction of internal gains.
        int_heat_flat: Value between 0 and 1 for latent fraction of internal gains.
        infil: Infiltration rate (ACH).
        vent: Ventilation rate (ACH).
        glazing_ratio: Value between 0 and 1 for glazing ratio.
        u_value: Window U-value including film coefficent [W/(m2-K)].
        shgc: Value between 0 and 1 for window Solar Heat Gain Coefficient (SHGC).
        condtype: Text string for cooling condensation system type. Choose from AIR or
            WATER.
        cop: Coefficient of Performance (COP) of cooling system (nominal).
        coolcap: Rated cooling system capacity [W/m2].
        heateff: Heating system capacity.
        initial_temp: Initial indoor air temperature [K].

    Properties:
        * floor_height -- floor height (m)
        * int_heat -- timestep internal heat gains (W m-2 bld) (sensible only)
        * int_heat_night -- nighttime internal heat gains (W m-2 floor)
        * int_heat_day -- daytime internal heat gains (W m-2 floor)
        * int_heat_frad -- radiant fraction of internal gains
        * int_heat_flat -- latent fraction of internal gains
        * infil -- Infiltration (ACH)
        * vent -- Ventilation (m^3/s/m^2)
        * glazing_ratio -- glazing ratio
        * u_value -- window U-value (W m-2 K-1) (including film coeff)
        * shgc -- window SHGC
        * condtype -- cooling condensation system type {'AIR', 'WATER'}
        * cop -- COP of the cooling system (nominal)
        * cop_adj -- adjusted COP per temperature
        * cool_setpoint_day -- daytime indoor cooling set-point (K)
        * cool_setpoint_night -- nighttime indoor cooling set-point (K)
        * heat_setpoint_day -- daytime indoor heating set-point (K)
        * heat_setpoint_night -- nighttime indoor heating set-point (K)
        * heat_cap -- rated heating system capacity (W m-2)
        * heateff -- heating system efficiency (-)
        * canyon_fraction -- fraction of waste heat released to canyon, default = 1
        * msys -- HVAC supply mass flowrate (kg s-1 m-2)
        * indoor_temp -- indoor air temperature (K)
        * indoor_hum -- indoor specific humidity (kg / kg)
        * Twb -- wetbulb temperature
        * Tdp -- dew point
        * indoorRhum -- indoor relative humidity
        * area_floor -- total floor space of the BEM
        * FanMax -- max fan flow rate (m^3/s) per DOE
        * nFloor -- number of floors
        * RadFOcc -- Radiant fraction of occupant
        * LatFOcc -- Latent fraction of occupant
        * RadFEquip -- Radiant fraction of equipment
        * RadFLight -- Radiant fraction of light
        * sensCoolDemand -- building sensible cooling demand (W m-2)
        * sensHeatDemand -- building sensible heating demand (W m-2)
        * copAdj -- adjusted COP per temperature
        * dehumDemand -- dehumidification energy (W m-2)
        * coolConsump -- cooling energy consumption (W m-2)
        * heatConsump -- heating energy consumption (W m-2)
        * sensWaste -- sensible waste heat (W m-2)
        * latWaste -- lat waste heat (W m-2)
        * fluxMass -- mass surface heat flux (W m-2) (mass to indoor air)
        * fluxWall -- wall surface heat flux (W m-2) (wall to inside)
        * fluxRoof -- roof surface heat flux (W m-2) (roof to inside)
        * fluxSolar -- solar heat gain (W m-2) through window (SHGC)
        * fluxWindow -- heat gain/loss from window (U-value)
        * fluxInterior -- internal heat gain adjusted for latent/LW heat (W m-2)
        * fluxInfil -- heat flux from infiltration (W m-2)
        * fluxVent -- heat flux from ventilation (W m-2)
        * ElecTotal -- total electricity consumption - (W/m^2) of floor
        * GasTotal -- total gas consumption - (W/m^2) of floor
        * Qhvac -- total heat removed (sensible + latent)
        * Qheat -- total heat added (sensible only)
    """

    TEMP_COEF_CONFLICT_MSG = \
        'FATAL ERROR! Try reducing the simulation timstep to fix this error.'

    def __init__(self, floor_height, int_heat_night, int_heat_day, int_heat_frad,
                 int_heat_flat, infil, vent, glazing_ratio, u_value, shgc, condtype, cop,
                 coolcap, heateff, initial_temp):

        self.floor_height = floor_height
        self.int_heat = int_heat_night
        self.int_heat_night = int_heat_night
        self.int_heat_day = int_heat_day
        self.int_heat_frad = int_heat_frad
        self.int_heat_flat = int_heat_flat
        self.infil = infil
        self.vent = vent
        self.glazing_ratio = glazing_ratio
        self.u_value = u_value
        self.shgc = shgc
        self.condtype = condtype
        self.cop = cop
        self.coolcap = coolcap
        self.heateff = heateff
        self.initial_temp = initial_temp
        self.indoor_temp = initial_temp
        self.cop_adj = cop

        # set properties
        self.indoor_hum = 0.012
        self.heat_cap = 999
        self.canyon_fraction = 1.0
        self.cool_setpoint_day = 297  # 24 C
        self.cool_setpoint_night = 297  # 24 C
        self.heat_setpoint_day = 293  # 24 C
        self.heat_setpoint_night = 293  # 24 C
        self.msys = \
            coolcap / 1004. / (min(self.cool_setpoint_day, self.cool_setpoint_night) -
                               14 - 273.15)

    @property
    def floor_height(self):
        """Get or set floor height in meters."""
        return self._floor_height

    @floor_height.setter
    def floor_height(self, value):
        self._floor_height = float_positive(value, 'floor_height')

    @property
    def int_heat_night(self):
        """Get or set nighttime internal sensible heat gain [W/m2]."""
        return self._int_heat_night

    @int_heat_night.setter
    def int_heat_night(self, value):
        self._int_heat_night = float_positive(value, 'int_heat_night')

    @property
    def int_heat_frad(self):
        """Get or set value between 0 and 1 radiant fraction of internal gains."""
        return self._int_heat_frad

    @int_heat_frad.setter
    def int_heat_frad(self, value):
        self._int_heat_frad = float_in_range(value, 0, 1, 'int_heat_frad')

    @property
    def int_heat_flat(self):
        """Get or set value between 0 and 1 for latent fraction of internal gains."""
        return self._int_heat_flat

    @int_heat_flat.setter
    def int_heat_flat(self, value):
        self._int_heat_flat = float_in_range(value, 0, 1, 'int_heat_flat')

    @property
    def infil(self):
        """Get or set infiltration rate (ACH)."""
        return self._infil

    @infil.setter
    def infil(self, value):
        self._infil = float_positive(value, 'infil')

    @property
    def vent(self):
        """Get or set ventilation rate (ACH)."""
        return self._vent

    @vent.setter
    def vent(self, value):
        self._vent = float_positive(value, 'vent')

    @property
    def glazing_ratio(self):
        """Get or set value between 0 and 1 for glazing ratio."""
        return self._glazing_ratio

    @glazing_ratio.setter
    def glazing_ratio(self, value):
        self._glazing_ratio = float_in_range(value, 0, 1, 'glazing_ratio')

    @property
    def u_value(self):
        """Get or set window U-value including film coefficent [W/(m2-K)]."""
        return self._u_value

    @u_value.setter
    def u_value(self, value):
        self._u_value = float_positive(value, 'u_value')

    @property
    def shgc(self):
        """Get or set building glazing Solar Heat Gain Coefficient."""
        return self._shgc

    @shgc.setter
    def shgc(self, value):
        self._shgc = float_in_range(value, 0, 1, 'shgc')

    @property
    def condtype(self):
        """Get or set text string for cooling condensation system type.

        Choose from:

        * "AIR"
        * "WATER".
        """
        return self._condtype

    @condtype.setter
    def condtype(self, value):
        value = value.upper()
        assert value in ('AIR', 'WATER'), 'condtype must be "AIR" or "WATER". ' \
            'Got: {}.'.format(value)
        self._condtype = value

    @property
    def cop(self):
        """Get or set the nominal Coefficient of Performance (COP) of cooling system."""
        return self._cop

    @cop.setter
    def cop(self, value):
        self._cop = float_positive(value, 'cop')

    @property
    def coolcap(self):
        """Get or set rated cooling system capacity [W/m2]."""
        return self._coolcap

    @coolcap.setter
    def coolcap(self, value):
        self._coolcap = float_positive(value, 'coolcap')

    @property
    def heateff(self):
        """Get or set heating system capacity."""
        return self._heateff

    @heateff.setter
    def heateff(self, value):
        self._heateff = float_positive(value, 'heateff')

    @property
    def initial_temp(self):
        """Get or set initial indoor air temperature [K]."""
        return self._initial_temp

    @initial_temp.setter
    def initial_temp(self, value):
        self._initial_temp = float_positive(value, 'initial_temp')

    @classmethod
    def from_dict(cls, data):
        """Create a Building object from a dictionary.

        Args:
            data: A Building dictionary following the format below.

        .. code-block:: python
            {
            "type": "Building"
            "floor_height": self.floor_height,
            "int_heat_night": self.int_heat_night,
            "int_heat_day": self.int_heat_day,
            "int_heat_frad": self.int_heat_frad,
            "int_heat_flat": self.int_heat_flat,
            "infil": self.infil,
            "vent": self.vent,
            "glazing_ratio": self.glazing_ratio,
            "u_value": self.u_value,
            "shgc": self.shgc,
            "condtype": self.condtype,
            "cop": self.cop,
            "coolcap": self.coolcap,
            "heateff": self.heateff,
            "initial_temp": self.initial_temp
            }
        """
        assert data['type'] == 'Building', 'Expected ' \
            'Building dictionary. Got {}.'.format(data['type'])

        return cls(data['floor_height'], data['int_heat_night'], data['int_heat_day'],
                   data['int_heat_frad'], data['int_heat_flat'], data['infil'],
                   data['vent'], data['glazing_ratio'], data['u_value'], data['shgc'],
                   data['condtype'], data['cop'], data['coolcap'], data['heateff'],
                   data['initial_temp'])

    def to_dict(self):
        """Building dictionary representation."""
        base = {'type': 'Building'}
        base['floor_height'] = self.floor_height
        base['int_heat_night'] = self.int_heat_night
        base['int_heat_day'] = self.int_heat_day
        base['int_heat_frad'] = self.int_heat_frad
        base['int_heat_flat'] = self.int_heat_flat
        base['infil'] = self.infil
        base['vent'] = self.vent
        base['glazing_ratio'] = self.glazing_ratio
        base['u_value'] = self.u_value
        base['shgc'] = self.shgc
        base['condtype'] = self.condtype
        base['cop'] = self.cop
        base['coolcap'] = self.coolcap
        base['heateff'] = self.heateff
        base['initial_temp'] = self.initial_temp
        return base

    def BEMCalc(self, UCM, BEM, forc, parameter, simTime):
        """Update BEM by a single timestep."""

        # total electricity consumption - (W/m^2) of floor
        self.ElecTotal = 0.0
        self.nFloor = max(UCM.bldHeight / float(self.floor_height), 1)
        self.Qheat = 0.0  # total sensible heat added
        self.sensCoolDemand = 0.0  # building sensible cooling demand (W m-2)
        self.sensHeatDemand = 0.0  # building sensible heating demand (W m-2)
        self.coolConsump = 0.0  # cooling energy consumption (W m-2)
        self.heatConsump = 0.0  # heating energy consumption (W m-2)
        self.sensWaste = 0.0  # Sensible waste heat (W m-2)
        self.dehumDemand = 0.0  # dehumidification energy (W m-2)
        self.Qhvac = 0.0  # Total heat removed (sensible + latent)
        Qdehum = 0.0
        # dens: Moist air density given dry bulb temp, humidity ratio, and pressure
        dens = moist_air_density(
            forc.pres, self.indoor_temp, self.indoor_hum)  # kgv/m3
        evapEff = 1.  # evaporation efficiency in the condenser
        # total vent volumetric flow [m3 s-1]
        volVent = self.vent * self.nFloor
        volInfil = self.infil * UCM.bldHeight / \
            3600.  # Change of units AC/H -> [m3 s-1]
        T_wall = BEM.wall.layerTemp[-1]  # Inner layer
        volSWH = BEM.swh * self.nFloor / \
            3600.  # Change of units l/hr per m^2 -> [L/s]
        T_ceil = BEM.roof.layerTemp[-1]  # Inner layer
        T_mass = BEM.mass.layerTemp[0]  # Outer layer
        T_indoor = self.indoor_temp  # Indoor temp (initial)
        T_can = UCM.canTemp  # Canyon temperature

        # Normalize areas to building foot print [m^2/m^2(bld)]
        facArea = UCM.verToHor / UCM.bldDensity  # [m2(facade)/m2(bld)]
        wallArea = facArea * (1. - self.glazing_ratio)  # [m2(wall)/m2(bld)]
        winArea = facArea * self.glazing_ratio  # [m2(window)/m2(bld)]
        massArea = 2 * self.nFloor - 1  # ceiling/floor (top & bottom)

        # Set temperature set points according to night/day setpoints in building
        # schedule & simTime hr
        isEqualNightStart = \
            is_near_zero((simTime.secDay / 3600.) - parameter.nightSetStart)
        if (simTime.secDay / 3600. < parameter.nightSetEnd) or \
           (simTime.secDay / 3600. > parameter.nightSetStart or isEqualNightStart):
            T_cool = self.cool_setpoint_night
            T_heat = self.heat_setpoint_night
            self.int_heat = self.int_heat_night * self.nFloor
        else:
            T_cool = self.cool_setpoint_day
            T_heat = self.heat_setpoint_day
            self.int_heat = self.int_heat_day * self.nFloor

        # Indoor convection heat transfer coefficients
        zac_in_wall = 3.076  # wall heat convection coefficeint
        zac_in_mass = 3.076  # mass heat convection coefficeint

        # Check that T_ceil and T_indoor within reasonable bounds
        converge_hi = 100.0 + 273.15
        converge_lo = -50.0 + 273.15

        try:
            chk_tin = converge_lo <= T_indoor <= converge_hi
            chk_tce = converge_lo <= T_ceil <= converge_hi

            if chk_tin is not True or chk_tce is not True:
                raise Exception("{}.\n Error at {}/{} {}s for bld {}.".format(
                    self.TEMP_COEF_CONFLICT_MSG, simTime.month, simTime.day,
                    simTime.secDay, BEM))

        except ValueError:
            raise Exception("{}.\n Error at {}/{} {}s for bld {}".format(
                self.TEMP_COEF_CONFLICT_MSG, simTime.month, simTime.day,
                simTime.secDay, BEM))

        # If temperature is reasonable assign coefficients
        if T_ceil > T_indoor:
            # set higher ceiling heat convection coefficient
            # based on heat is higher on ceiling
            zac_in_ceil = 0.948
        else:
            zac_in_ceil = 4.040

        # -------------------------------------------------------------
        # Heat fluxes (per m^2 of bld footprint)
        # -------------------------------------------------------------
        # Solar Heat Gain: solar radiation received (W m-2) * area * SHGC
        winTrans = BEM.wall.solRec * self.shgc * winArea

        # Latent heat infiltration & ventilation (W/m^2 of bld footprint) from
        # volInfil/volVent: [m3 s-1 m-2 (bld/facade#)]
        # parameter.lv: latent heat of evaporation [J kgv-1]
        # dens: kga m-3
        # UCM.canHum: canyon specific humidity (kgv kga-1)
        # indoorHum: indoor kv kga-1
        # QL = W m-2

        QLinfil = volInfil * dens * parameter.lv * \
            (UCM.canHum - self.indoor_hum)
        QLvent = volVent * dens * parameter.lv * (UCM.canHum - self.indoor_hum)
        # QLintload (Qlatent Internal load): timestep int gain * int gain latent frac
        QLintload = self.int_heat * self.int_heat_flat

        # Heat/Cooling load (W/m^2 of bld footprint), if any
        self.sensCoolDemand = max(
            wallArea * zac_in_wall * (T_wall - T_cool) +  # wall load
            massArea * zac_in_mass * (T_mass-T_cool) +  # mass load
            # window load due to temp delta
            winArea * self.u_value * (T_can - T_cool) +
            zac_in_ceil * (T_ceil-T_cool) +  # ceiling load
            self.int_heat +  # internal load
            volInfil * dens * parameter.cp * (T_can - T_cool) +  # m3 s-1
            volVent * dens * parameter.cp * (T_can-T_cool) +  # m3 s-1
            winTrans,  # solar load through window
            0.)

        self.sensHeatDemand = max(
            -(wallArea * zac_in_wall * (T_wall-T_heat) +  # wall load
              massArea * zac_in_mass * (T_mass-T_heat) +  # mass load
              winArea * self.u_value * (T_can - T_heat) +  # window load
              zac_in_ceil * (T_ceil-T_heat) +  # ceiling load
              self.int_heat +  # internal load
              volInfil * dens * parameter.cp * (T_can-T_heat) +  # m3 s-1
              volVent * dens * parameter.cp * (T_can-T_heat) +  # m3 s-1
              winTrans),  # solar load through window
            0.)

        # -------------------------------------------------------------
        # HVAC system (cooling demand = W/m^2 bld footprint)
        # -------------------------------------------------------------
        if self.sensCoolDemand > 0. and UCM.canTemp > 288.:
            # Cooling energy is the equivalent energy to bring a vol
            # where sensCoolDemand = dens * Cp * x * (T_indoor - 10C) &
            # given 7.8g/kg of air at 10C, assume 7g/kg of air
            # dehumDemand = x * dens * (self.indoorHum -
            # 0.9*0.0078)*parameter.lv
            VolCool = self.sensCoolDemand / \
                (dens * parameter.cp * (T_indoor - 283.15))
            self.dehumDemand = \
                max(VolCool * dens * (self.indoor_hum -
                                      0.9 * 0.0078) * parameter.lv, 0.)
            if (self.dehumDemand + self.sensCoolDemand) > (self.coolcap * self.nFloor):
                # if cooling demand greater then hvac cooling capacity
                self.Qhvac = self.coolcap * self.nFloor
                VolCool = (
                    VolCool / (self.dehumDemand + self.sensCoolDemand) *
                    (self.coolcap * self.nFloor))
                self.sensCoolDemand = (
                    self.sensCoolDemand * (self.coolcap * self.nFloor) /
                    (self.dehumDemand + self.sensCoolDemand))
                self.dehumDemand = (
                    self.dehumDemand * (self.coolcap * self.nFloor) /
                    (self.dehumDemand + self.sensCoolDemand))
            else:
                self.Qhvac = self.dehumDemand + self.sensCoolDemand

            Qdehum = VolCool * dens * parameter.lv * \
                (self.indoor_hum - 0.9 * 0.0078)
            self.coolConsump = \
                (max(self.sensCoolDemand + self.dehumDemand, 0.0)) / self.cop_adj

            # Waste heat from HVAC (per m^2 building foot print)
            if self.condtype == 'AIR':
                self.sensWaste = \
                    max(self.sensCoolDemand + self.dehumDemand, 0) + \
                    self.coolConsump
                self.latWaste = 0.0
            elif self.condtype == 'WATER':  # Not sure if this works well
                self.sensWaste = (
                    max(self.sensCoolDemand + self.dehumDemand, 0) + self.coolConsump *
                    (1. - evapEff))
                self.latWaste = (
                    max(self.sensCoolDemand + self.dehumDemand, 0) + self.coolConsump *
                    evapEff)
            self.sensHeatDemand = 0.

        # -------------------------------------------------------------
        # HVAC system (heating demand = W/m^2 bld footprint)
        # -------------------------------------------------------------
        elif self.sensHeatDemand > 0. and UCM.canTemp < 288.:
            # limit on heating capacity
            self.Qheat = min(self.sensHeatDemand, self.heat_cap * self.nFloor)

            self.heatConsump = self.Qheat / self.heateff
            self.sensWaste = self.heatConsump - self.Qheat  # waste per footprint
            self.heatConsump = self.heatConsump / self.nFloor  # adjust to per flr area
            self.sensHeatDemand = self.Qheat / self.nFloor  # adjust to per flr area
            Qdehum = 0.0
            self.sensCoolDemand = 0.0

        # -------------------------------------------------------------
        # Evolution of the internal temperature and humidity
        # -------------------------------------------------------------
        # wall, mass, roof, intload, infil, vent, hvac, heat, window

        Q = self.int_heat + winTrans + self.Qheat - self.sensCoolDemand

        H1 = (T_wall * wallArea * zac_in_wall +
              T_mass * massArea * zac_in_mass +
              T_ceil * zac_in_ceil +
              T_can * winArea * self.u_value +
              T_can * volInfil * dens * parameter.cp +
              T_can * volVent * dens * parameter.cp)

        H2 = (wallArea * zac_in_wall +
              massArea * zac_in_mass +
              zac_in_ceil +
              winArea * self.u_value +
              volInfil * dens * parameter.cp +
              volVent * dens * parameter.cp)

        # Assumes air temp of control volume is sum of surface boundary temps
        # weighted by area and heat transfer coefficient + generated heat
        self.indoor_temp = (H1 + Q) / H2

        self.indoor_hum = (
            self.indoor_hum + (simTime.dt / (dens * parameter.lv * UCM.bldHeight)) *
            (QLintload + QLinfil + QLvent - Qdehum))

        # Calculate relative hum (Pw/Pws*100) using pressurce, indoor temp, hum
        _Tdb, _w, _phi, _h, _Tdp, _v = \
            psychrometrics(self.indoor_temp, self.indoor_hum, forc.pres)
        self.indoorRhum = _phi

        # These are used for element calculation (per m^2 of element area)
        self.fluxWall = zac_in_wall * (T_indoor - T_wall)
        self.fluxRoof = zac_in_ceil * (T_indoor - T_ceil)
        self.fluxMass = (
            zac_in_mass * (T_indoor - T_mass) + self.int_heat * self.int_heat_f_rad /
            massArea)

        # These are for record keeping only, per m^2 of floor area (W m-2)
        self.fluxSolar = winTrans / self.nFloor
        self.fluxWindow = winArea * self.u_value * \
            (T_can - T_indoor) / self.nFloor
        self.fluxInterior = \
            self.int_heat * self.int_heat_f_rad * \
            (1. - self.int_heat_flat) / self.nFloor
        # volInfil: m3/s
        self.fluxInfil = \
            volInfil * dens * parameter.cp * (T_can - T_indoor) / self.nFloor
        # volVent: m3/s
        self.fluxVent = volVent * dens * parameter.cp * \
            (T_can - T_indoor) / self.nFloor
        self.coolConsump = self.coolConsump / self.nFloor
        self.sensCoolDemand = self.sensCoolDemand / self.nFloor

        # Total Electricity/building floor area (W/m^2)
        self.ElecTotal = self.coolConsump + BEM.elec + BEM.light

        # Waste heat to canyon, W/m^2 of building + water
        CpH20 = 4200.  # specific heat capacity of water J/(kg-C)
        T_hot = 49 + 273.15  # Service water temp (assume no storage)
        # N.B No L to kg conversion for water in this equation because
        # 1 L of water = 1 kg of water.
        self.sensWaste = (
            self.sensWaste + (1 / self.heateff - 1.) *
            (volSWH * CpH20 * (T_hot - forc.waterTemp)) + BEM.gas *
            (1 - self.heateff) * self.nFloor)

        # Gas equip per floor + water usage per floor + heating/floor
        self.GasTotal = (
            BEM.gas + volSWH * CpH20 * (T_hot - forc.waterTemp) / self.nFloor /
            self.heateff + self.heatConsump)

    def __repr__(self):
        return 'Building,\n floor_height: {}\n int_heat_night {}\n int_heat_day {}\n ' \
            'int_heat_frad {}\n int_heat_flat {}\n infil: {}\n vent: {}\n ' \
            'glazing_ratio: {}\n u_value: {}\n shgc: {}\n condtype: {}\n cop: {}\n ' \
            'cool_setpoint_day: {}\n cool_setpoint_night: {}\n heat_setpoint_day: {}\n' \
            'heat_setpoint_night: {}\n coolcap: {}\n heateff: {}\n ' \
            'initial_temp: {}'.format(
                self.floor_height, self.int_heat_night, self.int_heat_day,
                self.int_heat_frad, self.int_heat_flat, self.infil, self.vent,
                self.glazing_ratio, self.u_value, self.shgc, self.condtype,
                self.cop, self.cool_setpoint_day, self.cool_setpoint_night,
                self.heat_setpoint_day, self.heat_setpoint_night,
                self.coolcap, self.heateff, self.initial_temp)
