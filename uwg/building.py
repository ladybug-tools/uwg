"""Class of specified building characteristics."""
from __future__ import division

from .psychrometrics import psychrometrics, moist_air_density
from .utilities import is_near_zero


class Building(object):
    """Specified building characteristics.
    # TODO: fix docstring
    properties
        % Building parameters
        floorHeight         % floor height (m)
        intHeat;            % timestep internal heat gains (W m-2 bld) (sensible only)
        intHeatNight;       % nighttime internal heat gains (W m-2 floor)
        intHeatDay;         % daytime internal heat gains (W m-2 floor)
        intHeatFRad;        % radiant fraction of internal gains
        intHeatFLat;        % latent fraction of internal gains
        infil;              % Infiltration (ACH)
        vent;               % Ventilation (ACH)
        glazingRatio;       % glazing ratio
        uValue;             % window U-value (W m-2 K-1) (including film coeff)
        shgc;               % window SHGC
        condType;           % cooling condensation system type {'AIR', 'WATER'}
        cop;                % COP of the cooling system (nominal)
        coolSetpointDay;    % daytime indoor cooling set-point (K)
        coolSetpointNight;  % nighttime indoor cooling set-point (K)
        heatSetpointDay;    % daytime indoor heating set-point (K)
        heatSetpointNight;  % nighttime indoor heating set-point (K)
        coolCap;            % rated cooling system capacity (W m-2)
        heatCap;            % rated heating system capacity (W m-2)
        heatEff;            % heating system efficiency (-)
        canyon_fraction     # fraction of waste heat released to canyon, default = 1
        mSys;               % HVAC supply mass flowrate (kg s-1 m-2)
        indoorTemp;         % indoor air temperature (K)
        indoorHum;          % indoor specific humidity (kg / kg)
        Twb;                % wetbulb temperature
        Tdp;                % dew point
        indoorRhum;         % indoor relative humidity

        area_floor;         % total floor space of the BEM
        FanMax;             % max fan flow rate (m^3/s) per DOE
        nFloor;             % number of floors
        RadFOcc;            % Radiant fraction of occupant
        LatFOcc;            % Latent fraction of occupant
        RadFEquip;          % Radiant fraction of equipment
        RadFLight;          % Radiant fraction of light

        Type;               % DOE reference building type
        Era;                % PRE80, PST80, NEW
        Zone;               % Climate zone number

        % Calculated values
        sensCoolDemand;     % building sensible cooling demand (W m-2)
        sensHeatDemand;     % building sensible heating demand (W m-2)
        copAdj;             % adjusted COP per temperature
        dehumDemand;        % dehumidification energy (W m-2)
        coolConsump;        % cooling energy consumption (W m-2)
        heatConsump;        % heating energy consumption (W m-2)
        sensWaste;          % sensible waste heat (W m-2)
        latWaste;           % lat waste heat (W m-2)
        fluxMass;           % mass surface heat flux (W m-2) (mass to indoor air)
        fluxWall;           % wall surface heat flux (W m-2) (wall to inside)
        fluxRoof;           % roof surface heat flux (W m-2) (roof to inside)
        fluxSolar;          % solar heat gain (W m-2) through window (SHGC)
        fluxWindow;         % heat gain/loss from window (U-value)
        fluxInterior;       % internal heat gain adjusted for latent/LW heat (W m-2)
        fluxInfil;          % heat flux from infiltration (W m-2)
        fluxVent;           % heat flux from ventilation (W m-2)
        ElecTotal;          % total electricity consumption - (W/m^2) of floor
        GasTotal;           % total gas consumption - (W/m^2) of floor
        Qhvac;              % total heat removed (sensible + latent)
        Qheat;              % total heat added (sensible only)
    """

    TEMP_COEF_CONFLICT_MSG = 'FATAL ERROR! Try reducing the simulation timstep to fix ' \
        'this error.'

    def __init__(self, floorHeight, intHeatNight, intHeatDay, intHeatFRad, intHeatFLat,
                 infil, vent, glazingRatio, uValue, shgc, condType, cop, coolSetpointDay,
                 coolSetpointNight, heatSetpointDay, heatSetpointNight, coolCap, heatEff,
                 initialTemp):

        self.floorHeight = float(floorHeight)        # floor height
        self.intHeat = intHeatNight                 # timetep internal gains (W m-2 bld) (sensible only)
        self.intHeatNight = intHeatNight            # nighttime internal heat gains  (W m-2 floor)
        self.intHeatDay = intHeatDay                # daytime internal heat gains  (W m-2 floor)
        self.intHeatFRad = intHeatFRad              # internal gain radiant fraction
        self.intHeatFLat = intHeatFLat              # internal gain latent fraction
        self.infil = infil                          # Infiltration (ACH)
        self.vent = vent                            # Ventilation (ACH)
        self.glazingRatio = glazingRatio            # glazing ratio
        self.uValue = uValue                        # window U-value ( w m-2 K-1) including film coeff
        self.shgc = shgc                            # window SHGC
        self.condType = condType                    # cooling condensation system type: AIR, WATER
        self.cop = cop                              # COP of cooling system (nominal)
        self.coolSetpointDay = coolSetpointDay      # daytime indoor cooling setpoint [K]
        self.coolSetpointNight = coolSetpointNight  # nighttime indoor heating setpoint [K]
        self.heatSetpointDay = heatSetpointDay      # daytimge indoor heating setpoint [K]
        self.heatSetpointNight = heatSetpointNight  # nighttime indoor heating setpoint [K]
        self.coolCap = coolCap                      # rated cooling system capacity (W m-2)
        self.heatEff = heatEff                      # heating system capacity (-)
        self.mSys = coolCap/1004./(min(coolSetpointDay,coolSetpointNight)-14-273.15) # HVAC supply mass flowrate (kg s-1 m-2)
        self.initialTemp = initialTemp
        self.indoorTemp = initialTemp               # Indoor Air Temperature [K]
        self.indoorHum = 0.012                      # Indoor specific humidity [kgv/kga]
        self.heatCap = 999                          # Default heat capacity value
        self.copAdj = cop                           # adjusted COP per temperature
        self.canyon_fraction = 1.0                  # Default canyon fraction

    @classmethod
    def from_dict(cls, data):
        """Create a Building object from a dictionary.

        Args:
            data: A Building dictionary following the format below.

        .. code-block:: python
            {
            "type": "Building"
            "floorHeight": self.floorHeight,
            "intHeatNight": self.intHeatNight,
            "intHeatDay": self.intHeatDay,
            "intHeatFRad": self.intHeatFRad,
            "intHeatFLat": self.intHeatFlat,
            "infil": self.infil,
            "vent": self.vent,
            "glazingRatio": self.glazingRatio,
            "uValue": self.uValue,
            "shgc": self.shgc,
            "condType": self.condType,
            "cop": self.cop,
            "coolSetpointDay": self.coolSetpointDay,
            "coolSetpointNight": self.coolSetpointNight,
            "heatSetpointDay": self.heatSetpointDay,
            "heatSetpointNight": self.heatSetpointNight,
            "coolCap": self.coolCap,
            "heatEff": self.heatEff,
            "initialTemp": self.initialTemp
            }
        """
        assert data['type'] == 'Building', 'Expected ' \
            'Building dictionary. Got {}.'.format(data['type'])

        return cls(data['floorHeight'], data['intHeatNight'], data['intHeatDay'],
                   data['intHeatFRad'], data['intHeatFLat'], data['infil'], data['vent'],
                   data['glazingRatio'], data['uValue'], data['shgc'], data['condType'],
                   data['cop'], data['coolSetpointDay'], data['coolSetpointNight'],
                   data['heatSetpointDay'], data['heatSetpointNight'], data['coolCap'],
                   data['heatEff'], data['initialTemp'])

    def to_dict(self):
        """Building dictionary representation."""
        base = {'type': 'Building'}
        base['floorHeight'] = self.floorHeight
        base['intHeatNight'] = self.intHeatNight
        base['intHeatDay'] = self.intHeatDay
        base['intHeatFRad'] = self.intHeatFRad
        base['intHeatFLat'] = self.intHeatFLat
        base['infil'] = self.infil
        base['vent'] = self.vent
        base['glazingRatio'] = self.glazingRatio
        base['uValue'] = self.uValue
        base['shgc'] = self.shgc
        base['condType'] = self.condType
        base['cop'] = self.cop
        base['coolSetpointDay'] = self.coolSetpointDay
        base['coolSetpointNight'] = self.coolSetpointNight
        base['heatSetpointDay'] = self.heatSetpointDay
        base['heatSetpointNight'] = self.heatSetpointNight
        base['coolCap'] = self.coolCap
        base['heatEff'] = self.heatEff
        base['initialTemp'] = self.initialTemp
        return base

    def BEMCalc(self, UCM, BEM, forc, parameter, simTime):
        """Update BEM by a single timestep."""

        self.ElecTotal = 0.0                            # total electricity consumption - (W/m^2) of floor
        self.nFloor = max(UCM.bldHeight/float(self.floorHeight),1)   # At least one floor
        self.Qheat = 0.0                                # total sensible heat added
        self.sensCoolDemand = 0.0                       # building sensible cooling demand (W m-2)
        self.sensHeatDemand = 0.0                       # building sensible heating demand (W m-2)
        self.coolConsump  = 0.0                         # cooling energy consumption (W m-2)
        self.heatConsump  = 0.0                         # heating energy consumption (W m-2)
        self.sensWaste = 0.0                            # Sensible waste heat (W m-2)
        self.dehumDemand  = 0.0                         # dehumidification energy (W m-2)
        self.Qhvac = 0.0                                # Total heat removed (sensible + latent)
        Qdehum = 0.0
        dens =  moist_air_density(forc.pres,self.indoorTemp,self.indoorHum)# [kgv/ m-3] Moist air density given dry bulb temperature, humidity ratio, and pressure
        evapEff = 1.                                    # evaporation efficiency in the condenser
        volVent = self.vent * self.nFloor               # total vent volumetric flow [m3 s-1]
        volInfil = self.infil * UCM.bldHeight / 3600.   # Change of units AC/H -> [m3 s-1]
        T_wall = BEM.wall.layerTemp[-1]                 # Inner layer
        volSWH = BEM.swh * self.nFloor / 3600.            # Change of units l/hr per m^2 -> [L/s]
        T_ceil = BEM.roof.layerTemp[-1]                 # Inner layer
        T_mass = BEM.mass.layerTemp[0]                  # Outer layer
        T_indoor = self.indoorTemp                      # Indoor temp (initial)
        T_can = UCM.canTemp                             # Canyon temperature

        # Normalize areas to building foot print [m^2/m^2(bld)]
        facArea = UCM.verToHor/UCM.bldDensity           # [m2(facade)/m2(bld)]
        wallArea = facArea*(1.-self.glazingRatio)       # [m2(wall)/m2(bld)]
        winArea = facArea*self.glazingRatio             # [m2(window)/m2(bld)]
        massArea = 2*self.nFloor-1                      # ceiling/floor (top & bottom)

        # Set temperature set points according to night/day setpoints in building schedule & simTime hr
        isEqualNightStart = is_near_zero((simTime.secDay/3600.) - parameter.nightSetStart)
        if simTime.secDay/3600. < parameter.nightSetEnd or (simTime.secDay/3600. > parameter.nightSetStart or isEqualNightStart):

            T_cool = self.coolSetpointNight
            T_heat = self.heatSetpointNight
            self.intHeat = self.intHeatNight * self.nFloor
        else:
            T_cool = self.coolSetpointDay
            T_heat = self.heatSetpointDay
            self.intHeat = self.intHeatDay*self.nFloor

        # Indoor convection heat transfer coefficients
        zac_in_wall = 3.076                             # wall heat convection coefficeint
        zac_in_mass = 3.076                             # mass heat convection coefficeint

        # Check that T_ceil and T_indoor within reasonable bounds
        converge_hi = 100.0 + 273.15
        converge_lo = -50.0 + 273.15

        try:
            chk_tin = converge_lo <= T_indoor <= converge_hi
            chk_tce = converge_lo <= T_ceil <= converge_hi

            if chk_tin is not True or chk_tce is not True:
                raise Exception("{}.\n Error at {}/{} {}s for bld {}.".format(self.TEMP_COEF_CONFLICT_MSG, simTime.month, simTime.day, simTime.secDay, BEM))

        except ValueError:
            raise Exception("{}.\n Error at {}/{} {}s for bld {}".format(self.TEMP_COEF_CONFLICT_MSG, simTime.month, simTime.day, simTime.secDay, BEM))

        # If temperature is reasonable assign coefficients
        if T_ceil > T_indoor:                               # set higher ceiling heat convection coefficient
            zac_in_ceil  = 0.948                            # - based on heat is higher on ceiling
        else:
            zac_in_ceil  = 4.040

        # -------------------------------------------------------------
        # Heat fluxes (per m^2 of bld footprint)
        # -------------------------------------------------------------
        # Solar Heat Gain: solar radiation received (W m-2) * area * SHGC
        winTrans = (BEM.wall.solRec * self.shgc * winArea)

        # Latent heat infiltration & ventilation (W/m^2 of bld footprint) from
        # volInfil/volVent: [m3 s-1 m-2 (bld/facade#)]
        # parameter.lv: latent heat of evaporation [J kgv-1]
        # dens: kga m-3
        # UCM.canHum: canyon specific humidity (kgv kga-1)
        # indoorHum: indoor kv kga-1
        # QL = W m-2

        QLinfil = volInfil * dens * parameter.lv * (UCM.canHum - self.indoorHum)
        QLvent = volVent * dens * parameter.lv * (UCM.canHum - self.indoorHum)
        QLintload = self.intHeat * self.intHeatFLat # Qlatent Internal load = timestep internal gain * internal gain latent fraction

        # Heat/Cooling load (W/m^2 of bld footprint), if any
        self.sensCoolDemand = max(
            wallArea*zac_in_wall*(T_wall - T_cool) +            # wall load
            massArea*zac_in_mass*(T_mass-T_cool) +              # mass load
            winArea*self.uValue*(T_can-T_cool) +                # window load due to temp delta
            zac_in_ceil *(T_ceil-T_cool) +                      # ceiling load
            self.intHeat +                                      # internal load
            volInfil*dens*parameter.cp*(T_can-T_cool) +         # infiltration load (volInfil = m3 s-1)
            volVent*dens*parameter.cp*(T_can-T_cool) +          # ventilation load (volVent = m3 s-1)
            winTrans,                                           # solar load through window
            0.)

        self.sensHeatDemand = max(
            -(wallArea*zac_in_wall*(T_wall-T_heat) +            # wall load
            massArea*zac_in_mass*(T_mass-T_heat) +              # mass load
            winArea*self.uValue*(T_can-T_heat) +                # window load due to temp delta
            zac_in_ceil*(T_ceil-T_heat) +                       # ceiling load
            self.intHeat +                                      # internal load
            volInfil*dens*parameter.cp*(T_can-T_heat) +         # infiltration load (volInfil = m3 s-1)
            volVent*dens*parameter.cp*(T_can-T_heat) +          # ventilation load (volVent = m3 s-1)
            winTrans),                                          # solar load through window
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
            VolCool = self.sensCoolDemand / (dens*parameter.cp*(T_indoor-283.15)) # m3
            self.dehumDemand = max(VolCool * dens * (self.indoorHum - 0.9*0.0078)*parameter.lv, 0.)
            if (self.dehumDemand + self.sensCoolDemand) > (self.coolCap * self.nFloor): # if cooling demand greater then hvac cooling capacity
                self.Qhvac = self.coolCap * self.nFloor
                VolCool = VolCool / (self.dehumDemand + self.sensCoolDemand) * (self.coolCap * self.nFloor)
                self.sensCoolDemand = self.sensCoolDemand * (self.coolCap * self.nFloor) / (self.dehumDemand + self.sensCoolDemand)
                self.dehumDemand = self.dehumDemand * (self.coolCap * self.nFloor) / (self.dehumDemand + self.sensCoolDemand)
            else:
                self.Qhvac = self.dehumDemand + self.sensCoolDemand

            Qdehum = VolCool * dens * parameter.lv * (self.indoorHum - 0.9*0.0078)
            self.coolConsump = (max(self.sensCoolDemand+self.dehumDemand,0.0))/self.copAdj

            # Waste heat from HVAC (per m^2 building foot print)
            if (self.condType == 'AIR'):
                self.sensWaste = max(self.sensCoolDemand+self.dehumDemand,0)+self.coolConsump
                self.latWaste = 0.0
            elif (self.condType == 'WAT'): # Not sure if this works well
                self.sensWaste = max(self.sensCoolDemand+self.dehumDemand,0)+self.coolConsump*(1.-evapEff)
                self.latWaste = max(self.sensCoolDemand+self.dehumDemand,0)+self.coolConsump*evapEff

            self.sensHeatDemand = 0.

        # -------------------------------------------------------------
        # HVAC system (heating demand = W/m^2 bld footprint)
        # -------------------------------------------------------------
        elif self.sensHeatDemand > 0. and UCM.canTemp < 288.:
            # limit on heating capacity
            self.Qheat = min(self.sensHeatDemand, self.heatCap*self.nFloor)

            self.heatConsump  = self.Qheat / self.heatEff
            self.sensWaste = self.heatConsump - self.Qheat         # waste per footprint
            self.heatConsump = self.heatConsump/self.nFloor        # adjust to be per floor area
            self.sensHeatDemand = self.Qheat/self.nFloor           # adjust to be per floor area
            Qdehum = 0.0
            self.sensCoolDemand = 0.0


        # -------------------------------------------------------------
        # Evolution of the internal temperature and humidity
        # -------------------------------------------------------------
        # wall, mass, roof, intload, infil, vent, hvac, heat, window

        Q = self.intHeat + winTrans + self.Qheat - self.sensCoolDemand

        H1 = (T_wall*wallArea*zac_in_wall +
            T_mass*massArea*zac_in_mass +
            T_ceil*zac_in_ceil +
            T_can*winArea*self.uValue +
            T_can*volInfil * dens * parameter.cp +
            T_can*volVent * dens * parameter.cp)

        H2 = (wallArea*zac_in_wall +
            massArea*zac_in_mass +
            zac_in_ceil +
            winArea*self.uValue +
            volInfil * dens * parameter.cp +
            volVent * dens * parameter.cp)

        # Assumes air temperature of control volume is sum of surface boundary temperatures
        # weighted by area and heat transfer coefficient + generated heat
        self.indoorTemp = (H1 + Q)/H2

        self.indoorHum = self.indoorHum + (simTime.dt/(dens * parameter.lv * UCM.bldHeight)) * \
            (QLintload + QLinfil + QLvent - Qdehum)

        # Calculate relative humidity (Pw/Pws*100) using pressurce, indoor temperature, humidity
        _Tdb, _w, _phi, _h, _Tdp, _v = psychrometrics(self.indoorTemp, self.indoorHum, forc.pres)
        self.indoorRhum = _phi

        # These are used for element calculation (per m^2 of element area)
        self.fluxWall = zac_in_wall * (T_indoor - T_wall)
        self.fluxRoof = zac_in_ceil * (T_indoor - T_ceil)
        self.fluxMass = zac_in_mass * (T_indoor - T_mass) + self.intHeat * self.intHeatFRad/massArea

        # These are for record keeping only, per m^2 of floor area (W m-2)
        self.fluxSolar = winTrans/self.nFloor
        self.fluxWindow = winArea * self.uValue *(T_can - T_indoor)/self.nFloor
        self.fluxInterior = self.intHeat * self.intHeatFRad *(1.-self.intHeatFLat)/self.nFloor
        self.fluxInfil= volInfil * dens * parameter.cp *(T_can - T_indoor)/self.nFloor # volInfil = m3 s-1
        self.fluxVent = volVent * dens * parameter.cp *(T_can - T_indoor)/self.nFloor  # volVent = m3 s-1
        self.coolConsump = self.coolConsump/self.nFloor
        self.sensCoolDemand = self.sensCoolDemand/self.nFloor

        # Total Electricity/building floor area (W/m^2)
        self.ElecTotal = self.coolConsump + BEM.elec + BEM.light

        # Waste heat to canyon, W/m^2 of building + water
        CpH20 = 4200.           # heat capacity of water
        T_hot = 49 + 273.15     # Service water temp (assume no storage)
        self.sensWaste = self.sensWaste + (1/self.heatEff-1.) * (volSWH*CpH20*(T_hot - forc.waterTemp)) + BEM.gas * (1 - self.heatEff) * self.nFloor

        # Gas equip per floor + water usage per floor + heating/floor
        self.GasTotal = BEM.gas + volSWH * CpH20 * (T_hot - forc.waterTemp) / self.nFloor / self.heatEff + self.heatConsump

    def __repr__(self):
        return 'Building,\n floorHeight: {}\n glazingRatio: {}\n ' \
            'uValue: {}\n shgc: {}\n infil: {}\n vent: {}\n ' \
            'condType: {}\n'.format(
                self.floorHeight, self.glazingRatio, self.uValue, self.shgc,
                self.infil, self.vent, self.condType)

