
"""
Translated from: https://github.com/hansukyang/UWG_Matlab/blob/master/readDOE.m
Translated to Python by Chris Mackey (chris@mackeyarchitecture.com) and Saeran Vasanthakumar (saeranv@gmail.com) - August 2017
"""

class Building(object):
    """
    Building building class of specified building characteristics.

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

    def __init__(self,floorHeight,intHeatNight,intHeatDay,intHeatFRad,\
            intHeatFLat,infil,vent,glazingRatio,uValue,shgc,\
            condType,cop,coolSetpointDay,coolSetpointNight,\
            heatSetpointDay,heatSetpointNight,coolCap,heatEff,initialTemp):

            self.floorHeight =float(floorHeight)    # floor height
            self.intHeat = intHeatNight             # timetep internal gains
            self.intHeatNight = intHeatNight        # nighttime internal heat gains [W m-2]
            self.intHeatDay = intHeatDay            # daytime internal heat gains [W m-2]
            self.intHeatFRad = intHeatFRad          # internal gain radiant fraction
            self.intHeatFLat = intHeatFLat          # internal gain latent fraction
            self.infil = infil                      # Infiltration (ACH)
            self.vent = vent                        # Ventilation (ACH)
            self.glazingRatio = glazingRatio        # glazing ratio
            self.uValue = uValue                    # window U-value ( w m-2 K-1) including film coeff
            self.shgc = shgc                        # window SHGC
            self.condType = condType                # cooling condensation system type: AIR, WATER
            self.cop = cop                          # COP of cooling system (nominal)
            self.coolSetpointDay = coolSetpointDay  # daytime indoor cooling setpoint [K]
            self.coolSetpointNight = coolSetpointNight # nighttime indoor heating setpoint [K]
            self.heatSetpointDay = heatSetpointDay     # daytimge indoor heating setpoint [K]
            self.heatSetpointNight = heatSetpointNight # nighttime indoor heating setpoint [K]
            self.coolCap = coolCap                     # rated cooling system capacity (W m-2)
            self.heatEff = heatEff                     # heating system capacity (-)
            self.mSys = coolCap/1004./(min(coolSetpointDay,coolSetpointNight)-14-273.15) # HVAC supply mass flowrate (kg s-1 m-2)
            self.indoorTemp = initialTemp           # Indoor Air Temperature [K]
            self.indoorHum = 0.012                  # Indoor specific humidity [kgv/kga]
            self.heatCap = 999                      # Default heat capacity value
            self.copAdj = cop                       # adjusted COP per temperature
            self.canyon_fraction = 1.0              # Default canyon fraction

            self.Type = "null"                      # DOE reference building type
            self.Era = "null"                       # pre80, pst80, new
            self.Zone = "null"                      # Climate zone number

    def __repr__(self):
        return "Building: Type: {:s}, Era: {:s}, Zone: {:s}; @ Ti: {a}, WWR: {b}".format(
            self.Type,
            self.Era,
            self.Zone,
            a=self.indoorTemp-273.15,
            b=self.glazingRatio
            )

    def BEMCalc(self,UCM,BEM,forc,parameter,simTime):
        # Building Energy Model
        self.ElecTotal = 0.0                            # total electricity consumption - (W/m^2) of floor
        self.nFloor = max(UCM.bldHeight/float(self.floorHeight),1)   # At least one floor
        self.Qheat = 0.0                                # total sensible heat added
        self.sensCoolDemand = 0.0                       # building sensible cooling demand (W m-2)
        self.sensHeatDemand = 0.0                       # building sensible heating demand (W m-2)
        self.coolConsump  = 0.0                         # cooling energy consumption (W m-2)
        self.heatConsump  = 0.0                         # heating energy consumption (W m-2)
        self.sensWaste = 0.0                            # Sensible waste heat (W m-2)
        self.dehumDemand  = 0.0                         # dehumidification energy (W m-2)
        self.Qhvac = 0                                  # Total heat removed (sensible + latent)
        Qdehum = 0
        dens = map(lambda fP: fP/(1000*0.287042*self.indoorTemp*(1.+1.607858*self.indoorHum)), forc.pres) # [kgv/ m-3] Moist air density given dry bulb temperature, humidity ratio, and pressure
        evapEff = 1.                                    # evaporation efficiency in the condenser
        volVent = self.vent*self.nFloor                 # total vent volumetric flow for mass [m3 s-1 m-2 (bld/area)]
        volInfil = self.infil * UCM.bldHeight / 3600.   # Change of units AC/H -> [m3 s-1 m-2 (bld/facade#)]
        volSWH = BEM.SWH * self.nFloor/3600.            # Change of units l/hr per m^2 -> [L/s per m-2 (bld/area)]
        T_wall = BEM.wall.layerTemp[-1]                 # Inner layer
        T_ceil = BEM.roof.layerTemp[-1]                 # Inner layer
        T_mass = BEM.mass.layerTemp[0]                  # Outer layer
        T_indoor = self.indoorTemp                      # Indoor temp (initial)
        T_can = UCM.canTemp                             # Canyon temperature

        """
        % Normalize areas to building foot print [m^2/m^2(bld)]
        facArea = UCM.verToHor/UCM.bldDensity       % [m2/m2(bld)]
        wallArea = facArea*(1.-self.glazingRatio)    % [m2/m2(bld)]
        winArea = facArea*self.glazingRatio          % [m2/m2(bld)]
        massArea = 2*self.nFloor-1                   % ceiling/floor (top & bottom)

        % Temperature set points (updated per building schedule)
        if simTime.secDay/3600 < parameter.nightSetEnd || simTime.secDay/3600 >= parameter.nightSetStart
            T_cool = self.coolSetpointNight
            T_heat = self.heatSetpointNight
            self.intHeat = self.intHeatNight*self.nFloor
        else
            T_cool = self.coolSetpointDay
            T_heat = self.heatSetpointDay
            self.intHeat = self.intHeatDay*self.nFloor
        end

        % Indoor convection heat transfer coefficients
        zac_in_wall = 3.076
        zac_in_mass = 3.076
        if (T_ceil > T_indoor)
            zac_in_ceil  = 0.948
        elseif(T_ceil <= T_indoor)
            zac_in_ceil  = 4.040
        else
            disp('!!!!!FATAL ERROR!!!!!!')
            return
        end

        % -------------------------------------------------------------
        % Heat fluxes (per m^2 of bld footprint)
        % -------------------------------------------------------------
        % Solar Heat Gain
        winTrans = BEM.wall.solRec*self.shgc*winArea

        % Latent heat infiltration & ventilation (W/m^2 of bld footprint)
        QLinfil = volInfil * dens * parameter.lv *(UCM.canHum - self.indoorHum)
        QLvent = volVent * dens * parameter.lv *(UCM.canHum - self.indoorHum)
        QLintload = self.intHeat * self.intHeatFLat

        % Heat/Cooling load (W/m^2 of bld footprint), if any
        self.sensCoolDemand = max(wallArea*zac_in_wall*(T_wall-T_cool)+...
            massArea*zac_in_mass*(T_mass-T_cool)+...
            winArea*self.uValue*(T_can-T_cool)+...
            zac_in_ceil *(T_ceil-T_cool)+...
            self.intHeat+...
            volInfil * dens*parameter.cp*(T_can-T_cool)+...
            volVent * dens*parameter.cp*(T_can-T_cool) + ...
            winTrans,0)
        self.sensHeatDemand = max(-(wallArea*zac_in_wall*(T_wall-T_heat)+...
            massArea*zac_in_mass*(T_mass-T_heat)+...
            winArea*self.uValue*(T_can-T_heat)+...
            zac_in_ceil*(T_ceil-T_heat)+...
            self.intHeat+...
            volInfil*dens*parameter.cp*(T_can-T_heat)+...
            volVent*dens*parameter.cp*(T_can-T_heat) + ...
            winTrans),0)

        % -------------------------------------------------------------
        % HVAC system (cooling demand = W/m^2 bld footprint)
        % -------------------------------------------------------------
        if self.sensCoolDemand > 0 && UCM.canTemp > 288

            % Cooling energy is the equivalent energy to bring a vol
            % where sensCoolDemand = dens * Cp * x * (T_indoor - 10C) &
            % given 7.8g/kg of air at 10C, assume 7g/kg of air
            % dehumDemand = x * dens * (self.indoorHum -
            % 0.9*0.0078)*parameter.lv
            VolCool = self.sensCoolDemand / (dens*parameter.cp*(T_indoor-283.15))
            self.dehumDemand = max(VolCool * dens * (self.indoorHum - 0.9*0.0078)*parameter.lv,0)

            if (self.dehumDemand + self.sensCoolDemand) > (self.coolCap * self.nFloor)
                self.Qhvac = self.coolCap * self.nFloor
                VolCool = VolCool / (self.dehumDemand + self.sensCoolDemand) * (self.coolCap * self.nFloor)
                self.sensCoolDemand = self.sensCoolDemand * (self.coolCap * self.nFloor) / (self.dehumDemand + self.sensCoolDemand)
                self.dehumDemand = self.dehumDemand * (self.coolCap * self.nFloor) / (self.dehumDemand + self.sensCoolDemand)
            else
                self.Qhvac = self.dehumDemand + self.sensCoolDemand
            end
            Qdehum = VolCool * dens * parameter.lv * (self.indoorHum - 0.9*0.0078)
            self.coolConsump =(max(self.sensCoolDemand+self.dehumDemand,0))/self.copAdj

            % Waste heat from HVAC (per m^2 building foot print)
            if strcmp(self.condType,'AIR')
                self.sensWaste = max(self.sensCoolDemand+self.dehumDemand,0)+self.coolConsump
                self.latWaste = 0.0
            elseif strcmp(self.condType,'WAT') % Not sure if this works well
                self.sensWaste = max(self.sensCoolDemand+self.dehumDemand,0)+self.coolConsump*(1.-evapEff)
                self.latWaste = max(self.sensCoolDemand+self.dehumDemand,0)+self.coolConsump*evapEff
            end
            self.sensHeatDemand = 0

        % -------------------------------------------------------------
        % Heating system (heating demand = W/m^2 bld footprint)
        % -------------------------------------------------------------
        elseif self.sensHeatDemand > 0 && UCM.canTemp < 288

            % limit on heating capacity
            self.Qheat = min(self.sensHeatDemand,self.heatCap*self.nFloor)
            self.heatConsump  = self.Qheat / self.heatEff
            self.sensWaste = self.heatConsump - self.Qheat         % waste per footprint
            self.heatConsump = self.heatConsump/self.nFloor        % adjust to be per floor area
            self.sensHeatDemand = self.Qheat/self.nFloor           % adjust to be per floor area
            Qdehum = 0
            self.sensCoolDemand = 0
        end

        % -------------------------------------------------------------
        % Evolution of the internal temperature and humidity
        % -------------------------------------------------------------
        % wall, mass, roof, intload, infil, vent, hvac, heat, window
        Q = self.intHeat + winTrans + self.Qheat-self.sensCoolDemand

        H1 = T_wall*wallArea*zac_in_wall + ...
            T_mass*massArea*zac_in_mass + ...
            T_ceil*zac_in_ceil + ...
            T_can*winArea*self.uValue + ...
            T_can*volInfil * dens * parameter.cp + ...
            T_can*volVent * dens * parameter.cp

        H2 = wallArea*zac_in_wall + ...
            massArea*zac_in_mass + ...
            zac_in_ceil + ...
            winArea*self.uValue + ...
            volInfil * dens * parameter.cp + ...
            volVent * dens * parameter.cp

        self.indoorTemp = (H1 + Q)/H2
        self.indoorHum = self.indoorHum + simTime.dt/(dens * parameter.lv * UCM.bldHeight) * (...
            QLintload + QLinfil + QLvent - Qdehum)
        [~,~,self.indoorRhum,~,~,~] = Psychrometrics (self.indoorTemp, self.indoorHum, forc.pres)

        % These are used for element calculation (per m^2 of element area)
        self.fluxWall = zac_in_wall *(T_indoor - T_wall)
        self.fluxRoof = zac_in_ceil *(T_indoor - T_ceil)
        self.fluxMass = zac_in_mass *(T_indoor - T_mass) + self.intHeat * self.intHeatFRad/massArea

        % These are for record keeping only, per m^2 of floor area
        self.fluxSolar = winTrans/self.nFloor
        self.fluxWindow = winArea * self.uValue *(T_can - T_indoor)/self.nFloor
        self.fluxInterior = self.intHeat * self.intHeatFRad *(1.-self.intHeatFLat)/self.nFloor
        self.fluxInfil= volInfil * dens * parameter.cp *(T_can - T_indoor)/self.nFloor
        self.fluxVent = volVent * dens * parameter.cp *(T_can - T_indoor)/self.nFloor
        self.coolConsump = self.coolConsump/self.nFloor
        self.sensCoolDemand = self.sensCoolDemand/self.nFloor

        % Total Electricity/building floor area (W/m^2)
        self.ElecTotal = self.coolConsump + BEM.Elec + BEM.Light

        % Waste heat to canyon, W/m^2 of building + water
        CpH20 = 4200            % heat capacity of water
        T_hot = 49 + 273.15     % Service water temp (assume no storage)
        self.sensWaste = self.sensWaste + (1/self.heatEff-1)*(volSWH*CpH20*(T_hot - forc.waterTemp)) + BEM.Gas*(1-self.heatEff)*self.nFloor

        % Gas equip per floor + water usage per floor + heating/floor
        self.GasTotal = BEM.Gas + volSWH*CpH20*(T_hot - forc.waterTemp)/self.nFloor/self.heatEff + self.heatConsump
        """

"""
function r = root_finder(f,a,b)
    %root_finder calculates the roots of function f using the bisection search
    %method
    %INPUT function f, and interval a,b with the property that f(a) and f(b)
    %have opposite signs
    %OUTPUT r approximate value of root of f in interval [a,b]
    if (feval(f,a)*feval(f,b)) > 0
        disp('stop')
        error('Both endpoints have the same sign, please try again.')

    end

    while abs(b-a)>(10e-10)
        m = (a+b)/2
        x1 = feval(f,m)
        x2 = feval(f,a)
        if (x1 > 0 && x2 < 0) || (x1 < 0 && x2 > 0)
            b = m
        else
            a = m
        end
    end
    r = (a+b)/2
end

function w = mixing_ratio(T,p)
    %mixing_ratio finds the ratio of water vapor to the mass of dry air
    %INPUT Temperature and Pressure
    %OUTPUT MIXING RATIOs for inputting into wet_bulb.m
    p_a = 1013.246  % Standard sea-level atmospheric pressure in hPa
    T_a = 373.16  % Standard sea-level atmospheric temperature in Kelvin

    e1 = 11.344*(1-T/T_a)
    e2 = -3.49149*(T_a/T-1)
    f1 = -7.90298*(T_a/T-1)
    f2 = 5.02808*logn((T_a/T),10)
    f3 = -1.3816*((10^(e1)-1)/(1.e7))
    f4 = 8.1328*((10^(e2)-1)/(1.e3))
    f5 = logn(p_a,10)
    f = f1+f2+f3+f4+f5
    e = 10^(f)  % calculates vapor pressure in terms of T
    w = 0.62197*(e/(p-e))  % mass ratio g/kg
end

function [ z ] = logn(x,y)
    % logn
    %   Finds log base y of x
    z = log(x)/log(y);
end

"""
