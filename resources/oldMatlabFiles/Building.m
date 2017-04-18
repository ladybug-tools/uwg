classdef Building
    %   Building building class of specified building characteristics.
    
    properties
        % Building parameters
        floorHeight;        % floor height (m)
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
    end
    
    methods
        function obj = Building(floorHeight,intHeatNight,intHeatDay,intHeatFRad,...
                intHeatFLat,infil,vent,glazingRatio,uValue,shgc,...
                condType,cop,coolSetpointDay,coolSetpointNight,...
                heatSetpointDay,heatSetpointNight,coolCap,heatEff,initialTemp)
                
            % class constructor
                if (nargin > 0)
                    obj.floorHeight = floorHeight;
                    obj.intHeat = intHeatNight;
                    obj.intHeatNight = intHeatNight;
                    obj.intHeatDay = intHeatDay;
                    obj.intHeatFRad = intHeatFRad;    
                    obj.intHeatFLat = intHeatFLat;   
                    obj.infil = infil;      % ACH
                    obj.vent = vent;        
                    obj.glazingRatio = glazingRatio;
                    obj.uValue = uValue;
                    obj.shgc = shgc;
                    obj.condType = condType; 
                    obj.cop = cop;
                    obj.coolSetpointDay = coolSetpointDay;
                    obj.coolSetpointNight = coolSetpointNight;
                    obj.heatSetpointDay = heatSetpointDay; 
                    obj.heatSetpointNight = heatSetpointNight; 
                    obj.coolCap = coolCap;
                    obj.heatEff = heatEff;
                    obj.mSys = coolCap/1004./(min(coolSetpointDay,coolSetpointNight)-14-273.15);
                    obj.indoorTemp = initialTemp;
                    obj.indoorHum = 0.012;
                    obj.heatCap = 999;      % Default heat capacity value
                    obj.copAdj = cop;
                end
        end
        
        function obj = BEMCalc(obj,UCM,BEM,forc,parameter,simTime)
            
            % Building Energy Model (some of these can be moved up)
            obj.ElecTotal = 0;
            obj.nFloor = max(UCM.bldHeight/obj.floorHeight,1);  % At least one floor
             obj.Qheat = 0;
            obj.sensCoolDemand = 0.0;
            obj.sensHeatDemand = 0.0;
            obj.coolConsump  = 0.0;
            obj.heatConsump  = 0.0;
            obj.sensWaste = 0.0;
            obj.dehumDemand  = 0.0;
            obj.Qhvac = 0;
            Qdehum = 0;
            dens = forc.pres/(1000*0.287042*obj.indoorTemp*(1.+1.607858*obj.indoorHum));
            evapEff = 1.;                               % evaporation efficiency in the condenser
            volVent = obj.vent*obj.nFloor;              % [m3 s-1 m-2(bld)]
            volInfil = obj.infil*UCM.bldHeight/3600;    % Change of units AC/H -> [m3 s-1 m-2(bld)]
            volSWH = BEM.SWH * obj.nFloor/3600;
            T_wall = BEM.wall.layerTemp(end);           % Inner layer
            T_ceil = BEM.roof.layerTemp(end);           % Inner layer
            T_mass = BEM.mass.layerTemp(1);             % Outer layer
            T_indoor = obj.indoorTemp;                  % Indoor temp (initial)
            T_can = UCM.canTemp;                        % Canyon temperature
            
            % Normalize areas to building foot print [m^2/m^2(bld)]
            facArea = UCM.verToHor/UCM.bldDensity;      % [m2/m2(bld)]
            wallArea = facArea*(1.-obj.glazingRatio);   % [m2/m2(bld)]
            winArea = facArea*obj.glazingRatio;         % [m2/m2(bld)]
            massArea = 2*obj.nFloor-1;                  % ceiling/floor (top & bottom)

            % Temperature set points (updated per building schedule)
            if simTime.secDay/3600 < parameter.nightSetEnd || simTime.secDay/3600 >= parameter.nightSetStart
                T_cool = obj.coolSetpointNight;
                T_heat = obj.heatSetpointNight;
                obj.intHeat = obj.intHeatNight*obj.nFloor;
            else
                T_cool = obj.coolSetpointDay;
                T_heat = obj.heatSetpointDay;
                obj.intHeat = obj.intHeatDay*obj.nFloor;
            end

            % Indoor convection heat transfer coefficients
            zac_in_wall = 3.076;
            zac_in_mass = 3.076;
            if (T_ceil > T_indoor)
                zac_in_ceil  = 0.948;
            elseif(T_ceil <= T_indoor);
                zac_in_ceil  = 4.040;
            else
                disp('!!!!!FATAL ERROR!!!!!!');
                return;
            end
                        
            % -------------------------------------------------------------
            % Heat fluxes (per m^2 of bld footprint)
            % -------------------------------------------------------------
            % Solar Heat Gain
            winTrans = BEM.wall.solRec*obj.shgc*winArea;
            
            % Latent heat infiltration & ventilation (W/m^2 of bld footprint)
            QLinfil = volInfil * dens * parameter.lv *(UCM.canHum - obj.indoorHum);
            QLvent = volVent * dens * parameter.lv *(UCM.canHum - obj.indoorHum);
            QLintload = obj.intHeat * obj.intHeatFLat;
                        
            % Heat/Cooling load (W/m^2 of bld footprint), if any
            obj.sensCoolDemand = max(wallArea*zac_in_wall*(T_wall-T_cool)+...
                massArea*zac_in_mass*(T_mass-T_cool)+...
                winArea*obj.uValue*(T_can-T_cool)+...
                zac_in_ceil *(T_ceil-T_cool)+...
                obj.intHeat+...
                volInfil * dens*parameter.cp*(T_can-T_cool)+...
                volVent * dens*parameter.cp*(T_can-T_cool) + ...
                winTrans,0);
            obj.sensHeatDemand = max(-(wallArea*zac_in_wall*(T_wall-T_heat)+...
                massArea*zac_in_mass*(T_mass-T_heat)+...
                winArea*obj.uValue*(T_can-T_heat)+...
                zac_in_ceil*(T_ceil-T_heat)+...
                obj.intHeat+...
                volInfil*dens*parameter.cp*(T_can-T_heat)+...
                volVent*dens*parameter.cp*(T_can-T_heat) + ...
                winTrans),0);

            % -------------------------------------------------------------
            % HVAC system (cooling demand = W/m^2 bld footprint)
            % -------------------------------------------------------------
            if obj.sensCoolDemand > 0 && UCM.canTemp > 288
   
                % Cooling energy is the equivalent energy to bring a vol 
                % where sensCoolDemand = dens * Cp * x * (T_indoor - 10C) &
                % given 7.8g/kg of air at 10C, assume 7g/kg of air
                % dehumDemand = x * dens * (obj.indoorHum -
                % 0.9*0.0078)*parameter.lv
                VolCool = obj.sensCoolDemand / (dens*parameter.cp*(T_indoor-283.15));
                obj.dehumDemand = max(VolCool * dens * (obj.indoorHum - 0.9*0.0078)*parameter.lv,0);
                
                if (obj.dehumDemand + obj.sensCoolDemand) > (obj.coolCap * obj.nFloor)
                    obj.Qhvac = obj.coolCap * obj.nFloor;
                    VolCool = VolCool / (obj.dehumDemand + obj.sensCoolDemand) * (obj.coolCap * obj.nFloor);
                    obj.sensCoolDemand = obj.sensCoolDemand * (obj.coolCap * obj.nFloor) / (obj.dehumDemand + obj.sensCoolDemand);
                    obj.dehumDemand = obj.dehumDemand * (obj.coolCap * obj.nFloor) / (obj.dehumDemand + obj.sensCoolDemand);
                else
                    obj.Qhvac = obj.dehumDemand + obj.sensCoolDemand;
                end
                Qdehum = VolCool * dens * parameter.lv * (obj.indoorHum - 0.9*0.0078);
                obj.coolConsump =(max(obj.sensCoolDemand+obj.dehumDemand,0))/obj.copAdj;                
                
                % Waste heat from HVAC (per m^2 building foot print)
                if strcmp(obj.condType,'AIR')
                    obj.sensWaste = max(obj.sensCoolDemand+obj.dehumDemand,0)+obj.coolConsump;
                    obj.latWaste = 0.0;
                elseif strcmp(obj.condType,'WAT') % Not sure if this works well
                    obj.sensWaste = max(obj.sensCoolDemand+obj.dehumDemand,0)+obj.coolConsump*(1.-evapEff);
                    obj.latWaste = max(obj.sensCoolDemand+obj.dehumDemand,0)+obj.coolConsump*evapEff;
                end
                obj.sensHeatDemand = 0;

            % -------------------------------------------------------------
            % Heating system (heating demand = W/m^2 bld footprint)
            % -------------------------------------------------------------
            elseif obj.sensHeatDemand > 0 && UCM.canTemp < 288

                % limit on heating capacity
                obj.Qheat = min(obj.sensHeatDemand,obj.heatCap*obj.nFloor);
                obj.heatConsump  = obj.Qheat / obj.heatEff;
                obj.sensWaste = obj.heatConsump - obj.Qheat;        % waste per footprint 
                obj.heatConsump = obj.heatConsump/obj.nFloor;       % adjust to be per floor area
                obj.sensHeatDemand = obj.Qheat/obj.nFloor;          % adjust to be per floor area
                Qdehum = 0;
                obj.sensCoolDemand = 0;
            end

            % -------------------------------------------------------------
            % Evolution of the internal temperature and humidity
            % -------------------------------------------------------------
            % wall, mass, roof, intload, infil, vent, hvac, heat, window
            Q = obj.intHeat + winTrans + obj.Qheat-obj.sensCoolDemand;
            
            H1 = T_wall*wallArea*zac_in_wall + ...
                T_mass*massArea*zac_in_mass + ...
                T_ceil*zac_in_ceil + ...
                T_can*winArea*obj.uValue + ...
                T_can*volInfil * dens * parameter.cp + ...
                T_can*volVent * dens * parameter.cp;                
            
            H2 = wallArea*zac_in_wall + ...
                massArea*zac_in_mass + ...
                zac_in_ceil + ...
                winArea*obj.uValue + ...
                volInfil * dens * parameter.cp + ...
                volVent * dens * parameter.cp;

            obj.indoorTemp = (H1 + Q)/H2;
            obj.indoorHum = obj.indoorHum + simTime.dt/(dens * parameter.lv * UCM.bldHeight) * (...
                QLintload + QLinfil + QLvent - Qdehum);
            [~,~,obj.indoorRhum,~,~,~] = Psychrometrics (obj.indoorTemp, obj.indoorHum, forc.pres);

            % These are used for element calculation (per m^2 of element area)
            obj.fluxWall = zac_in_wall *(T_indoor - T_wall);
            obj.fluxRoof = zac_in_ceil *(T_indoor - T_ceil);
            obj.fluxMass = zac_in_mass *(T_indoor - T_mass) + obj.intHeat * obj.intHeatFRad/massArea;

            % These are for record keeping only, per m^2 of floor area
            obj.fluxSolar = winTrans/obj.nFloor;
            obj.fluxWindow = winArea * obj.uValue *(T_can - T_indoor)/obj.nFloor;
            obj.fluxInterior = obj.intHeat * obj.intHeatFRad *(1.-obj.intHeatFLat)/obj.nFloor;
            obj.fluxInfil= volInfil * dens * parameter.cp *(T_can - T_indoor)/obj.nFloor;
            obj.fluxVent = volVent * dens * parameter.cp *(T_can - T_indoor)/obj.nFloor; 
            obj.coolConsump = obj.coolConsump/obj.nFloor;
            obj.sensCoolDemand = obj.sensCoolDemand/obj.nFloor;

            % Total Electricity/building floor area (W/m^2)
            obj.ElecTotal = obj.coolConsump + BEM.Elec + BEM.Light;

            % Waste heat to canyon, W/m^2 of building + water
            CpH20 = 4200;           % heat capacity of water
            T_hot = 49 + 273.15;    % Service water temp (assume no storage)
            obj.sensWaste = obj.sensWaste + (1/obj.heatEff-1)*(volSWH*CpH20*(T_hot - forc.waterTemp)) + BEM.Gas*(1-obj.heatEff)*obj.nFloor;
            
            % Gas equip per floor + water usage per floor + heating/floor
            obj.GasTotal = BEM.Gas + volSWH*CpH20*(T_hot - forc.waterTemp)/obj.nFloor/obj.heatEff + obj.heatConsump;                        

        end
    end
end

function [Tdb, w, phi, h, Tdp, v] = Psychrometrics (Tdb_in, w_in, P)
    % Modified version of Psychometrics by Tea Zakula
    % MIT Building Technology Lab
    % Tdb (dry bulb temperature) and Tdp(dew point temperature) in C
    % w (humidity ratio) in kg/kg of dry air
    % phi (relative humidity) in %
    % h (enthalpy) in J/kg of dry air
    % v (specific volume) in m3/kg of dry air
    % P (Atmospheric Station Pressure) in Pa

    c_air = 1006;   %J/kg, value from ASHRAE Fundamentals
    hlg = 2501000;  %J/kg, value from ASHRAE Fundamentals
    cw  = 1860;     %J/kg, value from ASHRAE Fundamentals
    P = P/1000;     % convert from Pa to kPa

    Tdb = Tdb_in - 273.15;
    w = w_in;

    % phi calculation from Tdb and w
    Pw = w*P/(0.621945+w);              % Partial pressure of water wapor
    Pws = Saturation_pressure(Tdb);
    phi = Pw/Pws*100;

    h = c_air*Tdb+w*(hlg+cw*Tdb);       % Enthalpy
    v = 0.287042*(Tdb+273.15)*(1+1.607858*w)/P; % Specific volume 

    % Dew point 
    pw = (P*w)/(0.621945+w); % water vapor partial pressure in kPa
    alpha = log(pw);
    Tdp = 6.54 + 14.526*alpha+0.7389*(alpha^2)+0.09486*(alpha^3)+0.4569*(pw^0.1984); % valid for Tdp between 0 C and 93 C
end

function [Pws] = Saturation_pressure(Tdb)
    T = Tdb+273.15;
    Pws = exp(-(5.8002206e3)/T+1.3914993+-(4.8640239e-2)*T+(4.1764768e-5)*(T^2)-(1.4452093e-8)*(T^3)+6.5459673*log(T)); %in Pa
    Pws = Pws/1000; % in kPa
end

function psat = psat(temp,parameter)
    gamw  = (parameter.cl - parameter.cpv) / parameter.rv;
    betaw = (parameter.lvtt/parameter.rv) + (gamw * parameter.tt);
    alpw = log(parameter.estt) + (betaw /parameter.tt) + (gamw *log(parameter.tt));
    psat = zeros(size(temp));
    for jj=1:size(temp)
        psat = exp(alpw - betaw/temp - gamw*log(temp));
    end
end

% Not used for this release but saved for possible future use
function Twb = wet_bulb(Tdb,Tdp,pres)

    % Copyright (c) 2015, Rolf Henry Goodwin
    % All rights reserved.
    % 
    % Redistribution and use in source and binary forms, with or without
    % modification, are permitted provided that the following conditions are
    % met:
    % 
    %     * Redistributions of source code must retain the above copyright
    %       notice, this list of conditions and the following disclaimer.
    %     * Redistributions in binary form must reproduce the above copyright
    %       notice, this list of conditions and the following disclaimer in
    %       the documentation and/or other materials provided with the distribution
    % 
    % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    % AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    % IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    % ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
    % LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    % CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    % SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    % INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    % CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    % ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    % POSSIBILITY OF SUCH DAMAGE.
    
    % Code modified to merge into a single file - Joseph Yang, 2016


    % Tdb, Tdp, Twb in K
    % p in Pa (obtained function uses hPa, so /100 needed)    
    global T;
    global T_d;
    global p;
    T = Tdb;
    T_d = Tdp;
    p = pres/100;

    Twb = root_finder(@Delta_q,T_d,T);
end

function dQTw = Delta_q(T_w)
    %Delta_q finds the value of function dq(Tw)
    %INPUT wet bulb temperature T_w
    %OUTPUT dq(Tw)
    global T;
    global T_d;
    global p;

    Cp = 1005; % Heat capacity of water vapor in J/(kg*K)
    L = 2.501e6; % Latent heat of water vapor at 0 degC in J/kg
    w1 = mixing_ratio(T_d,p); % Mixing ratio corresponding to T_d and p
    w2 = mixing_ratio(T_w,p); % Mixing ratio corresponding to T_w and p

    dQTw = (L*(w2-w1))/(1+w2)-Cp*(T-T_w)*(1+0.8*w2); % Finds deltaq(Tw)

end

function r = root_finder(f,a,b)
    %root_finder calculates the roots of function f using the bisection search
    %method
    %INPUT function f, and interval a,b with the property that f(a) and f(b)
    %have opposite signs
    %OUTPUT r approximate value of root of f in interval [a,b]
    if (feval(f,a)*feval(f,b)) > 0
        disp('stop');
        error('Both endpoints have the same sign, please try again.')
        
    end

    while abs(b-a)>(10e-10)
        m = (a+b)/2;
        x1 = feval(f,m);
        x2 = feval(f,a);
        if (x1 > 0 && x2 < 0) || (x1 < 0 && x2 > 0)
            b = m;
        else
            a = m;
        end
    end
    r = (a+b)/2;
end

function w = mixing_ratio(T,p)
    %mixing_ratio finds the ratio of water vapor to the mass of dry air
    %INPUT Temperature and Pressure
    %OUTPUT MIXING RATIOs for inputting into wet_bulb.m
    p_a = 1013.246; % Standard sea-level atmospheric pressure in hPa
    T_a = 373.16; % Standard sea-level atmospheric temperature in Kelvin

    e1 = 11.344*(1-T/T_a);
    e2 = -3.49149*(T_a/T-1);
    f1 = -7.90298*(T_a/T-1);
    f2 = 5.02808*logn((T_a/T),10);
    f3 = -1.3816*((10^(e1)-1)/(1.e7));
    f4 = 8.1328*((10^(e2)-1)/(1.e3));
    f5 = logn(p_a,10);
    f = f1+f2+f3+f4+f5;
    e = 10^(f); % calculates vapor pressure in terms of T
    w = 0.62197*(e/(p-e)); % mass ratio g/kg
end

function [ z ] = logn(x,y)
    % logn
    %   Finds log base y of x
    z = log(x)/log(y);
end
