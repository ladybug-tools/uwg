classdef UBLDef
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
    
    methods
        function obj = UBLDef(location,charLength,initialTemp,maxdx,dayBLHeight,nightBLHeight)
            % class constructor
            if(nargin > 0)
                obj.location = location;
                obj.charLength = charLength;
                obj.perimeter = 4*charLength;
                obj.urbArea = charLength^2.;
                obj.orthLength = charLength;
                numdx = round(charLength/min(charLength,maxdx));
                obj.paralLength = charLength/numdx;
                obj.ublTemp = initialTemp;
                obj.ublTempdx = initialTemp*ones(1,numdx);
                obj.dayBLHeight = dayBLHeight;
                obj.nightBLHeight = nightBLHeight;
            end
        end  
        
        function obj = UBLModel(obj,UCM,RSM,rural,forc,parameter,simTime)
 
            % Note that only one urban canyon area is considered 
            obj.sensHeat = UCM.sensHeat;
            heatDif = max(obj.sensHeat - rural.sens,0);     
            Cp = parameter.cp;                  % Heat capacity of air (J/kg.K)
            k_w = parameter.circCoeff;          % k_w per Bueno 'the UWG', eq 8
            g = parameter.g;                    % Gravity
            v_wind = max(forc.wind,parameter.windMin);   % wind velocity

            % Air density
            refDens = 0;
            for iz=1:RSM.nzref
                refDens = refDens + RSM.densityProfC(iz)*RSM.dz(iz)/...
                    (RSM.z(RSM.nzref)+RSM.dz(RSM.nzref)/2);
            end
            forDens = 0;
            for iz=1:RSM.nzfor
                forDens = forDens + RSM.densityProfC(iz)*RSM.dz(iz)/...
                    (RSM.z(RSM.nzfor)+RSM.dz(RSM.nzfor)/2);
            end
            
            % ---------------------------------------------------------------------
            % Day
            % ---------------------------------------------------------------------
            time = simTime.secDay/3600;
            noon = 12;
            daylimit = parameter.dayThreshold;      % sunlight threshold for day (~150W/m^2)
            nightlimit = parameter.dayThreshold;    % sunlight threshold for night (~50W/m^2)
            sunlight = forc.dir+forc.dif;
            
            % If dir & dif light is greater than threshold, use day 
            if sunlight > daylimit && time <= noon ||...
                    sunlight > nightlimit && time > noon || obj.sensHeat > 150                
                
                % Circulation velocity per Bueno 'the UWG', eq 8
                h_UBL = obj.dayBLHeight;            % Day boundary layer height
                eqTemp = RSM.tempProf(RSM.nzref);
                eqWind = RSM.windProf(RSM.nzref);
                
                Csurf = UCM.Q_ubl*simTime.dt/(h_UBL*refDens*Cp);
                u_circ = k_w*(g*heatDif/Cp/refDens/eqTemp*h_UBL)^(1./3.);
                                
                if v_wind > u_circ    % Forced problem (usually this)
                    advCoef  = obj.orthLength*eqWind*simTime.dt/obj.urbArea*1.4;
                    obj.ublTemp = (Csurf+advCoef*eqTemp + obj.ublTemp)/(1 + advCoef);
                    obj.ublTempdx(:)= obj.ublTemp;

                else                  % Convective problem
                    advCoef  = obj.perimeter*u_circ*simTime.dt/obj.urbArea*1.4;
                    obj.ublTemp = (Csurf+advCoef*eqTemp + obj.ublTemp)/(1 + advCoef);
                    obj.ublTempdx(:)= obj.ublTemp;
                end
            % ---------------------------------------------------------------------
            % Night 
            % ---------------------------------------------------------------------
            else
                h_UBL = obj.nightBLHeight;      % Night boundary layer height
                Csurf = UCM.Q_ubl*simTime.dt/(h_UBL*refDens*Cp);
                 [obj.ublTemp,obj.ublTempdx] = NightForc(obj.ublTempdx,simTime.dt,...
                      h_UBL,obj.paralLength,obj.charLength,RSM,Csurf);
            end
        end
    end
end

function [ublTemp,ublTempdx] = NightForc(ublTempdx,dt,h_UBL,paralLength,charLength,RSM,Csurf)

    % Night forcing (RSM.nzfor = number of layers of forcing)
    % Average potential temperature & wind speed of the profile
    intAdv1 = 0;
    for iz=1:RSM.nzfor
        intAdv1 = intAdv1 + RSM.windProf(iz)*RSM.tempProf(iz)*RSM.dz(iz);
    end
    advCoef1 = 1.4*dt/paralLength/h_UBL*intAdv1;
    
    intAdv2 = 0;
    for iz=1:RSM.nzfor
        intAdv2 = intAdv2 + RSM.windProf(iz)*RSM.dz(iz);
    end
    advCoef2 = 1.4*dt/paralLength/h_UBL*intAdv2;
    
    ublTempdx(1) = (Csurf + advCoef1 + ublTempdx(1))/(1 + advCoef2);     
    ublTemp = ublTempdx(1);
    
    for i=2:(charLength/paralLength)
        eqTemp = ublTempdx(i-1);
        ublTempdx(i) = (Csurf + advCoef2*eqTemp + ublTempdx(i))/(1 + advCoef2); 
        ublTemp = ublTemp + ublTempdx(i);
    end
    ublTemp = ublTemp/charLength*paralLength;
end