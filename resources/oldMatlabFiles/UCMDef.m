classdef UCMDef
    % Definition of Urban Canopy - Building Energy Model Class
    
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
        roadConf;      % road-sky configuration factors
        alb_wall;      % average wall albedo
        wallConf;      % wall-sky configuration factors
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
    end
    
    methods
        function obj = UCMDef(bldHeight,bldDensity,verToHor,treeCoverage,sensAnthrop,latAnthrop,...
                initialTemp,initialHum,initialWind,parameter,r_glaze,SHGC,alb_wall,road)
            % class constructor
            if(nargin > 0)
                obj.road = road;
                obj.bldHeight = bldHeight;
                obj.verToHor = verToHor;            % for building only?
                obj.bldDensity = bldDensity;
                obj.treeCoverage = treeCoverage;
                obj.sensAnthrop = sensAnthrop;
                obj.latAnthrop = latAnthrop;
                obj.roadShad = min(treeCoverage/(1-bldDensity),1);
                obj.bldWidth = 4*bldHeight*bldDensity/verToHor; 
                d = obj.bldWidth/(bldDensity^0.5);
                obj.canWidth = d - obj.bldWidth;
                obj.canAspect = bldHeight/obj.canWidth;
                obj.roadConf = ((obj.canAspect)^2+1)^(1/2)-obj.canAspect;      
                obj.wallConf = 0.5*(obj.canAspect+1-...
                    ((obj.canAspect)^2+1)^(1/2))/(obj.canAspect); 
                obj.facArea = 4*obj.bldWidth*bldHeight;
                obj.roadArea = d^2 - (obj.bldWidth)^2;
                obj.roofArea = (obj.bldWidth)^2;
                obj.canTemp = initialTemp;
                obj.roadTemp = initialTemp;
                obj.canHum = initialHum;
                obj.ublWind = max(initialWind,parameter.windMin);
                obj.canWind = initialWind;
                obj.ustar = 0.1*initialWind;
                obj.ustarMod = 0.1*initialWind;
                frontDens = verToHor/4.;
                if frontDens < 0.15 
                  obj.z0u = frontDens*obj.bldHeight;
                else
                  obj.z0u = 0.15*obj.bldHeight;
                end
                if frontDens < 0.05 
                  obj.l_disp = 3*frontDens*obj.bldHeight;
                elseif frontDens < 0.15
                  obj.l_disp = (0.15+5.5*(frontDens-0.05))*obj.bldHeight;
                elseif frontDens < 1
                  obj.l_disp = (0.7+0.35*(frontDens-0.15))*obj.bldHeight;
                else 
                  obj.l_disp = 0.5*obj.bldHeight;
                end
                obj.alb_wall = alb_wall;
                obj.facAbsor = (1-r_glaze)*(1-alb_wall)+r_glaze*(1-0.75*SHGC);
                obj.roadAbsor = (1-road.vegCoverage)*(1-road.albedo);
                obj.sensHeat = 0.;
            end
        end
    
        function obj = UCModel(obj,BEM,T_ubl,forc,parameter)

            % Calculate the urban canyon temperature per The UWG (2012) Eq. 10 
            dens = forc.pres/(1000*0.287042*obj.canTemp*(1.+1.607858*obj.canHum));      % air density
            dens_ubl = forc.pres/(1000*0.287042*T_ubl*(1.+1.607858*forc.hum));      % air density
            Cp_air = parameter.cp;
            obj.Q_wall = 0;
            obj.Q_window = 0;
            obj.Q_road = 0;
            obj.Q_hvac = 0;
            obj.Q_traffic = 0;
            obj.Q_vent = 0;
            obj.Q_ubl = 0;
            obj.ElecTotal = 0;
            obj.GasTotal = 0;
            obj.roofTemp = 0;
            obj.wallTemp = 0;

            % Road to Canyon
            T_road = obj.road.layerTemp(1);
            h_conv = obj.road.aeroCond;
            H1 = T_road*h_conv*obj.roadArea;       % Heat (Sens) from road surface
            H2 = h_conv*obj.roadArea;              
            H1 = H1 + T_ubl*obj.roadArea*obj.uExch*Cp_air*dens_ubl; % Heat from UBL
            H2 = H2 + obj.roadArea*obj.uExch*Cp_air*dens_ubl;
            Q = (obj.roofArea+obj.roadArea)*(obj.sensAnthrop + obj.treeSensHeat*obj.treeCoverage);
            
            % Building energy output to canyon, in terms of absolute (total) values
            for j = 1:numel(BEM)

                % Re-naming variable for readability
                building = BEM(j).building;
                wall = BEM(j).wall;
                T_indoor = building.indoorTemp;
                T_wall = wall.layerTemp(1);
                R_glazing= building.glazingRatio;
                A_wall = (1-R_glazing)*obj.facArea;
                A_window = R_glazing*obj.facArea;
                U_window = building.uValue;
                
                H1 = H1 + BEM(j).frac*(...
                    T_indoor*A_window*U_window + ...    % window U
                    T_wall*A_wall*h_conv + ...          % Wall conv
                    T_indoor*obj.roofArea*BEM(j).building.vent*BEM(j).building.nFloor*Cp_air*dens + ...       % Vent
                    T_indoor*obj.roofArea*BEM(j).building.infil*obj.bldHeight/3600*Cp_air*dens);              % Infil
                
                H2 = H2 + BEM(j).frac*(...
                    A_window*U_window + ... 
                    A_wall*h_conv + ...
                    obj.roofArea*BEM(j).building.vent*BEM(j).building.nFloor*Cp_air*dens + ...    % Vent
                    obj.roofArea*BEM(j).building.infil*obj.bldHeight/3600*Cp_air*dens);           % Infil
                    
                Q = Q + BEM(j).frac*(...
                    obj.roofArea*building.sensWaste*obj.h_mix + ...         % HVAC waste heat
                    A_window*BEM(j).wall.solRec*(1-BEM(j).building.shgc));  % heat that didn't make it to inside
                
                obj.wallTemp = obj.wallTemp + BEM(j).frac*T_wall;
                obj.roofTemp = obj.roofTemp + BEM(j).frac*BEM(j).roof.layerTemp(1);
                obj.Q_ubl = obj.Q_ubl + BEM(j).frac*(BEM(j).roof.sens*obj.bldDensity+building.sensWaste*(1-obj.h_mix));
                    
            end
                
            % Solve for canyon temperature
            obj.canTemp = (H1 + Q)/H2;
            
            % Heat flux based per m^2 of urban area
            obj.Q_road = h_conv*(T_road-obj.canTemp)*(1-obj.bldDensity);  % Sensible heat from road (W/m^2 of urban area)
            obj.Q_ubl = obj.Q_ubl + obj.uExch*Cp_air*dens*(obj.canTemp-T_ubl)*(1-obj.bldDensity);
            obj.Q_wall = h_conv*(obj.wallTemp-obj.canTemp)*(obj.verToHor);
            obj.Q_traffic = obj.sensAnthrop;
            
            % Building energy output to canyon, per m^2 of urban area
            T_can = obj.canTemp;
            for j = 1:numel(BEM)
                V_vent = BEM(j).building.vent*BEM(j).building.nFloor; % ventilation volume per m^2 of building
                V_infil = BEM(j).building.infil*obj.bldHeight/3600;
                T_indoor = BEM(j).building.indoorTemp;
                R_glazing= building.glazingRatio;
                
                obj.Q_window = obj.Q_window + BEM(j).frac*obj.verToHor*R_glazing*U_window*(T_indoor-T_can);
                obj.Q_window = obj.Q_window + BEM(j).frac*obj.verToHor*R_glazing*BEM(j).wall.solRec*(1-BEM(j).building.shgc);
                obj.Q_vent = obj.Q_vent + BEM(j).frac*obj.bldDensity*Cp_air*dens*(V_vent + V_infil)*(T_indoor-T_can);
                obj.Q_hvac = obj.Q_hvac + BEM(j).frac*obj.bldDensity*BEM(j).building.sensWaste*obj.h_mix;

                obj.Q_roof = obj.Q_roof + BEM(j).frac*obj.bldDensity*BEM(j).roof.sens;
 
                % Total Electrical & Gas power in MW
                obj.ElecTotal = obj.ElecTotal + BEM(j).fl_area*BEM(j).building.ElecTotal/1e6;
                obj.GasTotal = obj.GasTotal + BEM(j).fl_area*BEM(j).building.GasTotal/1e6;
            end
            
            % Sensible Heat
            obj.sensHeat = obj.Q_wall + obj.Q_road + obj.Q_vent + obj.Q_window + obj.Q_hvac + obj.Q_traffic + obj.treeSensHeat + obj.Q_roof;
            
            % Error checking
            if obj.canTemp > 350 || obj.canTemp < 250
                disp('Something obviously went wrong (UCMDef.m)... ');
            end

        end
    end
end