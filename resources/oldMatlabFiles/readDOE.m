% Read Excel files of DOE buildings
% Sheet 1 = BuildingSummary
% Sheet 2 = ZoneSummary
% Sheet 3 = LocationSummary
% Sheet 4 = Schedules
% Note BLD8 & 10 = school

% DOE Building Types
bldType = {'FullServiceRestaurant';
    'Hospital';
    'LargeHotel';
    'LargeOffice';
    'MedOffice';
    'MidRiseApartment';
    'OutPatient';
    'PrimarySchool';
    'QuickServiceRestaurant';
    'SecondarySchool';
    'SmallHotel';
    'SmallOffice';
    'StandAloneRetail';
    'StripMall';
    'SuperMarket';
    'WareHouse'};

zoneType = {'1A (Miami)';
    '2A (Houston)';
    '2B (Phoenix)';
    '3A (Atlanta)';
    '3B-CA (Los Angeles)';
    '3B (Las Vegas)';
    '3C (San Francisco)';
    '4A (Baltimore)';
    '4B (Albuquerque)';
    '4C (Seattle)';
    '5A (Chicago)';
    '5B (Boulder)';
    '6A (Minneapolis)';
    '6B (Helena)';
    '7 (Duluth)';
    '8 (Fairbanks)'};

builtEra = {'Pre80';
    'Pst80';
    'New'};

% Building (Type,Era,Zone), Type = 1-16, Era = 1-3, Zone = 1-16
refDOE (16,3,16) = Building;
Schedule (16,3,16) = SchDef;
refBEM (16,3,16) = BEMDef;

for i = 1:16
    file = strcat('C:\Sim\UWG4.1\data\DOERefBuildings\BLD',num2str(i),'.xlsx');

    % Read building summary (Sheet 1)
    [num, ~, ~] = xlsread(file,1);
    nFloor = num(4,4:6);
    glazing = num(5,4:6);
    hCeiling = num(6,4:6);
    ver2hor = num(8,4:6);
    AreaRoof = num(9,4:6);

    % Read zone summary (Sheet 2)
    [num, ~, ~] = xlsread(file,2);
    AreaFloor = num(3:5,6);
    Volume = num(3:5,7);
    AreaWall = num(3:5,9);
    AreaWindow = num(3:5,10);
    Occupant = num(3:5,12);
    Light = num(3:5,13);
    Elec = num(3:5,14);
    Gas = num(3:5,15);
    SHW = num(3:5,16);      % litres/hour
    Infil = num(3:5,21);    % ACH
    Vent = num(3:5,18);     % L/s/m^2
    
    % Read location summary (Sheet 3)
    [num, text, ~] = xlsread(file,3);
    TypeWall = [text(4,5:20); text(15,5:20); text(26,5:20)];
    RvalWall = [num(5,5:20); num(16,5:20); num(27,5:20)];
    TypeRoof = [text(6,5:20); text(17,5:20); text(28,5:20)];
    RvalRoof = [num(7,5:20); num(18,5:20); num(29,5:20)];
    Uwindow = [num(8,5:20); num(19,5:20); num(30,5:20)];
    SHGC = [num(9,5:20); num(20,5:20); num(31,5:20)];
    HVAC = [num(10,5:20); num(21,5:20); num(32,5:20)];
    HEAT = [num(11,5:20); num(22,5:20); num(33,5:20)];
    COP = [num(12,5:20); num(23,5:20); num(34,5:20)];
    EffHeat = [num(13,5:20); num(24,5:20); num(35,5:20)];
    FanFlow = [num(14,5:20); num(25,5:20); num(36,5:20)];
    
    % Read schedule (Sheet 4)
    [num, text, ~] = xlsread(file,4);
    SchEquip = num(2:4,7:30);
    SchLight = num(5:7,7:30);
    SchOcc = num(8:10,7:30);
    SetCool = num(11:13,7:30);
    SetHeat = num(14:16,7:30);
    SchGas = num(17:19,7:30);
    SchSWH = num(20:22,7:30);
    
    for j = 1:3
        for k = 1:16
            refDOE (i,j,k) = Building(hCeiling(j),...  % floorHeight
                1,...                           % intHeatNight
                1,...                           % intHeatDay
                0.1,...                         % intHeatFRad
                0.1,...                         % intHeatFLat
                Infil(j),...                    % infil (ACH)
                Vent(j)/1000,...                % vent (m^3/s/m^2)
                glazing(j),...                  % glazing ratio
                Uwindow(j,k),...                % uValue
                SHGC(j,k),...                   % shgc
                'AIR',...                       % a/c type
                COP(j,k),...                    % cop
                297,...                         % coolSetpointDay
                297,...                         % coolSetpointNight
                293,...                         % heatSetpointDay
                293,...                         % heatSetpointNight
                HVAC(j,k)*1000/AreaFloor(j),... % coolCap (refDOE in kW)
                EffHeat(j,k),...                % heatEff
                293);                           % initialTemp
            
            refDOE(i,j,k).heatCap = HEAT(j,k)*1000/AreaFloor(j);
            refDOE(i,j,k).Type = bldType(i);
            refDOE(i,j,k).Era = builtEra(j);
            refDOE(i,j,k).Zone = zoneType(k);

            % Define wall, roof, road, and mass
            % Material (thermalCond, volHeat = specific heat * density);
            Stucco = Material (0.6918, 837.0 * 1858.0);
            Concrete = Material (1.311, 836.8 * 2240);
            Insulation = Material (0.049, 836.8 * 265.0);
            Gypsum = Material (0.16, 830.0 * 784.9);
            Wood = Material (0.11,1210.0*544.62);
                                    
            % Wall (1 in stucco, concrete, insulation, gypsum)
            if strcmp(TypeWall(j,k),'MassWall')
                % 1" stucco, 8" concrete, tbd insulation, 1/2" gypsum
                Rbase = 0.271087; % based on stucco, concrete, gypsum
                Rins = RvalWall(j,k) - Rbase;
                D_ins = Rins * Insulation.thermalCond;
                if D_ins > 0.01
                    thickness = [0.0254;0.0508;0.0508;0.0508;0.0508;D_ins;0.0127];
                    layers = [Stucco;Concrete;Concrete;Concrete;Concrete;Insulation;Gypsum];
                else
                    thickness = [0.0254;0.0508;0.0508;0.0508;0.0508;0.0127];
                    layers = [Stucco;Concrete;Concrete;Concrete;Concrete;Gypsum];
                end

                wall = Element(0.08,0.92,thickness,layers,0,293,0);
                
                % If mass wall, assume mass foor (4" concrete)
                % Mass (assume 4" concrete);
                alb = 0.2;
                emis = 0.9;
                thickness = [0.054;0.054];
                concrete = Material (1.31,2240.0*836.8);
                mass = Element(alb,emis,thickness,[concrete;concrete],0,293,1);
                                                
            elseif strcmp(TypeWall(j,k),'WoodFrame')
                % 0.01m wood siding, tbd insulation, 1/2" gypsum
                Rbase = 0.170284091;    % based on wood siding, gypsum
                Rins = RvalWall(j,k) - Rbase;
                D_ins = Rins * Insulation.thermalCond;
                
                if D_ins > 0.01
                    thickness = [0.01;D_ins;0.0127];
                    layers = [Wood;Insulation;Gypsum];
                else
                    thickness = [0.01;0.0127];
                    layers = [Wood;Gypsum];
                end
                wall = Element(0.22,0.92,thickness,layers,0,293,0);

                % If wood frame wall, assume wooden floor
                alb = 0.2;
                emis = 0.9;
                thickness = 0.05 * ones (2,1);
                wood = Material (1.31,2240.0*836.8);
                mass = Element(alb,emis,thickness,[wood;wood],0,293,1);
                
            elseif strcmp(TypeWall(j,k),'SteelFrame')
                % 1" stucco, 8" concrete, tbd insulation, 1/2" gypsum
                Rbase = 0.271087; % based on stucco, concrete, gypsum
                Rins = RvalWall(j,k) - Rbase;
                D_ins = Rins * Insulation.thermalCond;
                if D_ins > 0.01
                    thickness = [0.0254;0.0508;0.0508;0.0508;0.0508;D_ins;0.0127];
                    layers = [Stucco;Concrete;Concrete;Concrete;Concrete;Insulation;Gypsum];
                else    % If insulation is too thin, assume no insulation
                    thickness = [0.0254;0.0508;0.0508;0.0508;0.0508;0.0127];
                    layers = [Stucco;Concrete;Concrete;Concrete;Concrete;Gypsum];
                end
                wall = Element(0.15,0.92,thickness,layers,0,293,0);
                
%                 wall = Element(0.2,0.92,thickness,layers,0,293,0); % Singapore case

%                 % TOULOUS case
%                 thickness = [0.05;0.05;0.05;0.05;0.05;0.05;0.03;0.0127];
%                 layers = [Concrete;Concrete;Concrete;Concrete;Concrete;Concrete;Insulation;Gypsum];
%                 wall = Element(0.25,0.93,thickness,layers,0,293,0);
% 
%                 % BUBBLE case
%                 thickness = [0.05;0.05;0.05;0.05;0.03;0.0127];
%                 layers = [Concrete;Concrete;Concrete;Concrete;Insulation;Gypsum];
%                 wall = Element(0.15,0.93,thickness,layers,0,293,0);

                % If mass wall, assume mass foor
                % Mass (assume 4" concrete);
                alb = 0.2;      % Adjusted for Bubble/Toulouse case (0.08 per Energy Plus)
                emis = 0.93;
                thickness = 0.05 * ones (2,1); 
%                thickness = 0.102*ones(2,1);       % for Bubble 
                mass = Element(alb,emis,thickness,[Concrete;Concrete],0,293,1);
                
            elseif strcmp(TypeWall(j,k),'MetalWall')

                % metal siding, insulation, 1/2" gypsum
                alb = 0.2;
                emis = 0.9;
                D_ins = max((RvalWall(j,k) * Insulation.thermalCond)/2,0.01);                
                wall = Element(alb,emis,[D_ins;D_ins;0.0127],[Insulation;Insulation;Gypsum],0,293,0);

                % Mass (assume 4" concrete);
                alb = 0.2;
                emis = 0.9;
                thickness = 0.05 * ones (2,1);
                concrete = Material (1.31,2240.0*836.8);
                mass = Element(alb,emis,thickness,[concrete;concrete],0,293,1);

            end
            
            % Roof
            if strcmp(TypeRoof(j,k),'IEAD')
                % IEAD-> membrane, insulation, decking
                 alb = 0.2;
                 emis = 0.93;
                 D_ins = max(RvalRoof(j,k) * Insulation.thermalCond / 2,0.01);
                 roof = Element(alb,emis,[D_ins;D_ins],[Insulation;Insulation],0,293,0);
                
%                 % TOULOUS/BUBBLE Case
%                 alb = 0.15;
%                 emis = 0.93;
%                 thickness = [0.06;0.05;0.05;0.05;0.05;0.03];
%                 roof = Element(alb,emis,thickness,[Concrete;Wood;Wood;Wood;Wood;Insulation],0,293,0);
            
            elseif strcmp(TypeRoof(j,k),'Attic')
                % IEAD-> membrane, insulation, decking
                alb = 0.2;
                emis = 0.9;
                D_ins = max(RvalRoof(j,k) * Insulation.thermalCond/2,0.01);
                roof = Element(alb,emis,[D_ins;D_ins],[Insulation;Insulation],0,293,0);
                
            elseif strcmp(TypeRoof(j,k),'MetalRoof')
                % IEAD-> membrane, insulation, decking
                alb = 0.2;
                emis = 0.9;
                D_ins = max(RvalRoof(j,k) * Insulation.thermalCond/2,0.01);
                roof = Element(alb,emis,[D_ins;D_ins],[Insulation;Insulation],0,293,0);
            end
                        
            % Define bulding energy model, set fraction to zero
            refBEM(i,j,k) = BEMDef(refDOE(i,j,k),mass,wall,roof,0);
            refBEM(i,j,k).building.FanMax = FanFlow(j,k);
            
            Schedule(i,j,k).Elec = SchEquip;   % 3x24 matrix of schedule for electricity (WD,Sat,Sun)
            Schedule(i,j,k).Light = SchLight;  % 3x24 matrix of schedule for light (WD,Sat,Sun)
            Schedule(i,j,k).Gas = SchGas;      % 3x24 matrix of schedule for gas (WD,Sat,Sun)
            Schedule(i,j,k).Occ = SchOcc;      % 3x24 matrix of schedule for occupancy (WD,Sat,Sun)
            Schedule(i,j,k).Cool = SetCool;    % 3x24 matrix of schedule for cooling temp (WD,Sat,Sun) % off for BUBBLE case
            Schedule(i,j,k).Heat = SetHeat;    % 3x24 matrix of schedule for heating temp (WD,Sat,Sun)
            Schedule(i,j,k).SWH = SchSWH;      % 3x24 matrix of schedule for SWH (WD,Sat,Sun)
            
            Schedule(i,j,k).Qelec = Elec(j);         % W/m^2 (max) for electrical plug process
            Schedule(i,j,k).Qlight = Light(j);       % W/m^2 (max) for light
            Schedule(i,j,k).Nocc = Occupant(j)/AreaFloor(j); % Person/m^2
            Schedule(i,j,k).Qgas = Gas(j);           % W/m^2 (max) for gas
            Schedule(i,j,k).Vent = Vent(j)/1000;     % m^3/m^2 per person
            Schedule(i,j,k).Vswh = SHW(j)/AreaFloor(j);    % litres per hour per m^2 of floor
            
        end
    end
end

% % BUBBLE/TOULOUSE adjustment Case
% refBEM(6,2,5).building.glazingRatio = 0.3;

% % Singapore adjustment Case
% refBEM(6,2,1).building.glazingRatio = 0.3;

save ('RefDOE.mat','refDOE','refBEM','Schedule');

% wall & roof definition based on material
% Material,
%     1/2IN Gypsum,            !- Name
%     Smooth,                  !- Roughness
%     0.0127,                  !- Thickness {m}
%     0.1600,                  !- Conductivity {W/m-K}
%     784.9000,                !- Density {kg/m3}
%     830.0000,                !- Specific Heat {J/kg-K}
%     0.9000,                  !- Thermal Absorptance
%     0.9200,                  !- Solar Absorptance
%     0.9200;                  !- Visible Absorptance
% 
% Material,
%     1IN Stucco,              !- Name
%     Smooth,                  !- Roughness
%     0.0253,                  !- Thickness
%     0.6918,                  !- Conductivity
%     1858.0000,               !- Density
%     837.0000,                !- Specific Heat
%     0.9000,                  !- Thermal Absorptance
%     0.9200,                  !- Solar Absorptance
%     0.9200;                  !- Visible Absorptance
% 
% Material,
%     8IN CONCRETE HW,  !- Name
%     Rough,                   !- Roughness
%     0.2032,                  !- Thickness {m}
%     1.3110,                  !- Conductivity {W/m-K}
%     2240.0000,               !- Density {kg/m3}
%     836.8000,                !- Specific Heat {J/kg-K}
%     0.9000,                  !- Thermal Absorptance
%     0.7000,                  !- Solar Absorptance
%     0.7000;                  !- Visible Absorptance
%
% Material,
%     Mass NonRes Wall Insulation, !- Name
%     MediumRough,             !- Roughness
%     0.0484268844343858,      !- Thickness {m}
%     0.049,                   !- Conductivity {W/m-K}
%     265.0000,                !- Density {kg/m3}
%     836.8000,                !- Specific Heat {J/kg-K}
%     0.9000,                  !- Thermal Absorptance
%     0.7000,                  !- Solar Absorptance
%     0.7000;                  !- Visible Absorptance
%
% Material,
%     Std Wood 6inch,          !- Name
%     MediumSmooth,            !- Roughness
%     0.15,                    !- Thickness {m}
%     0.12,                    !- Conductivity {W/m-K}
%     540.0000,                !- Density {kg/m3}
%     1210,                    !- Specific Heat {J/kg-K}
%     0.9000000,               !- Thermal Absorptance
%     0.7000000,               !- Solar Absorptance
%     0.7000000;               !- Visible Absorptance! Common Materials
% 
% Material,
%     Wood Siding,             !- Name
%     MediumSmooth,            !- Roughness
%     0.0100,                  !- Thickness {m}
%     0.1100,                  !- Conductivity {W/m-K}
%     544.6200,                !- Density {kg/m3}
%     1210.0000,               !- Specific Heat {J/kg-K}
%     0.9000,                  !- Thermal Absorptance
%     0.7800,                  !- Solar Absorptance
%     0.7800;                  !- Visible Absorptance