classdef BEMDef
    %   Building Energy Model (BEM) class definition
    %   JY, March, 2016
    
    properties
        building;       % building class type
        mass;           % mass element class
        wall;           % wall element class
        roof;           % roof element class
        
        frac;           % fraction of the urban floor space of this typology
        fl_area;        % Total building floor area in the urban area

        ElecTotal;      % Actual total electricity (W/m^2)
        Elec;           % Actual electricity consumption(W/m^2)
        Gas;            % Actual gas consumption(W/m^2)
        Light;          % Actual light (W/m^2)
        Qocc;           % Actual heat load from occupant (W/m^2)
        SWH;            % Actual hot water usage
        Nocc;           % Actual number of people
        
        T_wallex;       % Wall surface temp (ext)
        T_wallin;       % Wall surface temp (int)
        T_roofex;       % Roof surface temp
        T_roofin;       % 
        Q_wallex;       % Wall surface flux (ext)
        Q_wallin;       % Wall surface flux (int)
        Q_roofex;
        Q_roofin;
    end
    
    methods
        function obj = BEMDef(building,mass,wall,roof,frac)
            % class constructor
            if(nargin > 0)
                obj.building = building;
                obj.mass = mass;
                obj.wall = wall;
                obj.roof = roof;
                obj.frac = frac;
            end
        end
    end
end