classdef SchDef
    %   Building Energy Model (BEM) class definition
    %   JY, March, 2016
    
    properties
        Elec;        % Schedule for electricity (WD,Sat,Sun)
        Gas;         % Schedule for gas (WD,Sat,Sun)
        Light;       % Schedule for light (WD,Sat,Sun)
        Occ;         % Schedule for occupant (WD,Sat,Sun)
        Cool;        % Temperature schedule for cooling (WD,Sat,Sun)
        Heat;        % Temperature schedule for heating (WD,Sat,Sun)
        SWH;         % Hot water schedule (tbd)
        
        % Internal Heat Load from DOE
        Qelec;       % W/m^2 (max) for electrical plug process
        Qgas;        % W/m^2 (max) for gas process
        Qlight;      % W/m^2 (max) for light process
        Nocc;        % #/m^2 (max) for occupancy
        Vswh;        % Hot water vol/hr (max)
        Vent;        % litres/s/person for ventilation
    end
    
    methods
        function obj = SchDef(Elec,Gas,Light,Occ,Cool,Heat,SWH)
            % class constructor
            if(nargin > 0)
                obj.Elec = Elec;
                obj.Gas = Gas;
                obj.Light = Light;
                obj.Occ = Occ;
                obj.Cool = Cool;
                obj.Heat = Heat;
                obj.SWH = SWH;
            end
        end
    end
end