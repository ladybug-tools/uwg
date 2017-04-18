classdef Material
    properties
        thermalCond; % thermal conductivity (W m-1 K-1)
        volHeat;     % volumetric heat capacity (J m-3 K-1)
    end

    methods
        function obj = Material(thermalCond,volHeat)
        % class constructor
            if(nargin > 0) % don't have to initialize the fields
                obj.thermalCond  = thermalCond;
                obj.volHeat      = volHeat;
            end
        end          
    end
end