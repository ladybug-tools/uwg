classdef Forcing
    %FORCING Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        deepTemp;   % deep soil temperature (K)
        infra;      % 
        wind;
        uDir;
        hum;
        pres;       % Pressure (Pa)
        temp;
        rHum;
        prec;       % Precipitation (m s-1)
        dir;
        dif;
        waterTemp;  % ground water temp, set to temp at 2m
    end
    
    methods
        function obj = Forcing(staTemp,weather)
            if(nargin > 0)
                obj.deepTemp = mean(staTemp);
                obj.waterTemp = mean(staTemp);
                obj.infra = [weather.staInfra];
                obj.uDir = [weather.staUdir];
                obj.hum  = [weather.staHum];
                obj.pres = [weather.staPres];
                obj.temp = [weather.staTemp];
                obj.rHum = [weather.staRhum];
                obj.dir = [weather.staDir];
                obj.dif = [weather.staDif];
                obj.prec = [weather.staRobs]/3.6e6;
                obj.wind = [weather.staUmod];                
            end
        end
    end
end