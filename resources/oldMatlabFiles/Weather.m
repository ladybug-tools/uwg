classdef Weather
    %FORCING Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        staTemp;   % air temperature (C)
        staRhum;   % air relative humidity (%)
        staPres;   % air pressure (Pa)
        staInfra;  % horizontal Infrared Radiation Intensity (W m-2)
        staHor;    % horizontal radiation
        staDir;    % normal solar direct radiation (W m-2)
        staDif;    % horizontal solar diffuse radiation (W m-2)
        staUdir;   % wind direction (º)
        staUmod;   % wind speed (m s-1)
        staRobs;   % Precipitation (mm h-1)
        staHum;    % specific humidty (kg kg-1)
    end
    
    methods
        function obj = Weather(climate_data,HI,HF)
             if(nargin > 0)
                obj.staTemp = csvread(climate_data,HI,6,[HI,6,HF,6]);
                obj.staRhum = csvread(climate_data,HI,8,[HI,8,HF,8]);
                obj.staPres = csvread(climate_data,HI,9,[HI,9,HF,9]);
                obj.staInfra = csvread(climate_data,HI,12,[HI,12,HF,12]);
                obj.staHor = csvread(climate_data,HI,13,[HI,13,HF,13]);
                obj.staDir = csvread(climate_data,HI,14,[HI,14,HF,14]);
                obj.staDif = csvread(climate_data,HI,15,[HI,15,HF,15]);
                obj.staUdir = csvread(climate_data,HI,20,[HI,20,HF,20]);
                obj.staUmod = csvread(climate_data,HI,21,[HI,21,HF,21]);
                obj.staRobs = csvread(climate_data,HI,33,[HI,33,HF,33]);
                obj.staHum = zeros(size(obj.staTemp,1),1);
                for i=1:size(obj.staTemp,1)
                  obj.staHum(i) = HumFromRHumTemp(obj.staRhum(i),...
                      obj.staTemp(i),obj.staPres(i));
                end
                obj.staTemp = obj.staTemp + 273.15;
             end
        end
    end
    
end

function W = HumFromRHumTemp(RH,T,P)

    % Saturation vapour pressure from ASHRAE
    C8=-5.8002206E3; C9=1.3914993; C10=-4.8640239E-2; C11=4.1764768E-5;
    C12=-1.4452093E-8; C13=6.5459673;
    T=T+273.15;
    PWS = exp(C8/T+C9+C10*T+C11*T^2+C12*T^3+C13*log(T));

    PW=RH*PWS/100;          % Vapour pressure
    W = 0.62198*PW/(P-PW);  % 4. Specific humidity

end
