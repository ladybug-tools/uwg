from utilities import read_csv

class Weather(object):
    """
    Weather
    Read a epw file

    properties
        staTemp   % air temperature (C)
        staRhum   % air relative humidity (%)
        staPres   % air pressure (Pa)
        staInfra  % horizontal Infrared Radiation Intensity (W m-2)
        staHor    % horizontal radiation
        staDir    % normal solar direct radiation (W m-2)
        staDif    % horizontal solar diffuse radiation (W m-2)
        staUdir   % wind direction ()
        staUmod   % wind speed (m s-1)
        staRobs   % Precipitation (mm h-1)
        staHum    % specific humidty (kg kg-1)
    """

    DIR_EPW_NAME = "data\\epw\\"
    def __init__(self,climate_file,HI,HF):
        self.climate_data = read_csv(self.DIR_EPW_NAME + climate_file)
        self.location = self.climate_data[0][1]
        """
        self.staTemp = csvread(climate_data,HI,6,[HI,6,HF,6])
        self.staRhum = csvread(climate_data,HI,8,[HI,8,HF,8])
        self.staPres = csvread(climate_data,HI,9,[HI,9,HF,9])
        self.staInfra = csvread(climate_data,HI,12,[HI,12,HF,12])
        self.staHor = csvread(climate_data,HI,13,[HI,13,HF,13])
        self.staDir = csvread(climate_data,HI,14,[HI,14,HF,14])
        self.staDif = csvread(climate_data,HI,15,[HI,15,HF,15])
        self.staUdir = csvread(climate_data,HI,20,[HI,20,HF,20])
        self.staUmod = csvread(climate_data,HI,21,[HI,21,HF,21])
        self.staRobs = csvread(climate_data,HI,33,[HI,33,HF,33])
        self.staHum = zeros(size(self.staTemp,1),1)
        for i=1:size(self.staTemp,1)
          self.staHum(i) = HumFromRHumTemp(self.staRhum(i),...
              self.staTemp(i),self.staPres(i))
        end
        self.staTemp = self.staTemp + 273.15
        """
    def __repr__(self):
        return "Weather @ {a}".format(a=self.location)
"""
function W = HumFromRHumTemp(RH,T,P)

    % Saturation vapour pressure from ASHRAE
    C8=-5.8002206E3 C9=1.3914993 C10=-4.8640239E-2 C11=4.1764768E-5
    C12=-1.4452093E-8 C13=6.5459673
    T=T+273.15
    PWS = exp(C8/T+C9+C10*T+C11*T^2+C12*T^3+C13*log(T))

    PW=RH*PWS/100          % Vapour pressure
    W = 0.62198*PW/(P-PW)  % 4. Specific humidity

 """
