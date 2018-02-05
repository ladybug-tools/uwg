from utilities import read_csv, str2fl
from math import pow, log, exp

class Weather(object):
    """
    Weather
    Read epw file
    http://bigladdersoftware.com/epx/docs/8-2/auxiliary-programs/epw-csv-format-inout.html
    properties
        location  # location name
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

    def __init__(self,climate_file,HI,HF):
        #HI: Julian start date
        #HF: Julian final date
        #H1 and HF define the row we want

        # Open .epw file and feed csv data to self.climate_data
        try:
            self.climate_data = read_csv(climate_file)
        except Exception as e:
            raise Exception("Failed to read .epw file! {}".format(e.message))

        self.location = self.climate_data[0][1]
        self.staTemp = str2fl(map(lambda r: r[6], self.climate_data[HI:HF+1]))      # drybulb [C]
        self.staRhum = str2fl(map(lambda r: r[8], self.climate_data[HI:HF+1]))      # air relative humidity (%)
        self.staPres = str2fl(map(lambda r: r[9], self.climate_data[HI:HF+1]))      # air pressure (Pa)
        self.staInfra = str2fl(map(lambda r: r[12], self.climate_data[HI:HF+1]))    # horizontal Infrared Radiation Intensity (W m-2)
        self.staHor = str2fl(map(lambda r: r[13], self.climate_data[HI:HF+1]))      # horizontal radiation [W m-2]
        self.staDir = str2fl(map(lambda r: r[14], self.climate_data[HI:HF+1]))      # normal solar direct radiation (W m-2)
        self.staDif = str2fl(map(lambda r: r[15], self.climate_data[HI:HF+1]))      # horizontal solar diffuse radiation (W m-2)
        self.staUdir = str2fl(map(lambda r: r[20], self.climate_data[HI:HF+1]))     # wind direction ()
        self.staUmod = str2fl(map(lambda r: r[21], self.climate_data[HI:HF+1]))     # wind speed (m s-1)
        self.staRobs = str2fl(map(lambda r: r[33], self.climate_data[HI:HF+1]))     # Precipitation (mm h-1)
        self.staHum = [0.0] * len(self.staTemp)                                     # specific humidty (kgH20 kgN202-1)
        for i in xrange(len(self.staTemp)):
            self.staHum[i] = self.HumFromRHumTemp(self.staRhum[i], self.staTemp[i], self.staPres[i])

        self.staTemp = map(lambda s: s+273.15, self.staTemp)                             # air temperature (K)

    def __repr__(self):
        return "Weather: {a}, HI Tdb:{b}, HF Tdb:{c}".format(
            a=self.location,
            b=self.staTemp[0]-273.15,
            c=self.staTemp[-1]-273.15
            )

    def HumFromRHumTemp(self,RH,T,P):
        # Derive Specific HUmidity [kgh20/kgn202] from RH, T and Pa
        # Saturation vapour pressure from ASHRAE
        C8 = -5.8002206e3
        C9 = 1.3914993
        C10 = -4.8640239e-2
        C11 = 4.1764768e-5
        C12 = -1.4452093e-8
        C13 = 6.5459673

        T += 273.15

        PWS = exp(C8/T + C9 + C10*T + C11 * pow(T,2) + C12 * pow(T,3) + C13 * log(T))
        PW = RH*PWS/100.0        # Vapour pressure
        W = 0.62198*PW/(P-PW)    # 4. Specific humidity
        return W
