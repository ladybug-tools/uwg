from utilities import read_csv, str2fl
from math import pow, log, exp
from psychrometrics import HumFromRHumTemp

class Weather(object):
    """
    Weather
    Read epw file
    http://bigladdersoftware.com/epx/docs/8-2/auxiliary-programs/epw-csv-format-inout.html
    properties
        location  # location name
        staTemp   % air temperature (C)
<<<<<<< HEAD
        staTdp    % dry bulb
=======
        staTdp    % dewpoint temperature (C)
>>>>>>> a923c1cef10d46cd50002b5a3994904ab34f528a
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

    #TODO: change to xrange
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
        self.staTemp = str2fl([r[6] for r in self.climate_data[HI:HF+1]])           # drybulb [C]
        self.staTdp = str2fl([r[7] for r in self.climate_data[HI:HF+1]])            # dewpoint [C]
        self.staRhum = str2fl([r[8] for r in self.climate_data[HI:HF+1]])           # air relative humidity (%)
        self.staPres = str2fl([r[9] for r in self.climate_data[HI:HF+1]])           # air pressure (Pa)
        self.staInfra = str2fl([r[12] for r in self.climate_data[HI:HF+1]])         # horizontal Infrared Radiation Intensity (W m-2)
        self.staHor = str2fl([r[13] for r in self.climate_data[HI:HF+1]])           # horizontal radiation [W m-2]
        self.staDir = str2fl([r[14] for r in self.climate_data[HI:HF+1]])           # normal solar direct radiation (W m-2)
        self.staDif = str2fl([r[15] for r in self.climate_data[HI:HF+1]])           # horizontal solar diffuse radiation (W m-2)
        self.staUdir = str2fl([r[20] for r in self.climate_data[HI:HF+1]])          # wind direction ()
        self.staUmod = str2fl([r[21] for r in self.climate_data[HI:HF+1]])          # wind speed (m s-1)
        self.staRobs = str2fl([r[33] for r in self.climate_data[HI:HF+1]])          # Precipitation (mm h-1)
        self.staHum = [0.0] * len(self.staTemp)                                     # specific humidty (kgH20 kgN202-1)
        for i in xrange(len(self.staTemp)):
            self.staHum[i] = HumFromRHumTemp(self.staRhum[i], self.staTemp[i], self.staPres[i])

        self.staTemp = [s+273.15 for s in self.staTemp]                             # air temperature (K)

    def __repr__(self):
        return "Weather: {a}, HI Tdb:{b}, HF Tdb:{c}".format(
            a=self.location,
            b=self.staTemp[0]-273.15,
            c=self.staTemp[-1]-273.15
            )
