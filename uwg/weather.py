from .utilities import read_csv, str2fl
from math import pow, log, exp
from .psychrometrics import HumFromRHumTemp

try:
    range = xrange
except NameError:
    pass


class Weather(object):
    """
    Weather
    Read epw file
    http://bigladdersoftware.com/epx/docs/8-2/auxiliary-programs/epw-csv-format-inout.html
    properties
        location  # location name
        staTemp   % air temperature (C)
        staTdp    % dewpoint temperature (C)
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
        cd = self.climate_data[HI:HF+1]
        self.staTemp = str2fl([cd[i][6] for i in range(len(cd))])           # drybulb [C]
        self.staTdp = str2fl([cd[i][7] for i in range(len(cd))])            # dewpoint [C]
        self.staRhum = str2fl([cd[i][8] for i in range(len(cd))])           # air relative humidity (%)
        self.staPres = str2fl([cd[i][9] for i in range(len(cd))])           # air pressure (Pa)
        self.staInfra = str2fl([cd[i][12] for i in range(len(cd))])         # horizontal Infrared Radiation Intensity (W m-2)
        self.staHor = str2fl([cd[i][13] for i in range(len(cd))])           # horizontal radiation [W m-2]
        self.staDir = str2fl([cd[i][14] for i in range(len(cd))])           # normal solar direct radiation (W m-2)
        self.staDif = str2fl([cd[i][15] for i in range(len(cd))])           # horizontal solar diffuse radiation (W m-2)
        self.staUdir = str2fl([cd[i][20] for i in range(len(cd))])          # wind direction ()
        self.staUmod = str2fl([cd[i][21] for i in range(len(cd))])          # wind speed (m s-1)
        self.staRobs = str2fl([cd[i][33] for i in range(len(cd))])          # Precipitation (mm h-1)
        self.staHum = [0.0] * len(self.staTemp)                                     # specific humidty (kgH20 kgN202-1)
        for i in range(len(self.staTemp)):
            self.staHum[i] = HumFromRHumTemp(self.staRhum[i], self.staTemp[i], self.staPres[i])

        self.staTemp = [s+273.15 for s in self.staTemp]                             # air temperature (K)

    def __repr__(self):
        return "Weather: City = {}, Max Tdb = {}C, Min Tdb = {}C".format(
            self.location,
            max(self.staTemp)-273.15,
            min(self.staTemp)-273.15
            )
