"""Class for Weather data from EPW."""
from .utilities import read_csv, str2fl
from .psychrometrics import HumFromRHumTemp

try:
    range = xrange
except NameError:
    pass


class Weather(object):
    """Read EPW file[1] and store weather timeseries data.

    Args:
        climate_file: Text string for the name of the rural epw file that will be
            morphed.
        HI: Number for initial EPW row based on final Julian data.
        HF: Number for final EPW row based on start Julian data.

    Properties:
        * location
        * staTemp
        * staTdp
        * staRhum
        * staPres
        * staInfra
        * staHor
        * staDir
        * staDif
        * staUdir
        * staUmod
        * staRobs
        * staHum

    Note:
        [1] EPW CSV Format (In/Out). Retrieved August 26, 2020,
        from http://bigladdersoftware.com/epx/docs/8-2/auxiliary-programs/epw-csv-format-inout.html

    """

    def __init__(self, epw_path, HI, HF):
        try:
            self._climate_data = read_csv(epw_path)
        except Exception as e:
            raise Exception("Failed to read .epw file! {}".format(e.message))

        self.location = self._climate_data[0][1]
        cd = self._climate_data[HI:HF + 1]
        self.staTemp = str2fl([cd[i][6] for i in range(len(cd))])  # drybulb [C]
        self.staTdp = str2fl([cd[i][7] for i in range(len(cd))])  # dewpoint [C]
        self.staRhum = str2fl([cd[i][8] for i in range(len(cd))])  # air RH (%)
        self.staPres = str2fl([cd[i][9] for i in range(len(cd))])  # air pressure (Pa)
        # staInfra: horizontal Infrared Radiation Intensity (W m-2)
        self.staInfra = str2fl([cd[i][12] for i in range(len(cd))])
        # staHor: horizontal radiation [W m-2]
        self.staHor = str2fl([cd[i][13] for i in range(len(cd))])
        # staDir: normal solar direct radiation (W m-2)
        self.staDir = str2fl([cd[i][14] for i in range(len(cd))])
        # staDif: horizontal solar diffuse radiation (W m-2)
        self.staDif = str2fl([cd[i][15] for i in range(len(cd))])
        self.staUdir = str2fl([cd[i][20] for i in range(len(cd))])  # wind direction
        self.staUmod = str2fl([cd[i][21] for i in range(len(cd))])  # wind speed (m s-1)
        self.staHum = [0.0] * len(self.staTemp)  # specific humidty (kgH20 kgN202-1)
        for i in range(len(self.staTemp)):
            self.staHum[i] = \
                HumFromRHumTemp(self.staRhum[i], self.staTemp[i], self.staPres[i])
        self.staTemp = [s + 273.15 for s in self.staTemp]  # air temperature (K)

        # Set precipitation to array of zeros. This is done to avoid errors
        # for EPW files with incomplete or missing precipitation data. The
        # current implementation of the UWG doesn't use precipitation data
        # due to the difficulty in validating latent heat impact, the poor
        # quality precipitation data from weather sensors, and the marginal
        # impact it has on UHI.
        self.staRobs = [0.0] * len(self.staTemp)  # Precipitation (mm h-1)

    def __repr__(self):
        return 'Weather,\n City: {}\n, Max Tdb: {}C\n, Min Tdb: {}C\n'.format(
            self.location, round(max(self.staTemp) - 273.15, 2),
            round(min(self.staTemp) - 273.15, 2))
