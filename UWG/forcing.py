class Forcing (object):
    """
    FORCING
    Force Restore method to estimate ground soil temperature?? -sv

    """

    def __init__(self,staTemp,weather):
        # weather: Weather obj
        # staTemp: list of hourly air temperature for simulation period
        self.deepTemp = sum(staTemp)/float(len(staTemp))    # deep soil temperature (K)
        self.waterTemp = sum(staTemp)/float(len(staTemp))   # ground water temp, set to temp at 2m
        self.infra = weather.staInfra     # horizontal Infrared Radiation Intensity (W m-2)
        self.uDir = weather.staUdir       # wind direction
        self.hum  = weather.staHum        # specific humidty (kg kg-1)
        self.pres = weather.staPres       # Pressure (Pa)
        self.temp = weather.staTemp       # air temperature (C)
        self.rHum = weather.staRhum       # Relative humidity (%)
        self.dir = weather.staDir         #  normal solar direct radiation (W m-2)
        self.dif = weather.staDif         # horizontal solar diffuse radiation (W m-2)
        self.prec = map(lambda p: p/3.6e6, weather.staRobs) #  Precipitation (mm h-1)
        self.wind = weather.staUmod       # wind speed (m s-1)

    def __repr__(self):
        return "Forcing: deepT={a}, waterT={b}".format(
            a=int(self.deepTemp),
            b=int(self.waterTemp)
            )
