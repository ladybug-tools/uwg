from __future__ import division


class Forcing (object):
    """
    FORCING
    Force method to estimate deep ground/water temperature

    self.deepTemp                # deep soil temperature (K)
    self.waterTemp               # ground water temp, set to temp at 2m
    self.infra                   # horizontal Infrared Radiation Intensity (W m-2)
    self.uDir                    # wind direction
    self.hum                     # specific humidty (kg kg-1)
    self.pres                    # Pressure (Pa)
    self.temp                    # air temperature (C)
    self.rHum                    # Relative humidity (%)
    self.dir                     # normal solar direct radiation (W m-2)
    self.dif                     # horizontal solar diffuse radiation (W m-2)
                                 # ...Amount of solar radiation received from the sky
                                 # ...(excluding the solar disk) on a horizontal surface
    self.prec                    # Precipitation (mm h-1)
    self.winde                   # wind speed (m s-1)

    """

    def __init__(self,staTemp=None,weather=None):
        # weather: Weather obj
        # staTemp: list of hourly air temperature for simulation period
        # Define default values for instance variables, when the type can be mutable:
        if staTemp==None and weather==None:
            self.deepTemp = None
            self.waterTemp = None
            self.infra = None
            self.uDir = None
            self.hum  = None
            self.pres = None
            self.temp = None
            self.rHum = None
            self.dir = None
            self.dif = None
            self.prec = None
            self.wind = None
        else:
            self.deepTemp = sum(staTemp)/float(len(staTemp))
            self.waterTemp = sum(staTemp)/float(len(staTemp))
            self.infra = weather.staInfra
            self.uDir = weather.staUdir
            self.hum  = weather.staHum
            self.pres = weather.staPres
            self.temp = weather.staTemp
            self.rHum = weather.staRhum
            self.dir = weather.staDir
            self.dif = weather.staDif
            self.prec = [p/3.6e6 for p in weather.staRobs]
            self.wind = weather.staUmod

    def __repr__(self):
        return "Forcing: deepT={a}, waterT={b}".format(
            a=int(self.deepTemp) if self.deepTemp else None,
            b=int(self.waterTemp) if self.waterTemp else None
            )
