class Forcing (object):
    """
    FORCING
    Force method to estimate deep ground/water temperature

    """

    def __init__(self,staTemp=None,weather=None):
        # weather: Weather obj
        # staTemp: list of hourly air temperature for simulation period
        
        # Define default values for instance variables, when the type can be mutable:
        if staTemp==None and weather==None:
            self.deepTemp = None                # deep soil temperature (K)
            self.waterTemp = None               # ground water temp, set to temp at 2m
            self.infra = None                   # horizontal Infrared Radiation Intensity (W m-2)
            self.uDir = None                    # wind direction
            self.hum  = None                    # specific humidty (kg kg-1)
            self.pres = None                    # Pressure (Pa)
            self.temp = None                    # air temperature (C)
            self.rHum = None                    # Relative humidity (%)
            self.dir = None                     # normal solar direct radiation (W m-2)
            self.dif = None                     # horizontal solar diffuse radiation (W m-2)
            self.prec = None                    # Precipitation (mm h-1)
            self.wind = None                    # wind speed (m s-1)
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
            self.prec = map(lambda p: p/3.6e6, weather.staRobs)
            self.wind = weather.staUmod

    def __repr__(self):
        return "Forcing: deepT={a}, waterT={b}".format(
            a=int(self.deepTemp) if self.deepTemp else None,
            b=int(self.waterTemp) if self.waterTemp else None
            )
