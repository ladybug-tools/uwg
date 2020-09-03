"""Class to estimate deep ground/water temperature."""
from __future__ import division


class Forcing (object):
    """Force method to estimate deep ground/water temperature.

    Args:
        staTemp: List of hourly air temperature for simulation period.
        weather: Weather object.

    Properties:
        * deepTemp
        * waterTemp
        * infra
        * uDir
        * hum
        * pres
        * temp
        * rHum
        * dir
        * dif
        * prec
        * wind
    """

    def __init__(self, staTemp=None, weather=None):

        # Define default values for instance variables

        # deepTemp: deep soil temperature (K)
        self.deepTemp = None
        # waterTemp: ground water temp, set to temp at 2m
        self.waterTemp = None
        # infra: horizontal Infrared Radiation Intensity (W m-2)
        self.infra = None
        # uDir: wind direction
        self.uDir = None
        # hum: specific humidty (kg kg-1)
        self.hum = None
        # pres: Pressure (Pa)
        self.pres = None
        # temp: air temperature (C)
        self.temp = None
        # rHum: Relative humidity (%)
        self.rHum = None
        # rHum: normal solar direct radiation (W m-2)
        self.dir = None
        # dif: horizontal solar diffuse radiation (W m-2)
        # ...Amount of solar radiation received from the sky
        # ...(excluding the solar disk) on a horizontal surface
        self.dif = None
        # prec: Precipitation (mm h-1)
        self.prec = None
        # wind: wind speed (m s-1)
        self.wind = None

        if staTemp and weather:
            self.deepTemp = sum(staTemp) / float(len(staTemp))
            self.waterTemp = sum(staTemp) / float(len(staTemp))
            self.infra = weather.staInfra
            self.uDir = weather.staUdir
            self.hum = weather.staHum
            self.pres = weather.staPres
            self.temp = weather.staTemp
            self.rHum = weather.staRhum
            self.dir = weather.staDir
            self.dif = weather.staDif
            self.prec = [p / 3.6e6 for p in weather.staRobs]
            self.wind = weather.staUmod

    def __repr__(self):
        deepT = int(self.deepTemp) if self.deepTemp else None
        waterT = int(self.waterTemp) if self.waterTemp else None
        return "Forcing: deepT={}, waterT={}".format(deepT, waterT)
