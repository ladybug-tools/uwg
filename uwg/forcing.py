"""Class to estimate deep ground/water temperature."""
from __future__ import division


class Forcing (object):
    """Force method to estimate deep ground/water temperature.

    Args:
        staTemp: List of hourly air temperature for simulation period.
        weather: Weather object.

    Properties:
        * deepTemp -- deep soil temperature (K)
        * waterTemp -- ground water temp, set to temp at 2m
        * infra -- horizontal Infrared Radiation Intensity (W m-2)
        * uDir -- wind direction
        * hum -- specific humidty (kg kg-1)
        * pres -- Pressure (Pa)
        * temp -- air temperature (C)
        * rHum -- Relative humidity (%)
        * dir -- normal solar direct radiation (W m-2)
        * dif -- horizontal solar diffuse radiation (W m-2)
        * prec -- precipitation (mm h-1)
        * wind -- wind speed (m s-1)
    """

    def __init__(self, staTemp=None, weather=None):

        if not (staTemp and weather):
            # Define default values for instance variables
            self.deepTemp = None
            self.waterTemp = None
            self.infra = None
            self.uDir = None
            self.hum = None
            self.pres = None
            self.temp = None
            self.rHum = None
            self.dir = None
            # dif: horizontal solar diffuse radiation (W m-2)
            # ...Amount of solar radiation received from the sky
            # ...(excluding the solar disk) on a horizontal surface
            self.dif = None
            self.prec = None
            self.wind = None
        else:
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
        return 'forcing,\n deepTemp: {}\n waterTemp: {}'.format(
            self.deepTemp, self.waterTemp)
