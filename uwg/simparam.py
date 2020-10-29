"""Class for simulation time parameters."""
from __future__ import division, print_function

try:
    range = xrange
except NameError:
    pass

import math
from .utilities import is_near_zero


class SimParam(object):
    """Calculates simulation time parameters.

    Simulation time parameters are calculated based on initial date and weather
    data time step, and simulation timestep.

    Args:
        dt: Number for simulation time step in seconds.
        timefor: Number for weather data time-step in seconds.
        M: Number between 1 and 12 indicating the simulation start month.
        DAY: Number between 1 and 31 indicating the simulation start day.
        days: Number indicating number of days to simulate.

    Properties:
        * dt -- uwg time simulation time step
        * timeForcing -- weather data timestep
        * month -- start month
        * day -- start day
        * days -- number of days in simulation
        * timePrint -- weather data timestep
        * timeDay -- number of times weather senses in a day
        * timeSim -- number of steps in weather data simulation
        * timeMax -- total seconds in simulation days
        * nt -- total number of timesteps
        * inobis -- list of julian dates for the start of the months
        * julian -- simulation start julian date
        * timeInitial -- epw sensor data for initial time based on julian day & timesteps
        * timeFinal -- epw sensor data for final time based on julian day & timesteps
        * secDay -- current seconds in day
        * hourDay -- current hours in day
    """
    TIMESTEP_CONFLICT_MSG = "TIMESTEP ERROR! Timestep must be a factor of 3600."

    def __init__(self, dt, timefor, M, DAY, days):
        self.dt = dt
        self.timeForcing = timefor
        self.month = int(M)
        self.day = DAY
        self.days = days
        self.timePrint = timefor
        self.timeDay = 24 * 3600 / timefor
        self.timeSim = self.timeDay * days
        self.timeMax = 24. * 3600. * days
        self.nt = int(round(self.timeMax / self.dt + 1))
        self.inobis = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]
        self.julian = self.inobis[self.month - 1] + DAY - 1
        # H1: (julian day * number of timesteps in a day) == sensor data index in epw
        H1 = int((self.inobis[self.month - 1] + DAY - 1) * self.timeDay)
        self.timeInitial = H1 + 8
        self.timeFinal = int(H1 + self.timeDay * self.days - 1 + 8)
        self.secDay = 0
        self.hourDay = 0

    def update_date(self):
        self.secDay = self.secDay + self.dt

        if is_near_zero(self.secDay - (3600 * 24)):
            self.day += 1
            self.julian = self.julian + 1
            self.secDay = 0.
            for j in range(12):
                if is_near_zero(self.julian - self.inobis[j]):
                    self.month = self.month + 1
                    self.day = 1

        if self.secDay > (3600 * 24):
            raise Exception("{}. CURRENTLY AT {}.".format(
                self.TIMESTEP_CONFLICT_MSG, self.dt))

        self.hourDay = int(math.floor(self.secDay / 3600.))  # 0 - 23hr

    def __repr__(self):
        return "SimParam: Start = {}/{}, days = {}, timestep = {}s".format(
            self.month, int(self.day), self.days, self.dt)
