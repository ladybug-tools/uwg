from __future__ import division, print_function

try:
    range = xrange
except NameError:
    pass

import math


class SimParam(object):
    """
    SimParam
    Calculates simulation time parameters based on initial date
    and weather data time step, and simulation timestep

    properties
        dt            % Simulation time-step
        timeForcing   % Weather data time-step
        month         % Begin month
        day           % Begin day of the month
        days          % Number of days of simulation
        timePrint     % time-step for printing outputs
        timeDay       % number of timesteps in a design-day
        timeSim       % number of timesteps of the simulation
        timeMax       % total simulation time (s)
        nt            % total number of timesteps
        timeFinal     % final timestep of simulation
        timeInitial   % initial timestep of simulation
        secDay        % seconds of one day (s)
        hourDay       % hour of the day (0 - 23hr)
        inobis        % julian day at the end of each month
        julian        % julian day
    """
    TIMESTEP_CONFLICT_MSG = "TIMESTEP ERROR! Timestep must be a factor of 3600."

    def __init__(self,dt,timefor,M,DAY,days):
        self.dt = dt                                                        # uwg time simulation time step
        self.timeForcing = timefor                                          # weather data timestep
        self.month = int(M)
        self.day = DAY
        self.days = days
        self.timePrint = timefor                                            # weather data timestep
        self.timeDay = 24*3600/timefor                                      # how many times weather senses in a day
        self.timeSim = self.timeDay*days                                    # how many steps in weather data simulation
        self.timeMax = 24.*3600.*days                                       # total seconds in simulation days
        self.nt = int(round(self.timeMax/self.dt+1))                        # total number of timesteps for uwg simuation
        self.inobis = [0,31,59,90,120,151,181,212,243,273,304,334]
        self.julian = self.inobis[self.month - 1] + DAY - 1
        #H1: (julian day * number of timesteps in a day) == sensor data index in epw
        H1 = int((self.inobis[self.month - 1] + DAY - 1) * self.timeDay)
        self.timeInitial = H1 + 8                                           # sensor data in epw for intial time based on julian day & timesteps
        self.timeFinal = int(H1 + self.timeDay * self.days - 1 + 8)         # sensor data in epw for final time based on julian day & timesteps
        self.secDay = 0                                                     # current seconds in day
        self.hourDay = 0                                                    # current hours in day

    def __repr__(self):
        return "SimParam: Start = {}/{}, days = {}, timestep = {}s".format(
            self.month,
            int(self.day),
            self.days,
            self.dt
            )

    def is_near_zero(self,num,eps=1e-10):
        return abs(float(num)) < eps

    def UpdateDate(self):
        self.secDay = self.secDay + self.dt

        if self.is_near_zero(self.secDay - 3600*24):
            self.day += 1
            self.julian = self.julian + 1
            self.secDay = 0.
            for j in range(12):
                if self.is_near_zero(self.julian - self.inobis[j]):
                    self.month = self.month + 1
                    self.day = 1

        if self.secDay > (3600*24):
            raise Exception("{}. CURRENTLY AT {}.".format(self.TIMESTEP_CONFLICT_MSG, self.dt))

        self.hourDay = int(math.floor(self.secDay/3600.))       # 0 - 23hr
