try:
    range = xrange
except NameError:
    pass

import pytest
import uwg
import os
import math
import pprint
from .test_base import TestBase


class TestSimParam(TestBase):

    def test_uwg_simparam_matlab_init(self):
        """ Matlab value comparison for simparam initialization in main uwg """

        # Initialize uwg from default initialize.uwg for 31 days simulation
        # Jan 1 00:00 - Feb 1 00:00

        # Month,1,        # starting month (1-12)
        # Day,1,          # starting day (1-31)
        # nDay,31,        # number of days to run simultion
        # dtSim,300,      # simulation time step (s)
        # dtWeather,3600, # weather time step (s)

        self.setup_uwg_integration(uwg_param_file="initialize_singapore.uwg")

        self.uwg.read_epw()
        self.uwg.set_input()
        self.uwg.init_BEM_obj()
        self.uwg.init_input_obj()
        # open matlab ref file
        uwg_matlab_val = self.setup_open_matlab_ref("matlab_simparam","matlab_ref_simparam_init.txt")

        uwg_python_val = [
        self.uwg.simTime.dt,            # uwg time simulation time step
        self.uwg.simTime.timeForcing,   # weather data timestep
        self.uwg.simTime.month,
        self.uwg.simTime.day,
        self.uwg.simTime.days,
        self.uwg.simTime.timePrint,     # weather data timestep
        self.uwg.simTime.timeDay,       # how many times weather senses in a day
        self.uwg.simTime.timeSim,       # how many steps in weather data simulation
        self.uwg.simTime.timeMax,       # total seconds in simulation days
        self.uwg.simTime.nt,            # total number of timesteps for uwg simuation
        self.uwg.simTime.julian,
        self.uwg.simTime.timeInitial,   # sensor data in epw for intial time based on julian day & timesteps
        self.uwg.simTime.timeFinal,     # sensor data in epw for final time based on julian day & timesteps
        self.uwg.simTime.secDay,        # current seconds in day
        self.uwg.simTime.hourDay
        ]

        # matlab ref checking
        assert len(uwg_matlab_val) == len(uwg_python_val)
        for i in range(len(uwg_matlab_val)):
            #print uwg_python_val[i], uwg_matlab_val[i]
            assert uwg_python_val[i] == pytest.approx(uwg_matlab_val[i], abs=1e-15), "error at index={}".format(i)

    def test_uwg_simparam_matlab_update_date(self):
        """ Matlab value comparison for simparam UpdateDate function in main uwg """

        # Initialize uwg from initialize_simparam.uwg for 150 day simulation
        # Mar 15 00:00 - Aug 12th, 00:00

        # Month,3,        # starting month (1-12)
        # Day,15,         # starting day (1-31)
        # nDay,150,       # number of days to run simultion
        # dtSim,300,      # simulation time step (s)
        # dtWeather,3600, # weather time step (s)

        self.setup_uwg_integration(uwg_param_file="initialize_simparam.uwg")

        self.uwg.read_epw()
        self.uwg.set_input()
        self.uwg.init_BEM_obj()
        self.uwg.init_input_obj()
        self.uwg.hvac_autosize()
        self.uwg.simulate()

        # open matlab ref file
        uwg_matlab_val = self.setup_open_matlab_ref("matlab_simparam","matlab_ref_simparam_update_date.txt")

        uwg_python_val = [
        self.uwg.simTime.secDay,
        self.uwg.simTime.day,
        self.uwg.simTime.julian,
        self.uwg.simTime.month,
        self.uwg.simTime.day,
        self.uwg.simTime.hourDay,
        self.uwg.ceil_time_step+1 # Add 1 to keep consistent with matlab list convenction
        ]

        # matlab ref checking
        assert len(uwg_matlab_val) == len(uwg_python_val)
        for i in range(len(uwg_matlab_val)):
            #print uwg_python_val[i], uwg_matlab_val[i]
            assert uwg_python_val[i] == pytest.approx(uwg_matlab_val[i], abs=1e-15), "error at index={}".format(i)


    def test_simparam(self):
        """ Tests simparam.py"""

        dtSim = 300             # Sim time step
        dtWeather = 3600        # Weather data time-step
        MONTH = 7               # Begin month
        DAY = 30                # Begin day of the month
        NUM_DAYS = 7            # Number of days of simulation

        simTime = uwg.SimParam(dtSim,dtWeather,MONTH,DAY,NUM_DAYS)

        # Simulation Parameters tests
        assert simTime.timeSim == pytest.approx(168, abs=1e-6)
        assert simTime.timeMax == pytest.approx(604800,abs=1e-6)
        assert simTime.nt == pytest.approx(2017,abs=1e-6)

        # Test UpdateDate() for < 1 hr
        for i in range(11): #11 * 300 = 3300 seconds = 55min
            simTime.UpdateDate()
        assert simTime.secDay == pytest.approx(3300., abs=1e-6)
        assert simTime.day == pytest.approx(30., abs=1e-6)
        assert simTime.hourDay == pytest.approx(0., abs=1e-6)
        # for == 1 hr
        simTime.UpdateDate()
        assert simTime.secDay == pytest.approx(3600., abs=1e-6)
        assert simTime.hourDay == pytest.approx(1., abs=1e-6)
        # for > 24hr
        for i in range(23 * 12):
            simTime.UpdateDate()
        assert simTime.secDay == pytest.approx(0., abs=1e-6)
        assert simTime.day == pytest.approx(31., abs=1e-6)
        assert simTime.hourDay == pytest.approx(0., abs=1e-6)

        # for == 1 month
        for i in range(24 * 12):
            simTime.UpdateDate()
        assert simTime.secDay == pytest.approx(0., abs=1e-6)
        assert simTime.day == pytest.approx(1., abs=1e-6)
        assert simTime.hourDay == pytest.approx(0., abs=1e-6)
        assert simTime.month == pytest.approx(8, abs=1e-6)

        # for + 1 month
        for i in range(24 * 12 * 31):
            simTime.UpdateDate()
        assert simTime.secDay == pytest.approx(0., abs=1e-6)
        assert simTime.day == pytest.approx(1., abs=1e-6)
        assert simTime.hourDay == pytest.approx(0., abs=1e-6)
        assert simTime.month == pytest.approx(9, abs=1e-6)




if __name__ == "__main__":

    tsp = TestSimParam()
    tsp.test_simparam()
    tsp.test_uwg_simparam_matlab_init()
    tsp.test_uwg_simparam_matlab_update_date()
