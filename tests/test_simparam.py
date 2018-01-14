import pytest
import UWG
import os
import math
import pprint

class TestSimParam(object):

    DIR_UP_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    DIR_EPW_PATH = os.path.join(DIR_UP_PATH,"resources/epw")
    DIR_MATLAB_PATH = os.path.join(os.path.abspath(os.path.dirname(__file__)), "matlab_ref","matlab_simparam")

    def test_uwg_simparam_matlab_init(self):
        """ Matlab value comparison for simparam initialization in main uwg """

        # Initialize UWG from default initialize.uwg for 31 days simulation
        # Jan 1 00:00 - Feb 1 00:00

        # Month,1,        # starting month (1-12)
        # Day,1,          # starting day (1-31)
        # nDay,31,        # number of days to run simultion
        # dtSim,300,      # simulation time step (s)
        # dtWeather,3600, # weather time step (s)

        epw_dir = self.DIR_EPW_PATH
        epw_file_name = "SGP_Singapore.486980_IWEC.epw"
        uwg_param_dir = os.path.join(self.DIR_UP_PATH,"resources")
        uwg_param_file_name = "initialize.uwg"

        self.uwg = UWG.UWG(epw_dir, epw_file_name, uwg_param_dir, uwg_param_file_name)
        self.uwg.read_epw()
        self.uwg.read_input()
        self.uwg.set_input()

        # Open matlab ref file
        matlab_path = os.path.join(self.DIR_MATLAB_PATH,"matlab_ref_simparam_init.txt")
        if not os.path.exists(matlab_path):
            raise Exception("Failed to open {}!".format(matlab_path))
        matlab_file = open(matlab_path,'r')
        uwg_matlab_val = [float(x) for x in matlab_file.readlines()]
        matlab_file.close()

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
        for i in xrange(len(uwg_matlab_val)):
            #print uwg_python_val[i], uwg_matlab_val[i]
            assert uwg_python_val[i] == pytest.approx(uwg_matlab_val[i], abs=1e-15), "error at index={}".format(i)

    def test_uwg_simparam_matlab_update_date(self):
        """ Matlab value comparison for simparam UpdateDate function in main uwg """

        # Initialize UWG from initialize_simparam.uwg for 150 day simulation
        # Mar 15 00:00 - Aug 12th, 00:00

        # Month,3,        # starting month (1-12)
        # Day,15,         # starting day (1-31)
        # nDay,150,       # number of days to run simultion
        # dtSim,300,      # simulation time step (s)
        # dtWeather,3600, # weather time step (s)

        epw_dir = self.DIR_EPW_PATH
        epw_file_name = "SGP_Singapore.486980_IWEC.epw"
        uwg_param_dir = self.DIR_MATLAB_PATH
        uwg_param_file_name = "initialize_simparam.uwg"

        self.uwg = UWG.UWG(epw_dir, epw_file_name, uwg_param_dir, uwg_param_file_name)
        self.uwg.read_epw()
        self.uwg.read_input()
        self.uwg.set_input()
        self.uwg.hvac_autosize()
        self.uwg.uwg_main()

        # open matlab ref file
        matlab_path = os.path.join(self.DIR_MATLAB_PATH,"matlab_ref_simparam_update_date.txt")
        if not os.path.exists(matlab_path):
            raise Exception("Failed to open {}!".format(matlab_path))
        matlab_file = open(matlab_path,'r')
        uwg_matlab_val = [float(x) for x in matlab_file.readlines()]
        matlab_file.close()

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
        for i in xrange(len(uwg_matlab_val)):
            #print uwg_python_val[i], uwg_matlab_val[i]
            assert uwg_python_val[i] == pytest.approx(uwg_matlab_val[i], abs=1e-15), "error at index={}".format(i)


    def test_simparam(self):
        """ Tests simparam.py"""

        dtSim = 300             # Sim time step
        dtWeather = 3600        # Weather data time-step
        MONTH = 7               # Begin month
        DAY = 30                # Begin day of the month
        NUM_DAYS = 7            # Number of days of simulation

        simTime = UWG.SimParam(dtSim,dtWeather,MONTH,DAY,NUM_DAYS)

        # Simulation Parameters tests
        assert simTime.timeSim == pytest.approx(168, abs=1e-6)
        assert simTime.timeMax == pytest.approx(604800,abs=1e-6)
        assert simTime.nt == pytest.approx(2017,abs=1e-6)

        # Test UpdateDate() for < 1 hr
        for i in xrange(11): #11 * 300 = 3300 seconds = 55min
            simTime.UpdateDate()
        assert simTime.secDay == pytest.approx(3300., abs=1e-6)
        assert simTime.day == pytest.approx(30., abs=1e-6)
        assert simTime.hourDay == pytest.approx(0., abs=1e-6)
        # for == 1 hr
        simTime.UpdateDate()
        assert simTime.secDay == pytest.approx(3600., abs=1e-6)
        assert simTime.hourDay == pytest.approx(1., abs=1e-6)
        # for > 24hr
        for i in xrange(23 * 12):
            simTime.UpdateDate()
        assert simTime.secDay == pytest.approx(0., abs=1e-6)
        assert simTime.day == pytest.approx(31., abs=1e-6)
        assert simTime.hourDay == pytest.approx(0., abs=1e-6)

        # for == 1 month
        for i in xrange(24 * 12):
            simTime.UpdateDate()
        assert simTime.secDay == pytest.approx(0., abs=1e-6)
        assert simTime.day == pytest.approx(1., abs=1e-6)
        assert simTime.hourDay == pytest.approx(0., abs=1e-6)
        assert simTime.month == pytest.approx(8, abs=1e-6)

        # for + 1 month
        for i in xrange(24 * 12 * 31):
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
