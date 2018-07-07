import pytest
import uwg
import os



class TestForcing(object):

    DIR_UP_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    DIR_EPW_PATH = os.path.join(DIR_UP_PATH,"resources/epw")

    def setup_forcing(self):
        """setup simparam and weather instance"""

        # From initialize.uwg referenece
        dtSim = 300             # Sim time step
        dtWeather = 3600        # Weather data time-step
        MONTH = 1               # Begin month
        DAY = 1                 # Begin day of the month
        NUM_DAYS = 31           # Number of days of simulation

        simTime = uwg.SimParam(dtSim,dtWeather,MONTH,DAY,NUM_DAYS)
        epw_name = "SGP_Singapore.486980_IWEC.epw"
        climate_file = os.path.join(self.DIR_EPW_PATH, epw_name)
        self.weather = uwg.Weather(climate_file, simTime.timeInitial, simTime.timeFinal)


    def test_forcing(self):
        """Test for forcing.py"""
        # setup
        self.setup_forcing()
        self.forcIP = uwg.Forcing(self.weather.staTemp, self.weather) # initialized Forcing class

        # Forcing tests
        assert self.forcIP.deepTemp == pytest.approx(299.8392473118278, abs=1e-12)
        assert self.forcIP.waterTemp == pytest.approx(299.8392473118278, abs=1e-12)
        assert self.forcIP.wind[0] == pytest.approx(3.2, abs=1e-10)
