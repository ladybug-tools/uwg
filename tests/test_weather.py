import pytest
import os

import uwg


class TestWeather(object):

    DIR_UP_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    DIR_EPW_PATH = os.path.join(DIR_UP_PATH,"resources/epw")

    def setup(self):
        """setup simparam instance"""

        dtSim = 300             # Sim time step
        dtWeather = 3600        # Weather data time-step
        MONTH = 7               # Begin month
        DAY = 30                # Begin day of the month
        NUM_DAYS = 7            # Number of days of simulation

        self.simTime = uwg.SimParam(dtSim,dtWeather,MONTH,DAY,NUM_DAYS)

    def test_weather(self):
        """Test for weather.py"""

        epw_name = "SGP_Singapore.486980_IWEC.epw"
        climate_file = os.path.join(self.DIR_EPW_PATH, epw_name)


        self.weather = uwg.Weather(climate_file,self.simTime.timeInitial,self.simTime.timeFinal)

        # Weather Tests
        assert len(self.weather.staDif) == pytest.approx(self.simTime.timeFinal - self.simTime.timeInitial + 1, abs=1e-6)
        assert len(self.weather.staHum) == pytest.approx(self.simTime.timeFinal - self.simTime.timeInitial + 1, abs=1e-6)
        assert len(self.weather.staTemp) == pytest.approx(self.simTime.timeFinal - self.simTime.timeInitial + 1, abs=1e-6)
        assert self.weather.staTemp[3] == pytest.approx(24.+273.15, abs=1e-6)
        assert self.weather.staTemp[-1] == pytest.approx(27.+273.15, abs=1e-6)
        assert self.weather.staUdir[2] == pytest.approx(270, abs=1e-1)  # 270 deg
        assert self.weather.staUmod[4] == pytest.approx(.5, abs=1e-6)   # 0.5 m/s
        assert self.weather.staPres[10] == pytest.approx(100600., abs=1e-1)
        assert self.weather.staInfra[13] == pytest.approx(428., abs=1e-1)
        assert self.weather.staDif[6] == pytest.approx(0., abs=1e-3)
        assert self.weather.staDif[8] == pytest.approx(95., abs=1e-6)
        assert self.weather.staRobs[8] == pytest.approx(0.0, abs=1e-3)  # 0. mm/hre
