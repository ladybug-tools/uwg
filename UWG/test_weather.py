import pytest
from simparam import SimParam
from weather import Weather


DIR_EPW_NAME = "data\\epw\\"

def setup_simparam():
    """setup simparam instance"""

    dtSim = 300             # Sim time step
    dtWeather = 3600        # Weather data time-step
    MONTH = 7               # Begin month
    DAY = 30                # Begin day of the month
    NUM_DAYS = 7            # Number of days of simulation

    return SimParam(dtSim,dtWeather,MONTH,DAY,NUM_DAYS)

def test_weather():
    """Test for weather.py"""

    simTime = setup_simparam()
    climate_file = DIR_EPW_NAME + "SGP_Singapore.486980_IWEC.epw"
    weather_ = Weather(climate_file,simTime.timeInitial,simTime.timeFinal)

    # Weather Tests
    assert len(weather_.staDif) == pytest.approx(simTime.timeFinal - simTime.timeInitial + 1, abs=1e-6)
    assert len(weather_.staHum) == pytest.approx(simTime.timeFinal - simTime.timeInitial + 1, abs=1e-6)
    assert len(weather_.staTemp) == pytest.approx(simTime.timeFinal - simTime.timeInitial + 1, abs=1e-6)
    assert weather_.staTemp[3] == pytest.approx(24.+273.15, abs=1e-6)
    assert weather_.staTemp[-1] == pytest.approx(27.+273.15, abs=1e-6)
    assert weather_.staUdir[2] == pytest.approx(270, abs=1e-1)  # 270 deg
    assert weather_.staUmod[4] == pytest.approx(.5, abs=1e-6)   # 0.5 m/s
    assert weather_.staPres[10] == pytest.approx(100600., abs=1e-1)
    assert weather_.staInfra[13] == pytest.approx(428., abs=1e-1)
    assert weather_.staDif[6] == pytest.approx(0., abs=1e-3)
    assert weather_.staDif[8] == pytest.approx(95., abs=1e-6)
    assert weather_.staRobs[8] == pytest.approx(0.0, abs=1e-3)  # 0. mm/hre

if __name__ == '__main__':
    test_weather()
