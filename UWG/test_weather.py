import pytest
from simparam import SimParam
from weather import Weather


DIR_EPW_NAME = "UWG\\data\\epw\\"

def setup_simparam():
    """setup by getting simparam instance"""

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
    pytest.approx(len(weather_.staDif),simTime.timeFinal - simTime.timeInitial + 1, abs=1e-6)
    pytest.approx(len(weather_.staHum),simTime.timeFinal - simTime.timeInitial + 1, abs=1e-6)
    pytest.approx(len(weather_.staTemp),simTime.timeFinal - simTime.timeInitial + 1, abs=1e-6)
    pytest.approx(weather_.staTemp[3],24.+273.15, abs=1e-6)
    pytest.approx(weather_.staTemp[-1],27.+273.15, abs=1e-6)
    pytest.approx(weather_.staUdir[2],270, abs=1e-1)  # 270 deg
    pytest.approx(weather_.staUmod[4],.5, abs=1e-6)   # 0.5 m/s
    pytest.approx(weather_.staPres[10],100600., abs=1e-1)
    pytest.approx(weather_.staInfra[13],428., abs=1e-1)
    pytest.approx(weather_.staDif[6],0., abs=1e-3)
    pytest.approx(weather_.staDif[8],95., abs=1e-6)
    pytest.approx(weather_.staRobs[8],0.0, abs=1e-3)  # 0. mm/hre

if __name__ == '__main__':
    test_weather()
