import pytest
from uwg import SimParam, Weather, Forcing
import os


EPW_PATH = \
    os.path.join(os.path.dirname(__file__), 'epw', 'SGP_Singapore.486980_IWEC.epw')


def test_forcing():
    """Test for forcing.py"""
    # setup_forcing
    dtSim = 300             # Sim time step
    dtWeather = 3600        # Weather data time-step
    MONTH = 1               # Begin month
    DAY = 1                 # Begin day of the month
    NUM_DAYS = 31           # Number of days of simulation

    simTime = SimParam(dtSim, dtWeather, MONTH, DAY, NUM_DAYS)
    print(EPW_PATH)
    weather = Weather(EPW_PATH, simTime.timeInitial, simTime.timeFinal)

    # initialized Forcing class
    forcIP = Forcing(weather.staTemp, weather)

    # Forcing tests
    assert forcIP.deepTemp == pytest.approx(299.8392473118278, abs=1e-12)
    assert forcIP.waterTemp == pytest.approx(299.8392473118278, abs=1e-12)
    assert forcIP.wind[0] == pytest.approx(3.2, abs=1e-10)
