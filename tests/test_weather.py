"""Test for Weather class."""

import pytest
from uwg import SimParam, Weather
from .test_base import DEFAULT_EPW_PATH


def test_weather():
    """Test for weather.py"""

    dtSim = 300             # Sim time step
    dtWeather = 3600        # Weather data time-step
    MONTH = 7               # Begin month
    DAY = 30                # Begin day of the month
    NUM_DAYS = 7            # Number of days of simulation

    simTime = SimParam(dtSim, dtWeather, MONTH, DAY, NUM_DAYS)
    weather = Weather(DEFAULT_EPW_PATH, simTime.timeInitial, simTime.timeFinal)

    # Weather Tests
    assert len(weather.staDif) == \
        pytest.approx(simTime.timeFinal - simTime.timeInitial + 1, abs=1e-6)
    assert len(weather.staHum) == \
        pytest.approx(simTime.timeFinal - simTime.timeInitial + 1, abs=1e-6)
    assert len(weather.staTemp) == \
        pytest.approx(simTime.timeFinal - simTime.timeInitial + 1, abs=1e-6)
    assert weather.staTemp[3] == pytest.approx(24.+273.15, abs=1e-6)
    assert weather.staTemp[-1] == pytest.approx(27.+273.15, abs=1e-6)
    assert weather.staUdir[2] == pytest.approx(270, abs=1e-1)  # 270 deg
    assert weather.staUmod[4] == pytest.approx(.5, abs=1e-6)   # 0.5 m/s
    assert weather.staPres[10] == pytest.approx(100600., abs=1e-1)
    assert weather.staInfra[13] == pytest.approx(428., abs=1e-1)
    assert weather.staDif[6] == pytest.approx(0., abs=1e-3)
    assert weather.staDif[8] == pytest.approx(95., abs=1e-6)
    assert weather.staRobs[8] == pytest.approx(0.0, abs=1e-3)  # 0. mm/hre
