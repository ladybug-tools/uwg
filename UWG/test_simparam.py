import pytest
from simparam import SimParam



def test_simparam():
    """Test for simparam.py"""

    dtSim = 300             # Sim time step
    dtWeather = 3600        # Weather data time-step
    MONTH = 7               # Begin month
    DAY = 30                # Begin day of the month
    NUM_DAYS = 7            # Number of days of simulation

    simTime = SimParam(dtSim,dtWeather,MONTH,DAY,NUM_DAYS)

    # Simulation Parameters tests
    assert pytest.approx(simTime.timeSim - 168, abs=1e-6)
    assert pytest.approx(simTime.timeMax - 604800,abs=1e-6)
    assert pytest.approx(simTime.timeMax - 2017,abs=1e-6)

if __name__ == '__main__':
    test_simparam()
