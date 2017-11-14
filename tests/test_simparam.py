import pytest
import UWG


def test_simparam():
    """Test for simparam.py"""

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
    test_simparam()
