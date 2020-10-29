"""Test SimParam class."""

import os
import pytest
from .test_base import auto_setup_uwg, setup_open_matlab_ref, TEST_DIR
from uwg import SimParam


def test_uwg_simparam_matlab_init():
    """Matlab value comparison for simparam initialization in main uwg."""

    # Initialize uwg from default initialize.uwg for 31 days simulation
    # Jan 1 00:00 - Feb 1 00:00

    # Month,1,        # starting month (1-12)
    # Day,1,          # starting day (1-31)
    # nDay,31,        # number of days to run simultion
    # dtSim,300,      # simulation time step (s)
    # dtWeather,3600, # weather time step (s)

    testuwg = auto_setup_uwg()
    testuwg.generate()

    # open matlab ref file
    uwg_matlab_val = \
        setup_open_matlab_ref('matlab_simparam', 'matlab_ref_simparam_init.txt')

    uwg_python_val = [
        testuwg.simTime.dt,            # uwg time simulation time step
        testuwg.simTime.timeForcing,   # weather data timestep
        testuwg.simTime.month,
        testuwg.simTime.day,
        testuwg.simTime.days,
        testuwg.simTime.timePrint,     # weather data timestep
        testuwg.simTime.timeDay,       # how many times weather senses in a day
        testuwg.simTime.timeSim,       # how many steps in weather data simulation
        testuwg.simTime.timeMax,       # total seconds in simulation days
        testuwg.simTime.nt,            # total number of timesteps for uwg simuation
        testuwg.simTime.julian,
        testuwg.simTime.timeInitial,   # sensor data in epw for intial time
        testuwg.simTime.timeFinal,     # sensor data in epw for final time
        testuwg.simTime.secDay,        # current seconds in day
        testuwg.simTime.hourDay]

    # matlab ref checking
    assert len(uwg_matlab_val) == len(uwg_python_val)
    for i in range(len(uwg_matlab_val)):
        assert uwg_python_val[i] == \
            pytest.approx(uwg_matlab_val[i], abs=1e-15), 'error at index={}'.format(i)


def test_uwg_simparam_matlab_update_date():
    """ Matlab value comparison for simparam UpdateDate function in main uwg """

    # Initialize uwg from initialize_simparam.uwg for 150 day simulation
    # Mar 15 00:00 - Aug 12th, 00:00

    # Month,3,        # starting month (1-12)
    # Day,15,         # starting day (1-31)
    # nDay,150,       # number of days to run simultion
    # dtSim,300,      # simulation time step (s)
    # dtWeather,3600, # weather time step (s)

    param_path = os.path.join(TEST_DIR, 'parameters', 'initialize_simparam.uwg')
    testuwg = auto_setup_uwg(param_path=param_path)

    testuwg.generate()
    testuwg.simulate()

    # open matlab ref file
    uwg_matlab_val = \
        setup_open_matlab_ref('matlab_simparam', 'matlab_ref_simparam_update_date.txt')

    uwg_python_val = [
        testuwg.simTime.secDay,
        testuwg.simTime.day,
        testuwg.simTime.julian,
        testuwg.simTime.month,
        testuwg.simTime.day,
        testuwg.simTime.hourDay,
        # Add 1 to keep consistent with matlab list convention
        testuwg.ceil_time_step + 1]

    # matlab ref checking
    assert len(uwg_matlab_val) == len(uwg_python_val)
    for i in range(len(uwg_matlab_val)):
        assert uwg_python_val[i] == \
            pytest.approx(uwg_matlab_val[i], abs=1e-15), 'error at index={}'.format(i)


def test_simparam():
    """Tests simparam.py."""

    dtSim = 300             # Sim time step
    dtWeather = 3600        # Weather data time-step
    MONTH = 7               # Begin month
    DAY = 30                # Begin day of the month
    NUM_DAYS = 7            # Number of days of simulation

    simTime = SimParam(dtSim, dtWeather, MONTH, DAY, NUM_DAYS)

    # Simulation Parameters tests
    assert simTime.timeSim == pytest.approx(168, abs=1e-6)
    assert simTime.timeMax == pytest.approx(604800, abs=1e-6)
    assert simTime.nt == pytest.approx(2017, abs=1e-6)

    # Test UpdateDate() for < 1 hr
    for i in range(11):  # 11 * 300 = 3300 seconds = 55min
        simTime.update_date()
    assert simTime.secDay == pytest.approx(3300., abs=1e-6)
    assert simTime.day == pytest.approx(30., abs=1e-6)
    assert simTime.hourDay == pytest.approx(0., abs=1e-6)

    # for == 1 hr
    simTime.update_date()
    assert simTime.secDay == pytest.approx(3600., abs=1e-6)
    assert simTime.hourDay == pytest.approx(1., abs=1e-6)

    # for > 24hr
    for i in range(23 * 12):
        simTime.update_date()
    assert simTime.secDay == pytest.approx(0., abs=1e-6)
    assert simTime.day == pytest.approx(31., abs=1e-6)
    assert simTime.hourDay == pytest.approx(0., abs=1e-6)

    # for == 1 month
    for i in range(24 * 12):
        simTime.update_date()
    assert simTime.secDay == pytest.approx(0., abs=1e-6)
    assert simTime.day == pytest.approx(1., abs=1e-6)
    assert simTime.hourDay == pytest.approx(0., abs=1e-6)
    assert simTime.month == pytest.approx(8, abs=1e-6)

    # for + 1 month
    for i in range(24 * 12 * 31):
        simTime.update_date()
    assert simTime.secDay == pytest.approx(0., abs=1e-6)
    assert simTime.day == pytest.approx(1., abs=1e-6)
    assert simTime.hourDay == pytest.approx(0., abs=1e-6)
    assert simTime.month == pytest.approx(9, abs=1e-6)