from simparam import SimParam



def test_simparam():
    """Test for simparam.py"""

    TOL = 1e-3
    dtSim = 300              # Sim time step
    dtWeather = 3600       # Weather data time-step
    MONTH = 7                # Begin month
    DAY = 30                 # Begin day of the month
    NUM_DAYS = 7             # Number of days of simulation

    simTime = SimParam(dtSim,dtWeather,MONTH,DAY,NUM_DAYS)

    # Simulation Parameters tests
    assert abs(simTime.timeSim-168) < TOL
    assert abs(simTime.timeMax-604800) < TOL
    assert abs(simTime.nt-2017) < TOL

if __name__ == '__main__':
    test_simparam()
