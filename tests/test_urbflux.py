"""Tests for urbflux functions."""

import pytest
from .test_base import auto_setup_uwg, setup_open_matlab_ref, calculate_tolerance


def test_urbflux_unit():
    """Test for urbflux."""

    testuwg = auto_setup_uwg()

    # Test Jan 1 (winter, no vegetation coverage)
    testuwg.month = 1
    testuwg.day = 1
    testuwg.nday = 1

    # set_input
    testuwg.generate()

    # In order to avoid integration effects. Test only first time step
    # Subtract timestep to stop at 300 sec
    testuwg.simTime.nt -= (23 * 12 + 11)

    # Run simulation
    testuwg.simulate()

    # check date
    assert testuwg.simTime.month == 1
    assert testuwg.simTime.day == 1
    assert testuwg.simTime.secDay == pytest.approx(300.0, abs=1e-15)

    # Calculated values
    uwg_python_val = [
        # building surface properties
        testuwg.BEM[1].roof.emissivity,        # roof emissivity
        testuwg.BEM[1].roof.layerTemp[0],      # roof surface temp (K)
        testuwg.BEM[1].roof.infra,             # roof net longwave radiation (W m-2)
        testuwg.BEM[1].mass.layerTemp[-1],
        testuwg.BEM[1].wall.emissivity,        # wall emissivity
        testuwg.BEM[1].wall.layerTemp[0],      # wall surface temp (K)
        testuwg.BEM[1].wall.infra,             # wall net longwave radiation (W m-2)
        testuwg.BEM[1].wall.layerTemp[-1],     # wall interior surface temp (K)
        # urban canyon properties
        testuwg.UCM.wallTemp,                  # canyon wall temp (K)
        testuwg.UCM.roofTemp,                  # canyon roof temp (K)
        testuwg.road.infra,                    # road net longwave radiation (W m-2)
        testuwg.UCM.roadTemp,                  # canyon road temperature (W m-2)
        testuwg.UCM.latHeat,                   # canyon latent heat (W m-2)
        testuwg.UBL.advHeat,                   # boundary layer advective heat (W m-2)
        testuwg.UCM.ustar,                     # canyon friction velocity (m s-1)
        testuwg.UCM.ustarMod,                  # modified friction velocity (m s-1)
        testuwg.UCM.uExch,                     # exchange velocity (m s-1)
        testuwg.UCM.canWind,                   # canyon wind velocity (m s-1)
        testuwg.UCM.turbU,                     # canyon turbulent velocity (m s-1)
        testuwg.UCM.turbV,                     # canyon turbulent velocity (m s-1)
        testuwg.UCM.turbW,                     # canyon turbulent velocity (m s-1)
        testuwg.UCM.windProf[-1]               # canyon wind profile
    ]

    # UBL.advHeat = Boundary layer advective heat rounding error
    # Past 10^-9 it's basically sero, but for accuracy's sake apply
    # Kahan's summation algorithm
    # -6.87022037029e-11  ==  -3.43511e-11 << w/o Kahan summation algo
    # -6.01144282401e-11  ==  -3.43511e-11 << w Kahan summation algo

    uwg_matlab_val = \
        setup_open_matlab_ref('matlab_urbflux', 'matlab_ref_urbflux_unit.txt')

    # matlab ref checking
    assert len(uwg_matlab_val) == len(uwg_python_val)

    for i in range(len(uwg_matlab_val)):
        if uwg_matlab_val[i] == '\n':
            uwg_matlab_val[i] = None
            assert uwg_python_val[i] == uwg_matlab_val[i], 'error at index={}'.format(i)
            continue

        elif i == 13:
            tol = 10.  # hardcoding tolerance b/c advheat rounding errors

        else:
            tol = calculate_tolerance(uwg_python_val[i], 15.0)

        assert uwg_python_val[i] == pytest.approx(uwg_matlab_val[i], abs=tol), \
            'error at index={}'.format(i)
