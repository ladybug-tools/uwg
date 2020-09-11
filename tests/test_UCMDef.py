"""Test UCMDef object."""

try:
    range = xrange
except NameError:
    pass

import pytest
import math
from .test_base import auto_setup_uwg, setup_open_matlab_ref, calculate_tolerance


def test_ucm_init():
    """Test ucm constructor."""

    testuwg = auto_setup_uwg()
    testuwg.generate()

    # Because we are not using latAnthrop, this value can be None in uwg.
    # So we will hardcode value for now.
    testuwg.UCM.latAnthrop = 2.0

    # Get uwg_python values
    uwg_python_val = [
        testuwg.UCM.h_mix,         # waste heat mix into canyon ratio
        testuwg.UCM.bldHeight,     # average building height (m)
        testuwg.UCM.verToHor,      # vertical-to-horizontal urban area ratio
        testuwg.UCM.bldDensity,    # horizontal building density (footprint)
        testuwg.UCM.treeCoverage,  # horizontal tree coverage (footprint)
        testuwg.UCM.sensAnthrop,   # sensible heat anthropogenic (W m-2)
        testuwg.UCM.latAnthrop,    # latent heat anthropogenic (W m-2)
        testuwg.UCM.roadShad,      # shadow fraction of road
        testuwg.UCM.bldWidth,      # building width (m)
        testuwg.UCM.canWidth,      # canyon width (m)
        testuwg.UCM.canAspect,     # canyon aspect ratio
        testuwg.UCM.roadConf,      # road sky view factor
        testuwg.UCM.wallConf,      # wall sky view factor
        testuwg.UCM.facArea,       # facade area (m2)
        testuwg.UCM.roadArea,      # road area (m2)
        testuwg.UCM.roofArea,      # roof area (m2)
        testuwg.UCM.canTemp,       # canyon air temperature (db) (K)
        testuwg.UCM.roadTemp,      # avg road temperature
        testuwg.UCM.canHum,        # canyon specific humidity (kgv kgda-1)
        testuwg.UCM.ublWind,       # urban boundary layer wind velocity (m s-1)
        testuwg.UCM.canWind,       # urban canyon wind velocity (m s-1)
        testuwg.UCM.ustar,         # friction velocity (m s-1)
        testuwg.UCM.ustarMod,      # modified friction velocity (m s-1)
        testuwg.UCM.z0u,           # urban roughness length (m)
        testuwg.UCM.l_disp,        # urban displacement length (m)
        testuwg.UCM.alb_wall,      # albedo of wall
        testuwg.UCM.facAbsor,      # absorptivity of wall
        testuwg.UCM.roadAbsor,     # average road absorptivity
        testuwg.UCM.sensHeat]      # urban sensible heat [W m-2]

    uwg_matlab_val = testuwg = setup_open_matlab_ref(
        'matlab_ucm', 'matlab_ref_ucm_init.txt')

    # matlab ref checking
    assert len(uwg_matlab_val) == len(uwg_python_val)
    for i in range(len(uwg_matlab_val)):
        tol = calculate_tolerance(uwg_python_val[i], 15.0)
        assert uwg_python_val[i] == pytest.approx(uwg_matlab_val[i], tol), \
            'error at index={}'.format(i)


def test_ucm_ucmodel():
    """Test ucm ucmodel."""

    testuwg = auto_setup_uwg()

    # Test Jan 1 (winter, no vegetation coverage)
    testuwg.month = 1
    testuwg.day = 1
    testuwg.nday = 1

    # set_input
    testuwg._read_epw()
    testuwg._compute_BEM()
    testuwg._compute_input()

    # In order to avoid integration effects. Test only first time step
    # Subtract timestep to stop at 300 sec
    testuwg.simTime.nt -= (23 * 12 + 11)

    # Run simulation
    testuwg._hvac_autosize()
    testuwg.simulate()

    # check date
    assert testuwg.simTime.month == 1
    assert testuwg.simTime.day == 1
    assert testuwg.simTime.secDay == pytest.approx(300.0, abs=1e-15)

    # Get uwg values
    uwg_python_val = [
        # heat load building
        testuwg.UCM.Q_wall,        # convective sensible heat flux from building roof
        testuwg.UCM.Q_window,      # sensible heat flux from window (via U-factor)
        testuwg.UCM.Q_hvac,        # sensible heat flux from HVAC waste
        testuwg.UCM.ElecTotal,     # total electricity consumption of urban area
        testuwg.UCM.GasTotal,      # total gas consumption of urban area
        # heat load road/canyon/ubl
        testuwg.UCM.Q_road,        # convective sensible heat flux from road
        testuwg.UCM.Q_traffic,     # net sensible heat flux from traffic
        testuwg.UCM.Q_vent,        # convective heat exchange from vent/infil
        testuwg.UCM.Q_ubl,         # convective heat exchange with UBL layer
        # Building temperature
        testuwg.UCM.roofTemp,      # average roof temp (K)
        testuwg.UCM.wallTemp,      # average wall temp (K)
        # Sensible heat
        testuwg.UCM.treeSensHeat,  # sensible heat from trees
        testuwg.UCM.sensHeat,      # urban sensible heat
        testuwg.UCM.canTemp        # canyon air temp (K)
    ]

    # Get uwg_matlab values
    uwg_matlab_val = setup_open_matlab_ref('matlab_ucm', 'matlab_ref_ucm_ucmodel.txt')

    # matlab ref checking
    assert len(uwg_matlab_val) == len(uwg_python_val)

    for i in range(len(uwg_matlab_val)):
        tol = calculate_tolerance(abs(uwg_python_val[i]), 15.0)
        tol = 10 ** (math.log10(tol) + 2)  # Lowering tolerance by two orders of magnitude
        assert uwg_python_val[i] == pytest.approx(uwg_matlab_val[i], abs=tol), \
            'error at index={}'.format(i)

