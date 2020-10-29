"""Tests for UBLDef object."""

import pytest
from .test_base import auto_setup_uwg, setup_open_matlab_ref, calculate_tolerance
from functools import reduce


def test_ubl_init():
    """Test ubl constructor."""

    testuwg = auto_setup_uwg()
    testuwg.generate()
    testuwg.UBL.__repr__()

    # Get uwg_python values
    uwg_python_val = [
        # characteristic length of the urban area (m)
        testuwg.UBL.charLength,
        testuwg.UBL.perimeter,
        # horizontal urban area (m2)
        testuwg.UBL.urbArea,
        # length of the side of the urban area orthogonal to wind dir (m)
        testuwg.UBL.orthLength,
        # length of the side of the urban area parallel to wind dir (m)
        testuwg.UBL.paralLength,
        # urban boundary layer temperature (K)
        testuwg.UBL.ublTemp,
        # urban boundary layer temperature discretization (K)
        reduce(lambda x, y: x + y, testuwg.UBL.ublTempdx),
        testuwg.UBL.dayBLHeight,  # daytime mixing height, orig = 700
        # Sing: 80, Bub-Cap: 50, nighttime boundary-layer height (m); orig 80
        testuwg.UBL.nightBLHeight]

    uwg_matlab_val = testuwg = setup_open_matlab_ref(
        'matlab_ubl', 'matlab_ref_ubl_init.txt')

    # matlab ref checking
    assert len(uwg_matlab_val) == len(uwg_python_val), \
        'matlab={}, python={}'.format(len(uwg_matlab_val), len(uwg_python_val))

    for i in range(len(uwg_matlab_val)):
        tol = calculate_tolerance(uwg_python_val[i], 15.0)
        assert uwg_python_val[i] == pytest.approx(uwg_matlab_val[i], tol), \
            'error at index={}'.format(i)


def test_ublmodel():
    """Test ubl constructor."""

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

    # Get uwg_python values
    uwg_python_val = [
        # urban boundary layer temperature (K)
        testuwg.UBL.ublTemp,
        # urban boundary layer temperature discretaization (K)
        reduce(lambda x, y: x + y, testuwg.UBL.ublTempdx),
        # night boundary layer height (m)
        testuwg.UBL.dayBLHeight,
        # night boundary layer height (m)
        testuwg.UBL.nightBLHeight]

    uwg_matlab_val = setup_open_matlab_ref(
        'matlab_ubl', 'matlab_ref_ublmodel.txt')

    # matlab ref checking
    assert len(uwg_matlab_val) == len(uwg_python_val)

    for i in range(len(uwg_matlab_val)):
        tol = calculate_tolerance(uwg_python_val[i], 15.0)
        assert uwg_python_val[i] == pytest.approx(uwg_matlab_val[i], tol), \
            'error at index={}'.format(i)
