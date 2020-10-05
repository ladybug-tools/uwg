"""Test the matlab generated epw file against the python generated epw files."""

import pytest
import os

from .test_base import auto_setup_uwg, calculate_tolerance, set_input_manually, \
    MATLAB_DIR, TEST_DIR
import uwg


def compare_epw(testuwg, matlab_fname, precision=10.0):
    # shorten some variable names
    ti = testuwg.simTime.timeInitial
    tf = testuwg.simTime.timeFinal

    # Make weather files for testing
    matlab_path_name = os.path.join(MATLAB_DIR, 'matlab_output', matlab_fname)

    pywtr = uwg.weather.Weather(testuwg.new_epw_path, ti, tf)
    matwtr = uwg.weather.Weather(matlab_path_name, ti, tf)

    assert len(pywtr.staTemp) == pytest.approx(len(matwtr.staTemp), abs=1e-15)

    for i in range(0, len(pywtr.staTemp)):

        # Check dry bulb [K]
        tol = calculate_tolerance(pywtr.staTemp[i], precision)
        assert pywtr.staTemp[i] == \
            pytest.approx(matwtr.staTemp[i], abs=tol), 'error at index={}'.format(i)

        # Check dew point [K]
        tol = calculate_tolerance(pywtr.staTdp[i], precision)
        assert pywtr.staTdp[i] == \
            pytest.approx(matwtr.staTdp[i], abs=tol), 'error at index={}'.format(i)

        # Check relative humidity [%]
        tol = calculate_tolerance(pywtr.staRhum[i], precision)
        assert pywtr.staRhum[i] == \
            pytest.approx(matwtr.staRhum[i], abs=tol), 'error at index={}'.format(i)

        # Check wind speed [m/s]
        tol = calculate_tolerance(pywtr.staUmod[i], precision)
        assert pywtr.staUmod[i] == \
            pytest.approx(matwtr.staUmod[i], abs=tol), 'error at index={}'.format(i)


def test_uwg_output_heatdemand_1_1_0000():
    """Initial conditions
        - night time
        - before vegstart
        - sensHeatDemand
    """
    epw_path = os.path.join(TEST_DIR, 'epw', 'CAN_ON_Toronto.716240_CWEC.epw')
    testuwg = auto_setup_uwg(epw_path=epw_path)

    # Test all year
    testuwg.month = 1
    testuwg.day = 1
    testuwg.nday = 365

    testuwg.generate()
    testuwg.simulate()
    testuwg.write_epw()

    compare_epw(testuwg, 'CAN_ON_Toronto.716240_CWEC_heatdemand_UWG_Matlab.epw')


def test_program_input():

    epw_path = os.path.join(TEST_DIR, 'epw', 'SGP_Singapore.486980_IWEC.epw')
    testuwg = auto_setup_uwg(epw_path=epw_path, param_path=None)

    # Assign manually
    testuwg = set_input_manually(testuwg)

    # main
    testuwg.generate()
    testuwg.simulate()
    testuwg.write_epw()

    # Check some of the inputs

    # Check building parameters
    assert testuwg.BEM[0].building.coolcap == \
           pytest.approx((3525.66904 * 1000.0) / 46320.0, abs=1e-3)
    assert testuwg.BEM[0].building.heat_cap == \
           pytest.approx((2875.97378 * 1000.0) / 46320.0, abs=1e-3)
    assert testuwg.BEM[1].building.coolcap == \
           pytest.approx((252.20895 * 1000.0) / 3135., abs=1e-2)
    assert testuwg.BEM[1].building.heat_cap == \
           pytest.approx((132.396 * 1000.0) / 3135., abs=1e-2)

    # Check that final day of timestep is at correct dayType
    assert testuwg.dayType == 1
    assert testuwg.schtraffic[testuwg.dayType - 1][testuwg.simTime.hourDay] == \
        pytest.approx(0.2, abs=1e-6)

    compare_epw(testuwg, 'SGP_Singapore.486980_IWEC_UWG_Matlab.epw')


def test_program_hybrid_input():
    """Testing inputting with api and with .uwg file."""

    testuwg = auto_setup_uwg()

    # Assign manually
    testuwg = set_input_manually(testuwg)

    # override (from intialize_singapore.uwg)
    testuwg.bldheight = 10
    testuwg.h_mix = 1
    testuwg.blddensity = 0.5
    testuwg.vertohor = 0.8
    testuwg.charlength = 1000
    testuwg.albroad = 0.1
    testuwg.d_road = 0.5
    testuwg.sensanth = 20

    # main
    testuwg._read_epw()
    testuwg._compute_BEM()
    testuwg._compute_input()
    testuwg._hvac_autosize()
    testuwg.simulate()
    testuwg.write_epw()

    # Check some of the inputs

    # Check building parameters
    assert testuwg.BEM[0].building.coolcap == \
           pytest.approx((3525.66904 * 1000.0) / 46320.0, abs=1e-3)
    assert testuwg.BEM[0].building.heat_cap == \
           pytest.approx((2875.97378 * 1000.0) / 46320.0, abs=1e-3)
    assert testuwg.BEM[1].building.coolcap \
        == pytest.approx((252.20895 * 1000.0) / 3135., abs=1e-2)
    assert testuwg.BEM[1].building.heat_cap == \
           pytest.approx((132.396 * 1000.0) / 3135., abs=1e-2)

    # Check that final day of timestep is at correct dayType
    assert testuwg.dayType == 1
    assert testuwg.schtraffic[testuwg.dayType - 1][testuwg.simTime.hourDay] == \
        pytest.approx(0.2, abs=1e-6)

    compare_epw(testuwg, 'SGP_Singapore.486980_IWEC_UWG_Matlab.epw')


def test_uwg_output_beijing():
    """Initial conditions:
        - day time
        - after vegstart
        - sensHeatDemand
    """

    epw_path = os.path.join(TEST_DIR, 'epw', 'CHN_Beijing.Beijing.545110_IWEC.epw')
    param_path = os.path.join(TEST_DIR, 'parameters', 'initialize_beijing.uwg')
    testuwg = auto_setup_uwg(epw_path=epw_path, param_path=param_path)

    # Test all year
    testuwg.month = 1
    testuwg.day = 1
    testuwg.nday = 365

    # main
    testuwg.generate()
    testuwg.simulate()
    testuwg.write_epw()

    compare_epw(testuwg, 'CHN_Beijing.Beijing.545110_IWEC_UWG_Matlab.epw')


def test_uwg_output_cooldemand_6_1_0000():
    """Initial conditions:
        - day time
        - after vegstart
        - sensHeatDemand
    """

    # set up the logger
    epw_path = os.path.join(TEST_DIR, 'epw', 'CAN_ON_Toronto.716240_CWEC.epw')
    testuwg = auto_setup_uwg(epw_path=epw_path)

    # Test 30 days in summer
    testuwg.month = 6
    testuwg.day = 1
    testuwg.nday = 30

    # main
    testuwg.generate()
    testuwg.simulate()
    testuwg.write_epw()

    testuwg = compare_epw(testuwg, 'CAN_ON_Toronto.716240_CWEC_cooldemand_UWG_Matlab.epw')


def test_uwg_output_cooldemand_1_1_0000():
    """Initial conditions:
        - night time
        - before vegstart
        - sensCoolDemand
    """
    testuwg = auto_setup_uwg()

    # Test all year
    testuwg.month = 1
    testuwg.day = 1
    testuwg.nday = 365

    # main
    testuwg.generate()
    testuwg.simulate()
    testuwg.write_epw()

    compare_epw(testuwg, 'SGP_Singapore.486980_IWEC_UWG_Matlab.epw')

