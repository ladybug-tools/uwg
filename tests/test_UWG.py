"""Test for uwg.py"""

import pytest
from .test_base import setup_uwg_integration, set_input_manually


def test_read_epw():
    """Test read epw"""
    testuwg = setup_uwg_integration()
    testuwg._read_epw()

    # test header
    assert testuwg._header[0][0] == 'LOCATION'
    assert testuwg._header[0][1] == 'SINGAPORE'
    assert testuwg.lat == pytest.approx(1.37, abs=1e-3)
    assert testuwg.lon == pytest.approx(103.98, abs=1e-3)
    assert testuwg.GMT == pytest.approx(8, abs=1e-3)
    # test soil data
    assert testuwg.nSoil == pytest.approx(3, abs=1e-2)
    # test soil depths
    assert testuwg.depth_soil[0][0] == pytest.approx(0.5, abs=1e-3)
    assert testuwg.depth_soil[1][0] == pytest.approx(2., abs=1e-3)
    assert testuwg.depth_soil[2][0] == pytest.approx(4., abs=1e-3)
    # test soil temps over 12 months
    assert testuwg.Tsoil[0][0] == pytest.approx(27.55+273.15, abs=1e-3)
    assert testuwg.Tsoil[1][2] == pytest.approx(28.01+273.15, abs=1e-3)
    assert testuwg.Tsoil[2][11] == pytest.approx(27.07+273.15, abs=1e-3)
    # test time step in weather file
    assert testuwg.epwinput[0][0] == '1989'
    assert float(testuwg.epwinput[3][6]) == pytest.approx(24.1, abs=1e-3)


def test_read_input():
    """Test read input."""
    testuwg = setup_uwg_integration()
    testuwg.read_input()
    testuwg.generate()

    # test uwg param dictionary first and last
    assert 'bldheight' in testuwg._init_param_dict
    assert 'h_obs' in testuwg._init_param_dict

    assert testuwg._init_param_dict['bldheight'] == pytest.approx(10., abs=1e-6)
    assert testuwg._init_param_dict['vegend'] == pytest.approx(10, abs=1e-6)
    assert testuwg._init_param_dict['albroof'] is None
    assert testuwg._init_param_dict['h_ubl1'] == pytest.approx(1000., abs=1e-6)
    assert testuwg._init_param_dict['h_ref'] == pytest.approx(150., abs=1e-6)

    # test SchTraffic schedule
    # first
    assert testuwg._init_param_dict['schtraffic'][0][0] == pytest.approx(0.2, abs=1e-6)
    # last
    assert testuwg._init_param_dict['schtraffic'][2][23] == pytest.approx(0.2, abs=1e-6)
    assert testuwg._init_param_dict['schtraffic'][0][19] == pytest.approx(0.8, abs=1e-6)
    assert testuwg._init_param_dict['schtraffic'][1][21] == pytest.approx(0.3, abs=1e-6)
    assert testuwg._init_param_dict['schtraffic'][2][6] == pytest.approx(0.4, abs=1e-6)

    # test bld fraction list
    assert testuwg._init_param_dict['bld'][0][0] == pytest.approx(0., abs=1e-6)
    assert testuwg._init_param_dict['bld'][3][1] == pytest.approx(0.4, abs=1e-6)
    assert testuwg._init_param_dict['bld'][5][1] == pytest.approx(0.6, abs=1e-6)
    assert testuwg._init_param_dict['bld'][15][2] == pytest.approx(0.0, abs=1e-6)

    # test BEMs
    assert len(testuwg.BEM) == pytest.approx(2., abs=1e-6)
    # test BEM office (BLD4 in DOE)
    assert testuwg.BEM[0].building.Type == 'LargeOffice'
    assert testuwg.BEM[0].building.Zone == '1A (Miami)'
    assert testuwg.BEM[0].building.Era == 'Pst80'
    assert testuwg.BEM[0].frac == 0.4

    # test BEM apartment
    assert testuwg.BEM[1].building.Type == 'MidRiseApartment'
    assert testuwg.BEM[1].building.Zone == '1A (Miami)'
    assert testuwg.BEM[1].building.Era == 'Pst80'
    assert testuwg.BEM[1].frac == 0.6

    # Check that schedules are called correctly
    # 9am on Weekday for Office
    assert testuwg.Sch[0].Light[0][8] == pytest.approx(0.9, abs=1e-6)
    # 9am on Weekday for Office
    assert testuwg.Sch[0].Light[0][7] == pytest.approx(0.3, abs=1e-6)
    # 12 noon on Weekend for apt
    assert testuwg.Sch[1].Occ[1][11] == pytest.approx(0.25, abs=1e-6)

    # Check that soil ground depth is set correctly
    assert testuwg.depth_soil[testuwg._soilindex1][0] == pytest.approx(0.5, abs=1e-6)
    assert testuwg.depth_soil[testuwg._soilindex2][0] == pytest.approx(0.5, abs=1e-6)

    # Check the road layer splitting
    assert len(testuwg.road.layerThickness) == pytest.approx(11., abs=1e-15)
    assert testuwg.road.layerThickness[0] == pytest.approx(0.05, abs=1e-15)

    # Check the road layer splitting for rural
    assert len(testuwg.rural.layerThickness) == pytest.approx(11., abs=1e-15)
    assert testuwg.rural.layerThickness[0] == pytest.approx(0.05, abs=1e-6)


def test_optional_blank_parameters():

    testuwg = setup_uwg_integration()
    testuwg.read_input()
    testuwg.generate()

    assert testuwg.BEM[0].building.glazingRatio == pytest.approx(0.38, abs=1e-15)
    assert testuwg.BEM[0].roof.albedo == pytest.approx(0.2, abs=1e-15)
    assert testuwg.BEM[0].roof.vegCoverage == pytest.approx(0.0, abs=1e-15)
    assert testuwg.BEM[1].roof.albedo == pytest.approx(0.2, abs=1e-15)
    assert testuwg.BEM[1].building.glazingRatio == pytest.approx(0.1499, abs=1e-15)
    assert testuwg.BEM[1].roof.vegCoverage == pytest.approx(0.0, abs=1e-15)


def test_optional_inputted_parameters():

    testuwg = setup_uwg_integration(param_path=None)
    testuwg = set_input_manually(testuwg)

    # Test setting values
    with pytest.raises(AssertionError):
        testuwg.albroad = 5

    # Test setting values
    with pytest.raises(AssertionError):
        testuwg.glzr = 5

    # Test setting values
    with pytest.raises(AssertionError):
        testuwg.month = 50

    # Set optional parameters
    testuwg.albroof = .5
    testuwg.vegroof = .1
    testuwg.glzr = .5
    testuwg.albwall = 0.91
    testuwg.shgc = 0.65
    testuwg.flr_h = 4.5

    # From blank inputs will be from DOE
    testuwg.generate()

    assert testuwg.BEM[0].building.glazingRatio == pytest.approx(0.5, abs=1e-15)
    assert testuwg.BEM[0].roof.albedo == pytest.approx(0.5, abs=1e-15)
    assert testuwg.BEM[0].roof.vegCoverage == pytest.approx(0.1, abs=1e-15)
    assert testuwg.BEM[1].building.glazingRatio == pytest.approx(0.5, abs=1e-15)
    assert testuwg.BEM[1].roof.albedo == pytest.approx(0.5, abs=1e-15)
    assert testuwg.BEM[1].roof.vegCoverage == pytest.approx(0.1, abs=1e-15)
    assert testuwg.BEM[0].wall.albedo == pytest.approx(0.91, abs=1e-15)
    assert testuwg.BEM[1].building.shgc == pytest.approx(0.65, abs=1e-15)
    assert testuwg.BEM[0].building.floorHeight == pytest.approx(4.5, abs=1e-15)


def test_procMat():
    """
    Test different max/min layer depths that generate different diffrent road layer
    thicknesses (to account for too deep elements with inaccurate heat transfer).
    """

    testuwg = setup_uwg_integration()
    testuwg.read_input()
    testuwg.generate()

    # test a 0.5m road split into 10 slices of 0.05m
    # base case; min=0.01, max=0.05, stays the same
    roadMat, newthickness = testuwg._procmat(testuwg.road, 0.05, 0.01)
    assert len(roadMat) == pytest.approx(11, abs=1e-6)
    assert len(newthickness) == pytest.approx(11, abs=1e-6)
    assert sum(newthickness) == pytest.approx(0.05*11, abs=1e-6)

    # modify to one layer for tests
    testuwg.road.layerThickness = [0.05]
    testuwg.road.layerThermalCond = testuwg.road.layerThermalCond[:1]
    testuwg.road.layerVolHeat = testuwg.road.layerVolHeat[:1]

    # 0.05 layer, will split in two
    roadMat, newthickness = testuwg._procmat(testuwg.road, 0.05, 0.01)
    assert len(roadMat) == pytest.approx(2, abs=1e-6)
    assert len(newthickness) == pytest.approx(2, abs=1e-6)
    assert sum(newthickness) == pytest.approx(0.025*2, abs=1e-6)

    # 0.015 layer, will split in min thickness in two
    testuwg.road.layerThickness = [0.015]
    roadMat, newthickness = testuwg._procmat(testuwg.road, 0.05, 0.01)
    assert len(roadMat) == pytest.approx(2, abs=1e-6)
    assert len(newthickness) == pytest.approx(2, abs=1e-6)
    assert sum(newthickness) == pytest.approx(0.005*2, abs=1e-6)

    # 0.12 layer, will split into 3 layers b/c > max_thickness
    testuwg.road.layerThickness = [0.12]
    roadMat, newthickness = testuwg._procmat(testuwg.road, 0.05, 0.01)
    assert len(roadMat) == pytest.approx(3, abs=1e-6)
    assert len(newthickness) == pytest.approx(3, abs=1e-6)
    assert sum(newthickness) == pytest.approx(0.04*3, abs=1e-6)


def test_hvac_autosize():
    """Test hvace autosize"""
    testuwg = setup_uwg_integration()
    testuwg.read_input()
    testuwg.generate()

    # Test setting values
    with pytest.raises(AssertionError):
        testuwg.autosize = [1, 2, 3]

    assert testuwg.autosize is False
    assert testuwg.autosize == 0
    assert len(testuwg.BEM) == pytest.approx(2, abs=1e-6)

    # coolCap and heatCap don't retain high accuracy when extracted from the
    # DOE reference csv, so we will reduce the tolerance here
    assert testuwg.BEM[0].building.coolCap == \
        pytest.approx((3525.66904 * 1000.0) / 46320.0, abs=1e-3)
    assert testuwg.BEM[0].building.heatCap == \
        pytest.approx((2875.97378 * 1000.0) / 46320.0, abs=1e-3)
    assert testuwg.BEM[1].building.coolCap \
        == pytest.approx((252.20895 * 1000.0) / 3135., abs=1e-2)
    assert testuwg.BEM[1].building.heatCap \
        == pytest.approx((132.396 * 1000.0) / 3135., abs=1e-2)

    testuwg.autosize = True
    assert testuwg.autosize == 1


def test_simulate():
    """Test UWG simulation."""
    testuwg = setup_uwg_integration()
    testuwg.read_input()
    testuwg.generate()
    testuwg.simulate()
    testuwg.write_epw()

    # Parameters from initialize.uwg
    # Month = 1;              % starting month (1-12)
    # Day = 1;                % starting day (1-31)
    # nDay = 31;              % number of days
    # dtSim = 300;            % simulation time step (s)
    # dtWeather = 3600;       % weather time step (s)

    assert testuwg.N == pytest.approx(744., abs=1e-6)       # total hours in simulation
    assert testuwg.ph == pytest.approx(0.083333, abs=1e-6)  # dt (sim time step) hours

    # test the weather data time series is equal to time step
    assert len(testuwg.forcIP.infra) == \
        pytest.approx((testuwg.simTime.nt - 1) / 12., abs=1e-3)
    # check that simulation time is happening every 5 minutes 8928
    assert testuwg.simTime.nt-1 == pytest.approx(31*24*3600/300., abs=1e-3)
    # check that weather step time is happening every 1 hour = 744
    assert len(testuwg.forcIP.dif) == pytest.approx(31 * 24, abs=1e-3)

    # check that final day of timestep is at correct dayType
    assert testuwg.dayType == pytest.approx(1., abs=1e-3)
    assert testuwg.schtraffic[testuwg.dayType - 1][testuwg.simTime.hourDay] == \
        pytest.approx(0.2, abs=1e-6)
