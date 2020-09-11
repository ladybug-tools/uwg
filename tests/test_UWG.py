"""Test for uwg.py"""

import pytest
from .test_base import auto_setup_uwg, set_input_manually
from uwg import SchDef
from uwg.readDOE import BLDTYPE, BUILTERA, ZONETYPE
from uwg.utilities import is_near_zero


def test_read_epw():
    """Test read epw"""
    testuwg = auto_setup_uwg()
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
    testuwg = auto_setup_uwg()
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
    assert len(testuwg.BEM) == 2
    # test BEM office (BLD4 in DOE)
    assert BLDTYPE[testuwg.BEM[0].bldtype] == 'LargeOffice'
    assert ZONETYPE[testuwg.BEM[0].zonetype] == '1A (Miami)'
    assert BUILTERA[testuwg.BEM[0].builtera] == 'Pst80'
    assert testuwg.BEM[0].frac == 0.4

    # test BEM apartment
    assert BLDTYPE[testuwg.BEM[1].bldtype] == 'MidRiseApartment'
    assert ZONETYPE[testuwg.BEM[1].zonetype] == '1A (Miami)'
    assert BUILTERA[testuwg.BEM[1].builtera] == 'Pst80'
    assert testuwg.BEM[1].frac == 0.6

    # Check that schedules are called correctly
    # 9am on Weekday for Office
    assert testuwg.Sch[0].light[0][8] == pytest.approx(0.9, abs=1e-6)
    # 9am on Weekday for Office
    assert testuwg.Sch[0].light[0][7] == pytest.approx(0.3, abs=1e-6)
    # 12 noon on Weekend for apt
    assert testuwg.Sch[1].occ[1][11] == pytest.approx(0.25, abs=1e-6)

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

    testuwg = auto_setup_uwg(param_path=None)
    testuwg = set_input_manually(testuwg)
    testuwg.generate()

    assert testuwg.BEM[0].building.glazingRatio == pytest.approx(0.38, abs=1e-15)
    assert testuwg.BEM[0].roof.albedo == pytest.approx(0.2, abs=1e-15)
    assert testuwg.BEM[0].roof.vegCoverage == pytest.approx(0.0, abs=1e-15)
    assert testuwg.BEM[1].roof.albedo == pytest.approx(0.2, abs=1e-15)
    assert testuwg.BEM[1].building.glazingRatio == pytest.approx(0.1499, abs=1e-15)
    assert testuwg.BEM[1].roof.vegCoverage == pytest.approx(0.0, abs=1e-15)


def test_optional_inputted_parameters():

    testuwg = auto_setup_uwg(param_path=None)
    testuwg = set_input_manually(testuwg)

    # test __repr__
    testuwg.generate()
    testuwg.__repr__()

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

    testuwg = auto_setup_uwg()
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
    testuwg = auto_setup_uwg()
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
    testuwg = auto_setup_uwg()
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


def test_dict():
    """Test uwg from dict method."""

    testuwg1 = auto_setup_uwg()

    # Set some optional values
    testuwg1.shgc = 0.3
    testuwg1.glzr = 0.5
    testuwg1.bld = [[0 for c in range(3)] for r in range(16)]
    testuwg1.bld[6][0] = 0.1
    testuwg1.bld[6][2] = 0.9

    # make dict
    uwgdict = testuwg1.to_dict()

    # test if dict and from_dict
    assert isinstance(uwgdict, dict)

    testuwg2 = testuwg1.from_dict(uwgdict)

    # Test attributes (including optional)
    testuwg1.flr_h is testuwg2.flr_h
    testuwg1.blddensity == pytest.approx(testuwg2.blddensity, abs=1e-10)
    testuwg1.c_circ == pytest.approx(testuwg2.c_circ, abs=1e-10)
    testuwg1.shgc == pytest.approx(testuwg2.shgc, abs=1e-10)
    testuwg1.glzr == pytest.approx(testuwg2.glzr, abs=1e-10)

    # bld matrix
    testuwg2.bld[6][0] == pytest.approx(0.1, abs=1e-10)
    testuwg2.bld[6][2] == pytest.approx(0.9, abs=1e-10)

    # Test error
    with pytest.raises(AssertionError):
        uwgdict['type'] = 'Error'
        testuwg1.from_dict(uwgdict)


def test_sch_refDOE():
    """Test uwg from dict method with refDOE override."""

    testuwg1 = auto_setup_uwg()

    # Set bld matrix and zone
    testuwg1.bld = [[0 for i in range(3)] for j in range(16)]
    testuwg1.bld[1][2] = 1
    testuwg1.zone = 1

    # add schedule to type=2, era=3
    testweek = [[0.1 for i in range(24)] for j in range(3)]
    testuwg1.refSchedule[1][2][0] = \
        SchDef(elec=testweek, gas=testweek, light=testweek, occ=testweek, cool=testweek,
               heat=testweek, swh=testweek, bldtype=1, builtera=2)

    testuwg1.generate()  # initialize BEM, Sch objects

    # make dict
    uwgdict = testuwg1.to_dict(include_refDOE=True)
    assert 'ref_sch_vector' in uwgdict
    assert len(uwgdict['ref_sch_vector']) == 1

    testuwg2 = testuwg1.from_dict(uwgdict)

    # Check values
    assert testuwg2.bld[1][2] == pytest.approx(1, abs=1e-10)
    testsch = testuwg2.refSchedule[1][2][0]
    for i in range(3):
        for j in range(24):
            assert testsch.elec[i][j] == pytest.approx(0.1, abs=1e-10)
            assert testsch.swh[i][j] == pytest.approx(0.1, abs=1e-10)

    # Test adding 1 extra type to bld on second row
    testuwg1 = auto_setup_uwg()

    # Set bld matrix and zone
    testuwg1.bld = [[0 for i in range(3)] for j in range(18)]
    testuwg1.bld[17][2] = 1
    testuwg1.zone = 1

    testweek = [[0.2 for i in range(24)] for j in range(3)]
    newsch = SchDef(elec=testweek, gas=testweek, light=testweek, occ=testweek,
                    cool=testweek, heat=testweek, swh=testweek, bldtype=17, builtera=2)

    # extend reference matrices
    testuwg1.refBEM.extend([[None for i in range(3)] for j in range(2)])
    testuwg1.refSchedule.extend([[[None for k in range(16)]
                                 for i in range(3)] for j in range(2)])
    testuwg1.refSchedule[17][2][0] = newsch

    uwgdict = testuwg1.to_dict(include_refDOE=True)
    testuwg2 = testuwg1.from_dict(uwgdict)

    # check lengths
    assert len(uwgdict['ref_sch_vector']) == 1
    assert len(testuwg2.refSchedule) == 18
    assert len(testuwg2.bld) == 18

    # Check values
    testsch = testuwg2.refSchedule[17][2][0]
    for i in range(3):
        for j in range(24):
            assert testsch.elec[i][j] == pytest.approx(0.2, abs=1e-10)
            assert testsch.swh[i][j] == pytest.approx(0.2, abs=1e-10)


def test_customize_reference_data():
    """Test adding reference DOE data to UWG."""

    testuwg = auto_setup_uwg()

    # set bld matrix and zone
    testuwg.bld = [[0 for i in range(3)] for j in range(20)]
    testuwg.bld[5][0] = 0.5  # test insertion
    testuwg.bld[19][2] = 0.5  # test extention
    testuwg.zone = 15
    zi = testuwg.zone - 1

    # make new sched
    testweek = [[2000.0 for i in range(24)] for j in range(3)]
    newsch1 = SchDef(elec=testweek, gas=testweek, light=testweek, occ=testweek,
                     cool=testweek, heat=testweek, swh=testweek,
                     bldtype=5, builtera=0)
    testweek = [[1000.0 for i in range(24)] for j in range(3)]
    newsch2 = SchDef(elec=testweek, gas=testweek, light=testweek, occ=testweek,
                     cool=testweek, heat=testweek, swh=testweek,
                     bldtype=19, builtera=2)

    # Test default lengths
    assert len(testuwg.refSchedule) == 16
    for day in testuwg.refSchedule[5][0][zi].heat:
        for hr in day:
            assert not is_near_zero(hr - 2000.0, 1e-10)

    # TODO add bem
    ref_sch_vec = [newsch1, newsch2]
    ref_bem_vec = [None, None]
    testuwg.customize_reference_data(ref_bem_vec, ref_sch_vec, zi)

    # Test lengths
    assert len(testuwg.refSchedule) == 20
    for day in testuwg.refSchedule[5][0][zi].heat:
        for hr in day:
            assert is_near_zero(hr - 2000.0, 1e-10)




