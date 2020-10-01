"""Test for uwg.py"""

import os
import pytest
from copy import deepcopy
from .test_base import auto_setup_uwg, set_input_manually
from uwg import SchDef, BEMDef, Building, Element, Material, UWG
from uwg.readDOE import BLDTYPE, BUILTERA, ZONETYPE
from uwg.utilities import is_near_zero


def test_init():
    """Test initialization methods."""

    test_dir = os.path.abspath(os.path.dirname(__file__))
    param_path = os.path.join(test_dir, 'parameters', 'initialize_singapore.uwg')
    epw_path = os.path.join(test_dir, 'epw', 'SGP_Singapore.486980_IWEC.epw')

    refBEM, refSch = UWG.load_refDOE()
    refBEM[0][2][0].frac = 0.9
    ref_bem_vec = [refBEM[0][2][0], refBEM[2][2][0]]
    ref_sch_vec = [refSch[0][2][0], refSch[2][2][0]]

    # base init
    UWG(epw_path)

    # from param_file
    UWG.from_param_file(epw_path, param_path)

    # from args
    UWG.from_param_args(epw_path, 10.0, 0.5, 0.5, 1)
    model = UWG.from_param_args(epw_path, 10.0, 0.5, 0.5, 1, ref_bem_vector=[],
                                ref_sch_vector=[])
    model.generate()
    assert model.ref_bem_vector is None
    assert model.ref_sch_vector is None

    UWG.from_param_args(epw_path, 10.0, 0.5, 0.5, 1, ref_bem_vector=ref_bem_vec,
                        ref_sch_vector=ref_sch_vec)
    with pytest.raises(AssertionError):
        UWG.from_param_args(epw_path, 10.0, 0.5, 0.5, 1, ref_bem_vector=ref_bem_vec,
                            ref_sch_vector=ref_sch_vec[:1])
    with pytest.raises(AssertionError):
        UWG.from_param_args(epw_path, 10.0, 0.5, 0.5, 1, ref_bem_vector=None,
                            ref_sch_vector=ref_sch_vec)

    # from dict
    data = UWG.from_param_args(epw_path, 10.0, 0.5, 0.5, 1).to_dict(include_refDOE=False)
    UWG.from_dict(data)
    model1 = UWG.from_param_args(
        epw_path, 10.0, 0.5, 0.5, 1, ref_bem_vector=ref_bem_vec,
        ref_sch_vector=ref_sch_vec)
    data = model1.to_dict(include_refDOE=True)
    model2 = UWG.from_dict(data)
    model2.generate()
    assert model2.ref_bem_vector[0].frac == pytest.approx(0.9, abs=1e-10)
    assert model2.refBEM[0][2][0].frac == pytest.approx(0.9, abs=1e-10)


def test_dict():
    """Test uwg to/from dict method."""

    model1 = auto_setup_uwg()

    # Set some optional values
    model1.shgc = 0.3
    model1.glzr = 0.5
    model1.bld = [[0 for c in range(3)] for r in range(16)]
    model1.bld[6][0] = 0.1
    model1.bld[6][2] = 0.9

    # make dict
    uwgdict = model1.to_dict()

    # test if dict and from_dict
    assert isinstance(uwgdict, dict)

    model2 = model1.from_dict(uwgdict)

    # Test attributes (including optional)
    model1.flr_h is model2.flr_h
    model1.blddensity == pytest.approx(model2.blddensity, abs=1e-10)
    model1.c_circ == pytest.approx(model2.c_circ, abs=1e-10)
    model1.shgc == pytest.approx(model2.shgc, abs=1e-10)
    model1.glzr == pytest.approx(model2.glzr, abs=1e-10)

    # bld matrix
    model2.bld[6][0] == pytest.approx(0.1, abs=1e-10)
    model2.bld[6][2] == pytest.approx(0.9, abs=1e-10)

    # Test error
    with pytest.raises(AssertionError):
        uwgdict['type'] = 'Error'
        model1.from_dict(uwgdict)


def test_sch_refDOE():
    """Test uwg from dict method with refDOE override for Schs."""

    model1 = auto_setup_uwg()

    # Set bld matrix and zone
    model1.bld = [[0 for i in range(3)] for j in range(16)]
    model1.bld[1][2] = 1
    model1.zone = 1

    # add schedule to type=2, era=3
    testweek = [[0.1 for i in range(24)] for j in range(3)]
    model1._ref_bem_vector = [model1.refBEM[0][0][0]]
    model1._ref_sch_vector = \
        [SchDef(elec=testweek, gas=testweek, light=testweek, occ=testweek, cool=testweek,
                heat=testweek, swh=testweek, q_elec=18.9, q_gas=3.2, q_light=18.9,
                n_occ=0.12, vent=0.0013, v_swh=0.2846, bldtype=1, builtera=2)]

    model1.generate()  # initialize BEM, Sch objects

    # make dict
    uwgdict = model1.to_dict(include_refDOE=True)
    assert 'ref_sch_vector' in uwgdict
    assert len(uwgdict['ref_sch_vector']) == 1

    model2 = model1.from_dict(uwgdict)
    model2.generate()

    # Check values
    assert model2.bld[1][2] == pytest.approx(1, abs=1e-10)
    testsch = model2.refSchedule[1][2][0]
    for i in range(3):
        for j in range(24):
            assert testsch.elec[i][j] == pytest.approx(0.1, abs=1e-10)
            assert testsch.swh[i][j] == pytest.approx(0.1, abs=1e-10)

    # Test adding 1 extra type to bld on second row
    model1 = auto_setup_uwg()

    # Set bld matrix and zone
    model1.bld = [[0 for i in range(3)] for j in range(18)]
    model1.bld[17][2] = 1
    model1.zone = 1

    testweek = [[0.2 for i in range(24)] for j in range(3)]
    newsch = SchDef(elec=testweek, gas=testweek, light=testweek, occ=testweek,
                    cool=testweek, heat=testweek, swh=testweek, q_elec=18.9,
                    q_gas=3.2, q_light=18.9, n_occ=0.12, vent=0.0013, v_swh=0.2846,
                    bldtype=17, builtera=2)

    newbem = _generate_bemdef()
    model1.ref_bem_vector = [newbem]
    model1.ref_sch_vector = [newsch]

    uwgdict = model1.to_dict(include_refDOE=True)
    model2 = model1.from_dict(uwgdict)
    model2.generate()

    # check lengths
    assert len(uwgdict['ref_sch_vector']) == 1
    assert len(model2.refSchedule) == 18
    assert len(model2.bld) == 18

    # Check values
    testsch = model2.refSchedule[17][2][0]
    for i in range(3):
        for j in range(24):
            assert testsch.elec[i][j] == pytest.approx(0.2, abs=1e-10)
            assert testsch.swh[i][j] == pytest.approx(0.2, abs=1e-10)


def test_bem_refDOE():
    """Test uwg from dict method with refDOE override for BEMs."""

    model1 = auto_setup_uwg()

    # Set bld matrix and zone
    model1.bld = [[0 for i in range(3)] for j in range(16)]
    model1.bld[1][2] = 0.814
    model1.zone = 1

    # add schedule to type=1, era=2
    bem = _generate_bemdef()
    bem.bldtype = 1
    bem.builtera = 2
    bem.building.cop = 4000.0
    bem.roof.emissivity = 0.001
    model1.ref_bem_vector = [bem]
    model1.ref_sch_vector = [model1.refSchedule[0][0][0]]

    # make dict
    uwgdict = model1.to_dict(include_refDOE=True)
    assert 'ref_bem_vector' in uwgdict
    assert len(uwgdict['ref_bem_vector']) == 1

    model2 = model1.from_dict(uwgdict)

    # Test default values being overwritten with compute_BEM
    assert model2.refBEM[1][2][0].frac == pytest.approx(0.0, abs=1e-10)
    model2.generate()
    # Object will be linked therefore modified
    assert model2.refBEM[1][2][0].frac == pytest.approx(0.814, abs=1e-10)

    # Check values
    assert len(model2.BEM) == 1
    testbem = model2.refBEM[1][2][0]

    assert testbem.building.cop == pytest.approx(4000.0, 1e-10)
    assert testbem.roof.emissivity == pytest.approx(0.001, abs=1e-10)


def test_customize_reference_data():
    """Test adding reference DOE data to UWG."""

    model = auto_setup_uwg()
    model.zone = 15
    zi = model.zone - 1

    # make new sched and unrealistic values
    testweek = [[2000.0 for i in range(24)] for j in range(3)]
    newsch1 = SchDef(elec=testweek, gas=testweek, light=testweek, occ=testweek,
                     cool=testweek, heat=testweek, swh=testweek, q_elec=18.9,
                     q_gas=3.2, q_light=18.9, n_occ=0.12, vent=0.0013, v_swh=0.2846,
                     bldtype=5, builtera=0)
    testweek = [[1000.0 for i in range(24)] for j in range(3)]
    newsch2 = SchDef(elec=testweek, gas=testweek, light=testweek, occ=testweek,
                     cool=testweek, heat=testweek, swh=testweek, q_elec=18.9,
                     q_gas=3.2, q_light=18.9, n_occ=0.12, vent=0.0013, v_swh=0.2846,
                     bldtype=19, builtera=2)

    # make new blds and add unrealistic values
    bem1 = _generate_bemdef()
    bem1.bldtype = 5
    bem1.builtera = 0
    bem1.frac = 0.314
    bem1.building.cop = 3000.0
    bem1.roof.emissivity = 0.0

    bem2 = deepcopy(_generate_bemdef())
    bem2.bldtype = 19
    bem2.builtera = 2
    bem2.frac = 0.714
    bem2.building.cop = 4000.0
    bem2.roof.emissivity = 0.001

    # test default lengths
    assert len(model.refSchedule) == 16
    assert len(model.refBEM) == 16

    for day in model.refSchedule[5][0][zi].heat:
        for hr in day:
            assert not is_near_zero(hr - 2000.0, 1e-10)

    assert not is_near_zero(model.refBEM[5][0][zi].frac - 0.314, 1e-10)
    assert not is_near_zero(model.refBEM[5][0][zi].building.cop - 3000.0, 1e-10)
    assert not is_near_zero(model.refBEM[5][0][zi].roof.emissivity - 0.0, 1e-10)

    # run method
    ref_sch_vec = [newsch1, newsch2]
    ref_bem_vec = [bem1, bem2]

    # test bld matrix error
    with pytest.raises(AssertionError):
        model._check_reference_data(ref_bem_vec, ref_sch_vec)

    # set bld matrix and zone
    model.bld = [[0 for i in range(3)] for j in range(20)]
    model.bld[5][0] = 0.5  # test insertion
    model.bld[19][2] = 0.5  # test extention

    model.ref_bem_vector, model.ref_sch_vector = \
        model._check_reference_data(ref_bem_vec, ref_sch_vec)

    with pytest.raises(Exception):
        model.ref_bem_vector = ref_bem_vec
    with pytest.raises(Exception):
        model.ref_sch_vector = ref_sch_vec

    model._customize_reference_data()

    # Test customized schedules
    assert len(model.refSchedule) == 20
    for day in model.refSchedule[5][0][zi].heat:
        for hr in day:
            assert is_near_zero(hr - 2000.0, 1e-10)
    for day in model.refSchedule[19][2][zi].heat:
        for hr in day:
            assert is_near_zero(hr - 1000.0, 1e-10)

    # Test customised bemdefs
    assert len(model.refBEM) == 20
    assert is_near_zero(model.refBEM[5][0][zi].frac - 0.314, 1e-10)
    assert is_near_zero(model.refBEM[5][0][zi].building.cop - 3000.0, 1e-10)
    assert is_near_zero(model.refBEM[5][0][zi].roof.emissivity - 0.0, 1e-10)

    assert is_near_zero(model.refBEM[19][2][zi].frac - 0.714, 1e-10)
    assert is_near_zero(model.refBEM[19][2][zi].building.cop - 4000.0, 1e-10)
    assert is_near_zero(model.refBEM[19][2][zi].roof.emissivity - 0.001, 1e-10)


def test_read_epw():
    """Test read epw"""
    model = auto_setup_uwg()
    model._read_epw()

    # test header
    assert model._header[0][0] == 'LOCATION'
    assert model._header[0][1] == 'SINGAPORE'
    assert model.lat == pytest.approx(1.37, abs=1e-3)
    assert model.lon == pytest.approx(103.98, abs=1e-3)
    assert model.gmt == pytest.approx(8, abs=1e-3)
    # test soil data
    assert model.nSoil == pytest.approx(3, abs=1e-2)
    # test soil depths
    assert model.depth_soil[0][0] == pytest.approx(0.5, abs=1e-3)
    assert model.depth_soil[1][0] == pytest.approx(2., abs=1e-3)
    assert model.depth_soil[2][0] == pytest.approx(4., abs=1e-3)
    # test soil temps over 12 months
    assert model.Tsoil[0][0] == pytest.approx(27.55+273.15, abs=1e-3)
    assert model.Tsoil[1][2] == pytest.approx(28.01+273.15, abs=1e-3)
    assert model.Tsoil[2][11] == pytest.approx(27.07+273.15, abs=1e-3)
    # test time step in weather file
    assert model.epwinput[0][0] == '1989'
    assert float(model.epwinput[3][6]) == pytest.approx(24.1, abs=1e-3)


def test_read_input():
    """Test read input."""
    model = auto_setup_uwg()
    model.generate()

    # test uwg param dictionary first and last
    assert 'bldheight' in model._init_param_dict
    assert 'h_obs' in model._init_param_dict

    assert model._init_param_dict['bldheight'] == pytest.approx(10., abs=1e-6)
    assert model._init_param_dict['vegend'] == pytest.approx(10, abs=1e-6)
    assert model._init_param_dict['albroof'] is None
    assert model._init_param_dict['h_ubl1'] == pytest.approx(1000., abs=1e-6)
    assert model._init_param_dict['h_ref'] == pytest.approx(150., abs=1e-6)

    # test SchTraffic schedule
    # first
    assert model._init_param_dict['schtraffic'][0][0] == pytest.approx(0.2, abs=1e-6)
    # last
    assert model._init_param_dict['schtraffic'][2][23] == pytest.approx(0.2, abs=1e-6)
    assert model._init_param_dict['schtraffic'][0][19] == pytest.approx(0.8, abs=1e-6)
    assert model._init_param_dict['schtraffic'][1][21] == pytest.approx(0.3, abs=1e-6)
    assert model._init_param_dict['schtraffic'][2][6] == pytest.approx(0.4, abs=1e-6)

    # test bld fraction list
    assert model._init_param_dict['bld'][0][0] == pytest.approx(0., abs=1e-6)
    assert model._init_param_dict['bld'][3][1] == pytest.approx(0.4, abs=1e-6)
    assert model._init_param_dict['bld'][5][1] == pytest.approx(0.6, abs=1e-6)
    assert model._init_param_dict['bld'][15][2] == pytest.approx(0.0, abs=1e-6)

    # test BEMs
    assert len(model.BEM) == 2
    # test BEM office (BLD4 in DOE)
    assert BLDTYPE[model.BEM[0].bldtype] == 'LargeOffice'
    assert ZONETYPE[model.BEM[0].zonetype] == '1A (Miami)'
    assert BUILTERA[model.BEM[0].builtera] == 'Pst80'
    assert model.BEM[0].frac == 0.4

    # test BEM apartment
    assert BLDTYPE[model.BEM[1].bldtype] == 'MidRiseApartment'
    assert ZONETYPE[model.BEM[1].zonetype] == '1A (Miami)'
    assert BUILTERA[model.BEM[1].builtera] == 'Pst80'
    assert model.BEM[1].frac == 0.6

    # Check that schedules are called correctly
    # 9am on Weekday for Office
    assert model.Sch[0].light[0][8] == pytest.approx(0.9, abs=1e-6)
    # 9am on Weekday for Office
    assert model.Sch[0].light[0][7] == pytest.approx(0.3, abs=1e-6)
    # 12 noon on Weekend for apt
    assert model.Sch[1].occ[1][11] == pytest.approx(0.25, abs=1e-6)

    # Check that soil ground depth is set correctly
    assert model.depth_soil[model._soilindex1][0] == pytest.approx(0.5, abs=1e-6)
    assert model.depth_soil[model._soilindex2][0] == pytest.approx(0.5, abs=1e-6)

    # Check the road layer splitting
    assert len(model.road.layer_thickness_lst) == pytest.approx(11., abs=1e-15)
    assert model.road.layer_thickness_lst[0] == pytest.approx(0.05, abs=1e-15)

    # Check the road layer splitting for rural
    assert len(model.rural.layer_thickness_lst) == pytest.approx(11., abs=1e-15)
    assert model.rural.layer_thickness_lst[0] == pytest.approx(0.05, abs=1e-6)


def test_optional_blank_parameters():

    model = auto_setup_uwg(param_path=None)
    model = set_input_manually(model)
    model.generate()

    assert model.BEM[0].building.glazing_ratio == pytest.approx(0.38, abs=1e-15)
    assert model.BEM[0].roof.albedo == pytest.approx(0.2, abs=1e-15)
    assert model.BEM[0].roof.vegcoverage == pytest.approx(0.0, abs=1e-15)
    assert model.BEM[1].roof.albedo == pytest.approx(0.2, abs=1e-15)
    assert model.BEM[1].building.glazing_ratio == pytest.approx(0.1499, abs=1e-15)
    assert model.BEM[1].roof.vegcoverage == pytest.approx(0.0, abs=1e-15)


def test_optional_inputted_parameters():

    model = auto_setup_uwg(param_path=None)
    model = set_input_manually(model)

    # test __repr__
    model.generate()
    model.__repr__()

    # Test setting values
    with pytest.raises(AssertionError):
        model.albroad = 5

    # Test setting values
    with pytest.raises(AssertionError):
        model.glzr = 5

    # Test setting values
    with pytest.raises(AssertionError):
        model.month = 50

    # Set optional parameters
    model.albroof = .5
    model.vegroof = .1
    model.glzr = .5
    model.albwall = 0.91
    model.shgc = 0.65
    model.flr_h = 4.5

    # From blank inputs will be from DOE
    model.generate()

    assert model.BEM[0].building.glazing_ratio == pytest.approx(0.5, abs=1e-15)
    assert model.BEM[0].roof.albedo == pytest.approx(0.5, abs=1e-15)
    assert model.BEM[0].roof.vegcoverage == pytest.approx(0.1, abs=1e-15)
    assert model.BEM[1].building.glazing_ratio == pytest.approx(0.5, abs=1e-15)
    assert model.BEM[1].roof.albedo == pytest.approx(0.5, abs=1e-15)
    assert model.BEM[1].roof.vegcoverage == pytest.approx(0.1, abs=1e-15)
    assert model.BEM[0].wall.albedo == pytest.approx(0.91, abs=1e-15)
    assert model.BEM[1].building.shgc == pytest.approx(0.65, abs=1e-15)
    assert model.BEM[0].building.floor_height == pytest.approx(4.5, abs=1e-15)


def test_procMat():
    """
    Test different max/min layer depths that generate different diffrent road layer
    thicknesses (to account for too deep elements with inaccurate heat transfer).
    """

    model = auto_setup_uwg()
    model.generate()

    # test a 0.5m road split into 10 slices of 0.05m
    # base case; min=0.01, max=0.05, stays the same
    roadMat, newthickness = model._procmat(model.road, 0.05, 0.01)
    assert len(roadMat) == pytest.approx(11, abs=1e-6)
    assert len(newthickness) == pytest.approx(11, abs=1e-6)
    assert sum(newthickness) == pytest.approx(0.05*11, abs=1e-6)

    # modify to one layer for tests
    model.road.layer_thickness_lst = [0.05]
    model.road.layerThermalCond = model.road.layerThermalCond[:1]
    model.road.layerVolHeat = model.road.layerVolHeat[:1]

    # 0.05 layer, will split in two
    roadMat, newthickness = model._procmat(model.road, 0.05, 0.01)
    assert len(roadMat) == pytest.approx(2, abs=1e-6)
    assert len(newthickness) == pytest.approx(2, abs=1e-6)
    assert sum(newthickness) == pytest.approx(0.025*2, abs=1e-6)

    # 0.015 layer, will split in min thickness in two
    model.road.layer_thickness_lst = [0.015]
    roadMat, newthickness = model._procmat(model.road, 0.05, 0.01)
    assert len(roadMat) == pytest.approx(2, abs=1e-6)
    assert len(newthickness) == pytest.approx(2, abs=1e-6)
    assert sum(newthickness) == pytest.approx(0.005*2, abs=1e-6)

    # 0.12 layer, will split into 3 layers b/c > max_thickness
    model.road.layer_thickness_lst = [0.12]
    roadMat, newthickness = model._procmat(model.road, 0.05, 0.01)
    assert len(roadMat) == pytest.approx(3, abs=1e-6)
    assert len(newthickness) == pytest.approx(3, abs=1e-6)
    assert sum(newthickness) == pytest.approx(0.04*3, abs=1e-6)


def test_hvac_autosize():
    """Test hvace autosize"""
    model = auto_setup_uwg()
    model.generate()

    # Test setting values
    with pytest.raises(AssertionError):
        model.autosize = [1, 2, 3]

    assert model.autosize is False
    assert model.autosize == 0
    assert len(model.BEM) == pytest.approx(2, abs=1e-6)

    # coolCap and heatCap don't retain high accuracy when extracted from the
    # DOE reference csv, so we will reduce the tolerance here
    assert model.BEM[0].building.coolcap == \
           pytest.approx((3525.66904 * 1000.0) / 46320.0, abs=1e-3)
    assert model.BEM[0].building.heat_cap == \
           pytest.approx((2875.97378 * 1000.0) / 46320.0, abs=1e-3)
    assert model.BEM[1].building.coolcap \
        == pytest.approx((252.20895 * 1000.0) / 3135., abs=1e-2)
    assert model.BEM[1].building.heat_cap \
        == pytest.approx((132.396 * 1000.0) / 3135., abs=1e-2)

    model.autosize = True
    assert model.autosize == 1


def test_simulate():
    """Test UWG simulation."""
    model = auto_setup_uwg()
    model.generate()

    model.simulate()
    model.write_epw()

    # Parameters from initialize.uwg
    # Month = 1;              % starting month (1-12)
    # Day = 1;                % starting day (1-31)
    # nDay = 31;              % number of days
    # dtSim = 300;            % simulation time step (s)
    # dtWeather = 3600;       % weather time step (s)

    assert model.N == pytest.approx(744., abs=1e-6)       # total hours in simulation
    assert model.ph == pytest.approx(0.083333, abs=1e-6)  # dt (sim time step) hours

    # test the weather data time series is equal to time step
    assert len(model.forcIP.infra) == \
        pytest.approx((model.simTime.nt - 1) / 12., abs=1e-3)
    # check that simulation time is happening every 5 minutes 8928
    assert model.simTime.nt-1 == pytest.approx(31*24*3600/300., abs=1e-3)
    # check that weather step time is happening every 1 hour = 744
    assert len(model.forcIP.dif) == pytest.approx(31 * 24, abs=1e-3)

    # check that final day of timestep is at correct dayType
    assert model.dayType == pytest.approx(1., abs=1e-3)
    assert model.schtraffic[model.dayType - 1][model.simTime.hourDay] == \
        pytest.approx(0.2, abs=1e-6)


def _generate_bemdef():
    """Create BEMDef: LargeOffce, Pst80, Zone 1A (Miami)."""

    # Material: (thermalCond, volHeat = specific heat * density)
    concrete = Material(1.311, 836.8 * 2240, 'Concrete')
    gypsum = Material(0.16, 830.0 * 784.9, 'Gypsum')
    stucco = Material(0.6918, 837.0 * 1858.0, 'Stucco')
    insulation = Material(0.049, 836.8 * 265.0, "Insulation")

    # Mass wall for LargeOffce, Pst80, Zone 1A (Miami)
    thicknessLst = [0.0254, 0.0508, 0.0508, 0.0508, 0.0508, 0.0127]
    materialLst = [stucco, concrete, concrete, concrete, concrete, gypsum]
    wall = Element(albedo=0.08, emissivity=0.92, layer_thickness_lst=thicknessLst,
                   material_lst=materialLst, vegcoverage=0, t_init=293,
                   horizontal=False, name='MassWall')

    # IEAD roof
    thicknessLst = [0.058, 0.058]
    materialLst = [insulation, insulation]
    roof = Element(albedo=0.2, emissivity=0.93, layer_thickness_lst=thicknessLst,
                   material_lst=materialLst, vegcoverage=0.5, t_init=293,
                   horizontal=True, name='IEAD')

    # Mass floor
    thicknessLst = [0.054, 0.054]
    materialLst = [concrete, concrete]
    floor = Element(albedo=0.2, emissivity=0.9, layer_thickness_lst=thicknessLst,
                    material_lst=materialLst, vegcoverage=0.0, t_init=293,
                    horizontal=True, name='MassFloor')

    bld = Building(floor_height=3.5, int_heat_night=1, int_heat_day=1, int_heat_frad=0.1,
                   int_heat_flat=0.1, infil=0.26, vent=0.0005, glazing_ratio=0.4,
                   u_value=5.8, shgc=0.2, condtype='AIR', cop=5.2, cool_setpoint_day=297,
                   cool_setpoint_night=297, heat_setpoint_day=293,
                   heat_setpoint_night=293, coolcap=76, heateff=0.7, initial_temp=293)

    return BEMDef(bld, floor, wall, roof, frac=0.1, bldtype=0, builtera=1)