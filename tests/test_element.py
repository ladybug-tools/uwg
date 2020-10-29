"""Test Element class."""

import pytest
from .test_base import auto_setup_uwg, setup_open_matlab_ref, calculate_tolerance
from uwg import Element, Material


def test_element_init():
    """Test Element init method."""
    # Material: (thermalCond, volHeat = specific heat * density)
    concrete = Material(1.311, 836.8 * 2240, 'Concrete')
    gypsum = Material(0.16, 830.0 * 784.9, 'Gypsum')
    stucco = Material(0.6918, 837.0 * 1858.0, 'Stucco')

    # Mass wall for LargeOffce, Pst80, Zone 1A (Miami)
    thicknessLst = [0.0254, 0.0508, 0.0508, 0.0508, 0.0508, 0.0127]
    materialLst = [stucco, concrete, concrete, concrete, concrete, gypsum]
    wall = Element(albedo=0.08, emissivity=0.92, layer_thickness_lst=thicknessLst,
                   material_lst=materialLst, vegcoverage=0, t_init=293,
                   horizontal=False, name='MassWall')

    # test repr
    wall.__repr__()


def test_element_dict():
    """Test Element dict method."""

    # Init
    # Material: (thermalCond, volHeat = specific heat * density)
    concrete = Material(1.311, 836.8 * 2240, 'Concrete')
    gypsum = Material(0.16, 830.0 * 784.9, 'Gypsum')
    stucco = Material(0.6918, 837.0 * 1858.0, 'Stucco')

    # Mass wall for LargeOffce, Pst80, Zone 1A (Miami)
    thicknessLst = [0.0254, 0.0508, 0.0508, 0.0508, 0.0508, 0.0127]
    materialLst = [stucco, concrete, concrete, concrete, concrete, gypsum]
    wall = Element(albedo=0.08, emissivity=0.92, layer_thickness_lst=thicknessLst,
                   material_lst=materialLst, vegcoverage=0, t_init=293,
                   horizontal=False, name='MassWall')

    # make dict
    eldict = wall.to_dict()
    wall2 = Element.from_dict(eldict)

    assert wall.albedo == pytest.approx(wall2.albedo, abs=1e-10)
    assert wall.material_lst[0].thermalcond == \
           pytest.approx(wall2.material_lst[0].thermalcond, abs=1e-10)
    assert wall.emissivity == pytest.approx(wall2.emissivity, abs=1e-10)
    assert wall.layer_thickness_lst[1] == \
           pytest.approx(wall2.layer_thickness_lst[1], abs=1e-10)

    with pytest.raises(AssertionError):
        eldict['type'] = 'BemDef'
        Element.from_dict(eldict)


def test_SurfFlux_with_waterStorage_start():
    """Edge case: test element SurfFlux against matlab reference

    When waterStorage > 0.0. This has to be hardcoded b/c doesn't get used
    otherwise.
    """

    testuwg = auto_setup_uwg()
    testuwg.generate()

    # We subtract 30 days and 11 hours
    # New time: Jan 1, 1:00
    testuwg.simTime.nt -= (30 * 24 * 12 + 23 * 12 + 11)

    # turn rural road waterStorage to 1.
    testuwg.rural.waterStorage = 0.005  # .5cm thick film (from wgmax constant)

    # run simulation
    testuwg._hvac_autosize()
    testuwg.simulate()

    # check date
    assert testuwg.simTime.month == 1
    assert testuwg.simTime.day == 1
    assert testuwg.simTime.secDay == pytest.approx(300., abs=1e-15)

    # Check waterStorage
    assert 0.005 == pytest.approx(testuwg.rural.waterStorage, 1e-15)


def test_SurfFlux_with_waterStorage_middle():
    """ Edge case: test element SurfFlux against matlab reference.

    When waterStorage > 0.0. This has to be hardcoded b/c doesn't get used otherwise.
    """
    testuwg = auto_setup_uwg()
    testuwg.generate()

    # We subtract 30 days and 11 hours
    # New time: Jan 1, 1:00
    testuwg.simTime.nt -= (30 * 24 * 12 + 11 * 12)

    # turn rural road waterStorage to 1.
    testuwg.rural.waterStorage = 0.005  # .5cm thick film (from wgmax constant)

    # run simulation
    testuwg._hvac_autosize()
    testuwg.simulate()

    # check date
    assert testuwg.simTime.month == 1
    assert testuwg.simTime.day == 1
    assert testuwg.simTime.secDay/3600. == pytest.approx(13.0, abs=1e-15)

    # Check waterStorage
    assert 0.004643553960210 == pytest.approx(testuwg.rural.waterStorage, 1e-15)


def test_SurfFlux_unit():
    """Test element SurfFlux against matlab references at the start of timestep."""

    testuwg = auto_setup_uwg()

    testuwg._read_epw()
    testuwg.generate()

    # We subtract 23 hours and 55 minutes so we can test
    # initial timestep (1, 1, 300). New time: Jan 1, 5min
    testuwg.simTime.nt -= (23 * 12 + 24 * 12 * 30 + 11)

    # run simulation
    testuwg._hvac_autosize()
    testuwg.simulate()

    # check date is Jan 1st, 300s
    assert testuwg.simTime.month == 1
    assert testuwg.simTime.day == 1
    assert testuwg.simTime.secDay == pytest.approx(300.0, abs=1e-15)

    uwg_python_val = [
        testuwg.rural.aeroCond,        # Convection coef (refL uwg, eq.12)
        testuwg.rural.waterStorage,    # thickness of water film (m) (only for hor surf)
        testuwg.rural.solAbs,          # solar radiation absorbed (W m-2)
        testuwg.rural.lat,             # surface latent heat flux (W m-2)
        testuwg.rural.sens,            # surface sensible heat flux (W m-2)
        testuwg.rural.flux,            # external surface heat flux (W m-2)
        testuwg.rural.T_ext,           # external surface temperature (K)
        testuwg.rural.T_int]           # interior surface temperature (K)

    # Matlab Checking for rural road
    uwg_matlab_val = setup_open_matlab_ref(
        'matlab_element', 'matlab_ref_element_surfflux_winter.txt')

    # matlab ref checking
    assert len(uwg_matlab_val) == len(uwg_python_val)

    for i in range(len(uwg_matlab_val)):
        tol = calculate_tolerance(uwg_matlab_val[i], 15.0)
        assert uwg_python_val[i] == pytest.approx(uwg_matlab_val[i], abs=tol), \
            'error at index={}'.format(i)


def test_SurfFlux_integration():
    """Test element SurfFlux against matlab references at the end of timestep.

    Integration test as it requires other uwg classes to be working before functioning
    correctly.
    """

    testuwg = auto_setup_uwg()

    # Change time and vegCoverage parameters so we can get
    # effect of vegetation on surface heat flux
    testuwg.vegstart = 2    # February
    testuwg.nday = 31 + 15  # February, 15

    testuwg.generate()

    # We subtract 11 hours from total timestep so still have sun. New time: 1300
    testuwg.simTime.nt -= 12 * 11

    testuwg.simulate()

    # check date is February 15th, 13:00
    assert testuwg.simTime.month == 2
    assert testuwg.simTime.day == 15
    assert testuwg.simTime.secDay/3600. == pytest.approx(13.0, abs=1e-15)

    uwg_python_val = [
        testuwg.rural.aeroCond,        # Convection coef (refL uwg, eq.12)
        testuwg.rural.waterStorage,    # thickness of water film (m) (only for hor surf)
        testuwg.rural.solAbs,          # solar radiation absorbed (W m-2)
        testuwg.rural.lat,             # surface latent heat flux (W m-2)
        testuwg.rural.sens,            # surface sensible heat flux (W m-2)
        testuwg.rural.flux,            # external surface heat flux (W m-2)
        testuwg.rural.T_ext,           # external surface temperature (K)
        testuwg.rural.T_int]           # interior surface temperature (K)

    # Matlab Checking for rural road
    uwg_matlab_val = \
        setup_open_matlab_ref('matlab_element', 'matlab_ref_element_surfflux.txt')

    # matlab ref checking
    assert len(uwg_matlab_val) == len(uwg_python_val)

    for i in range(len(uwg_matlab_val)):
        tol = calculate_tolerance(uwg_matlab_val[i], 15.0)
        assert uwg_python_val[i] == pytest.approx(uwg_matlab_val[i], abs=tol), \
            'error at index={}'.format(i)
