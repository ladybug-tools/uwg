"""Tests for BEM object."""

# Fix for Python 2.X
try:
    unicode
except NameError:
    unicode = str

import pytest
from .test_base import auto_setup_uwg, setup_open_matlab_ref, calculate_tolerance
from uwg import BEMDef, Building, Element, Material


def _bemdef():
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
                   u_value=5.8, shgc=0.2, condtype='AIR', cop=5.2,
                   coolcap=76, heateff=0.7, initial_temp=293)

    return BEMDef(bld, floor, wall, roof, bldtype='fullservicerestaurant',
                  builtera='pst80')


def test_bem_init_manual():
    """Test BEMDef init methods."""

    # test init
    bem = _bemdef()

    # test repr
    bem.__repr__()


def test_bemdef_dict():
    """Test BEMDef dict methods."""

    bem = _bemdef()

    # make dict
    bemdict = bem.to_dict()

    # test if dict and from_dict
    assert isinstance(bemdict, dict)

    bem2 = BEMDef.from_dict(bemdict)

    assert bem.wall.albedo == pytest.approx(bem2.wall.albedo, abs=1e-10)
    assert bem.frac == pytest.approx(bem2.frac, abs=1e-10)
    assert bem.building.cop == pytest.approx(bem2.building.cop, abs=1e-10)
    assert bem.mass.layerThermalCond[0] == \
        pytest.approx(bem2.mass.layerThermalCond[0], abs=1e-10)

    with pytest.raises(AssertionError):
        bemdict['type'] = 'BemDef'
        BEMDef.from_dict(bemdict)


def test_bem_building_init_largeoffice_readDOE():
    """Test for BEM.building init"""

    testuwg = auto_setup_uwg()
    testuwg.generate()

    uwg_python_val = [
        testuwg.BEM[0].building.floor_height,
        testuwg.BEM[0].building.int_heat,
        testuwg.BEM[0].building.int_heat_night,
        testuwg.BEM[0].building.int_heat_day,
        testuwg.BEM[0].building.int_heat_frad,
        testuwg.BEM[0].building.int_heat_flat,
        testuwg.BEM[0].building.infil,
        testuwg.BEM[0].building.vent,
        testuwg.BEM[0].building.glazing_ratio,
        testuwg.BEM[0].building.u_value,
        testuwg.BEM[0].building.shgc,
        testuwg.BEM[0].building.condtype,
        testuwg.BEM[0].building.cop,
        testuwg.BEM[0].building.cool_setpoint_day,
        testuwg.BEM[0].building.cool_setpoint_night,
        testuwg.BEM[0].building.heat_setpoint_day,
        testuwg.BEM[0].building.heat_setpoint_night,
        testuwg.BEM[0].building.coolcap,
        testuwg.BEM[0].building.heat_cap,
        testuwg.BEM[0].building.heateff,
        testuwg.BEM[0].building.msys,
        testuwg.BEM[0].building.indoor_temp,
        testuwg.BEM[0].building.indoor_hum,
        testuwg.BEM[0].building.FanMax,
        testuwg.BEM[0].bldtype,
        testuwg.BEM[0].builtera,
        testuwg.BEM[0].zonetype]

    uwg_matlab_val = setup_open_matlab_ref(
        'matlab_bem', 'matlab_ref_bem_building_init_largeoffice.txt')

    # matlab ref checking
    assert len(uwg_matlab_val) == len(uwg_python_val)
    for i in range(len(uwg_matlab_val)):
        if isinstance(uwg_python_val[i], (str, unicode)):
            assert ''.join(uwg_python_val[i].lower().split()) == \
                ''.join(uwg_matlab_val[i].lower().split()), \
                'error at index={}'.format(i)
        else:
            tol = calculate_tolerance(uwg_python_val[i], 14.0)
            assert uwg_python_val[i] == \
                pytest.approx(
                    uwg_matlab_val[i], abs=tol), 'error at index={}'.format(i)


def test_bem_building_bemcalc_largeoffice_cooling_readDOE():
    """Test for bem.building bemcalc during cooling period."""

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
        testuwg.BEM[0].building.int_heat,
        testuwg.BEM[0].building.int_heat_night,
        testuwg.BEM[0].building.int_heat_day,
        testuwg.BEM[0].building.indoor_temp,
        testuwg.BEM[0].building.indoor_hum,
        testuwg.BEM[0].building.indoorRhum,
        testuwg.BEM[0].building.nFloor,
        testuwg.BEM[0].building.sensCoolDemand,
        testuwg.BEM[0].building.sensHeatDemand,
        testuwg.BEM[0].building.cop_adj,
        testuwg.BEM[0].building.dehumDemand,
        testuwg.BEM[0].building.coolConsump,
        testuwg.BEM[0].building.heatConsump,
        testuwg.BEM[0].building.sensWaste,
        testuwg.BEM[0].building.fluxMass,
        testuwg.BEM[0].building.fluxWall,
        testuwg.BEM[0].building.fluxRoof,
        testuwg.BEM[0].building.fluxSolar,
        testuwg.BEM[0].building.fluxWindow,
        testuwg.BEM[0].building.fluxInterior,
        testuwg.BEM[0].building.fluxInfil,
        testuwg.BEM[0].building.fluxVent,
        testuwg.BEM[0].building.ElecTotal,
        testuwg.BEM[0].building.GasTotal,
        testuwg.BEM[0].building.Qhvac,
        testuwg.BEM[0].building.Qheat]

    # Note: the following properties are instantiated but never set in UWG_Matlab
    # testuwg.BEM[0].building.latWaste
    # testuwg.BEM[0].building.RadFOcc
    # testuwg.BEM[0].building.LatFOcc
    # testuwg.BEM[0].building.RadFEquip
    # testuwg.BEM[0].building.RadFLight
    # testuwg.BEM[0].building.Tdp

    uwg_matlab_val = setup_open_matlab_ref(
        'matlab_bem', 'matlab_ref_bem_building_bemcalc_largeoffice_cooling.txt')

    # matlab ref checking
    assert len(uwg_matlab_val) == len(uwg_python_val)

    for i in range(len(uwg_matlab_val)):
        tol = calculate_tolerance(uwg_python_val[i], 15.0)
        assert uwg_python_val[i] == \
            pytest.approx(uwg_matlab_val[i],
                          abs=tol), 'error at index={}'.format(i)
