"""Tests for BEM object."""

import pytest
from .test_base import setup_uwg_integration, setup_open_matlab_ref, calculate_tolerance


def test_bem_building_init_largeoffice():
    """Test for BEM.building init"""

    testuwg = setup_uwg_integration()
    testuwg.read_input()
    testuwg.generate()

    uwg_python_val = [
        testuwg.BEM[0].building.floorHeight,
        testuwg.BEM[0].building.intHeat,
        testuwg.BEM[0].building.intHeatNight,
        testuwg.BEM[0].building.intHeatDay,
        testuwg.BEM[0].building.intHeatFRad,
        testuwg.BEM[0].building.intHeatFLat,
        testuwg.BEM[0].building.infil,
        testuwg.BEM[0].building.vent,
        testuwg.BEM[0].building.glazingRatio,
        testuwg.BEM[0].building.uValue,
        testuwg.BEM[0].building.shgc,
        testuwg.BEM[0].building.condType,
        testuwg.BEM[0].building.cop,
        testuwg.BEM[0].building.coolSetpointDay,
        testuwg.BEM[0].building.coolSetpointNight,
        testuwg.BEM[0].building.heatSetpointDay,
        testuwg.BEM[0].building.heatSetpointNight,
        testuwg.BEM[0].building.coolCap,
        testuwg.BEM[0].building.heatCap,
        testuwg.BEM[0].building.heatEff,
        testuwg.BEM[0].building.mSys,
        testuwg.BEM[0].building.indoorTemp,
        testuwg.BEM[0].building.indoorHum,
        testuwg.BEM[0].building.FanMax,
        testuwg.BEM[0].building.Type,
        testuwg.BEM[0].building.Era,
        testuwg.BEM[0].building.Zone]

    uwg_matlab_val = setup_open_matlab_ref(
        'matlab_bem', 'matlab_ref_bem_building_init_largeoffice.txt')

    # matlab ref checking
    assert len(uwg_matlab_val) == len(uwg_python_val)

    for i in range(len(uwg_matlab_val)):
        if isinstance(uwg_python_val[i], str):
            assert ''.join(uwg_python_val[i].split()) == \
                ''.join(uwg_matlab_val[i].split()), 'error at index={}'.format(i)
        else:
            tol = calculate_tolerance(uwg_python_val[i], 14.0)
            assert uwg_python_val[i] == \
                pytest.approx(uwg_matlab_val[i], abs=tol), 'error at index={}'.format(i)


def test_bem_building_bemcalc_largeoffice_cooling():
    """Test for bem.building bemcalc during cooling period."""

    testuwg = setup_uwg_integration()
    testuwg._read_epw()
    testuwg.read_input()

    # Test Jan 1 (winter, no vegetation coverage)
    testuwg.month = 1
    testuwg.day = 1
    testuwg.nday = 1

    # set_input
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

    # Calculated values
    # commented out properties are instantiated building variables that are never
    # actually set in UWG_Matlab
    uwg_python_val = [
        # changes from constructor
        testuwg.BEM[0].building.intHeat,
        testuwg.BEM[0].building.intHeatNight,
        testuwg.BEM[0].building.intHeatDay,
        testuwg.BEM[0].building.indoorTemp,
        testuwg.BEM[0].building.indoorHum,
        # new values
        # testuwg.BEM[0].building.Tdp,
        testuwg.BEM[0].building.indoorRhum,
        testuwg.BEM[0].building.nFloor,
        # testuwg.BEM[0].building.RadFOcc,
        # testuwg.BEM[0].building.LatFOcc,
        # testuwg.BEM[0].building.RadFEquip,
        # testuwg.BEM[0].building.RadFLight,
        testuwg.BEM[0].building.sensCoolDemand,
        testuwg.BEM[0].building.sensHeatDemand,
        testuwg.BEM[0].building.copAdj,
        testuwg.BEM[0].building.dehumDemand,
        testuwg.BEM[0].building.coolConsump,
        testuwg.BEM[0].building.heatConsump,
        testuwg.BEM[0].building.sensWaste,
        # testuwg.BEM[0].building.latWaste,
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

    uwg_matlab_val = setup_open_matlab_ref(
        'matlab_bem', 'matlab_ref_bem_building_bemcalc_largeoffice_cooling.txt')

    # matlab ref checking
    assert len(uwg_matlab_val) == len(uwg_python_val)

    for i in range(len(uwg_matlab_val)):
        tol = calculate_tolerance(uwg_python_val[i], 15.0)
        assert uwg_python_val[i] == \
            pytest.approx(uwg_matlab_val[i], abs=tol), 'error at index={}'.format(i)
