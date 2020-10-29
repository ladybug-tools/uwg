"""Tests for SolarCalcs object."""

import pytest
from .test_base import auto_setup_uwg, setup_open_matlab_ref
from uwg import SolarCalcs


def test_solarangles():
    """Test solar angles."""

    model = auto_setup_uwg()
    model.generate()

    solar = SolarCalcs(model.UCM, model.BEM, model.simTime, model.RSM,
                       model.forc, model.geoParam, model.rural)

    # Timestep every 5 minutes (300s)
    for i in range(int(12 * 24 * 1.5)):  # 1.5 days, Jan 2nd, 12 noon
        solar.simTime.update_date()
    solar.solarangles()

    # test simtime for solar after 1.5 days
    assert solar.ut == pytest.approx(12.0, abs=1e-15)
    assert solar.simTime.day == pytest.approx(2.0, abs=1e-15)
    assert solar.ad == pytest.approx(0.197963373, abs=1e-8)

    # Recalculate solar angles with new time simulation, add to previous time increment
    # Increment time to 22 days, 13hrs, 30min = Jan 22 at 1330
    for i in range(12 * 24 * 20 + 12 * 1 + 6):
        solar.simTime.update_date()

    # Run simulation
    solar.solarangles()

    uwg_python_val = [
        solar.simTime.month,
        solar.simTime.secDay,
        solar.simTime.day,
        solar.simTime.hourDay,
        solar.ut,
        solar.ad,
        solar.eqtime,
        solar.decsol,
        solar.zenith,
        solar.tanzen,
        solar.critOrient]

    # open matlab ref file
    uwg_matlab_val = \
        setup_open_matlab_ref('matlab_solarcalcs', 'matlab_ref_solarangles.txt')

    # matlab ref checking
    assert len(uwg_matlab_val) == len(uwg_python_val)
    for i in range(len(uwg_matlab_val)):
        assert uwg_python_val[i] == \
            pytest.approx(uwg_matlab_val[i], abs=1e-15), 'error at index={}'.format(i)


def test_solarcalcs():
    """Test solar calculation."""

    model = auto_setup_uwg()
    model.generate()

    # We subtract 11 hours from total timestep so
    # we can stop simulation while we still have sun!
    # New time: Jan 31, 1300
    model.simTime.nt -= 12 * 11

    model.simulate()

    # check date
    assert model.simTime.month == 1
    assert model.simTime.day == 31
    assert model.simTime.secDay / 3600. == pytest.approx(13.0, abs=1e-15)

    # Manual checking

    # Total Radiation (W m-2) at 13:00 on Jan 31 w/o reflected sunlight
    assert model.solar.bldSol == pytest.approx(170.4189961148073, abs=1e-12)
    assert model.solar.roadSol == pytest.approx(307.4032644249752, abs=1e-12)

    # Matlab Checking
    uwg_matlab_val = setup_open_matlab_ref(
        'matlab_solarcalcs', 'matlab_ref_solarcalcs.txt')

    # Modify treeSensHeat based on changes to veg/tree heat balance 10/20
    # grassSensHeat = (
    #     (1 - model.geoParam.vegAlbedo) * (1 - model.geoParam.grassFLat) *
    #     model.UCM.SolRecRoad * (model.grasscover / (1 - model.blddensity)))
    # vegSensHeat = \
    #     (model.solar.UCM.treeSensHeat - grassSensHeat) / model.UCM.treeCoverage
    treeSensHeat = (
        (1 - model.solar.parameter.vegAlbedo) * (1 - model.solar.parameter.treeFLat) *
        model.solar.UCM.SolRecRoad)
    treeLatHeat = (
        (1 - model.solar.parameter.vegAlbedo) * model.solar.parameter.treeFLat *
        model.solar.UCM.SolRecRoad)

    uwg_python_val = [
        model.solar.horSol,  # 0
        model.solar.Kw_term,  # 1
        model.solar.Kr_term,  # 2
        model.solar.mr,  # 3
        model.solar.mw,  # 4
        model.solar.BEM[0].roof.solRec,  # 5
        model.solar.BEM[0].wall.solRec,  # 6
        model.solar.BEM[1].roof.solRec,  # 7
        model.solar.BEM[1].wall.solRec,  # 8
        model.solar.rural.solRec,  # 9
        model.solar.UCM.SolRecRoof,  # 10
        model.solar.UCM.SolRecRoad,  # 11
        model.solar.UCM.SolRecWall,  # 12
        treeSensHeat,  # 13
        treeLatHeat]  # 14

    # matlab ref checking
    assert len(uwg_matlab_val) == len(uwg_python_val)
    for i in range(len(uwg_matlab_val)):
        assert uwg_python_val[i] == \
            pytest.approx(uwg_matlab_val[i], abs=1e-15), 'error at index={}'.format(i)
