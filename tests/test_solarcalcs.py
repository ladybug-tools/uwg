"""Tests for SolarCalcs object."""

import json
import os

import pytest
from .test_base import auto_setup_uwg
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
    assert solar.ad == pytest.approx(2 * 3.141592653589793 / 365, abs=1e-12)

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

    # The historical MATLAB angles encode the equations corrected against NOAA.
    # Keep the fixed-time geometry regression; the external equation checks are
    # in test_solar_noaa_regression.py.
    with open(os.path.join(os.path.dirname(__file__), 'solar_corrected_reference.json')) as stream:
        solar_reference = json.load(stream)['angles']

    # corrected solar reference checking
    assert len(solar_reference) == len(uwg_python_val)
    for i in range(len(solar_reference)):
        assert uwg_python_val[i] == \
            pytest.approx(solar_reference[i], abs=1e-15), 'error at index={}'.format(i)


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
    assert model.solar.bldSol == pytest.approx(157.56521874370145, abs=1e-12)
    assert model.solar.roadSol == pytest.approx(351.3780296908984, abs=1e-12)

    # Only directly affected solar flux point expectations are updated here.
    # Other building/model MATLAB reference files are unchanged.
    with open(os.path.join(os.path.dirname(__file__), 'solar_corrected_reference.json')) as stream:
        solar_reference = json.load(stream)['solar']

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

    # corrected solar reference checking
    assert len(solar_reference) == len(uwg_python_val)
    for i in range(len(solar_reference)):
        assert uwg_python_val[i] == \
            pytest.approx(solar_reference[i], abs=1e-15), 'error at index={}'.format(i)
