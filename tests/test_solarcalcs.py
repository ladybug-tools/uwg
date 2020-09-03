"""Tests for SolarCalcs object."""

import pytest
from .test_base import setup_uwg_integration, setup_open_matlab_ref
from uwg import SolarCalcs


def test_solarangles():
    """Test solar angles."""

    testuwg = setup_uwg_integration()
    testuwg.read_input()
    testuwg.generate()

    solar = SolarCalcs(testuwg.UCM, testuwg.BEM, testuwg.simTime, testuwg.RSM,
                       testuwg.forc, testuwg.geoParam, testuwg.rural)

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

    testuwg = setup_uwg_integration()
    testuwg.read_input()
    testuwg.generate()

    # We subtract 11 hours from total timestep so
    # we can stop simulation while we still have sun!
    # New time: Jan 31, 1300
    testuwg.simTime.nt -= 12 * 11

    testuwg._hvac_autosize()
    testuwg.simulate()

    # check date
    assert testuwg.simTime.month == 1
    assert testuwg.simTime.day == 31
    assert testuwg.simTime.secDay / 3600. == pytest.approx(13.0, abs=1e-15)

    # Manual checking

    # Total Radiation (W m-2) at 13:00 on Jan 31 w/o reflected sunlight
    assert testuwg.solar.bldSol == pytest.approx(170.4189961148073, abs=1e-12)
    assert testuwg.solar.roadSol == pytest.approx(307.4032644249752, abs=1e-12)

    # Matlab Checking
    uwg_matlab_val = setup_open_matlab_ref(
        'matlab_solarcalcs', 'matlab_ref_solarcalcs.txt')

    uwg_python_val = [
        testuwg.solar.horSol,
        testuwg.solar.Kw_term,
        testuwg.solar.Kr_term,
        testuwg.solar.mr,
        testuwg.solar.mw,
        testuwg.solar.BEM[0].roof.solRec,
        testuwg.solar.BEM[0].wall.solRec,
        testuwg.solar.BEM[1].roof.solRec,
        testuwg.solar.BEM[1].wall.solRec,
        testuwg.solar.rural.solRec,
        testuwg.solar.UCM.SolRecRoof,
        testuwg.solar.UCM.SolRecRoad,
        testuwg.solar.UCM.SolRecWall,
        testuwg.solar.UCM.treeSensHeat,
        testuwg.solar.UCM.treeLatHeat]

    # matlab ref checking
    assert len(uwg_matlab_val) == len(uwg_python_val)
    for i in range(len(uwg_matlab_val)):
        assert uwg_python_val[i] == \
            pytest.approx(uwg_matlab_val[i], abs=1e-15), 'error at index={}'.format(i)
