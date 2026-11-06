"""Unavailable cooling must not remove sensible heat from the room."""
from types import SimpleNamespace as NS

import pytest

from uwg import Building


def cutoff_case(canyon=285., capacity=1000., floors=1, internal=200., surface_temp=297.):
    building = Building(3., internal, internal, 0., 0., 0., 0., 0., 1., .5,
                        'AIR', 4., capacity, .8, 300.)
    building.int_heat_f_rad = 0.  # Normally initialized by UWG.simulate.
    building.indoor_hum = .00702  # No latent cooling in this sensible balance.
    surface = NS(layerTemp=[surface_temp], solRec=0.)
    bem = NS(wall=surface, roof=surface, mass=surface, swh=0., elec=0., light=0., gas=0.)
    ucm = NS(bldHeight=3. * floors, canTemp=canyon, canHum=.00702,
             verToHor=1., bldDensity=.5)
    forcing = NS(pres=101325., waterTemp=285.)
    parameter = NS(nightSetStart=22., nightSetEnd=6., lv=2500000., cp=1004.)
    time = NS(secDay=43200., dt=300., month=7, day=1)
    return building, ucm, bem, forcing, parameter, time


def step(case):
    building, ucm, bem, forcing, parameter, time = case
    building.BEMCalc(ucm, bem, forcing, parameter, time)
    return building


def conductance(floors):
    # Per footprint: wall area 2, mass area 2*N-1, roof area 1.
    # Fixed-surface, no-window/no-airflow case; ceiling is cooler than initial air.
    return 2. * 3.076 + (2. * floors - 1.) * 3.076 + 4.040


@pytest.mark.parametrize('canyon', [285., 287.999999, 288.])
@pytest.mark.parametrize('capacity', [0., 20., 1000.])
@pytest.mark.parametrize('floors', [1, 3])
def test_unavailable_cooling_keeps_requested_load_without_removing_heat(canyon, capacity, floors):
    building = step(cutoff_case(canyon, capacity, floors))
    assert building.nFloor == floors
    assert building.sensCoolDemand == pytest.approx(200.)  # Requested per floor.
    assert building.Qhvac == 0.
    assert building.coolConsump == 0.
    assert building.sensWaste == 0.
    assert building.Qheat == 0.
    assert building.indoor_hum == pytest.approx(.00702, abs=1e-14)
    expected = 297. + 200. * floors / conductance(floors)
    assert building.indoor_temp == pytest.approx(expected, abs=1e-10)
    assert conductance(floors) * (building.indoor_temp - 297.) == pytest.approx(200. * floors, abs=1e-9)


@pytest.mark.parametrize('canyon', [288.000001, 300.])
@pytest.mark.parametrize('capacity', [20., 1000.])
@pytest.mark.parametrize('floors', [1, 3])
def test_active_cooling_retains_capacity_power_waste_and_room_balance(canyon, capacity, floors):
    building = step(cutoff_case(canyon, capacity, floors))
    delivered = min(200., capacity) * floors
    assert building.Qhvac == pytest.approx(delivered)
    assert building.sensCoolDemand * floors == pytest.approx(delivered)
    assert building.dehumDemand == pytest.approx(0.)
    assert building.coolConsump * floors == pytest.approx(delivered / 4.)
    assert building.sensWaste == pytest.approx(delivered * 1.25)
    assert building.indoor_temp == pytest.approx(297. + (200. * floors - delivered) / conductance(floors), abs=1e-10)


@pytest.mark.parametrize('canyon', [285., 288., 300.])
def test_zero_cooling_demand_control(canyon):
    building = step(cutoff_case(canyon, internal=0.))
    assert building.sensCoolDemand == 0.
    assert building.Qhvac == 0.
    assert building.indoor_temp == pytest.approx(297., abs=1e-10)


def test_active_to_unavailable_to_active_does_not_reuse_delivered_cooling():
    case = cutoff_case(300.)
    building = step(case)
    assert building.indoor_temp == pytest.approx(297.)
    case[1].canTemp = 288.
    for _ in range(2):
        building = step(case)
        assert building.Qhvac == 0.
        assert building.sensCoolDemand == pytest.approx(200.)
        assert building.indoor_temp == pytest.approx(297. + 200. / conductance(1))
    case[1].canTemp = 300.
    building = step(case)
    assert building.Qhvac == pytest.approx(200.)
    assert building.indoor_temp == pytest.approx(297.)


@pytest.mark.parametrize('canyon, heats', [(285., True), (288., False), (300., False)])
def test_existing_heating_cutoff_is_preserved(canyon, heats):
    building = step(cutoff_case(canyon, internal=0., surface_temp=285.))
    assert building.sensCoolDemand == 0.
    assert building.Qhvac == 0.
    assert building.indoor_temp == pytest.approx(293. if heats else 285.)
    assert (building.Qheat > 0.) == heats
