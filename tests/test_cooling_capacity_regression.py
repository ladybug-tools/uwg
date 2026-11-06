"""Capacity and coil-energy conservation with synthetic boundary conditions."""
from types import SimpleNamespace

import pytest
from uwg import Building
from uwg.psychrometrics import moist_air_density


@pytest.mark.parametrize('capacity', [5, 20, 100, 1000])
@pytest.mark.parametrize('floors', [1, 3])
def test_cooling_capacity_conserves_sensible_and_latent_energy(capacity, floors):
    building = Building(3, 200, 200, 0, 0, 0, 0, 0, 1, .5, 'AIR', 4, capacity, .8, 300)
    building.int_heat_f_rad = 0  # field set by UWG's simulation orchestrator
    initial_humidity = building.indoor_hum
    density = moist_air_density(101325, 300, initial_humidity)
    surface = SimpleNamespace(layerTemp=[297.], solRec=0.)
    bem = SimpleNamespace(wall=surface, roof=surface, mass=surface,
                          swh=0., elec=0., light=0., gas=0.)
    canyon = SimpleNamespace(bldHeight=3*floors, canTemp=300.,
                             canHum=initial_humidity, verToHor=1., bldDensity=.5)
    parameters = SimpleNamespace(nightSetStart=22, nightSetEnd=6, lv=2500000., cp=1004.)
    forcing = SimpleNamespace(pres=101325., waterTemp=285.)
    time = SimpleNamespace(secDay=43200., dt=300., month=7, day=1)
    building.BEMCalc(canyon, bem, forcing, parameters, time)
    sensible = building.sensCoolDemand * building.nFloor
    total = sensible + building.dehumDemand
    # All quantities below are per footprint, unlike coolConsump/sensCoolDemand.
    assert total == pytest.approx(building.Qhvac)
    assert total <= capacity * floors + 1e-10
    assert building.coolConsump * floors * building.cop_adj == pytest.approx(total)
    assert building.sensWaste == pytest.approx(total * (1 + 1 / building.cop_adj))
    latent_from_moisture = ((initial_humidity - building.indoor_hum) * density *
                           parameters.lv * canyon.bldHeight / time.dt)
    assert building.dehumDemand == pytest.approx(latent_from_moisture)
