"""Tests for building.py"""

import pytest
from uwg import Building


def test_building_init():
    """Test BEMDef init method."""

    # LargeOffce, Pst80, Zone 1A (Miami)
    # coolCap for UTM ~ 200 W/m2, U ~ 2.0
    bld = Building(floorHeight=3.5, intHeatNight=1, intHeatDay=1, intHeatFRad=0.1,
                   intHeatFLat=0.1, infil=0.26, vent=0.0005, glazingRatio=0.4,
                   uValue=5.8, shgc=0.2, condType='AIR', cop=5.2, coolSetpointDay=297,
                   coolSetpointNight=297, heatSetpointDay=293, heatSetpointNight=293,
                   coolCap=76, heatEff=0.7, initialTemp=293)
    # test repr
    bld.__repr__()


def test_building_dict():
    """Test building dict methods."""

    # init
    bld1 = Building(floorHeight=3.5, intHeatNight=1, intHeatDay=1, intHeatFRad=0.1,
                    intHeatFLat=0.1, infil=0.26, vent=0.0005, glazingRatio=0.4,
                    uValue=5.8, shgc=0.2, condType='AIR', cop=5.2, coolSetpointDay=297,
                    coolSetpointNight=297, heatSetpointDay=293, heatSetpointNight=293,
                    coolCap=76, heatEff=0.7, initialTemp=293)

    # make dict
    blddict = bld1.to_dict()

    # test if dict and from_dict
    assert isinstance(blddict, dict)

    bld2 = Building.from_dict(blddict)

    assert bld1.floorHeight == pytest.approx(bld2.floorHeight, abs=1e-10)
    assert bld1.uValue == pytest.approx(bld2.uValue, abs=1e-10)
    assert bld1.coolCap == pytest.approx(bld2.coolCap, abs=1e-10)
    assert bld1.vent == pytest.approx(bld2.vent, abs=1e-10)

    with pytest.raises(AssertionError):
        blddict['type'] = 'BemDef'
        Building.from_dict(blddict)

