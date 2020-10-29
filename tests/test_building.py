"""Tests for building.py"""

import pytest
from uwg import Building


def test_building_init():
    """Test BEMDef init method."""

    # LargeOffce, Pst80, Zone 1A (Miami)
    # coolCap for UTM ~ 200 W/m2, U ~ 2.0
    bld = Building(floor_height=3.5, int_heat_night=1, int_heat_day=1, int_heat_frad=0.1,
                   int_heat_flat=0.1, infil=0.26, vent=0.0005, glazing_ratio=0.4,
                   u_value=5.8, shgc=0.2, condtype='AIR', cop=5.2,
                   coolcap=76, heateff=0.7, initial_temp=293)
    # test repr
    bld.__repr__()


def test_building_dict():
    """Test building dict methods."""

    # init
    bld1 = Building(floor_height=3.5, int_heat_night=1, int_heat_day=1, int_heat_frad=0.1,
                    int_heat_flat=0.1, infil=0.26, vent=0.0005, glazing_ratio=0.4,
                    u_value=5.8, shgc=0.2, condtype='AIR', cop=5.2,
                    coolcap=76, heateff=0.7, initial_temp=293)

    # make dict
    blddict = bld1.to_dict()

    # test if dict and from_dict
    assert isinstance(blddict, dict)

    bld2 = Building.from_dict(blddict)

    assert bld1.floor_height == pytest.approx(bld2.floor_height, abs=1e-10)
    assert bld1.u_value == pytest.approx(bld2.u_value, abs=1e-10)
    assert bld1.coolcap == pytest.approx(bld2.coolcap, abs=1e-10)
    assert bld1.vent == pytest.approx(bld2.vent, abs=1e-10)

    with pytest.raises(AssertionError):
        blddict['type'] = 'BemDef'
        Building.from_dict(blddict)

