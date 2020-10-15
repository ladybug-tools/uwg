"""Test for schdef.py"""

import pytest
from uwg import SchDef


DEFAULT_WEEK = [
    [0.2, 0.2, 0.2, 0.2, 0.2, 0.4, 0.7, 0.9, 0.9, 0.6, 0.6, 0.6, 0.6, 0.6,
     0.7, 0.8, 0.9, 0.9, 0.8, 0.8, 0.7, 0.3, 0.2, 0.2],  # Weekday
    [0.2, 0.2, 0.2, 0.2, 0.2, 0.3, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
        0.6, 0.7, 0.7, 0.7, 0.7, 0.5, 0.4, 0.3, 0.2, 0.2],  # Saturday
    [0.2, 0.2, 0.2, 0.2, 0.2, 0.3, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
        0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.3, 0.3, 0.2, 0.2]]  # Sunday


def test_schdef_init():
    """Test SchDef init method."""

    with pytest.raises(AssertionError):
        SchDef(elec=DEFAULT_WEEK, gas=DEFAULT_WEEK, light=DEFAULT_WEEK, occ=DEFAULT_WEEK,
               cool=DEFAULT_WEEK, heat=DEFAULT_WEEK, swh={}, q_elec=18.9, q_gas=3.2,
               q_light=18.9, n_occ=0.12, vent=0.0013, v_swh=0.2846,
               bldtype='fullservicerestaurant', builtera='pre80')

    with pytest.raises(AssertionError):
        testweek = [[0.1] * 24] * 2
        SchDef(elec=DEFAULT_WEEK, gas=DEFAULT_WEEK, light=DEFAULT_WEEK, occ=DEFAULT_WEEK,
               cool=DEFAULT_WEEK, heat=testweek, swh=DEFAULT_WEEK, q_elec=18.9,
               q_gas=3.2, q_light=18.9, n_occ=0.12, vent=0.0013, v_swh=0.2846,
               bldtype='fullservicerestaurant', builtera='pre80')

    with pytest.raises(AssertionError):
        testweek = [[0.1] * 24] * 2
        testweek += [0.1] * 23
        SchDef(elec=DEFAULT_WEEK, gas=DEFAULT_WEEK, light=DEFAULT_WEEK, occ=DEFAULT_WEEK,
               cool=DEFAULT_WEEK, heat=testweek, swh=DEFAULT_WEEK, q_elec=18.9,
               q_gas=3.2, q_light=18.9, n_occ=0.12, vent=0.0013, v_swh=0.2846,
               bldtype='fullservicerestaurant', builtera='pre80')

    with pytest.raises(AssertionError):
        testweek = [[0.1] * 24] * 3
        testweek[2][23] = 'a'
        SchDef(elec=DEFAULT_WEEK, gas=DEFAULT_WEEK, light=DEFAULT_WEEK, occ=DEFAULT_WEEK,
               cool=DEFAULT_WEEK, heat=testweek, swh=DEFAULT_WEEK, q_elec=18.9,
               q_gas=3.2, q_light=18.9, n_occ=0.12, vent=0.0013, v_swh=0.2846,
               bldtype='fullservicerestaurant', builtera='pre80')

    # init
    sch = SchDef(elec=DEFAULT_WEEK, gas=DEFAULT_WEEK, light=DEFAULT_WEEK,
                 occ=DEFAULT_WEEK, cool=DEFAULT_WEEK, heat=DEFAULT_WEEK,
                 swh=DEFAULT_WEEK, q_elec=18.9, q_gas=3.2, q_light=18.9, n_occ=0.12,
                 vent=0.0013, v_swh=0.2846, bldtype='fullservicerestaurant',
                 builtera='pre80')

    # test __repr__
    sch.__repr__()


def test_schdef_dict():
    """Test SchDef dict methods."""

    # init
    testsch = SchDef(
        elec=DEFAULT_WEEK, gas=DEFAULT_WEEK, light=DEFAULT_WEEK, occ=DEFAULT_WEEK,
        cool=DEFAULT_WEEK, heat=DEFAULT_WEEK, swh=DEFAULT_WEEK, q_elec=18.9,
        q_gas=3.2, q_light=18.9, n_occ=0.12, vent=0.0013, v_swh=0.2846,
        bldtype='fullservicerestaurant', builtera='pre80')

    # make dict
    schdict = testsch.to_dict()

    # test if dict and from_dict
    assert isinstance(schdict, dict)

    testsch2 = SchDef.from_dict(schdict)

    assert testsch.elec == testsch2.elec
    assert testsch.gas == testsch2.gas
    assert testsch.light == testsch2.light
    assert testsch.occ == testsch2.occ
    assert testsch.cool == testsch2.cool
    assert testsch.heat == testsch2.heat
    assert testsch.swh == testsch2.swh

    with pytest.raises(AssertionError):
        schdict['type'] = 'BemDef'
        SchDef.from_dict(schdict)
