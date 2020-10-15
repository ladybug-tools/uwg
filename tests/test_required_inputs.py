"""Test for uwg.py"""
import os
import pytest
from .test_base import auto_setup_uwg

DIR_CURR = os.path.abspath(os.path.dirname(__file__))
DEFAULT_EPW_PATH = os.path.join(
    DIR_CURR, 'epw', 'SGP_Singapore.486980_IWEC.epw')
DEFAULT_PARAM_PATH = os.path.join(
    DIR_CURR, 'parameters', 'initialize_singapore.uwg')


def test_required_inputs_wrong_type():
    """Test that we can catch wrong types, list lengths etc."""

    testuwg = auto_setup_uwg(DEFAULT_EPW_PATH, DEFAULT_PARAM_PATH)

    with pytest.raises(AssertionError):
        testuwg.h_temp = -6
    testuwg.h_temp = 2.0

    with pytest.raises(AssertionError):
        testuwg.windmin = -100
    testuwg.windmin = 2.0

    # assert type(self.Day) == float or type(self.Day) == int
    with pytest.raises(ValueError):
        testuwg.day = 'a'

    testuwg.day = 2
    testuwg.day = 4.0

    # assert isinstance(self.SchTraffic, list)
    with pytest.raises(AssertionError):
        testuwg.schtraffic = 76

    testuwg.schtraffic = [[0.1 for i in range(24)]] * 3  # length 3 list

    with pytest.raises(AssertionError):
        testuwg.bld = 4543

    testuwg.bld = [['test', 'pre80', 0.1], ['test2', 'pre80', 0.9]]
