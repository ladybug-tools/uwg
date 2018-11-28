import os
import pytest
import uwg
import math
from .test_base import TestBase

from pprint import pprint
from decimal import Decimal
pp = pprint
dd = Decimal.from_float

class TestRequiredParams(TestBase):
    """Test for uwg.py
    """

    def test_required_inputs_from_file(self):
        # From a .uwg file
        self.setup_uwg_integration("SGP_Singapore.486980_IWEC.epw", "initialize_singapore.uwg")
        self.uwg.set_input()
        self.uwg.check_required_inputs()

    def test_required_inputs_wrong_type(self):
        # Test that we can catch wrong types, list lengths etc
        self.setup_uwg_integration("SGP_Singapore.486980_IWEC.epw", "initialize_singapore.uwg")
        self.uwg.set_input()

        with pytest.raises(AssertionError):
            self.uwg.h_temp = 6
            self.uwg.check_required_inputs()
        self.uwg.h_temp = 2.0
        self.uwg.check_required_inputs()

        with pytest.raises(AssertionError):
            self.uwg.windMin = 100
            self.uwg.check_required_inputs()
        self.uwg.windMin = 2.0
        self.uwg.check_required_inputs()

        # assert type(self.Day) == float or type(self.Day) == int
        with pytest.raises(AssertionError):
            self.uwg.Day = "a"
            self.uwg.check_required_inputs()
        self.uwg.Day = 2
        self.uwg.check_required_inputs()
        self.uwg.Day = 4.0
        self.uwg.check_required_inputs()

        # assert isinstance(self.SchTraffic, list)
        with pytest.raises(AssertionError):
            self.uwg.SchTraffic = 76
            self.uwg.check_required_inputs()
        self.uwg.SchTraffic = [0] * 3 # length 3 list
        self.uwg.check_required_inputs()

        # assert len(self.bld) == 16
        with pytest.raises(AssertionError):
            self.uwg.bld = 4543
            self.uwg.check_required_inputs()
        self.uwg.bld = [0] * 16
        self.uwg.check_required_inputs()



if __name__ == "__main__":
    pass
    #testreq = TestRequiredParams()
    #testreq.test_required_inputs_from_file()