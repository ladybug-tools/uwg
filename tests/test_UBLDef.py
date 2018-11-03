try:
    range = xrange
except NameError:
    pass

from functools import reduce

import pytest
import uwg
import os
import math
import pprint
import decimal
from .test_base import TestBase

dd = decimal.Decimal.from_float

class TestUBLDef(TestBase):

    def test_ubl_init(self):
        """ test ubl constructor """

        self.setup_uwg_integration()
        self.uwg.read_epw()
        self.uwg.set_input()
        self.uwg.init_BEM_obj()
        self.uwg.init_input_obj()

        # Get uwg_python values
        uwg_python_val = [
            self.uwg.UBL.charLength,          # characteristic length of the urban area (m)
            self.uwg.UBL.perimeter,
            self.uwg.UBL.urbArea,             # horizontal urban area (m2)
            self.uwg.UBL.orthLength,          # length of the side of the urban area orthogonal to wind dir (m)
            self.uwg.UBL.paralLength,         # length of the side of the urban area parallel to wind dir (m)
            self.uwg.UBL.ublTemp,             # urban boundary layer temperature (K)
            reduce(lambda x,y:x+y, self.uwg.UBL.ublTempdx), # urban boundary layer temperature discretization (K)
            self.uwg.UBL.dayBLHeight,         # daytime mixing height, orig = 700
            self.uwg.UBL.nightBLHeight        # Sing: 80, Bub-Cap: 50, nighttime boundary-layer height (m); orig 80
        ]

        uwg_matlab_val = self.setup_open_matlab_ref("matlab_ubl","matlab_ref_ubl_init.txt")

        # matlab ref checking
        assert len(uwg_matlab_val) == len(uwg_python_val), "matlab={}, python={}".format(len(uwg_matlab_val), len(uwg_python_val))

        for i in range(len(uwg_matlab_val)):
            #print uwg_python_val[i], uwg_matlab_val[i]
            tol = self.CALCULATE_TOLERANCE(uwg_python_val[i],15.0)
            assert uwg_python_val[i] == pytest.approx(uwg_matlab_val[i], tol), "error at index={}".format(i)

    def test_ublmodel(self):
        """ test ubl constructor """

        self.setup_uwg_integration()
        self.uwg.read_epw()
        self.uwg.set_input()

        # Test Jan 1 (winter, no vegetation coverage)
        self.uwg.Month = 1
        self.uwg.Day = 1
        self.uwg.nDay = 1

        # set_input
        self.uwg.init_BEM_obj()
        self.uwg.init_input_obj()

        # In order to avoid integration effects. Test only first time step
        # Subtract timestep to stop at 300 sec
        self.uwg.simTime.nt -= (23*12 + 11)

        # Run simulation
        self.uwg.hvac_autosize()
        self.uwg.simulate()

        # check date
        #print self.uwg.simTime
        assert self.uwg.simTime.month == 1
        assert self.uwg.simTime.day == 1
        assert self.uwg.simTime.secDay == pytest.approx(300.0,abs=1e-15)

        # Get uwg_python values
        uwg_python_val = [
            self.uwg.UBL.ublTemp,               # urban boundary layer temperature (K)
            reduce(lambda x,y:x+y,self.uwg.UBL.ublTempdx), # urban boundary layer temperature discretaization (K)
            self.uwg.UBL.dayBLHeight,         # night boundary layer height (m)
            self.uwg.UBL.nightBLHeight         # night boundary layer height (m)
        ]

        uwg_matlab_val = self.setup_open_matlab_ref("matlab_ubl","matlab_ref_ublmodel.txt")

        # matlab ref checking
        assert len(uwg_matlab_val) == len(uwg_python_val)

        for i in range(len(uwg_matlab_val)):
            #print uwg_python_val[i], uwg_matlab_val[i]
            tol = self.CALCULATE_TOLERANCE(uwg_python_val[i],15.0)
            assert uwg_python_val[i] == pytest.approx(uwg_matlab_val[i], tol), "error at index={}".format(i)


if __name__ == "__main__":
    ubl = TestUBLDef()
    ubl.test_ubl_init()
    ubl.test_ublmodel()
