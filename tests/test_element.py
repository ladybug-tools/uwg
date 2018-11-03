try:
    range = xrange
except NameError:
    pass

import pytest
import uwg
import os
import math
import pprint
from .test_base import TestBase


class TestElement(TestBase):

    def test_SurfFlux_with_waterStorage_start(self):
        """ Edge case: test element SurfFlux against matlab reference
        when waterStorage > 0.0.
        This has to be hardcoded b/c doesn't get used otherwise

        """
        self.setup_uwg_integration()
        self.uwg.read_epw()
        self.uwg.set_input()
        self.uwg.init_BEM_obj()
        self.uwg.init_input_obj()
        # We subtract 30 days and 11 hours
        # New time: Jan 1, 1:00
        self.uwg.simTime.nt -= (30*24*12 + 23*12 + 11)

        # turn rural road waterStorage to 1.
        self.uwg.rural.waterStorage = 0.005 # .5cm thick film (from wgmax constant)

        # run simulation
        self.uwg.hvac_autosize()
        self.uwg.simulate()

        # check date
        #print self.uwg.simTime
        assert self.uwg.simTime.month == 1
        assert self.uwg.simTime.day == 1
        assert self.uwg.simTime.secDay == pytest.approx(300.,abs=1e-15)

        # Check waterStorage
        assert 0.005 == pytest.approx(self.uwg.rural.waterStorage, 1e-15)

    def test_SurfFlux_with_waterStorage_middle(self):
        """ Edge case: test element SurfFlux against matlab reference
        when waterStorage > 0.0.
        This has to be hardcoded b/c doesn't get used otherwise

        """
        self.setup_uwg_integration()
        self.uwg.read_epw()
        self.uwg.set_input()
        self.uwg.init_BEM_obj()
        self.uwg.init_input_obj()
        # We subtract 30 days and 11 hours
        # New time: Jan 1, 1:00
        self.uwg.simTime.nt -= (30*24*12 + 11*12)

        # turn rural road waterStorage to 1.
        self.uwg.rural.waterStorage = 0.005 # .5cm thick film (from wgmax constant)

        # run simulation
        self.uwg.hvac_autosize()
        self.uwg.simulate()

        # check date
        #print self.uwg.simTime
        assert self.uwg.simTime.month == 1
        assert self.uwg.simTime.day == 1
        assert self.uwg.simTime.secDay/3600. == pytest.approx(13.0,abs=1e-15)

        # Check waterStorage
        #print self.uwg.rural.waterStorage
        assert 0.004643553960210 == pytest.approx(self.uwg.rural.waterStorage, 1e-15)

    def test_SurfFlux_unit(self):
        """ test element SurfFlux against matlab references at the
        start of timestep.

        """

        self.setup_uwg_integration()

        self.uwg.read_epw()
        self.uwg.set_input()
        self.uwg.init_BEM_obj()
        self.uwg.init_input_obj()

        # We subtract 23 hours and 55 minutes so we can test
        # initial timestep (1, 1, 300). New time: Jan 1, 5min
        self.uwg.simTime.nt -= (23*12 + 24*12*30 + 11)

        # run simulation
        self.uwg.hvac_autosize()
        self.uwg.simulate()

        # check date is Jan 1st, 300s
        #print self.uwg.simTime
        assert self.uwg.simTime.month == 1
        assert self.uwg.simTime.day == 1
        assert self.uwg.simTime.secDay == pytest.approx(300.0,abs=1e-15)

        uwg_python_val = [
        self.uwg.rural.aeroCond,        # Convection coef (refL uwg, eq.12)
        self.uwg.rural.waterStorage,    # thickness of water film (m) (only for horizontal surfaces)
        self.uwg.rural.solAbs,          # solar radiation absorbed (W m-2)
        self.uwg.rural.lat,             # surface latent heat flux (W m-2)
        self.uwg.rural.sens,            # surface sensible heat flux (W m-2)
        self.uwg.rural.flux,            # external surface heat flux (W m-2)
        self.uwg.rural.T_ext,           # external surface temperature (K)
        self.uwg.rural.T_int            # interior surface temperature (K)
        ]

        # Matlab Checking for rural road
        uwg_matlab_val = self.setup_open_matlab_ref("matlab_element","matlab_ref_element_surfflux_winter.txt")

        # matlab ref checking
        assert len(uwg_matlab_val) == len(uwg_python_val)

        for i in range(len(uwg_matlab_val)):
            #print uwg_python_val[i], uwg_matlab_val[i]
            tol = self.CALCULATE_TOLERANCE(uwg_matlab_val[i],15.0)
            assert uwg_python_val[i] == pytest.approx(uwg_matlab_val[i], abs=tol), "error at index={}".format(i)


    def test_SurfFlux_integration(self):
        """ test element SurfFlux against matlab references at the
        end of timestep. Integration test as it requires other uwg
        classes to be working before functioning correctly.

        """
        self.setup_uwg_integration()

        self.uwg.read_epw()
        self.uwg.set_input()
        # Change time and vegCoverage parameters so we can get
        # effect of vegetation on surface heat flux
        self.uwg.vegStart = 2   # February
        self.uwg.nDay = 31 + 15 # February, 15

        self.uwg.init_BEM_obj()
        self.uwg.init_input_obj()

        # We subtract 11 hours from total timestep so still have sun. New time: 1300
        self.uwg.simTime.nt -= 12*11

        self.uwg.hvac_autosize()
        self.uwg.simulate()

        # check date is February 15th, 13:00
        assert self.uwg.simTime.month == 2
        assert self.uwg.simTime.day == 15
        assert self.uwg.simTime.secDay/3600. == pytest.approx(13.0,abs=1e-15)

        uwg_python_val = [
        self.uwg.rural.aeroCond,        # Convection coef (refL uwg, eq.12)
        self.uwg.rural.waterStorage,    # thickness of water film (m) (only for horizontal surfaces)
        self.uwg.rural.solAbs,          # solar radiation absorbed (W m-2)
        self.uwg.rural.lat,             # surface latent heat flux (W m-2)
        self.uwg.rural.sens,            # surface sensible heat flux (W m-2)
        self.uwg.rural.flux,             # external surface heat flux (W m-2)
        self.uwg.rural.T_ext,           # external surface temperature (K)
        self.uwg.rural.T_int            # interior surface temperature (K)
        ]

        # Matlab Checking for rural road
        uwg_matlab_val = self.setup_open_matlab_ref("matlab_element","matlab_ref_element_surfflux.txt")

        # matlab ref checking
        assert len(uwg_matlab_val) == len(uwg_python_val)

        for i in range(len(uwg_matlab_val)):
            #print uwg_python_val[i], uwg_matlab_val[i]
            tol = self.CALCULATE_TOLERANCE(uwg_matlab_val[i],15.0)
            assert uwg_python_val[i] == pytest.approx(uwg_matlab_val[i], abs=tol), "error at index={}".format(i)


if __name__ == "__main__":
    te = TestElement()
    #te.test_SurfFlux_with_waterStorage_start()
    #te.test_SurfFlux_with_waterStorage_middle()
    #te.test_SurfFlux_unit()
    te.test_SurfFlux_integration()
