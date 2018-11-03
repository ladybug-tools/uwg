try:
    range = xrange
except NameError:
    pass

import pytest
import uwg
import os
import math
import pprint
import decimal

import logging

from .test_base import TestBase

dd = decimal.Decimal.from_float
pp = pprint.pprint


class TestOutput(TestBase):
    """
    Compare the matlab generated epw file against the python generated epw files using the same .m/.uwg inputs.
    """
    def compare_epw(self,matlab_fname,precision = 10.0):

        # shorten some variable names
        ti = self.uwg.simTime.timeInitial
        tf = self.uwg.simTime.timeFinal

        # Make weather files for testing
        matlab_path_name = os.path.join(self.DIR_CURR,"..","tests","matlab_ref","matlab_output",
            matlab_fname)

        pywtr = uwg.weather.Weather(self.uwg.newPathName, ti, tf)
        matwtr = uwg.weather.Weather(matlab_path_name, ti, tf)
        #epwwtr = self.uwg.weather # for reference only

        assert len(pywtr.staTemp) == pytest.approx(len(matwtr.staTemp), abs=1e-15)

        for i in range(0,len(pywtr.staTemp)):

            # Check dry bulb [K]
            tol = self.CALCULATE_TOLERANCE(pywtr.staTemp[i],precision)
            assert pywtr.staTemp[i] == pytest.approx(matwtr.staTemp[i], abs=tol), "error at index={}".format(i)

            # Check dew point [K]
            tol = self.CALCULATE_TOLERANCE(pywtr.staTdp[i],precision)
            assert pywtr.staTdp[i] == pytest.approx(matwtr.staTdp[i], abs=tol), "error at index={}".format(i)

            # Check relative humidity [%]
            tol = self.CALCULATE_TOLERANCE(pywtr.staRhum[i],precision)
            assert pywtr.staRhum[i] == pytest.approx(matwtr.staRhum[i], abs=tol), "error at index={}".format(i)

            # Check wind speed [m/s]
            tol = self.CALCULATE_TOLERANCE(pywtr.staUmod[i],precision)
            assert pywtr.staUmod[i] == pytest.approx(matwtr.staUmod[i], abs=tol), "error at index={}".format(i)

    def test_uwg_output_heatdemand_1_1_0000(self):
        """
        Initial conditions:
            - night time
            - before vegstart
            - sensHeatDemand
        """

        self.setup_uwg_integration(epw_file="CAN_ON_Toronto.716240_CWEC.epw")
        self.uwg.read_epw()
        self.uwg.set_input()

        # Test all year
        self.uwg.Month = 1
        self.uwg.Day = 1
        self.uwg.nDay = 365

        self.uwg.init_BEM_obj()
        self.uwg.init_input_obj()
        self.uwg.hvac_autosize()
        self.uwg.simulate()
        self.uwg.write_epw()

        self.compare_epw("CAN_ON_Toronto.716240_CWEC_heatdemand_UWG_Matlab.epw")#,precision=3.0)

    def test_program_input(self):
        self.setup_uwg_integration("SGP_Singapore.486980_IWEC.epw", None)

        # Check some random variables
        assert self.uwg.vegCover == None
        assert self.uwg.treeCoverage == None
        assert self.uwg.vegStart == None

        # Assign manually
        self.set_input_manually()

        # main
        self.uwg.read_epw()
        self.uwg.set_input()
        self.uwg.init_BEM_obj()
        self.uwg.init_input_obj()
        self.uwg.hvac_autosize()
        self.uwg.simulate()
        self.uwg.write_epw()

        # Check some of the inputs

        # Check building parameters
        assert self.uwg.BEM[0].building.coolCap == pytest.approx((3525.66904*1000.0)/46320.0, abs=1e-3)
        assert self.uwg.BEM[0].building.heatCap == pytest.approx((2875.97378*1000.0)/46320.0, abs=1e-3)
        assert self.uwg.BEM[1].building.coolCap == pytest.approx((252.20895*1000.0)/3135., abs=1e-2)
        assert self.uwg.BEM[1].building.heatCap == pytest.approx((132.396*1000.0)/3135., abs=1e-2)

        # Check that final day of timestep is at correct dayType
        assert self.uwg.dayType == pytest.approx(1., abs=1e-3)
        assert self.uwg.SchTraffic[self.uwg.dayType-1][self.uwg.simTime.hourDay] == pytest.approx(0.2, abs=1e-6)

        self.compare_epw("SGP_Singapore.486980_IWEC_UWG_Matlab.epw")


    def test_program_hybrid_input(self):
        """ Testing inputting with api and with .uwg file"""

        self.setup_uwg_integration("SGP_Singapore.486980_IWEC.epw", "initialize_singapore.uwg")

        # Check some random variables
        assert self.uwg.vegCover == None
        assert self.uwg.treeCoverage == None
        assert self.uwg.vegStart == None

        # Assign manually
        self.set_input_manually()

        # Empty some parameters so that it can be defined by uwg file
        # Urban characteristics
        self.uwg.bldHeight = None
        self.uwg.h_mix = None
        self.uwg.bldDensity = None
        self.uwg.verToHor = None
        self.uwg.charLength = None
        self.uwg.alb_road = None
        self.uwg.d_road = None
        self.uwg.sensAnth = None
        self.uwg.latAnth = None

        # main
        self.uwg.read_epw()
        self.uwg.set_input()
        self.uwg.init_BEM_obj()
        self.uwg.init_input_obj()
        self.uwg.hvac_autosize()
        self.uwg.simulate()
        self.uwg.write_epw()

        # Check some of the inputs

        # Check building parameters
        assert self.uwg.BEM[0].building.coolCap == pytest.approx((3525.66904*1000.0)/46320.0, abs=1e-3)
        assert self.uwg.BEM[0].building.heatCap == pytest.approx((2875.97378*1000.0)/46320.0, abs=1e-3)
        assert self.uwg.BEM[1].building.coolCap == pytest.approx((252.20895*1000.0)/3135., abs=1e-2)
        assert self.uwg.BEM[1].building.heatCap == pytest.approx((132.396*1000.0)/3135., abs=1e-2)

        # Check that final day of timestep is at correct dayType
        assert self.uwg.dayType == pytest.approx(1., abs=1e-3)
        assert self.uwg.SchTraffic[self.uwg.dayType-1][self.uwg.simTime.hourDay] == pytest.approx(0.2, abs=1e-6)

        self.compare_epw("SGP_Singapore.486980_IWEC_UWG_Matlab.epw")

    def test_uwg_output_beijing(self):
        """
        Initial conditions:
            - day time
            - after vegstart
            - sensHeatDemand
        """

        self.setup_uwg_integration(
            epw_file="CHN_Beijing.Beijing.545110_IWEC.epw",
            uwg_param_file="initialize_beijing.uwg"
            )

        self.uwg.read_epw()
        self.uwg.set_input()

        # Test all year
        self.uwg.Month = 1
        self.uwg.Day = 1
        self.uwg.nDay = 365


        # main
        self.uwg.init_BEM_obj()
        self.uwg.init_input_obj()
        self.uwg.hvac_autosize()
        self.uwg.simulate()
        self.uwg.write_epw()

        self.compare_epw("CHN_Beijing.Beijing.545110_IWEC_UWG_Matlab.epw")

    def test_uwg_output_cooldemand_6_1_0000(self):
        """
        Initial conditions:
            - day time
            - after vegstart
            - sensHeatDemand
        """
        # set up the logger

        self.setup_uwg_integration(epw_file="CAN_ON_Toronto.716240_CWEC.epw")#,log_file_name="test_uwg_output_cooldemand_6_1_1300.log",log_level=logging.INFO)

        #self.uwg.logger.critical("Cool demand output")

        self.uwg.read_epw()
        self.uwg.set_input()

        # Test 30 days in summer
        self.uwg.Month = 6
        self.uwg.Day = 1
        self.uwg.nDay = 30

        # main
        self.uwg.init_BEM_obj()
        self.uwg.init_input_obj()
        self.uwg.hvac_autosize()
        self.uwg.simulate()
        self.uwg.write_epw()

        self.compare_epw("CAN_ON_Toronto.716240_CWEC_cooldemand_UWG_Matlab.epw")

    def test_uwg_output_cooldemand_1_1_0000(self):
        """
        Initial conditions:
            - night time
            - before vegstart
            - sensCoolDemand
        """
        self.setup_uwg_integration()
        self.uwg.read_epw()
        self.uwg.set_input()

        # Test all year
        self.uwg.Month = 1
        self.uwg.Day = 1
        self.uwg.nDay = 365

        # main
        self.uwg.init_BEM_obj()
        self.uwg.init_input_obj()
        self.uwg.hvac_autosize()
        self.uwg.simulate()
        self.uwg.write_epw()

        self.compare_epw("SGP_Singapore.486980_IWEC_UWG_Matlab.epw")


if __name__ == "__main__":
    test = TestOutput()
    #test.test_uwg_output_heatdemand_1_1_0000()
    #test.test_uwg_output_beijing()
    #test.test_uwg_output_cooldemand_6_1_0000()
    #test.test_uwg_output_cooldemand_1_1_0000()
