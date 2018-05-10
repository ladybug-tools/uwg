import pytest
import UWG
import os
import math
import pprint
import decimal

import logging

from test_base import TestBase

dd = decimal.Decimal.from_float
pp = pprint.pprint


class TestOutput(TestBase):
    """
    Three checks
    1 check at time < parameter.nightSetStart && time > parameter.nightSetEnd || check is day || check dir + dif > 0
    2 check month > vegStart and month < vegEnd
    3 check if self.sensCoolDemand > 0. and UCM.canTemp > 288. #288~15 || check if self.sensHeatDemand > 0. and UCM.canTemp < 288.

    """
    """
    def test_uwg_output_heatdemand_1_1_0000(self):
        #
        #Initial conditions:
        #    - night time
        #    - before vegstart
        #    - sensHeatDemand
        #

        self.setup_uwg_integration(epw_file="CAN_ON_Toronto.716240_CWEC.epw")
        self.uwg.read_epw()
        self.uwg.read_input()

        # Test all year
        self.uwg.Month = 1
        self.uwg.Day = 1
        self.uwg.nDay = 1

        self.uwg.set_input()

        self.uwg.hvac_autosize()

        #for i in xrange(12*13):
        #    self.uwg.simTime.UpdateDate()

        #print 'mth', self.uwg.simTime.month
        #print 'day', self.uwg.simTime.day
        #print 'hr', self.uwg.simTime.secDay/3600.
        #print '------------------'

        self.uwg.simulate()
        self.uwg.write_epw()

        # shorten some variable names
        ti = self.uwg.simTime.timeInitial
        tf = self.uwg.simTime.timeFinal
        matlab_fname = "CAN_ON_Toronto.716240_CWEC_heatdemand_UWG_Matlab.epw"

        # Get matlab data
        matlab_path_name = os.path.join(self.DIR_CURR,"..","tests","matlab_ref","matlab_output",matlab_fname)

        # Get Matlab EPW file
        if not os.path.exists(matlab_path_name):
            raise Exception("Param file: '{}' does not exist.".format(matlab_path_name))

        # Open .epw file and feed csv data to initializeDataFile
        try:
            new_epw = UWG.utilities.read_csv(matlab_path_name)
        except Exception as e:
            raise Exception("Failed to read .uwg file! {}".format(e.message))

        matlab_weather = UWG.weather.Weather(matlab_path_name,ti,tf)
        print self.uwg.newPathName
        # Make weather files for testing
        pywtr = UWG.weather.Weather(self.uwg.newPathName,ti,tf)
        matwtr = UWG.weather.Weather(matlab_path_name, ti, tf)
        epwwtr = self.uwg.weather
        #print '\n'

        assert len(pywtr.staTemp) == pytest.approx(len(matwtr.staTemp), abs=1e-15)

        # compare per hours
        lstlen = len(pywtr.staTemp)
        for i in xrange(lstlen-1,lstlen-2,-1):
            print 'hr:', i
            # dry bulb temperature  [?C]
            print 'Tdb'
            print 'p-', pywtr.staTemp[i]-273.15
            print 'm-', matwtr.staTemp[i]-273.15
            print 'e-', epwwtr.staTemp[i]-273.15,'\n--'
            # dew point temperature [?C]
            #print 'Tdp'
            #print 'p-', pywtr.staTdp[i]
            #print 'm-', matwtr.staTdp[i]
            #print 'e-', epwwtr.staTdp[i],'\n--'
            # relative humidity     [%]
            print 'RH'
            print 'p-', pywtr.staRhum[i]
            print 'm-', matwtr.staRhum[i]
            print 'e-', epwwtr.staRhum[i],'\n--'
            # wind speed [m/s]
            #print pywtr.staUmod[i]
            #print matwtr.staUmod[i]
            #print epwwtr.staUmod[i],'\n--'
            print '----#'

    def test_uwg_output_work(self):
        #"
        #Initial conditions:
        #    - day time
        #    - after vegstart
        #    - sensHeatDemand
        #"

        # Set up logger
        #log_path = os.path.join(self.DIR_UP_PATH,"tests","log","test_uwg_output_cooldemand_6_1_1300.log")
        #logging.basicConfig(filename=log_path, filemode="w",level=logging.DEBUG)

        self.setup_uwg_integration(
            epw_file="CHN_Beijing.Beijing.545110_IWEC.epw",#USA_MA_Boston-Logan.Intl.AP.725090_TMY3.epw",
            uwg_param_file="initialize_beijing.uwg"
            )#"SGP_Singapore.486980_IWEC.epw") #"CAN_ON_Toronto.716240_CWEC.epw")

        self.uwg.read_epw()
        self.uwg.read_input()

        # Test all year
        self.uwg.Month = 8
        self.uwg.Day = 7
        self.uwg.nDay = 14

        self.uwg.set_input()
        self.uwg.hvac_autosize()

        # 0 - 23
        self.uwg.simulate(0,0)
        self.uwg.write_epw()
    """
    def test_uwg_output_cooldemand_6_1_1300(self):
        """
        Initial conditions:
            - day time
            - after vegstart
            - sensHeatDemand
        """
        # set up the logger
        self.set_log_file("test_uwg_output_cooldemand_6_1_1300.log",log_level=logging.INFO)

        self.setup_uwg_integration(epw_file="CAN_ON_Toronto.716240_CWEC.epw")

        self.uwg.logger.critical("Cool demand output")

        self.uwg.read_epw()
        self.uwg.read_input()

        # Test all year
        self.uwg.Month = 1
        self.uwg.Day = 1
        self.uwg.nDay = 1#365

        self.uwg.set_input()
        self.uwg.hvac_autosize()

        # 0 - 23
        self.uwg.simulate(0,0)
        self.uwg.write_epw()

        # shorten some variable names
        ti = self.uwg.simTime.timeInitial
        tf = self.uwg.simTime.timeFinal

        """
        #TODO: need to generate matlab cooldemand
        matlab_fname = "CAN_ON_Toronto.716240_CWEC_cooldemand_UWG_Matlab.epw"
        # Get matlab data
        matlab_path_name = os.path.join(self.DIR_CURR,"..","tests","matlab_ref","matlab_output",matlab_fname)
        # Get Matlab EPW file
        if not os.path.exists(matlab_path_name):
            raise Exception("Param file: '{}' does not exist.".format(matlab_path_name))
        # Open .uwg file and feed csv data to initializeDataFile
        try:
            new_epw = UWG.utilities.read_csv(matlab_path_name)
        except Exception as e:
            raise Exception("Failed to read .uwg file! {}".format(e.message))
        matlab_weather = UWG.weather.Weather(matlab_path_name,ti,tf)

        # Make weather files for testing
        pywtr = UWG.weather.Weather(self.uwg.newPathName,ti,tf)
        matwtr = UWG.weather.Weather(matlab_path_name, ti, tf)

        epwwtr = self.uwg.weather
        #print '\n'

        assert len(pywtr.staTemp) == pytest.approx(len(matwtr.staTemp), abs=1e-15)

        # compare per hours
        lstlen = len(pywtr.staTemp)
        for i in xrange(0,lstlen,1):
            #print 'hr:', i
            # dry bulb temperature  [?C]
            #print 'Tdb'
            print (pywtr.staTemp[i]-273.15)-(matwtr.staTemp[i]-273.15)
            #print epwwtr.staTemp[i]-273.15,'\n--'
            # dew point temperature [?C]
            #print 'Tdp'
            #print pywtr.staTdp[i]
            #print matwtr.staTdp[i]
            #print epwwtr.staTdp[i],'\n--'
            # relative humidity     [%]
            #print 'RH'
            #print pywtr.staRhum[i]-matwtr.staRhum[i]
            #print epwwtr.staRhum[i],'\n--'
            # wind speed [m/s]
            #print pywtr.staUmod[i]
            #print matwtr.staUmod[i]
            #print epwwtr.staUmod[i],'\n--'
            #print '----#'
        """
    """
    def test_uwg_output_cooldemand_1_1_0000(self):
        "
        Initial conditions:
            - night time
            - before vegstart
            - sensCoolDemand
        "
        self.setup_uwg_integration()
        self.uwg.read_epw()
        self.uwg.read_input()

        # Test all year
        self.uwg.Month = 1
        self.uwg.Day = 1
        self.uwg.nDay = 1

        self.uwg.set_input()
        self.uwg.hvac_autosize()
        self.uwg.simulate()
        self.uwg.write_epw()

        # shorten some variable names
        ti = self.uwg.simTime.timeInitial
        tf = self.uwg.simTime.timeFinal
        matlab_fname = "SGP_Singapore.486980_IWEC_UWG_Matlab.epw"

        # Get matlab data
        matlab_path_name = os.path.join(self.DIR_CURR,"..","tests","matlab_ref","matlab_output",matlab_fname)

        # Get Matlab EPW file
        if not os.path.exists(matlab_path_name):
            raise Exception("Param file: '{}' does not exist.".format(matlab_path_name))

        # Open .uwg file and feed csv data to initializeDataFile
        try:
            new_epw = UWG.utilities.read_csv(matlab_path_name)
        except Exception as e:
            raise Exception("Failed to read .uwg file! {}".format(e.message))

        matlab_weather = UWG.weather.Weather(matlab_path_name,ti,tf)

        # Make weather files for testing
        pywtr = UWG.weather.Weather(self.uwg.newPathName,ti,tf)
        matwtr = UWG.weather.Weather(matlab_path_name, ti, tf)
        epwwtr = self.uwg.weather
        #print '\n'

        assert len(pywtr.staTemp) == pytest.approx(len(matwtr.staTemp), abs=1e-15)

        # compare per hours
        lstlen = len(pywtr.staTemp)
        for i in xrange(lstlen-1,lstlen-2,-1):
            print 'hr:', i
            # dry bulb temperature  [?C]
            print 'Tdb'
            print pywtr.staTemp[i]-273.15
            print matwtr.staTemp[i]-273.15
            print epwwtr.staTemp[i]-273.15,'\n--'
            # dew point temperature [?C]
            #print 'Tdp'
            #print pywtr.staTdp[i]
            #print matwtr.staTdp[i]
            #print epwwtr.staTdp[i],'\n--'
            # relative humidity     [%]
            print 'RH'
            print pywtr.staRhum[i]
            print matwtr.staRhum[i]
            print epwwtr.staRhum[i],'\n--'
            # wind speed [m/s]
            #print pywtr.staUmod[i]
            #print matwtr.staUmod[i]
            #print epwwtr.staUmod[i],'\n--'
            print '----#'

    """

if __name__ == "__main__":
    test = TestOutput()
    #print "coolemand 1 1 0000\n"
    #test.test_uwg_output_cooldemand_1_1_0000()
    #print "cooldemand 1 1 13000\n"

    #$ warble bug
    test.test_uwg_output_cooldemand_6_1_1300()

    #print "heatdemand 1 1 0000\n"
    #test.test_uwg_output_heatdemand_1_1_0000()

    #
    #test.test_uwg_output_work()
