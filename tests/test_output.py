import pytest
import UWG
import os
import math

import pprint
import decimal

dd = lambda x: decimal.Decimal.from_float(x)
pp = lambda x: pprint.pprint(x)


class TestOutput(object):

    DIR_UP_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    DIR_EPW_PATH = os.path.join(DIR_UP_PATH,"resources/epw")
    DIR_MATLAB_PATH = os.path.join(os.path.abspath(os.path.dirname(__file__)), "matlab_ref","matlab_ubl")
    CALCULATE_TOLERANCE = lambda s,x,p: 1*10**-(p - 1 - int(math.log10(x))) if abs(float(x)) > 1e-15 else 1e-15

    def setup_uwg_integration(self):
        """ set up uwg object from initialize.uwg """

        epw_dir = self.DIR_EPW_PATH
        epw_file_name = "SGP_Singapore.486980_IWEC.epw"
        uwg_param_dir = os.path.join(self.DIR_UP_PATH,"resources")
        uwg_param_file_name = "initialize.uwg"

        self.uwg = UWG.UWG(epw_dir, epw_file_name, uwg_param_dir, uwg_param_file_name)

    def setup_open_matlab_ref(self,matlab_ref_file_path):
        """ open the matlab reference file """

        matlab_path = os.path.join(self.DIR_MATLAB_PATH,matlab_ref_file_path)
        if not os.path.exists(matlab_path):
            raise Exception("Failed to open {}!".format(matlab_path))
        matlab_file = open(matlab_path,'r')
        uwg_matlab_val_ = [float(x) for x in matlab_file.readlines()]
        matlab_file.close()
        return uwg_matlab_val_


    def test_uwg_output(self):

        self.setup_uwg_integration()
        self.uwg.read_epw()
        self.uwg.read_input()

        # Test all year
        self.uwg.Month = 6
        self.uwg.Day = 21
        self.uwg.nDay = 1

        self.uwg.set_input()
        self.uwg.hvac_autosize()
        self.uwg.uwg_main()
        self.uwg.write_epw()

        # shorten some variable names
        ti = self.uwg.simTime.timeInitial
        tf = self.uwg.simTime.timeFinal
        matlab_fname = "SGP_Singapore.486980_IWEC_UWG_Matlab.epw"

        # Get matlab data
        matlab_path_name = os.path.join(self.DIR_UP_PATH,"tests","matlab_ref","matlab_output",matlab_fname)

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
        for i in xrange(len(pywtr.staTemp)):
            print 'hr:', i
            # dry bulb temperature  [?C]
            print 'Tdb'
            print pywtr.staTemp[i]-273.15
            print matwtr.staTemp[i]-273.15
            print epwwtr.staTemp[i]-273.15,'\n--'
            # dew point temperature [?C]
            print 'Tdp'
            print pywtr.staTdp[i]
            print matwtr.staTdp[i]
            print epwwtr.staTdp[i],'\n--'
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


if __name__ == "__main__":
    test = TestOutput()
    test.test_uwg_output()
