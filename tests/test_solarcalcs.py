import pytest
import UWG
import os
import math
import pprint

class TestSolarCalcs(object):

    DIR_UP_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    DIR_EPW_PATH = os.path.join(DIR_UP_PATH,"resources/epw")
    DIR_MATLAB_PATH = os.path.join(os.path.abspath(os.path.dirname(__file__)), "matlab_ref","matlab_solarcalcs")

    def setup_solarcalcs(self):
        """ set up solarcalcs object from initialize.uwg """

        epw_dir = self.DIR_EPW_PATH
        epw_file_name = "SGP_Singapore.486980_IWEC.epw"
        uwg_param_dir = os.path.join(self.DIR_UP_PATH,"resources")
        uwg_param_file_name = "initialize.uwg"

        self.uwg = UWG.UWG(epw_dir, epw_file_name, uwg_param_dir, uwg_param_file_name)

    def test_solarangles(self):
        """ test solar angles """

        self.setup_solarcalcs()
        self.uwg.read_epw()
        self.uwg.read_input()

        solar = UWG.SolarCalcs(self.uwg.UCM, self.uwg.BEM, self.uwg.simTime, self.uwg.RSM, \
            self.uwg.forc, self.uwg.geoParam, self.uwg.rural)

        #timestep every 5 minutes (300s)
        for i in xrange(int(12*24*1.5)): #1.5 days, Jan 2nd, 12 noon
            solar.simTime.UpdateDate()
        solar.solarangles()

        # test simtime for solar after 1.5 days
        assert solar.ut == pytest.approx(12.0,abs=1e-15)
        assert solar.simTime.day == pytest.approx(2.0,abs=1e-15)
        assert solar.ad == pytest.approx(0.197963373,abs=1e-8)

        # recalculate solar angles with new time simulation, add to previous time increment
        for i in xrange(12*24*20 + 12*1 + 6): # Increment time to 22 days, 13hrs, 30min = Jan 22 at 1330
            solar.simTime.UpdateDate()

        # Run simulation
        solar.solarangles()

        # open matlab ref file
        matlab_path = os.path.join(self.DIR_MATLAB_PATH,"matlab_ref_solarangles.txt")
        if not os.path.exists(matlab_path):
            raise Exception("Failed to open {}!".format(matlab_path))
        matlab_file = open(matlab_path,'r')
        uwg_matlab_val = [float(x) for x in matlab_file.readlines()]
        matlab_file.close()

        uwg_python_val = [
        solar.simTime.month,
        solar.simTime.secDay,
        solar.simTime.day,
        solar.simTime.hourDay,
        solar.ut,
        solar.ad,
        solar.eqtime,
        solar.decsol,
        solar.zenith,
        solar.tanzen,
        solar.critOrient
        ]

        # matlab ref checking
        assert len(uwg_matlab_val) == len(uwg_python_val)
        for i in xrange(len(uwg_matlab_val)):
            #print uwg_python_val[i], uwg_matlab_val[i]
            assert uwg_python_val[i] == pytest.approx(uwg_matlab_val[i], abs=1e-15), "error at index={}".format(i)

    def test_solarcalcs(self):
        """ test solar calculation """

        self.setup_solarcalcs()
        self.uwg.read_epw()
        self.uwg.read_input()

        # We subtract 11 hours from total timestep so
        # we can stop simulation while we still have sun!
        # New time: Jan 31, 1300
        self.uwg.simTime.nt -= 12*11

        self.uwg.hvac_autosize()
        self.uwg.uwg_main()

        # check date
        assert self.uwg.simTime.month == 1
        assert self.uwg.simTime.day == 31
        assert self.uwg.simTime.secDay/3600. == pytest.approx(13.0,abs=1e-15)

        # open matlab ref file
        matlab_path = os.path.join(self.DIR_MATLAB_PATH,"matlab_ref_solarcalcs.txt")
        if not os.path.exists(matlab_path):
            raise Exception("Failed to open {}!".format(matlab_path))
        matlab_file = open(matlab_path,'r')
        uwg_matlab_val = [float(x) for x in matlab_file.readlines()]
        matlab_file.close()

        # redo matlab sim for new date
        uwg_python_val = [
        self.uwg.solar.horSol,
        self.uwg.solar.Kw_term,
        self.uwg.solar.Kr_term
        ]
        
        # matlab ref checking
        assert len(uwg_matlab_val) == len(uwg_python_val)
        for i in xrange(len(uwg_matlab_val)):
            print uwg_python_val[i], uwg_matlab_val[i]
            #assert uwg_python_val[i] == pytest.approx(uwg_matlab_val[i], abs=1e-15), "error at index={}".format(i)


if __name__ == "__main__":
    tsc = TestSolarCalcs()
    tsc.test_solarangles()
    tsc.test_solarcalcs()
