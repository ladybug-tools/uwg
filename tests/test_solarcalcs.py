import pytest
import UWG
import os
import pprint

class TestSolarCalcs(object):

    DIR_UP_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    DIR_EPW_PATH = os.path.join(DIR_UP_PATH,"resources/epw")
    DIR_MATLAB_PATH = os.path.join(os.path.abspath(os.path.dirname(__file__)), \
        "matlab_ref","matlab_solarcalcs")

    def setup_solarcalcs(self):
        """ set up solarcalcs object from initialize.uwg """

        epw_dir = self.DIR_EPW_PATH
        epw_file_name = "SGP_Singapore.486980_IWEC.epw"
        uwg_param_dir = os.path.join(self.DIR_UP_PATH,"resources")
        uwg_param_file_name = "initialize.uwg"

        self.uwg = UWG.UWG(epw_dir, epw_file_name, uwg_param_dir, uwg_param_file_name)
        self.uwg.read_epw()
        self.uwg.read_input()

        self.UCM = self.uwg.UCM              # Urban Canopy - Building Energy Model object
        self.BEM = self.uwg.BEM              # Building Energy Model object
        #time step is 5 min
        self.simTime = self.uwg.simTime      # Simulation time bbject
        self.RSM = self.uwg.RSM              # Rural Site & Vertical Diffusion Model Object
        self.forc = self.uwg.forc            # Forcing object
        self.geoParam = self.uwg.geoParam    # Geo Param Object
        self.rural = self.uwg.rural          # Rural road Element object

        self.solarcalcs = UWG.SolarCalcs(self.UCM,self.BEM,self.simTime,self.RSM,self.forc,self.geoParam,self.rural)

    def test_solarangles(self):
        """ test solar angles """

        self.setup_solarcalcs()
        #timestep every 5 minutes (300s)
        for i in xrange(int(12*24*1.5)): #1.5 days, Jan 2nd, 12 noon
            self.solarcalcs.simTime.UpdateDate()
        zenith, tanzen, theta0 = self.solarcalcs.solarangles()

        # test simtime for solar after 1.5 days
        assert self.solarcalcs.ut == pytest.approx(12.0,abs=1e-15)
        assert self.solarcalcs.simTime.day == pytest.approx(2.0,abs=1e-15)
        assert self.solarcalcs.ad == pytest.approx(0.197963373,abs=1e-8)

        # matlab ref checking
        for i in xrange(12*24*20 + 12*1 + 6): # Incremen time to 22 days, 13hrs, 30min = Jan 22 at 1330
            self.solarcalcs.simTime.UpdateDate()

        # open matlab ref file
        matlab_path = os.path.join(self.DIR_MATLAB_PATH,"matlab_ref_solarangles.txt")
        if not os.path.exists(matlab_path):
            raise Exception("Failed to open {}!".format(matlab_path))
        matlab_file = open(matlab_path,'r')
        uwg_matlab_val = [float(x) for x in matlab_file.readlines()]
        matlab_file.close()

        uwg_python_val = [
        self.solarcalcs.simTime.month,
        self.solarcalcs.simTime.secDay,
        self.solarcalcs.simTime.day,
        self.solarcalcs.simTime.hourDay,
        self.solarcalcs.ut,
        self.solarcalcs.ad,
        self.solarcalcs.eqtime,
        self.solarcalcs.decsol,
        zenith,
        tanzen,
        theta0
        ]

        #   print len(uwg_matlab_val), len(uwg_python_val)
        assert len(uwg_matlab_val) == len(uwg_python_val)

        print self.solarcalcs.ut, 13.500000000000000
        #for i in xrange(len(uwg_matlab_val)):
            #print uwg_python_val[i], uwg_matlab_val[i]

            #assert uwg_python_val[i] == pytest.approx(uwg_matlab_val[i], abs=1e-15), "error at index={}".format(i)





if __name__ == "__main__":
    tsc = TestSolarCalcs()
    tsc.test_solarangles()
