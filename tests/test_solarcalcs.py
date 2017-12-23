import pytest
import UWG
import os

class TestSolarCalcs(object):

    DIR_UP_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    DIR_EPW_PATH = os.path.join(DIR_UP_PATH,"resources/epw")

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
        for i in xrange(12*36): #1.5 days
            self.simTime.UpdateDate()

        zenith, tanzern, theta0 = self.solarcalcs.solarangles()

        # test simtime for solar
        assert self.solarcalcs.ut == pytest.approx(12.0,abs=1e-15)
        assert self.solarcalcs.date == pytest.approx(1.0,abs=1e-15)
        assert self.solarcalcs.ad == pytest.approx(0.197963373,abs=1e-8)
        #matlab checking
        #self.eqtime
        #self.decsol




if __name__ == "__main__":
    tsc = TestSolarCalcs()
    tsc.test_solarangles()
