#import pytest
import UWG


class TestSolarCalcs(object):

    DIR_UP_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    DIR_EPW_PATH = os.path.join(DIR_UP_PATH,"resources/epw")

    def setup_solarcalcs(self):
        """ set up solarcalcs object from initialize.uwg """

        epw_dir = self.DIR_EPW_PATH
        epw_file_name = "SGP_Singapore.486980_IWEC.epw"
        uwg_param_dir = os.path.join(self.DIR_UP_PATH,"resources")
        uwg_param_file_name = "initialize.uwg"

        UCM = self.uwg.UCM              # Urban Canopy - Building Energy Model object
        BEM = self.uwg.BEM              # Building Energy Model object
        simTime = self.uwg.simTime      # Simulation time bbject
        RSM = self.uwg.RSM              # Rural Site & Vertical Diffusion Model Object
        forc = self.uwg.forc            # Forcing object
        geoParam = self.uwg.geoParam    # Geo Param Object
        rural = self.wug.rural          # Rural road Element object

        self.uwg = UWG.UWG(epw_dir, epw_file_name, uwg_param_dir, uwg_param_file_name)
        self.solarcalcs = UWG.SolarCalcs(UCM,BEM,simTime,RSM,forc,parameter,rural)

    def test_solarangles(self):
        """ test solar angles """
        print 'test solar'
        self.setup_SolarCalcs()
        zenith, tanzern, theta0 = self.solarcalcs.solarangles(UCM.canAspect,simTime,RSM.lon,RSM.lat,RSM.GMT)
