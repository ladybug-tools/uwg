import os
import pytest
import UWG

class TestRSMDef(object):
    """Test for RSMDef.py - Rural Site & Vertical Diffusion Model (RSM & VDM)

    Naming: Test prefixed test classes (without an __init__ method)
    for test autodetection by pytest
    """
    DIR_UP_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    DIR_EPW_PATH = os.path.join(DIR_UP_PATH,"resources/epw")

    def setup_init(self):
        """ Initialize rsm class """

        # Make UWG instance
        epw_dir = self.DIR_EPW_PATH
        epw_file_name = "SGP_Singapore.486980_IWEC.epw"
        uwg_param_dir = os.path.join(self.DIR_UP_PATH,"resources")
        uwg_param_file_name = "initialize.uwg"
        self.uwg = UWG.UWG(epw_dir, epw_file_name, uwg_param_dir, uwg_param_file_name)
        self.uwg.read_epw()
        self.uwg.read_input()

        # Make RSM instance for rural site parameters
        lat = 1.37
        lon = 103.98
        GMT = 8.0
        rural_height = 0.1            # average obstacle height from initialize.uwg
        urban_height = 1.0            # average building height / 10.
        T_init = 297.85               # initial dry bulb
        P_init = 100900.0             # initial pressure

        geo_param = self.uwg.geoParam # geographic parameters
        self.rural_rsm = UWG.RSMDef(lat,lon,GMT,rural_height,T_init,P_init,geo_param)
        self.urban_rsm = UWG.RSMDef(lat,lon,GMT,urban_height,T_init,P_init,geo_param)

    def test_rsm(self):
        self.setup_init()
        print self.rural_rsm

if __name__ == "__main__":
    test_rsm = TestRSMDef()
    test_rsm.test_rsm()
