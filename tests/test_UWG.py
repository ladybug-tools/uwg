import os
import pytest
import UWG



class TestUWG(object):
    """Test for UWG.py
    Naming: Test prefixed test classes (without an __init__ method)
    for test autodetection by pytest
    """
    DIR_UP_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    DIR_EPW_PATH = os.path.join(DIR_UP_PATH,"resources/epw")

    def setup_init_uwg(self):
        epw_dir = self.DIR_EPW_PATH
        epw_file_name = "SGP_Singapore.486980_IWEC.epw"
        uwg_param_dir = os.path.join(self.DIR_UP_PATH,"resources")
        uwg_param_file_name = "initialize.uwg"

        self.uwg = UWG.UWG(epw_dir, epw_file_name, uwg_param_dir, uwg_param_file_name)

    def test_read_epw(self):
        self.setup_init_uwg()
        self.uwg.read_epw()

        # test header
        assert self.uwg.header[0][0] == "LOCATION"
        assert self.uwg.header[0][1] == "SINGAPORE"
        assert self.uwg.lat == pytest.approx(1.37, abs=1e-3)
        assert self.uwg.lon == pytest.approx(103.98, abs=1e-3)
        assert self.uwg.GMT == pytest.approx(8, abs=1e-3)
        # test soil data
        assert self.uwg.nSoil == pytest.approx(3, abs=1e-2)
        # test soil depths
        assert self.uwg.depth[0][0] == pytest.approx(0.5, abs=1e-3)
        assert self.uwg.depth[1][0] == pytest.approx(2., abs=1e-3)
        assert self.uwg.depth[2][0] == pytest.approx(4., abs=1e-3)
        # test soil temps over 12 months
        assert self.uwg.Tsoil[0][0] == pytest.approx(27.55+273.15, abs=1e-3)
        assert self.uwg.Tsoil[1][2] == pytest.approx(28.01+273.15, abs=1e-3)
        assert self.uwg.Tsoil[2][11] == pytest.approx(27.07+273.15, abs=1e-3)
        # test time step in weather file
        assert self.uwg.epwinput[0][0] == "1989"
        assert float(self.uwg.epwinput[3][6]) == pytest.approx(24.1,abs=1e-3)

    def test_read_input(self):
        self.setup_init_uwg()
        self.uwg.read_epw()
        self.uwg.read_input()

        #test uwg param dictionary first and last
        assert self.uwg.init_param_dict.has_key('bldHeight') == True
        assert self.uwg.init_param_dict.has_key('h_obs') == True
        #test values
        assert self.uwg.init_param_dict['bldHeight'] == pytest.approx(10., abs=1e-6)
        assert self.uwg.init_param_dict['vegEnd'] == pytest.approx(10, abs=1e-6)
        assert self.uwg.init_param_dict['albRoof'] == pytest.approx(0.5, abs=1e-6)
        assert self.uwg.init_param_dict['h_ubl1'] == pytest.approx(1000., abs=1e-6)
        assert self.uwg.init_param_dict['h_ref'] == pytest.approx(150., abs=1e-6)

        # test SchTraffic schedule
        assert self.uwg.init_param_dict['SchTraffic'][0][0] == pytest.approx(0.2, abs=1e-6) # first
        assert self.uwg.init_param_dict['SchTraffic'][2][23] == pytest.approx(0.2, abs=1e-6) # last
        assert self.uwg.init_param_dict['SchTraffic'][0][19] == pytest.approx(0.8, abs=1e-6)
        assert self.uwg.init_param_dict['SchTraffic'][1][21] == pytest.approx(0.3, abs=1e-6)
        assert self.uwg.init_param_dict['SchTraffic'][2][6] == pytest.approx(0.4, abs=1e-6)

        # test bld fraction list
        assert self.uwg.init_param_dict['bld'][0][0] == pytest.approx(0., abs=1e-6)
        assert self.uwg.init_param_dict['bld'][3][1] == pytest.approx(0.4, abs=1e-6)
        assert self.uwg.init_param_dict['bld'][5][1] == pytest.approx(0.6, abs=1e-6)
        assert self.uwg.init_param_dict['bld'][15][2] == pytest.approx(0.0, abs=1e-6)

        # test BEMs
        assert len(self.uwg.BEM) == pytest.approx(2.,abs=1e-6)
        # test BEM office
        assert self.uwg.BEM[0].building.Type == "LargeOffice"
        assert self.uwg.BEM[0].building.Zone == "1A (Miami)"
        assert self.uwg.BEM[0].building.Era == "Pst80"
        assert self.uwg.BEM[0].frac == 0.4

        # test BEM apartment
        assert self.uwg.BEM[1].building.Type == "MidRiseApartment"
        assert self.uwg.BEM[1].building.Zone == "1A (Miami)"
        assert self.uwg.BEM[1].building.Era == "Pst80"
        assert self.uwg.BEM[1].frac == 0.6

        # Check that schedules are called correctly
        assert self.uwg.Sch[0].Light[0][8] == pytest.approx(0.9, abs=1e-6)   #9am on Weekday for Office
        assert self.uwg.Sch[0].Light[0][7] == pytest.approx(0.3, abs=1e-6)   #9am on Weekday for Office
        assert self.uwg.Sch[1].Occ[1][11] == pytest.approx(0.25, abs=1e-6)     #12 noon on Weekend for apt

if __name__ == "__main__":
    test = TestUWG()
    test.test_read_epw()
    test.test_read_input()
