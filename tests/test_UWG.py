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
        uwg_param_dir = None
        uwg_param_file_name = None

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

if __name__ == "__main__":
    test = TestUWG()
    test.test_read_epw()
