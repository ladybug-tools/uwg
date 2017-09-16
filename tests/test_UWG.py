import os
import pytest
import UWG

class TestUWG(object):
    """Test for UWG.py
    Naming: Test prefixed test classes (without an __init__ method)
    for test autodetection by pytest
    """
    DIR_UP_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    DIR_EPW_PATH = os.path.join(DIR_UP_PATH,"UWG\\data\\epw\\")

    def setup(self):
        epw_dir = self.DIR_EPW_PATH
        epw_file_name = "SGP_Singapore.486980_IWEC.epw"
        uwg_param_dir = None
        uwg_param_file_name = None

        self.uwg = UWG.UWG(epw_dir, epw_file_name, uwg_param_dir, uwg_param_file_name)


    def test_read_epw(self):
        self.uwg.read_epw()
        assert self.uwg.climateDataFile[0][0] == "LOCATION"
        assert self.uwg.climateDataFile[0][1] == "SINGAPORE"
