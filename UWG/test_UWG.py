import os
import pytest
from UWG import UWG


DIR_EPW_PATH = os.path.join(os.path.dirname(__file__),"data/epw/")

def test_read_epw(uwg_):
    uwg_.read_epw()

def test_UWG():
    """Test for UWG.py"""
    epw_dir = DIR_EPW_PATH
    epw_file_name = "SGP_Singapore.486980_IWEC.epdw"
    uwg_param_dir = None
    uwg_param_file_name = None

    uwg = UWG(epw_dir, epw_file_name, uwg_param_dir, uwg_param_file_name)
    test_read_epw(uwg)


if __name__ == '__main__':
    test_UWG()
