import pytest
import UWG
import os
import math
import pprint

class TestUCMDef(object):

    DIR_UP_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    DIR_EPW_PATH = os.path.join(DIR_UP_PATH,"resources/epw")
    DIR_MATLAB_PATH = os.path.join(os.path.abspath(os.path.dirname(__file__)), "matlab_ref","matlab_ucm")

    def setup_uwg(self):
        """ set up uwg object from initialize.uwg """

        epw_dir = self.DIR_EPW_PATH
        epw_file_name = "SGP_Singapore.486980_IWEC.epw"
        uwg_param_dir = os.path.join(self.DIR_UP_PATH,"resources")
        uwg_param_file_name = "initialize.uwg"

        self.uwg = UWG.UWG(epw_dir, epw_file_name, uwg_param_dir, uwg_param_file_name)

    def test_ucm_init(self):
        """ test ucm constructor """

        self.setup_uwg()
        self.uwg.read_epw()
        self.uwg.read_input()


        # open matlab ref file
        matlab_path = os.path.join(self.DIR_MATLAB_PATH,"matlab_ref_ucm_init.txt")
        if not os.path.exists(matlab_path):
            raise Exception("Failed to open {}!".format(matlab_path))
        matlab_file = open(matlab_path,'r')
        uwg_matlab_val = [float(x) for x in matlab_file.readlines()]
        matlab_file.close()

        uwg_python_val = [
        self.uwg.UCM.h_mix,
        self.uwg.UCM.bldHeight,
        self.uwg.UCM.verToHor,
        self.uwg.UCM.bldDensity,
        self.uwg.UCM.treeCoverage,
        self.uwg.UCM.sensAnthrop,
        self.uwg.UCM.latAnthrop,
        self.uwg.UCM.roadShad,
        self.uwg.UCM.bldWidth,
        self.uwg.UCM.canWidth,
        self.uwg.UCM.canAspect,
        self.uwg.UCM.roadConf,
        self.uwg.UCM.wallConf,
        self.uwg.UCM.facArea,
        self.uwg.UCM.roadArea,
        self.uwg.UCM.roofArea,
        self.uwg.UCM.canTemp,
        self.uwg.UCM.roadTemp,
        self.uwg.UCM.canHum,
        self.uwg.UCM.ublWind,
        self.uwg.UCM.canWind,
        self.uwg.UCM.ustar,
        self.uwg.UCM.ustarMod,
        self.uwg.UCM.z0u,
        self.uwg.UCM.l_disp,
        self.uwg.UCM.alb_wall,
        self.uwg.UCM.facAbsor,
        self.uwg.UCM.roadAbsor,
        self.uwg.UCM.sensHeat
        ]

        # matlab ref checking
        assert len(uwg_matlab_val) == len(uwg_python_val)
        for i in xrange(len(uwg_matlab_val)):
            #print uwg_python_val[i], uwg_matlab_val[i]
            assert uwg_python_val[i] == pytest.approx(uwg_matlab_val[i], abs=1e-15), "error at index={}".format(i)

if __name__ == "__main__":
    tucm = TestUCMDef()
    tucm.test_ucm_init()
