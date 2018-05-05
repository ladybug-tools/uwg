import os
import UWG
import math

class TestBase(object):
    """
    Test base class. This is inherited by all other tests.
    """

    DIR_CURR = os.path.abspath(os.path.dirname(__file__))
    DIR_EPW_PATH = os.path.join(DIR_CURR,"epw")
    DIR_UWGPARAM_PATH = os.path.join(DIR_CURR,"parameters")
    DIR_DESTINATION_PATH = os.path.join(DIR_CURR,"epw_uwg")
    DIR_MATLAB_PATH = os.path.join(DIR_CURR, "matlab_ref")

    def CALCULATE_TOLERANCE(self,x,p):
        if abs(float(x)) > 1e-15:
            tol = 1*10**-(p - 1 - int(math.log10(abs(x))))
        else:
            tol = 1e-15
        return tol

    def setup_uwg_integration(self, epw_file="SGP_Singapore.486980_IWEC.epw", uwg_param_file="initialize_singapore.uwg", uwg_param_dir=None):
        """ set up uwg object from initialize.uwg """

        epw_dir = self.DIR_EPW_PATH
        uwg_param_dir_ = self.DIR_UWGPARAM_PATH if uwg_param_dir == None else uwg_param_dir
        destination_dir = self.DIR_DESTINATION_PATH

        self.uwg = UWG.UWG(epw_dir, epw_file, uwg_param_dir_, uwg_param_file, destinationDir=destination_dir)

    def setup_open_matlab_ref(self, matlab_class_dir, matlab_ref_file_path):
        """ open the matlab reference file """
        def convert_type(x):
            try:
                return float(x)
            except ValueError:
                return x

        matlab_path = os.path.join(self.DIR_MATLAB_PATH, matlab_class_dir, matlab_ref_file_path)
        if not os.path.exists(matlab_path):
            raise Exception("Failed to open {}!".format(matlab_path))
        matlab_file = open(matlab_path,'r')
        uwg_matlab_val_ = [convert_type(x) for x in matlab_file.readlines()]
        matlab_file.close()
        return uwg_matlab_val_
