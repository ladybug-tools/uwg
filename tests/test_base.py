import os
from uwg import UWG
import math
import logging

TEST_DIR = os.path.dirname(__file__)
DEFAULT_EPW_PATH = os.path.join(
    TEST_DIR, 'epw', 'SGP_Singapore.486980_IWEC.epw')
DEFAULT_PARAM_PATH = os.path.join(
    TEST_DIR, 'parameters', 'initialize_singapore.uwg')
NEW_DIR = os.path.join(TEST_DIR, 'epw_uwg')
MATLAB_DIR = os.path.join(TEST_DIR, 'matlab_ref')


def calculate_tolerance(x, p):
    if abs(float(x)) > 1e-15:
        tol = 1 * 10 ** -(p - 1 - int(math.log10(abs(x))))
    else:
        tol = 1e-15
    return tol


def auto_setup_uwg(epw_path=DEFAULT_EPW_PATH, param_path=DEFAULT_PARAM_PATH,
                   log_file_name=None, log_level=None, uwg_param_dir=None):
    """Set up uwg object from initialize.uwg."""
    if log_file_name and log_level:
        setup_log_file(log_file_name, log_level)

    if param_path:
        model = UWG.from_param_file(
            param_path, epw_path=epw_path, new_epw_dir=NEW_DIR)
    else:
        model = UWG(epw_path, new_epw_dir=NEW_DIR)

    # Increase precision for testing
    model.epw_precision = 16

    return model


def set_input_manually(model):
    """Assign everything manually from ../resources/initialize_singapore.uwg"""

    # Define Simulation and Weather parameters
    model.month = 1.
    model.day = 1.
    model.nday = 31.
    model.dtsim = 300.
    model.dtweather = 3600.

    # HVAC system and internal laod
    model.autosize = 0.
    model.sensocc = 100.
    model.latfocc = 0.3
    model.radfocc = 0.2
    model.radfequip = 0.5
    model.radflight = 0.7

    # Define Urban microclimate parameters
    model.h_ubl1 = 1000.
    model.h_ubl2 = 80.
    model.h_ref = 150.
    model.h_temp = 2.
    model.h_wind = 10.
    model.c_circ = 1.2
    model.c_exch = 1.
    model.maxday = 150.
    model.maxnight = 20.
    model.windmin = 1.
    model.h_obs = 0.1

    # Urban characteristics
    model.bldheight = 10.
    model.h_mix = 1.
    model.blddensity = 0.5
    model.vertohor = 0.8
    model.charlength = 1000.
    model.albroad = 0.1
    model.droad = 0.5
    model.sensanth = 20.

    # Define optional Building characteristics
    model.bld = [('largeoffice', 'pst80', 0.4),
                 ('midriseapartment', 'pst80', 0.6)]

    # climate Zone
    model.zone = '1A'

    # Vegetation parameters
    model.vegstart = 4.0
    model.vegend = 10.0
    model.grasscover = 0.1
    model.treecover = 0.1
    model.albveg = 0.25
    model.rurvegcover = 0.9
    model.latgrss = 0.4
    model.lattree = 0.6

    # Define Traffic schedule
    model.schtraffic = [
        [0.2, 0.2, 0.2, 0.2, 0.2, 0.4, 0.7, 0.9, 0.9, 0.6, 0.6, 0.6, 0.6, 0.6, 0.7, 0.8,
         0.9, 0.9, 0.8, 0.8, 0.7, 0.3, 0.2, 0.2],  # Weekday
        [0.2, 0.2, 0.2, 0.2, 0.2, 0.3, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.6, 0.7,
         0.7, 0.7, 0.7, 0.5, 0.4, 0.3, 0.2, 0.2],  # Saturday
        [0.2, 0.2, 0.2, 0.2, 0.2, 0.3, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
         0.4, 0.4, 0.4, 0.4, 0.3, 0.3, 0.2, 0.2]]  # Sunday

    # Define Road (Assume 0.5m of asphalt)
    model.kroad = 1.
    model.croad = 1600000.

    # Optional parameters are not defined as they are already
    # initialized as None
    return model


def setup_open_matlab_ref(matlab_class_dir, matlab_ref_file_path):
    """ Open the matlab reference file """
    def convert_type(x):
        try:
            return float(x)
        except ValueError:
            return x

    matlab_path = os.path.join(
        MATLAB_DIR, matlab_class_dir, matlab_ref_file_path)
    if not os.path.exists(matlab_path):
        raise Exception("Failed to open {}!".format(matlab_path))

    # Open matlab file, read lines and return values as python simple datatype
    matlab_file = open(matlab_path, 'r')
    uwg_matlab_val_ = [convert_type(x) for x in matlab_file.readlines()]
    matlab_file.close()
    return uwg_matlab_val_


def setup_log_file(log_file_name, log_level=None):
    """
    https://fangpenlin.com/posts/2012/08/26/good-logging-practice-in-python/

    # LOGGER LEVELS
    DEBUG                               # Detailed information
    INFO                                # Confirm things are working as expected
    WARNING                             # Indication something unexpected has happened
    ERROR                               # Error occurs
    CRITICAL                            # Program may not be able to run
    """

    log_file_path = os.path.join(TEST_DIR, "..", "logs", log_file_name)

    if log_level is None:
        log_level = logging.DEBUG  # Default is set to debug for everything

    logging.basicConfig(level=log_level,
                        filename=log_file_path,
                        filemode="w")

    # To stream to console
    # logging.basicConfig(level=logging.DEBUG)
