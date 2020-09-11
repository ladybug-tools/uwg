import os
from uwg import UWG
import math
import logging

TEST_DIR = os.path.abspath(os.path.dirname(__file__))
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
        testuwg = UWG.from_param_file(epw_path, param_path, new_epw_dir=NEW_DIR)
    else:
        testuwg = UWG(epw_path, new_epw_dir=NEW_DIR)

    # Increase precision for testing
    testuwg.epw_precision = 16

    return testuwg


def set_input_manually(testuwg):

    """Assign everything manually from ../resources/initialize_singapore.uwg"""

    # Define Simulation and Weather parameters
    testuwg.month = 1.
    testuwg.day = 1.
    testuwg.nday = 31.
    testuwg.dtsim = 300.
    testuwg.dtweather = 3600.

    # HVAC system and internal laod
    testuwg.autosize = 0.
    testuwg.sensocc = 100.
    testuwg.latfocc = 0.3
    testuwg.radfocc = 0.2
    testuwg.radfequip = 0.5
    testuwg.radflight = 0.7

    # Define Urban microclimate parameters
    testuwg.h_ubl1 = 1000.
    testuwg.h_ubl2 = 80.
    testuwg.h_ref = 150.
    testuwg.h_temp = 2.
    testuwg.h_wind = 10.
    testuwg.c_circ = 1.2
    testuwg.c_exch = 1.
    testuwg.maxday = 150.
    testuwg.maxnight = 20.
    testuwg.windmin = 1.
    testuwg.h_obs = 0.1

    # Urban characteristics
    testuwg.bldheight = 10.
    testuwg.h_mix = 1.
    testuwg.blddensity = 0.5
    testuwg.vertohor = 0.8
    testuwg.charlength = 1000.
    testuwg.albroad = 0.1
    testuwg.droad = 0.5
    testuwg.sensanth = 20.

    # Define optional Building characteristics
    testuwg.bld = [
        [0, 0, 0],    # FullServiceRestaurant
        [0, 0, 0],    # Hospital
        [0, 0, 0],    # LargeHotel
        [0, .4, 0],   # LargeOffice
        [0, 0, 0],    # MediumOffice
        [0, .6, 0],   # MidRiseApartment
        [0, 0, 0],    # OutPatient
        [0, 0, 0],    # PrimarySchool
        [0, 0, 0],    # QuickServiceRestaurant
        [0, 0, 0],    # SecondarySchool
        [0, 0, 0],    # SmallHotel
        [0, 0, 0],    # SmallOffice
        [0, 0, 0],    # Stand-aloneRetail
        [0, 0, 0],    # StripMall
        [0, 0, 0],    # SuperMarket
        [0, 0, 0]     # Warehouse
    ]

    # climate Zone
    testuwg.zone = 1.0

    # Vegetation parameters
    testuwg.vegcover = 0.2
    testuwg.vegstart = 4.0
    testuwg.treecoverage = 0.1
    testuwg.vegend = 10.0
    testuwg.albveg = 0.25
    testuwg.rurvegcover = 0.9
    testuwg.latgrss = 0.4
    testuwg.lattree = 0.6

    # Define Traffic schedule
    testuwg.schtraffic = [
        [0.2, 0.2, 0.2, 0.2, 0.2, 0.4, 0.7, 0.9, 0.9, 0.6, 0.6, 0.6, 0.6, 0.6, 0.7, 0.8,
         0.9, 0.9, 0.8, 0.8, 0.7, 0.3, 0.2, 0.2],  # Weekday
        [0.2, 0.2, 0.2, 0.2, 0.2, 0.3, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.6, 0.7,
         0.7, 0.7, 0.7, 0.5, 0.4, 0.3, 0.2, 0.2],  # Saturday
        [0.2, 0.2, 0.2, 0.2, 0.2, 0.3, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
         0.4, 0.4, 0.4, 0.4, 0.3, 0.3, 0.2, 0.2]]  # Sunday

    # Define Road (Assume 0.5m of asphalt)
    testuwg.kroad = 1.
    testuwg.croad = 1600000.

    # Optional parameters are not defined as they are already
    # initialized as None
    return testuwg


def setup_open_matlab_ref(matlab_class_dir, matlab_ref_file_path):
    """ Open the matlab reference file """
    def convert_type(x):
        try:
            return float(x)
        except ValueError:
            return x

    matlab_path = os.path.join(MATLAB_DIR, matlab_class_dir, matlab_ref_file_path)
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
