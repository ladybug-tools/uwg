import os
import uwg
import math
import logging


class TestBase(object):
    """
    Test base class. This is inherited by all other tests.
    """

    DIR_CURR = os.path.abspath(os.path.dirname(__file__))
    DIR_EPW_PATH = os.path.join(DIR_CURR,"epw")
    DIR_UWGPARAM_PATH = os.path.join(DIR_CURR,"parameters")
    DIR_DESTINATION_PATH = os.path.join(DIR_CURR,"epw_uwg")
    DIR_MATLAB_PATH = os.path.join(DIR_CURR, "matlab_ref")

    RUN_LOG = True

    def CALCULATE_TOLERANCE(self,x,p):
        if abs(float(x)) > 1e-15:
            tol = 1*10**-(p - 1 - int(math.log10(abs(x))))
        else:
            tol = 1e-15
        return tol

    def setup_uwg_integration(self, epw_file="SGP_Singapore.486980_IWEC.epw", uwg_param_file="initialize_singapore.uwg",
        log_file_name=None,log_level=None, uwg_param_dir=None):

        """ set up uwg object from initialize.uwg """

        epw_dir = self.DIR_EPW_PATH
        uwg_param_dir_ = self.DIR_UWGPARAM_PATH if uwg_param_dir == None else uwg_param_dir
        destination_dir = self.DIR_DESTINATION_PATH

        if log_file_name and log_level:
            self.setup_log_file(log_file_name, log_level)

        self.uwg = uwg.uwg(epw_file, uwg_param_file, epwDir=epw_dir, uwgParamDir=uwg_param_dir_, destinationDir=destination_dir)

        # Increase precision for testing
        self.uwg.epw_precision = 16

    def set_input_manually(self):

        """Assign everything manually from ../resources/initialize_singapore.uwg"""

        # Define Simulation and Weather parameters
        self.uwg.Month = 1.
        self.uwg.Day = 1.
        self.uwg.nDay = 31.
        self.uwg.dtSim = 300.
        self.uwg.dtWeather = 3600.

        # HVAC system and internal laod
        self.uwg.autosize = 0.
        self.uwg.sensOcc = 100.
        self.uwg.LatFOcc = 0.3
        self.uwg.RadFOcc = 0.2
        self.uwg.RadFEquip = 0.5
        self.uwg.RadFLight = 0.7

        # Define Urban microclimate parameters
        self.uwg.h_ubl1 = 1000.
        self.uwg.h_ubl2 = 80.
        self.uwg.h_ref = 150.
        self.uwg.h_temp = 2.
        self.uwg.h_wind = 10.
        self.uwg.c_circ = 1.2
        self.uwg.c_exch = 1.
        self.uwg.maxDay = 150.
        self.uwg.maxNight = 20.
        self.uwg.windMin = 1.
        self.uwg.h_obs = 0.1

        # Urban characteristics
        self.uwg.bldHeight = 10.
        self.uwg.h_mix = 1.
        self.uwg.bldDensity = 0.5
        self.uwg.verToHor = 0.8
        self.uwg.charLength = 1000.
        self.uwg.alb_road = 0.1
        self.uwg.d_road = 0.5
        self.uwg.sensAnth = 20.
        self.uwg.latAnth = 2.

        # Define optional Building characteristics
        self.uwg.bld = [
            [0,0,0],    # FullServiceRestaurant
            [0,0,0],    # Hospital
            [0,0,0],    # LargeHotel
            [0,.4,0],   # LargeOffice
            [0,0,0],    # MediumOffice
            [0,.6,0],   # MidRiseApartment
            [0,0,0],    # OutPatient
            [0,0,0],    # PrimarySchool
            [0,0,0],    # QuickServiceRestaurant
            [0,0,0],    # SecondarySchool
            [0,0,0],    # SmallHotel
            [0,0,0],    # SmallOffice
            [0,0,0],    # Stand-aloneRetail
            [0,0,0],    # StripMall
            [0,0,0],    # SuperMarket
            [0,0,0]     # Warehouse
        ]

        # climate Zone
        self.uwg.zone = 1.0

        # Vegetation parameters
        self.uwg.vegCover = 0.2
        self.uwg.vegStart = 4.0
        self.uwg.treeCoverage = 0.1
        self.uwg.vegEnd = 10.0
        self.uwg.albVeg = 0.25
        self.uwg.rurVegCover = 0.9
        self.uwg.latGrss = 0.4
        self.uwg.latTree = 0.6

        # Define Traffic schedule
        self.uwg.SchTraffic = [
            [0.2,0.2,0.2,0.2,0.2,0.4,0.7,0.9,0.9,0.6,0.6,0.6,0.6,0.6,0.7,0.8,0.9,0.9,0.8,0.8,0.7,0.3,0.2,0.2], # Weekday
            [0.2,0.2,0.2,0.2,0.2,0.3,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.6,0.7,0.7,0.7,0.7,0.5,0.4,0.3,0.2,0.2], # Saturday
            [0.2,0.2,0.2,0.2,0.2,0.3,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.3,0.3,0.2,0.2]  # Sunday
        ]

        # Define Road (Assume 0.5m of asphalt)
        self.uwg.kRoad = 1.
        self.uwg.cRoad = 1600000.

        # Optional parameters are not defined as they are already
        # initialized as None

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

        # Open matlab file, read lines and return values as python simple datatype
        matlab_file = open(matlab_path,'r')
        uwg_matlab_val_ = [convert_type(x) for x in matlab_file.readlines()]
        matlab_file.close()
        return uwg_matlab_val_


    def setup_log_file(self, log_file_name, log_level=None):
        """
        https://fangpenlin.com/posts/2012/08/26/good-logging-practice-in-python/

        # LOGGER LEVELS
        DEBUG                               # Detailed information
        INFO                                # Confirm things are working as expected
        WARNING                             # Indication something unexpected has happened
        ERROR                               # Error occurs
        CRITICAL                            # Program may not be able to run
        """

        log_file_path = os.path.join(self.DIR_CURR, "..", "logs", log_file_name)

        if log_level==None:
            log_level = logging.DEBUG # Default is set to debug for everything

        logging.basicConfig(level=log_level,
                            filename=log_file_path,
                            filemode="w")
        # To stream to console
        #logging.basicConfig(level=logging.DEBUG)
