"""
=========================================================================
 THE URBAN WEATHER GENERATOR (uwg)
=========================================================================
Version 4.2

Original Author: B. Bueno
Edited by A. Nakano & Lingfu Zhang
Modified by Joseph Yang (joeyang@mit.edu) - May, 2016
Translated to Python by Saeran Vasanthakumar - February, 2018

Original Pulbication on the uwg's Methods:
Bueno, Bruno; Norford, Leslie; Hidalgo, Julia; Pigeon, Gregoire (2013).
The urban weather generator, Journal of Building Performance Simulation. 6:4,269-281.
doi: 10.1080/19401493.2012.718797
=========================================================================
"""
from __future__ import division, print_function
from functools import reduce

try:
    range = xrange
except NameError:
    pass

import os
import math
import copy
import logging

try:
    import cPickle as pickle
except ImportError:
    import pickle

from .simparam import SimParam
from .weather import Weather
from .building import Building
from .material import Material
from .element import Element
from .BEMDef import BEMDef
from .schdef import SchDef
from .param import Param
from .UCMDef import UCMDef
from .forcing import Forcing
from .UBLDef import UBLDef
from .RSMDef import RSMDef
from .solarcalcs import SolarCalcs
from .psychrometrics import psychrometrics
from .readDOE import readDOE
from .urbflux import urbflux
from . import utilities

# For debugging only
#from pprint import pprint
#from decimal import Decimal
#pp = pprint
#dd = Decimal.from_float


class uwg(object):
    """Morph a rural EPW file to urban conditions using a file with a list of urban parameters.

    args:
        epwDir: The directory in which the rural EPW file sits.
        epwFileName: The name of the rural epw file that will be morphed.
        uwgParamDir: The directory in which the uwg Parameter File (.uwg) sits.
        uwgParamFileName: The name of the uwg Parameter File (.uwg).
        destinationDir: Optional destination directory for the morphed EPW file.
            If left blank, the morphed file will be written into the same directory
            as the rural EPW file (the epwDir).
        destinationFileName: Optional destination file name for the morphed EPW file.
            If left blank, the morphed file will append "_UWG" to the original file name.
    returns:
        newClimateFile: the path to a new EPW file that has been morphed to account
            for uban conditions.
    """

    """ Section 1 - Definitions for constants / other parameters """
    MINTHICKNESS = 0.01    # Minimum layer thickness (to prevent crashing) (m)
    MAXTHICKNESS = 0.05    # Maximum layer thickness (m)
    # http://web.mit.edu/parmstr/Public/NRCan/nrcc29118.pdf (Figly & Snodgrass)
    SOILTCOND = 1
    # http://www.europment.org/library/2013/venice/bypaper/MFHEEF/MFHEEF-21.pdf (average taken from Table 1)
    SOILVOLHEAT = 2e6
    # Soil material used for soil-depth padding
    SOIL = Material(SOILTCOND, SOILVOLHEAT, name="soil")

    # Physical constants
    G = 9.81               # gravity (m s-2)
    CP = 1004.             # heat capacity for air (J/kg K)
    VK = 0.40              # von karman constant (dimensionless)
    R = 287.               # gas constant dry air (J/kg K)
    RV = 461.5             # gas constant water vapor (J/kg K)
    LV = 2.26e6            # latent heat of evaporation (J/kg)
    SIGMA = 5.67e-08       # Stefan Boltzmann constant (W m-2 K-4)
    WATERDENS = 1000.      # water density (kg m-3)
    LVTT = 2.5008e6        #
    TT = 273.16            #
    ESTT = 611.14          #
    CL = 4.218e3           #
    CPV = 1846.1           #
    B = 9.4                # Coefficients derived by Louis (1979)
    CM = 7.4               #
    # (Pr/Sc)^(2/3) for Colburn analogy in water evaporation
    COLBURN = math.pow((0.713/0.621), (2/3.))

    # Site-specific parameters
    WGMAX = 0.005  # maximum film water depth on horizontal surfaces (m)

    # File path parameter
    RESOURCE_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "resources"))
    CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))

    def __init__(self, epwFileName, uwgParamFileName=None, epwDir=None, uwgParamDir=None, destinationDir=None, destinationFileName=None):

        # Logger will be disabled by default unless explicitly called in tests
        self.logger = logging.getLogger(__name__)

        # User defined
        self.epwFileName = epwFileName if epwFileName.lower().endswith('.epw') else epwFileName + \
            '.epw'  # Revise epw file name if not end with epw
        # If file name is entered then will uwg will set input from .uwg file
        self.uwgParamFileName = uwgParamFileName

        # If user does not overload
        self.destinationFileName = destinationFileName if destinationFileName else self.epwFileName.strip(
            '.epw') + '_UWG.epw'
        self.epwDir = epwDir if epwDir else os.path.join(self.RESOURCE_PATH, "epw")
        self.uwgParamDir = uwgParamDir if uwgParamDir else os.path.join(
            self.RESOURCE_PATH, "parameters")
        self.destinationDir = destinationDir if destinationDir else os.path.join(
            self.RESOURCE_PATH, "epw_uwg")

        # refdata: Serialized DOE reference data, z_meso height data
        self.readDOE_file_path = os.path.join(self.CURRENT_PATH, "refdata", "readDOE.pkl")
        self.z_meso_dir_path = os.path.join(self.CURRENT_PATH, "refdata")

        # EPW precision
        self.epw_precision = 1

        # init uwg variables
        self._init_param_dict = None

        # Define Simulation and Weather parameters
        self.Month = None        # starting month (1-12)
        self.Day = None          # starting day (1-31)
        self.nDay = None         # number of days
        self.dtSim = None        # simulation time step (s)
        self.dtWeather = None    # seconds (s)

        # HVAC system and internal laod
        self.autosize = None     # autosize HVAC (1 or 0)
        self.sensOcc = None      # Sensible heat from occupant
        self.LatFOcc = None      # Latent heat fraction from occupant (normally 0.3)
        self.RadFOcc = None      # Radiant heat fraction from occupant (normally 0.2)
        self.RadFEquip = None    # Radiant heat fraction from equipment (normally 0.5)
        self.RadFLight = None    # Radiant heat fraction from light (normally 0.7)

        # Define Urban microclimate parameters
        self.h_ubl1 = None       # ubl height - day (m)
        self.h_ubl2 = None       # ubl height - night (m)
        self.h_ref = None        # inversion height
        self.h_temp = None       # temperature height
        self.h_wind = None       # wind height
        self.c_circ = None       # circulation coefficient
        self.c_exch = None       # exchange coefficient
        self.maxDay = None       # max day threshhold
        self.maxNight = None     # max night threshhold
        self.windMin = None      # min wind speed (m/s)
        self.h_obs = None        # rural average obstacle height

        # Urban characteristics
        self.bldHeight = None    # average building height (m)
        self.h_mix = None        # mixing height (m)
        self.bldDensity = None   # building density (0-1)
        self.verToHor = None     # building aspect ratio
        # radius defining the urban area of study [aka. characteristic length] (m)
        self.charLength = None
        self.alb_road = None     # road albedo
        self.d_road = None       # road pavement thickness
        self.sensAnth = None     # non-building sensible heat (W/m^2)
        self.latAnth = None      # non-building latent heat heat (W/m^2). Not used, taken out by JH.


        # Fraction of building typology stock
        self.bld = None         # 16x3 matrix of fraction of building type by era

        # climate Zone
        self.zone = None

        # Vegetation parameters
        self.vegCover = None     # urban area veg coverage ratio
        self.treeCoverage = None  # urban area tree coverage ratio
        self.vegStart = None     # vegetation start month
        self.vegEnd = None       # vegetation end month
        self.albVeg = None       # Vegetation albedo
        self.rurVegCover = None  # rural vegetation cover
        self.latGrss = None      # latent fraction of grass
        self.latTree = None      # latent fraction of tree

        # Define Traffic schedule
        self.SchTraffic = None

        # Define Road (Assume 0.5m of asphalt)
        self.kRoad = None       # road pavement conductivity (W/m K)
        self.cRoad = None       # road volumetric heat capacity (J/m^3 K)

        # Define optional Building characteristics
        self.flr_h = None       # floor-to-floor height
        self.albRoof = None     # roof albedo (0 - 1)
        self.vegRoof = None     # Fraction of the roofs covered in grass/shrubs (0-1)
        self.glzR = None        # Glazing Ratio
        self.SHGC = None       # Solar Heat Gain Coefficient
        self.albWall = None    # Wall albedo

    def ToString(self):
        """Overwrite .NET ToString method."""
        return self.__repr__()

    def __repr__(self):
        def _split_string(s):
            return s[0] + ":\n  " + s[1].replace(",", "\n  ")

        def _tabbed(s):
            return _split_string(s.__repr__().split(":"))

        def _list_2_tabbed(b):
            return reduce(lambda a, b: a+"\n"+b, [_tabbed(_b) for _b in b])

        return "uwg for {}:\n\n{}{}{}{}{}{}{}{}".format(
            self.epwFileName,
            _tabbed(self.simTime)+"\n" if hasattr(self, "simTime") else "No simTime attr.\n",
            _tabbed(self.weather)+"\n" if hasattr(self, "weather") else "No weather attr.\n",
            _tabbed(self.geoParam)+"\n" if hasattr(self, "geoParam") else "No geoParam attr.\n",
            _tabbed(self.UBL)+"\n" if hasattr(self, "UBL") else "No UBL attr.\n",
            "Rural "+_tabbed(self.RSM)+"\n" if hasattr(self, "RSM") else "No Rural RSM attr.\n",
            "Urban "+_tabbed(self.USM)+"\n" if hasattr(self, "USM") else "No Urban RSM attr.\n",
            _tabbed(self.UCM)+"\n" if hasattr(self, "UCM") else "No UCM attr.\n",
            _list_2_tabbed(self.BEM) if hasattr(self, "BEM") else "No BEM attr."
            )

    def is_near_zero(self,num,eps=1e-10):
        return abs(float(num)) < eps

    def read_epw(self):
        """Section 2 - Read EPW file
        properties:
            self.climateDataPath
            self.newPathName
            self._header    # header data
            self.epwinput   # timestep data for weather
            self.lat        # latitude
            self.lon        # longitude
            self.GMT        # GMT
            self.nSoil      # Number of soil depths
            self.Tsoil      # nSoil x 12 matrix for soil temperture (K)
            self.depth_soil # nSoil x 1 matrix for soil depth (m)
        """

        # Make dir path to epw file
        self.climateDataPath = os.path.join(self.epwDir, self.epwFileName)

        # Open epw file and feed csv data to climate_data
        try:
            climate_data = utilities.read_csv(self.climateDataPath)
        except Exception as e:
            raise Exception("Failed to read epw file! {}".format(e.message))

        # Read header lines (1 to 8) from EPW and ensure TMY2 format.
        self._header = climate_data[0:8]

        # Read weather data from EPW for each time step in weather file. (lines 8 - end)
        self.epwinput = climate_data[8:]

        # Read Lat, Long (line 1 of EPW)
        self.lat = float(self._header[0][6])
        self.lon = float(self._header[0][7])
        self.GMT = float(self._header[0][8])

        # Read in soil temperature data (assumes this is always there)
        # ref: http://bigladdersoftware.com/epx/docs/8-2/auxiliary-programs/epw-csv-format-inout.html
        soilData = self._header[3]
        self.nSoil = int(soilData[1])           # Number of ground temperature depths
        self.Tsoil = utilities.zeros(self.nSoil, 12)  # nSoil x 12 matrix for soil temperture (K)
        self.depth_soil = utilities.zeros(self.nSoil, 1)   # nSoil x 1 matrix for soil depth (m)

        # Read monthly data for each layer of soil from EPW file
        for i in range(self.nSoil):
            self.depth_soil[i][0] = float(soilData[2 + (i*16)])  # get soil depth for each nSoil
            # Monthly data
            for j in range(12):
                # 12 months of soil T for specific depth
                self.Tsoil[i][j] = float(soilData[6 + (i*16) + j]) + 273.15

        # Set new directory path for the moprhed EPW file
        self.newPathName = os.path.join(self.destinationDir, self.destinationFileName)

    def read_input(self):
        """Section 3 - Read Input File (.m, file)
        Note: UWG_Matlab input files are xlsm, XML, .m, file.
        properties:
            self._init_param_dict   # dictionary of simulation initialization parameters

            self.sensAnth           # non-building sensible heat (W/m^2)
            self.SchTraffic         # Traffice schedule

            self.BEM                # list of BEMDef objects extracted from readDOE
            self.Sch                # list of Schedule objects extracted from readDOE

        """

        uwg_param_file_path = os.path.join(self.uwgParamDir, self.uwgParamFileName)

        if not os.path.exists(uwg_param_file_path):
            raise Exception("Param file: '{}' does not exist.".format(uwg_param_file_path))

        # Open .uwg file and feed csv data to initializeDataFile
        try:
            uwg_param_data = utilities.read_csv(uwg_param_file_path)
        except Exception as e:
            raise Exception("Failed to read .uwg file! {}".format(e.message))

        # The initialize.uwg is read with a dictionary so that users changing
        # line endings or line numbers doesn't make reading input incorrect
        self._init_param_dict = {}
        count = 0
        while count < len(uwg_param_data):
            row = uwg_param_data[count]
            row = [row[i].replace(" ", "") for i in range(len(row))]  # strip white spaces

            # Optional parameters might be empty so handle separately
            is_optional_parameter = (
                row != [] and
                (
                    row[0] == "albRoof" or
                    row[0] == "vegRoof" or
                    row[0] == "glzR" or
                    row[0] == "hvac" or
                    row[0] == "albWall" or
                    row[0] == "SHGC"
                )
            )
            try:
                if row == [] or "#" in row[0]:
                    count += 1
                    continue
                elif row[0] == "SchTraffic":
                    # SchTraffic: 3 x 24 matrix
                    trafficrows = uwg_param_data[count+1:count+4]
                    self._init_param_dict[row[0]] = [utilities.str2fl(r[:24]) for r in trafficrows]
                    count += 4
                elif row[0] == "bld":
                    # bld: 17 x 3 matrix
                    bldrows = uwg_param_data[count+1:count+17]
                    self._init_param_dict[row[0]] = [utilities.str2fl(r[:3]) for r in bldrows]
                    count += 17
                elif is_optional_parameter:
                    self._init_param_dict[row[0]] = float(row[1]) if row[1] != "" else None
                    count += 1
                else:
                    self._init_param_dict[row[0]] = float(row[1])
                    count += 1
            except ValueError:
                print("Error while reading parameter at {} {}".format(count, row))

        ipd = self._init_param_dict

        # Define Simulation and Weather parameters
        if self.Month is None: self.Month = ipd['Month']
        if self.Day is None: self.Day = ipd['Day']
        if self.nDay is None: self.nDay = ipd['nDay']
        if self.dtSim is None: self.dtSim = ipd['dtSim']
        if self.dtWeather is None: self.dtWeather = ipd['dtWeather']

        # HVAC system and internal laod
        if self.autosize is None: self.autosize = ipd['autosize']
        if self.sensOcc is None: self.sensOcc = ipd['sensOcc']
        if self.LatFOcc is None: self.LatFOcc = ipd['LatFOcc']
        if self.RadFOcc is None: self.RadFOcc = ipd['RadFOcc']
        if self.RadFEquip is None: self.RadFEquip = ipd['RadFEquip']
        if self.RadFLight is None: self.RadFLight = ipd['RadFLight']

        # Define Urban microclimate parameters
        if self.h_ubl1 is None: self.h_ubl1 = ipd['h_ubl1']
        if self.h_ubl2 is None: self.h_ubl2 = ipd['h_ubl2']
        if self.h_ref is None: self.h_ref = ipd['h_ref']
        if self.h_temp is None: self.h_temp = ipd['h_temp']
        if self.h_wind is None: self.h_wind = ipd['h_wind']
        if self.c_circ is None: self.c_circ = ipd['c_circ']
        if self.c_exch is None: self.c_exch = ipd['c_exch']
        if self.maxDay is None: self.maxDay = ipd['maxDay']
        if self.maxNight is None: self.maxNight = ipd['maxNight']
        if self.windMin is None: self.windMin = ipd['windMin']
        if self.h_obs is None: self.h_obs = ipd['h_obs']

        # Urban characteristics
        if self.bldHeight is None: self.bldHeight = ipd['bldHeight']
        if self.h_mix is None: self.h_mix = ipd['h_mix']
        if self.bldDensity is None: self.bldDensity = ipd['bldDensity']
        if self.verToHor is None: self.verToHor = ipd['verToHor']
        if self.charLength is None: self.charLength = ipd['charLength']
        if self.alb_road is None: self.alb_road = ipd['albRoad']
        if self.d_road is None: self.d_road = ipd['dRoad']
        if self.sensAnth is None: self.sensAnth = ipd['sensAnth']
        # if self.latAnth is None: self.latAnth = ipd['latAnth'] # Not used, taken out by JH.

        # climate Zone
        if self.zone is None: self.zone = ipd['zone']

        # Vegetation parameters
        if self.vegCover is None: self.vegCover = ipd['vegCover']
        if self.treeCoverage is None: self.treeCoverage = ipd['treeCoverage']
        if self.vegStart is None: self.vegStart = ipd['vegStart']
        if self.vegEnd is None: self.vegEnd = ipd['vegEnd']
        if self.albVeg is None: self.albVeg = ipd['albVeg']
        if self.rurVegCover is None: self.rurVegCover = ipd['rurVegCover']
        if self.latGrss is None: self.latGrss = ipd['latGrss']
        if self.latTree is None: self.latTree = ipd['latTree']

        # Define Traffic schedule
        if self.SchTraffic is None: self.SchTraffic = ipd['SchTraffic']

        # Define Road (Assume 0.5m of asphalt)
        if self.kRoad is None: self.kRoad = ipd['kRoad']
        if self.cRoad is None: self.cRoad = ipd['cRoad']

        # Building stock fraction
        if self.bld is None: self.bld = ipd['bld']

        # Optional parameters
        if self.albRoof is None: self.albRoof = ipd['albRoof']
        if self.vegRoof is None: self.vegRoof = ipd['vegRoof']
        if self.glzR is None: self.glzR = ipd['glzR']
        if self.albWall is None: self.albWall = ipd['albWall']
        if self.SHGC is None: self.SHGC = ipd['SHGC']

    def check_required_inputs(self):
        # Fail if required parameters aren't correct
        assert isinstance(self.Month, (float, int)), \
            'Month must be a number. Got {}'.format(type(self.Month))
        assert isinstance(self.Day, (float, int)), \
            'Day must be a number. Got {}'.format(type(self.Day))
        assert isinstance(self.nDay, (float, int)), \
            'nDay must be a number. Got {}'.format(type(self.nDay))
        assert isinstance(self.dtSim, float), \
            'dtSim must be a float. Got {}'.format(type(self.dtSim))
        assert isinstance(self.dtWeather, float), \
            'dtWeather must be a float. Got {}'.format(type(self.dtWeather))
        assert isinstance(self.autosize, (float, int)), \
            'autosize must be a number. Got {}'.format(type(self.autosize))
        assert isinstance(self.sensOcc, float), \
            'sensOcc must be a float. Got {}'.format(type(self.sensOcc))
        assert isinstance(self.LatFOcc, float), \
            'LatFOcc must be a float. Got {}'.format(type(self.LatFOcc))
        assert isinstance(self.RadFOcc, float), \
            'RadFOcc must be a float. Got {}'.format(type(self.RadFOcc))
        assert isinstance(self.RadFEquip, float), \
            'RadFEquip must be a float. Got {}'.format(type(self.RadFEquip))
        assert isinstance(self.RadFLight, float), \
            'RadFLight must be a float. Got {}'.format(type(self.RadFLight))
        assert isinstance(self.h_ubl1, float), \
            'h_ubl1 must be a float. Got {}'.format(type(self.h_ubl1))
        assert isinstance(self.h_ubl2, float), \
            'h_ubl2 must be a float. Got {}'.format(type(self.h_ubl2))
        assert isinstance(self.h_ref, float), \
            'h_ref must be a float. Got {}'.format(type(self.h_ref))
        assert isinstance(self.h_temp, float), \
            'h_temp must be a float. Got {}'.format(type(self.h_temp))
        assert isinstance(self.h_wind, float), \
            'h_wind must be a float. Got {}'.format(type(self.h_wind))
        assert isinstance(self.c_circ, float), \
            'c_circ must be a float. Got {}'.format(type(self.c_circ))
        assert isinstance(self.c_exch, float), \
            'c_exch must be a float. Got {}'.format(type(self.c_exch))
        assert isinstance(self.maxDay, float), \
            'maxDay must be a float. Got {}'.format(type(self.maxDay))
        assert isinstance(self.maxNight, float), \
            'maxNight must be a float. Got {}'.format(type(self.maxNight))
        assert isinstance(self.windMin, float), \
            'windMin must be a float. Got {}'.format(type(self.windMin))
        assert isinstance(self.h_obs, float), \
            'h_obs must be a float. Got {}'.format(type(self.h_obs))
        assert isinstance(self.bldHeight, float), \
            'bldHeight must be a float. Got {}'.format(type(self.bldHeight))
        assert isinstance(self.h_mix, float), \
            'h_mix must be a float. Got {}'.format(type(self.h_mix))
        assert isinstance(self.bldDensity, float), \
            'bldDensity must be a float. Got {}'.format(type(self.bldDensity))
        assert isinstance(self.verToHor, float), \
            'verToHor must be a float. Got {}'.format(type(self.verToHor))
        assert isinstance(self.charLength, float), \
            'charLength must be a float. Got {}'.format(type(self.charLength))
        assert isinstance(self.alb_road, float), \
            'alb_road must be a float. Got {}'.format(type(self.alb_road))
        assert isinstance(self.d_road, float), \
            'd_road must be a float. Got {}'.format(type(self.d_road))
        assert isinstance(self.sensAnth, float), \
            'sensAnth must be a float. Got {}'.format(type(self.sensAnth))
        # assert isinstance(self.latAnth, float) # Take this out as isn't being used
        assert isinstance(self.bld, list), \
            'bld must be a list. Got {}'.format(type(self.bld))
        assert len(self.bld) == 16, \
            'length of bld must be 16. Got {}'.format(len(self.bld))
        assert isinstance(self.latTree, float), \
            'latTree must be a float. Got {}'.format(type(self.latTree))
        assert isinstance(self.latGrss, float), \
            'latGrss must be a float. Got {}'.format(type(self.latGrss))
        assert isinstance(self.zone, (float, int)), \
            'zone must be a number. Got {}'.format(type(self.zone))
        assert isinstance(self.vegStart, (float, int)), \
            'vegStart must be a number. Got {}'.format(type(self.vegStart))
        assert isinstance(self.vegEnd, (float, int)), \
            'vegEnd must be a number. Got {}'.format(type(self.vegEnd))
        assert isinstance(self.vegCover, float), \
            'vegCover must be a float. Got {}'.format(type(self.vegCover))
        assert isinstance(self.treeCoverage, float), \
            'treeCoverage must be a float. Got {}'.format(type(self.treeCoverage))
        assert isinstance(self.albVeg, float), \
            'albVeg must be a float. Got {}'.format(type(self.albVeg))
        assert isinstance(self.rurVegCover, float), \
            'rurVegCover must be a float. Got {}'.format(type(self.rurVegCover))
        assert isinstance(self.kRoad, float), \
            'kRoad must be a float. Got {}'.format(type(self.kRoad))
        assert isinstance(self.cRoad, float), \
            'cRoad must be a float. Got {}'.format(type(self.cRoad))
        assert isinstance(self.SchTraffic, list), \
            'SchTraffic must be a list. Got {}'.format(type(self.SchTraffic))
        assert len(self.SchTraffic) == 3, \
            'length of SchTraffic must be 3. Got {}'.format(len(self.SchTraffic))

    def set_input(self):
        """ Set inputs from .uwg input file if not already defined, the check if all
        the required input parameters are there.
        """

        # If a uwgParamFileName is set, then read inputs from .uwg file.
        # User-defined class properties will override the inputs from the .uwg file.
        if self.uwgParamFileName is not None:
            print("\nReading uwg file input.")
            self.read_input()
        else:
            print("\nNo .uwg file input.")

        self.check_required_inputs()

        # Modify zone to be used as python index
        self.zone = int(self.zone)-1

    def init_BEM_obj(self):
        """
        Define BEM for each DOE type (read the fraction)
        self.BEM                # list of BEMDef objects
        self.r_glaze            # Glazing ratio for total building stock
        self.SHGC               # SHGC addition for total building stock
        self.alb_wall           # albedo wall addition for total building stock
        """

        if not os.path.exists(self.readDOE_file_path):
            raise Exception("readDOE.pkl file: '{}' does not exist.".format(readDOE_file_path))

        readDOE_file = open(self.readDOE_file_path, 'rb')  # open pickle file in binary form
        refDOE = pickle.load(readDOE_file)
        refBEM = pickle.load(readDOE_file)
        refSchedule = pickle.load(readDOE_file)
        readDOE_file.close()

        # Define building energy models
        k = 0
        self.r_glaze_total = 0.             # Glazing ratio for total building stock
        self.SHGC_total = 0.                # SHGC addition for total building stock
        self.alb_wall_total = 0.            # albedo wall addition for total building stock
        h_floor = self.flr_h or 3.05  # average floor height

        total_urban_bld_area = math.pow(self.charLength, 2)*self.bldDensity * \
            self.bldHeight/h_floor  # total building floor area
        area_matrix = utilities.zeros(16, 3)

        self.BEM = []           # list of BEMDef objects
        self.Sch = []           # list of Schedule objects

        for i in range(16):    # 16 building types
            for j in range(3):  # 3 built eras
                if self.bld[i][j] > 0.:
                    # Add to BEM list
                    self.BEM.append(refBEM[i][j][self.zone])
                    self.BEM[k].frac = self.bld[i][j]
                    self.BEM[k].fl_area = self.bld[i][j] * total_urban_bld_area

                    # Overwrite with optional parameters if provided
                    if self.glzR:
                        self.BEM[k].building.glazingRatio = self.glzR
                    if self.albRoof:
                        self.BEM[k].roof.albedo = self.albRoof
                    if self.vegRoof:
                        self.BEM[k].roof.vegCoverage = self.vegRoof
                    if self.SHGC:
                        self.BEM[k].building.shgc = self.SHGC
                    if self.albWall:
                        self.BEM[k].wall.albedo = self.albWall
                    if self.flr_h:
                        self.BEM[k].building.floorHeight = self.flr_h

                    # Keep track of total urban r_glaze, SHGC, and alb_wall for UCM model
                    self.r_glaze_total += self.BEM[k].frac * self.BEM[k].building.glazingRatio
                    self.SHGC_total += self.BEM[k].frac * self.BEM[k].building.shgc
                    self.alb_wall_total += self.BEM[k].frac * self.BEM[k].wall.albedo
                    # Add to schedule list
                    self.Sch.append(refSchedule[i][j][self.zone])
                    k += 1

    def init_input_obj(self):
        """Section 4 - Create uwg objects from input parameters

            self.simTime            # simulation time parameter obj
            self.weather            # weather obj for simulation time period
            self.forcIP             # Forcing obj
            self.forc               # Empty forcing obj
            self.geoParam           # geographic parameters obj
            self.RSM                # Rural site & vertical diffusion model obj
            self.USM                # Urban site & vertical diffusion model obj
            self.UCM                # Urban canopy model obj
            self.UBL                # Urban boundary layer model

            self.road               # urban road element
            self.rural              # rural road element

            self.soilindex1         # soil index for urban rsoad depth
            self.soilindex2         # soil index for rural road depth
            self.Sch                # list of Schedule objects
        """

        climate_file_path = os.path.join(self.epwDir, self.epwFileName)

        self.simTime = SimParam(self.dtSim, self.dtWeather, self.Month,
                                self.Day, self.nDay)  # simulation time parametrs
        # weather file data for simulation time period
        self.weather = Weather(climate_file_path, self.simTime.timeInitial, self.simTime.timeFinal)
        self.forcIP = Forcing(self.weather.staTemp, self.weather)  # initialized Forcing class
        self.forc = Forcing()  # empty forcing class

        # Initialize geographic Param and Urban Boundary Layer Objects
        nightStart = 18.        # arbitrary values for begin/end hour for night setpoint
        nightEnd = 8.
        maxdx = 250.            # max dx (m)

        self.geoParam = Param(self.h_ubl1, self.h_ubl2, self.h_ref, self.h_temp, self.h_wind, self.c_circ,
                              self.maxDay, self.maxNight, self.latTree, self.latGrss, self.albVeg, self.vegStart, self.vegEnd,
                              nightStart, nightEnd, self.windMin, self.WGMAX, self.c_exch, maxdx, self.G, self.CP, self.VK, self.R,
                              self.RV, self.LV, math.pi, self.SIGMA, self.WATERDENS, self.LVTT, self.TT, self.ESTT, self.CL,
                              self.CPV, self.B, self.CM, self.COLBURN)

        self.UBL = UBLDef(
            'C', self.charLength, self.weather.staTemp[0], maxdx, self.geoParam.dayBLHeight, self.geoParam.nightBLHeight)

        # Defining road
        emis = 0.93
        asphalt = Material(self.kRoad, self.cRoad, 'asphalt')
        road_T_init = 293.
        road_horizontal = 1
        # fraction of surface vegetation coverage
        road_veg_coverage = min(self.vegCover/(1-self.bldDensity), 1.)

        # define road layers
        road_layer_num = int(math.ceil(self.d_road/0.05))
        # 0.5/0.05 ~ 10 x 1 matrix of 0.05 thickness
        thickness_vector = [0.05 for r in range(road_layer_num)]
        material_vector = [asphalt for r in range(road_layer_num)]

        self.road = Element(self.alb_road, emis, thickness_vector, material_vector, road_veg_coverage,
                            road_T_init, road_horizontal, name="urban_road")

        self.rural = copy.deepcopy(self.road)
        self.rural.vegCoverage = self.rurVegCover
        self.rural._name = "rural_road"

        # Reference site class (also include VDM)
        self.RSM = RSMDef(self.lat, self.lon, self.GMT, self.h_obs,
                          self.weather.staTemp[0], self.weather.staPres[0], self.geoParam, self.z_meso_dir_path)
        self.USM = RSMDef(self.lat, self.lon, self.GMT, self.bldHeight/10.,
                          self.weather.staTemp[0], self.weather.staPres[0], self.geoParam, self.z_meso_dir_path)

        T_init = self.weather.staTemp[0]
        H_init = self.weather.staHum[0]

        self.UCM = UCMDef(self.bldHeight, self.bldDensity, self.verToHor, self.treeCoverage, self.sensAnth, self.latAnth, T_init, H_init,
                          self.weather.staUmod[0], self.geoParam, self.r_glaze_total, self.SHGC_total, self.alb_wall_total, self.road)
        self.UCM.h_mix = self.h_mix

        # Define Road Element & buffer to match ground temperature depth
        roadMat, newthickness = procMat(self.road, self.MAXTHICKNESS, self.MINTHICKNESS)

        for i in range(self.nSoil):
            # if soil depth is greater then the thickness of the road
            # we add new slices of soil at max thickness until road is greater or equal

            is_soildepth_equal = self.is_near_zero(self.depth_soil[i][0] - sum(newthickness), 1e-15)

            if is_soildepth_equal or (self.depth_soil[i][0] > sum(newthickness)):
                while self.depth_soil[i][0] > sum(newthickness):
                    newthickness.append(self.MAXTHICKNESS)
                    roadMat.append(self.SOIL)
                self.soilindex1 = i
                break

        self.road = Element(self.road.albedo, self.road.emissivity, newthickness, roadMat,
                            self.road.vegCoverage, self.road.layerTemp[0], self.road.horizontal, self.road._name)

        # Define Rural Element
        ruralMat, newthickness = procMat(self.rural, self.MAXTHICKNESS, self.MINTHICKNESS)

        for i in range(self.nSoil):
            # if soil depth is greater then the thickness of the road
            # we add new slices of soil at max thickness until road is greater or equal

            is_soildepth_equal = self.is_near_zero(self.depth_soil[i][0] - sum(newthickness), 1e-15)

            if is_soildepth_equal or (self.depth_soil[i][0] > sum(newthickness)):
                while self.depth_soil[i][0] > sum(newthickness):
                    newthickness.append(self.MAXTHICKNESS)
                    ruralMat.append(self.SOIL)

                self.soilindex2 = i
                break

        self.rural = Element(self.rural.albedo, self.rural.emissivity, newthickness,
                             ruralMat, self.rural.vegCoverage, self.rural.layerTemp[0], self.rural.horizontal, self.rural._name)

    def hvac_autosize(self):
        """ Section 6 - HVAC Autosizing (unlimited cooling & heating) """

        for i in range(len(self.BEM)):
            if self.is_near_zero(self.autosize) == False:
                self.BEM[i].building.coolCap = 9999.
                self.BEM[i].building.heatCap = 9999.

    def simulate(self):
        """ Section 7 - uwg main section

            self.N                  # Total hours in simulation
            self.ph                 # per hour
            self.dayType            # 3=Sun, 2=Sat, 1=Weekday
            self.ceil_time_step     # simulation timestep (dt) fitted to weather file timestep

            # Output of object instance vector
            self.WeatherData        # Nx1 vector of forc instance
            self.UCMData            # Nx1 vector of UCM instance
            self.UBLData            # Nx1 vector of UBL instance
            self.RSMData            # Nx1 vector of RSM instance
            self.USMData            # Nx1 vector of USM instance
        """

        self.N = int(self.simTime.days * 24)       # total number of hours in simulation
        n = 0                                      # weather time step counter
        self.ph = self.simTime.dt/3600.            # dt (simulation time step) in hours

        # Data dump variables
        time = range(self.N)

        self.WeatherData = [None for x in range(self.N)]
        self.UCMData = [None for x in range(self.N)]
        self.UBLData = [None for x in range(self.N)]
        self.RSMData = [None for x in range(self.N)]
        self.USMData = [None for x in range(self.N)]

        print('\nSimulating new temperature and humidity values for {} days from {}/{}.\n'.format(
            int(self.nDay), int(self.Month), int(self.Day)))
        self.logger.info("Start simulation")

        for it in range(1, self.simTime.nt, 1):  # for every simulation time-step (i.e 5 min) defined by uwg
            # Update water temperature (estimated)
            if self.is_near_zero(self.nSoil):
                # for BUBBLE/CAPITOUL/Singapore only
                self.forc.deepTemp = sum(self.forcIP.temp)/float(len(self.forcIP.temp))
                self.forc.waterTemp = sum(
                    self.forcIP.temp)/float(len(self.forcIP.temp)) - 10.      # for BUBBLE/CAPITOUL/Singapore only
            else:
                # soil temperature by depth, by month
                self.forc.deepTemp = self.Tsoil[self.soilindex1][self.simTime.month-1]
                self.forc.waterTemp = self.Tsoil[2][self.simTime.month-1]

            # There's probably a better way to update the weather...
            self.simTime.UpdateDate()

            self.logger.info("\n{0} m={1}, d={2}, h={3}, s={4}".format(
                __name__, self.simTime.month, self.simTime.day, self.simTime.secDay/3600., self.simTime.secDay))

            # simulation time increment raised to weather time step
            self.ceil_time_step = int(math.ceil(it * self.ph))-1
            # minus one to be consistent with forcIP list index
            # Updating forcing instance
            # horizontal Infrared Radiation Intensity (W m-2)
            self.forc.infra = self.forcIP.infra[self.ceil_time_step]
            # wind speed (m s-1)
            self.forc.wind = max(self.forcIP.wind[self.ceil_time_step], self.geoParam.windMin)
            self.forc.uDir = self.forcIP.uDir[self.ceil_time_step]          # wind direction
            # specific humidty (kg kg-1)
            self.forc.hum = self.forcIP.hum[self.ceil_time_step]
            self.forc.pres = self.forcIP.pres[self.ceil_time_step]          # Pressure (Pa)
            self.forc.temp = self.forcIP.temp[self.ceil_time_step]          # air temperature (C)
            self.forc.rHum = self.forcIP.rHum[self.ceil_time_step]          # Relative humidity (%)
            self.forc.prec = self.forcIP.prec[self.ceil_time_step]          # Precipitation (mm h-1)
            # horizontal solar diffuse radiation (W m-2)
            self.forc.dif = self.forcIP.dif[self.ceil_time_step]
            # normal solar direct radiation (W m-2)
            self.forc.dir = self.forcIP.dir[self.ceil_time_step]
            # Canyon humidity (absolute) same as rural
            self.UCM.canHum = copy.copy(self.forc.hum)

            # Update solar flux
            self.solar = SolarCalcs(self.UCM, self.BEM, self.simTime,
                                    self.RSM, self.forc, self.geoParam, self.rural)
            self.rural, self.UCM, self.BEM = self.solar.solarcalcs()

            # Update building & traffic schedule
            # Assign day type (1 = weekday, 2 = sat, 3 = sun/other)
            if self.is_near_zero(self.simTime.julian % 7):
                self.dayType = 3                                        # Sunday
            elif self.is_near_zero(self.simTime.julian % 7 - 6.):
                self.dayType = 2                                        # Saturday
            else:
                self.dayType = 1                                        # Weekday

            # Update anthropogenic heat load for each hour (building & UCM)
            self.UCM.sensAnthrop = self.sensAnth * (self.SchTraffic[self.dayType-1][self.simTime.hourDay])

            # Update the energy components for building types defined in initialize.uwg
            for i in range(len(self.BEM)):
                # Set temperature
                self.BEM[i].building.coolSetpointDay = self.Sch[i].Cool[self.dayType -
                                                                        1][self.simTime.hourDay] + 273.15  # add from temperature schedule for cooling
                self.BEM[i].building.coolSetpointNight = self.BEM[i].building.coolSetpointDay
                self.BEM[i].building.heatSetpointDay = self.Sch[i].Heat[self.dayType -
                                                                        1][self.simTime.hourDay] + 273.15  # add from temperature schedule for heating
                self.BEM[i].building.heatSetpointNight = self.BEM[i].building.heatSetpointDay

                # Internal Heat Load Schedule (W/m^2 of floor area for Q)
                self.BEM[i].Elec = self.Sch[i].Qelec * self.Sch[i].Elec[self.dayType -
                                                                        1][self.simTime.hourDay]      # Qelec x elec fraction for day
                self.BEM[i].Light = self.Sch[i].Qlight * self.Sch[i].Light[self.dayType -
                                                                           1][self.simTime.hourDay]    # Qlight x light fraction for day
                self.BEM[i].Nocc = self.Sch[i].Nocc * self.Sch[i].Occ[self.dayType -
                                                                      1][self.simTime.hourDay]        # Number of occupants x occ fraction for day
                # Sensible Q occupant * fraction occupant sensible Q * number of occupants
                self.BEM[i].Qocc = self.sensOcc * (1 - self.LatFOcc) * self.BEM[i].Nocc

                # SWH and ventilation schedule
                self.BEM[i].SWH = self.Sch[i].Vswh * self.Sch[i].SWH[self.dayType -
                                                                     1][self.simTime.hourDay]          # litres per hour x SWH fraction for day
                # m^3/s/m^2 of floor
                self.BEM[i].building.vent = self.Sch[i].Vent
                self.BEM[i].Gas = self.Sch[i].Qgas * self.Sch[i].Gas[self.dayType -
                                                                     1][self.simTime.hourDay]          # Gas Equip Schedule, per m^2 of floor

                # This is quite messy, should update
                # Update internal heat and corresponding fractional loads
                intHeat = self.BEM[i].Light + self.BEM[i].Elec + self.BEM[i].Qocc
                # W/m2 from light, electricity, occupants
                self.BEM[i].building.intHeatDay = intHeat
                self.BEM[i].building.intHeatNight = intHeat
                # fraction of radiant heat from light and equipment of whole internal heat
                self.BEM[i].building.intHeatFRad = (
                    self.RadFLight * self.BEM[i].Light + self.RadFEquip * self.BEM[i].Elec) / intHeat
                # fraction of latent heat (from occupants) of whole internal heat
                self.BEM[i].building.intHeatFLat = self.LatFOcc * \
                    self.sensOcc * self.BEM[i].Nocc/intHeat

                # Update envelope temperature layers
                self.BEM[i].T_wallex = self.BEM[i].wall.layerTemp[0]
                self.BEM[i].T_wallin = self.BEM[i].wall.layerTemp[-1]
                self.BEM[i].T_roofex = self.BEM[i].roof.layerTemp[0]
                self.BEM[i].T_roofin = self.BEM[i].roof.layerTemp[-1]

            # Update rural heat fluxes & update vertical diffusion model (VDM)
            self.rural.infra = self.forc.infra - self.rural.emissivity * self.SIGMA * \
                self.rural.layerTemp[0]**4.    # Infrared radiation from rural road

            self.rural.SurfFlux(self.forc, self.geoParam, self.simTime,
                                self.forc.hum, self.forc.temp, self.forc.wind, 2., 0.)
            self.RSM.VDM(self.forc, self.rural, self.geoParam, self.simTime)

            # Calculate urban heat fluxes, update UCM & UBL
            self.UCM, self.UBL, self.BEM = urbflux(
                self.UCM, self.UBL, self.BEM, self.forc, self.geoParam, self.simTime, self.RSM)
            self.UCM.UCModel(self.BEM, self.UBL.ublTemp, self.forc, self.geoParam)
            self.UBL.UBLModel(self.UCM, self.RSM, self.rural,
                              self.forc, self.geoParam, self.simTime)

            """
            # Experimental code to run diffusion model in the urban area
            # N.B Commented out in python uwg because computed wind speed in
            # urban VDM: y = =0.84*ln((2-x/20)/0.51) results in negative log
            # for building heights >= 40m.

            Uroad = copy.copy(self.UCM.road)
            Uroad.sens = copy.copy(self.UCM.sensHeat)
            Uforc = copy.copy(self.forc)
            Uforc.wind = copy.copy(self.UCM.canWind)
            Uforc.temp = copy.copy(self.UCM.canTemp)
            self.USM.VDM(Uforc,Uroad,self.geoParam,self.simTime)
            """

            self.logger.info("dbT = {}".format(self.UCM.canTemp-273.15))
            if n > 0:
                logging.info("dpT = {}".format(self.UCM.Tdp))
                logging.info("RH  = {}".format(self.UCM.canRHum))

            if self.is_near_zero(self.simTime.secDay % self.simTime.timePrint) and n < self.N:

                self.logger.info("{0} ----sim time step = {1}----\n\n".format(__name__, n))

                self.WeatherData[n] = copy.copy(self.forc)
                _Tdb, _w, self.UCM.canRHum, _h, self.UCM.Tdp, _v = psychrometrics(
                    self.UCM.canTemp, self.UCM.canHum, self.forc.pres)

                self.UBLData[n] = copy.copy(self.UBL)
                self.UCMData[n] = copy.copy(self.UCM)
                self.RSMData[n] = copy.copy(self.RSM)

                self.logger.info("dbT = {}".format(self.UCMData[n].canTemp-273.15))
                self.logger.info("dpT = {}".format(self.UCMData[n].Tdp))
                self.logger.info("RH  = {}".format(self.UCMData[n].canRHum))

                n += 1

    def write_epw(self):
        """ Section 8 - Writing new EPW file
        """
        epw_prec = self.epw_precision  # precision of epw file input

        for iJ in range(len(self.UCMData)):
            # [iJ+self.simTime.timeInitial-8] = increments along every weather timestep in epw
            # [6 to 21]                       = column data of epw
            self.epwinput[iJ+self.simTime.timeInitial-8][6] = "{0:.{1}f}".format(
                self.UCMData[iJ].canTemp - 273.15, epw_prec)  # dry bulb temperature  [?C]
            # dew point temperature [?C]
            self.epwinput[iJ+self.simTime.timeInitial -
                          8][7] = "{0:.{1}f}".format(self.UCMData[iJ].Tdp, epw_prec)
            # relative humidity     [%]
            self.epwinput[iJ+self.simTime.timeInitial -
                          8][8] = "{0:.{1}f}".format(self.UCMData[iJ].canRHum, epw_prec)
            self.epwinput[iJ+self.simTime.timeInitial-8][21] = "{0:.{1}f}".format(
                self.WeatherData[iJ].wind, epw_prec)        # wind speed [m/s]

        # Writing new EPW file
        epw_new_id = open(self.newPathName, "w")

        for i in range(8):
            new_epw_line = '{}\n'.format(reduce(lambda x, y: x+","+y, self._header[i]))
            epw_new_id.write(new_epw_line)

        for i in range(len(self.epwinput)):
            printme = ""
            for ei in range(34):
                printme += "{}".format(self.epwinput[i][ei]) + ','
            printme = printme + "{}".format(self.epwinput[i][ei])
            new_epw_line = "{0}\n".format(printme)
            epw_new_id.write(new_epw_line)

        epw_new_id.close()

        print("New climate file '{}' is generated at {}.".format(
            self.destinationFileName, self.destinationDir))

    def run(self):

        # run main class methods
        self.read_epw()
        self.set_input()
        self.init_BEM_obj()
        self.init_input_obj()
        self.hvac_autosize()
        self.simulate()
        self.write_epw()


def procMat(materials, max_thickness, min_thickness):
    """ Processes material layer so that a material with single
    layer thickness is divided into two and material layer that is too
    thick is subdivided
    """
    newmat = []
    newthickness = []
    k = materials.layerThermalCond
    Vhc = materials.layerVolHeat

    if len(materials.layerThickness) > 1:

        for j in range(len(materials.layerThickness)):
            # Break up each layer that's more than max thickness (0.05m)
            if materials.layerThickness[j] > max_thickness:
                nlayers = math.ceil(materials.layerThickness[j]/float(max_thickness))
                for i in range(int(nlayers)):
                    newmat.append(Material(k[j], Vhc[j], name=materials._name))
                    newthickness.append(materials.layerThickness[j]/float(nlayers))
            # Material that's less then min_thickness is not added.
            elif materials.layerThickness[j] < min_thickness:
                print("WARNING: Material '{}' layer found too thin (<{:.2f}cm), ignored.").format(
                    materials._name, min_thickness*100)
            else:
                newmat.append(Material(k[j], Vhc[j], name=materials._name))
                newthickness.append(materials.layerThickness[j])

    else:

        # Divide single layer into two (uwg assumes at least 2 layers)
        if materials.layerThickness[0] > max_thickness:
            nlayers = math.ceil(materials.layerThickness[0]/float(max_thickness))
            for i in range(int(nlayers)):
                newmat.append(Material(k[0], Vhc[0], name=materials._name))
                newthickness.append(materials.layerThickness[0]/float(nlayers))
        # Material should be at least 1cm thick, so if we're here,
        # should give warning and stop. Only warning given for now.
        elif materials.layerThickness[0] < min_thickness*2:
            newthickness = [min_thickness/2., min_thickness/2.]
            newmat = [Material(k[0], Vhc[0], name=materials._name),
                      Material(k[0], Vhc[0], name=materials._name)]
            print("WARNING: a thin (<2cm) single material '{}' layer found. May cause error.".format(
                materials._name))
        else:
            newthickness = [materials.layerThickness[0]/2., materials.layerThickness[0]/2.]
            newmat = [Material(k[0], Vhc[0], name=materials._name),
                      Material(k[0], Vhc[0], name=materials._name)]
    return newmat, newthickness
