"""Urban Weather Generator (UWG) Version 4.2

Original Author: B. Bueno[1]
Edited by A. Nakano & Lingfu Zhang
Modified by Joseph Yang (joeyang@mit.edu) - May, 2016
Translated to Python by Saeran Vasanthakumar - February, 2018

Note:
    [1] Bueno, Bruno; Norford, Leslie; Hidalgo, Julia; Pigeon, Gregoire (2012a).
    The urban weather generator, Journal of Building Performance Simulation. 6:4,269-281.
    doi: 10.1080/19401493.2012.718797
"""

from __future__ import division, print_function
from functools import reduce

try:
    range = xrange
except NameError:
    pass

try:
    str = basestring
except NameError:
    pass

import os
import math
import copy

try:
    import cPickle as pickle
except ImportError:
    import pickle

from .simparam import SimParam
from .weather import Weather
from .material import Material
from .element import Element
from .param import Param
from .UCMDef import UCMDef
from .forcing import Forcing
from .UBLDef import UBLDef
from .RSMDef import RSMDef
from .solarcalcs import SolarCalcs
from .psychrometrics import psychrometrics
from .urbflux import urbflux
from .schdef import SchDef
from .BEMDef import BEMDef
from . import utilities
from .utilities import int_in_range, float_in_range, int_positive, float_positive
from .utilities import REF_BLDTYPE, REF_BUILTERA, REF_ZONETYPE, REF_BLDTYPE_SET, \
    REF_BUILTERA_SET, REF_ZONETYPE_SET


class UWG(object):
    """Morph a rural EPW file to urban conditions based on defined urban parameters.

    Args:
        epw_path: Text string for full path of the rural .epw file that will be
            morphed. If set to None, other input parameters can be assigned but the UWG
            model cannot be generated from the inputs, which is useful in cases where a
            UWG model needs to be serialized but the file path structure is not known.
            (Default: None).
        new_epw_dir: Optional text string for the destination directory into which the
            morphed .epw file is written. If None the morphed file will be written into
            the same directory as the rural .epw file. (Default: None).
        new_epw_name: Optional text string for the destination file name of the morphed
            .epw file. If None the morphed file will append '_UWG' to the original file
            name. (Default: None).

    Properties:
        * epw_path -- Full path of the rural .epw file that will be morphed.
        * new_epw_path -- Full path of the file name of the morphed .epw file.
        * refBEM -- Reference BEMDef matrix defined by built type, era, and zone.
        * refSchedule -- Reference SchDef matrix defined by built type, era, and zone.
        * month -- Number (1-12) representing simulation start month.
        * day -- Number (1-31) representing simulation start day.
        * nday -- Number of days to simulate.
        * dtsim -- Simlation time step in seconds.
        * dtweather -- Number for weather data time-step in seconds.
        * autosize -- Boolean to set HVAC autosize.
        * sensocc -- Sensible heat from occupant [W].
        * latfocc -- Latent heat fraction from occupant.
        * radfocc -- Radiant heat fraction from occupant.
        * radfequip -- Radiant heat fraction from equipment.
        * radflight -- Radiant heat fraction from electric light.
        * h_ubl1 -- Daytime urban boundary layer height in meters.
        * h_ubl2 -- Nighttime urban boundary layer height in meters.
        * h_ref -- Inversion height in meters.
        * h_temp -- Temperature height in meters.
        * h_wind -- Wind height in meters.
        * c_circ -- Wind scaling coefficient.
        * c_exch -- Exchange velocity coefficient.
        * maxday -- Maximum heat flux threshold for daytime conditions [W/m2].
        * maxnight -- Maximum heat flux threshold for nighttime conditions [W/m2].
        * windmin -- Minimum wind speed in m/s.
        * h_obs -- Rural average obstacle height in meters.
        * bldheight -- Urban building height in meters.
        * h_mix -- Fraction of HVAC waste heat released to street canyon.
        * blddensity -- Building footprint density as fraction of urban area.
        * vertohor -- Vertical-to-horizontal urban area ratio.
        * charlength -- Urban characteristic length in meters.
        * albroad -- Urban road albedo.
        * droad -- Thickness of urban road pavement thickness in meters.
        * sensanth -- Street level anthropogenic sensible heat [W/m2].
        * zone -- Index representing an ASHRAE climate zone.
        * grasscover -- Fraction of urban ground covered in grass only.
        * treecover -- Fraction of urban ground covered in trees.
        * vegstart -- Month in which vegetation starts to evapotranspire.
        * vegend -- Month in which vegetation stops evapotranspiration.
        * albveg -- Vegetation albedo.
        * rurvegcover -- Fraction of rural ground covered by vegetation.
        * latgrss -- Fraction of latent heat absorbed by urban grass.
        * lattree -- Fraction latent heat absorbed by urban trees.
        * schtraffic -- Schedule of fractional anthropogenic heat load.
        * kroad -- Road pavement conductivity [W/m-K].
        * croad -- Road pavement volumetric heat capacity [J/m^3K].
        * bld -- Matrix of numbers representing fraction of urban building stock.
        * albroof -- Average building roof albedo.
        * vegroof -- Fraction of roof covered in grass/shrubs.
        * glzr -- Building glazing ratio.
        * albwall -- Building albedo.
        * shgc -- Building glazing Solar Heat Gain Coefficient (SHGC).
        * flr_h -- Building floor height in meters.
        * ref_bem_vector -- List of custom BEMDef objects to override the refBEM.
        * ref_sch_vector -- List of custom SchDef objects to override the refSchedule.
    """

    # Definitions for constants / other parameters
    MINTHICKNESS = 0.01    # Minimum layer thickness (to prevent crashing) (m)
    MAXTHICKNESS = 0.05    # Maximum layer thickness (m)
    # http://web.mit.edu/parmstr/Public/NRCan/nrcc29118.pdf (Figly & Snodgrass)
    SOILTCOND = 1
    # http://www.europment.org/library/2013/venice/bypaper/MFHEEF/MFHEEF-21.pdf
    # (average taken from Table 1)
    SOILVOLHEAT = 2e6
    # Soil material used for soil-depth padding
    SOIL = Material(SOILTCOND, SOILVOLHEAT, name='soil')
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
    COLBURN = math.pow((0.713 / 0.621), (2 / 3.))
    # Site-specific parameters
    WGMAX = 0.005  # maximum film water depth on horizontal surfaces (m)

    # UWG object constants
    PARAMETER_LIST = ('month', 'day', 'nday', 'dtsim', 'dtweather', 'autosize',
                      'sensocc', 'latfocc', 'radfocc', 'radfequip', 'radflight',
                      'h_ubl1', 'h_ubl2', 'h_ref', 'h_temp', 'h_wind', 'c_circ',
                      'c_exch', 'maxday', 'maxnight', 'windmin', 'h_obs', 'bldheight',
                      'h_mix', 'blddensity', 'vertohor', 'charlength', 'albroad',
                      'droad', 'sensanth', 'zone', 'grasscover', 'treecover',
                      'vegstart', 'vegend', 'albveg', 'rurvegcover', 'latgrss',
                      'lattree', 'schtraffic', 'kroad', 'croad', 'bld', 'shgc',
                      'albroof', 'glzr', 'vegroof', 'albwall', 'flr_h')
    OPTIONAL_PARAMETER_SET = {'shgc', 'albroof',
                              'glzr', 'vegroof', 'albwall', 'flr_h'}
    DEFAULT_BLD = (('largeoffice', 'pst80', 0.4),
                   ('midriseapartment', 'pst80', 0.6))
    DEFAULT_SCHTRAFFIC = (
        (0.2, 0.2, 0.2, 0.2, 0.2, 0.4, 0.7, 0.9, 0.9, 0.6, 0.6, 0.6, 0.6, 0.6, 0.7, 0.8,
         0.9, 0.9, 0.8, 0.8, 0.7, 0.3, 0.2, 0.2),  # Weekday
        (0.2, 0.2, 0.2, 0.2, 0.2, 0.3, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.6, 0.7,
         0.7, 0.7, 0.7, 0.5, 0.4, 0.3, 0.2, 0.2),  # Saturday
        (0.2, 0.2, 0.2, 0.2, 0.2, 0.3, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
         0.4, 0.4, 0.4, 0.4, 0.3, 0.3, 0.2, 0.2))  # Sunday

    # Constant file paths
    CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
    Z_MESO_PATH = os.path.join(CURRENT_PATH, 'refdata', 'z_meso.txt')
    REFDOE_PATH = os.path.join(CURRENT_PATH, 'refdata', 'readDOE.pkl')

    def __init__(self, epw_path=None, new_epw_dir=None, new_epw_name=None):

        self.epw_path = epw_path
        self._new_epw_dir, self._new_epw_name = new_epw_dir, new_epw_name
        self._new_epw_path = None

        # set defaults for reference data
        self._refBEM = None
        self._refSchedule = None
        self._ref_bem_vector = None
        self._ref_sch_vector = None

        # Parameters for UWG computation
        self.epw_precision = 1

        # Set default for optional parameters as None, if user doesn't modify
        # use DOE reference equivalents.
        self._shgc = None
        self._flr_h = None
        self._albroof = None
        self._albwall = None
        self._glzr = None
        self._vegroof = None

        # Empty latanth not used.
        self.latanth = None

    @classmethod
    def from_param_file(cls, param_path, epw_path=None, new_epw_dir=None,
                        new_epw_name=None):
        """Create a UWG object from the .uwg parameter file.

        Note: this method of initializing the UWG object doesn't permit adding custom
        reference data.

        Args:
            param_path: Optional text string for full path of the the .uwg parameter file
                path.
            epw_path: Text string for full path of the rural .epw file that will be
                morphed. If set to None, other input parameters can be assigned but the
                UWG model cannot be generated from the inputs, which is useful in cases
                where a UWG model needs to be serialized but the file path structure is
                not known. (Default: None).
            new_epw_dir: Optional text string destination directory for the morphed
                .epw file. If None the morphed file will be written into the same
                directory as the rural .epw file. (Default: None).
            new_epw_name: Optional destination file name for the morphed .epw file.
                If None the morphed file will append '_UWG' to the original file name.
                (Default: None).
        """

        assert os.path.exists(param_path), 'Parameter file "{}" does not ' \
            'exist.'.format(param_path)

        assert param_path.endswith('.uwg'), 'Parameter file must be a ".uwg" ' \
            'filetype. Got: {}.'.format(param_path)

        model = UWG(epw_path, new_epw_dir, new_epw_name)
        model._read_input(param_path)
        return model

    @classmethod
    def from_param_args(cls, bldheight, blddensity, vertohor, grasscover,
                        treecover, zone, month=1, day=1, nday=31, dtsim=300,
                        dtweather=3600, bld=DEFAULT_BLD, autosize=False,
                        h_mix=1, sensocc=100, latfocc=0.3, radfocc=0.2, radfequip=0.5,
                        radflight=0.7, charlength=1000, albroad=0.1,
                        droad=0.5, kroad=1, croad=1600000, rurvegcover=0.9, vegstart=4, vegend=10,
                        albveg=0.25, latgrss=0.4, lattree=0.6, sensanth=20,
                        schtraffic=DEFAULT_SCHTRAFFIC, h_ubl1=1000, h_ubl2=80, h_ref=150,
                        h_temp=2, h_wind=10, c_circ=1.2, c_exch=1, maxday=150,
                        maxnight=20, windmin=1, h_obs=0.1, epw_path=None,
                        new_epw_dir=None, new_epw_name=None, ref_bem_vector=None,
                        ref_sch_vector=None):
        """Create an UWG object based on default method arguments.

        The default parameters are set from example parameters defined by Bueno et al.
        (2012a)[1] for Singapore. The original file can be accessed here:
        https://github.com/ladybug-tools/uwg/blob/master/resources/initialize_singapore.uwg.

        Note:
            [1] Bueno, Bruno; Norford, Leslie; Hidalgo, Julia; Pigeon, Gregoire (2012a).
            The urban weather generator, Journal of Building Performance Simulation. 6:4,
            269-281. doi: 10.1080/19401493.2012.718797
        """

        model = UWG(epw_path, new_epw_dir, new_epw_name)

        model.bldheight = bldheight
        model.blddensity = blddensity
        model.vertohor = vertohor
        model.zone = zone
        model.month = month
        model.day = day
        model.nday = nday
        model.dtsim = dtsim
        model.dtweather = dtweather
        model.autosize = autosize
        model.h_mix = h_mix
        model.sensocc = sensocc
        model.latfocc = latfocc
        model.radfocc = radfocc
        model.radfequip = radfequip
        model.radflight = radflight
        model.bld = bld
        model.charlength = charlength
        model.albroad = albroad
        model.droad = droad
        model.sensanth = sensanth
        model.kroad = kroad
        model.croad = croad
        model.treecover = treecover
        model.grasscover = grasscover
        model.vegstart = vegstart
        model.vegend = vegend
        model.albveg = albveg
        model.rurvegcover = rurvegcover
        model.latgrss = latgrss
        model.lattree = lattree
        model.schtraffic = schtraffic
        model.h_ubl1 = h_ubl1
        model.h_ubl2 = h_ubl2
        model.h_ref = h_ref
        model.h_temp = h_temp
        model.h_wind = h_wind
        model.c_circ = c_circ
        model.c_exch = c_exch
        model.maxday = maxday
        model.maxnight = maxnight
        model.windmin = windmin
        model.h_obs = h_obs

        # check and add reference data
        refcheck = int(ref_sch_vector is not None) + \
            int(ref_bem_vector is not None)
        assert refcheck != 1, 'The ref_sch_vector and ref_bem_vector must both be ' \
            'defined in order to modify the UWG reference data. Only {} is ' \
            'defined.'.format(
                'ref_sch_vector' if ref_bem_vector is None else 'ref_bem_vector')

        if refcheck == 2:
            model.ref_bem_vector, model.ref_sch_vector = \
                model._check_reference_data(ref_bem_vector, ref_sch_vector)

        return model

    @classmethod
    def from_dict(cls, data, epw_path=None, new_epw_dir=None, new_epw_name=None):
        """Create an UWG object from a dictionary.

        Args:
            data: An UWG dictionary following the format below. Note that
                this example has been truncated for the sake of brevity. For
                the full list of required properties in the UWG, see the
                initialization docstrings.
            epw_path: Text string for full path of the rural .epw file that will be
                morphed. If set to None, other input parameters can be assigned but the
                UWG model cannot be generated from the inputs, which is useful in cases
                where a UWG model needs to be serialized but the file path structure is
                not known. (Default: None).
            new_epw_dir: Optional text string for the destination directory into which
                the morphed .epw file is written. If None the morphed file will be
                written into the same directory as the rural .epw file. (Default: None).
            new_epw_name: Optional text string for the destination file name of the
                morphed .epw file. If None the morphed file will append '_UWG' to the
                original file name. (Default: None).

        .. code-block:: python

            {
            "type": "UWG",
            "bldheight": 10,
            "blddensity": 0.5,
            "vertohor": 0.8,
            ...
            "h_obs": 0.1,
            "flr_h": 3.5,
            "shgc": None,
            # Optional vector of SchDef dictionary.
            "ref_sch_vector": [sch.to_dict()]
            # Optional vector of BEMDef dictionary.
            "ref_bem_vector": [bem.to_dict()]
            }
        """
        assert data['type'] == 'UWG', \
            'Expected UWG dictionary. Got {}.'.format(data['type'])

        model = UWG(epw_path, new_epw_dir, new_epw_name)

        # set UWG parameters
        for attr in cls.PARAMETER_LIST:
            setattr(model, attr, data[attr])

        # check and add reference data
        check_sch = 'ref_sch_vector' in data and data['ref_sch_vector'] is not None
        check_bem = 'ref_bem_vector' in data and data['ref_bem_vector'] is not None
        refcheck = int(check_sch) + int(check_bem)
        assert refcheck != 1, 'The ref_sch_vector and ref_bem_vector ' \
            'properties must both be defined in order to modify the UWG reference ' \
            'data. Only {} is defined.'.format(
                'ref_sch_vector' if 'ref_sch_vector' in data else 'ref_bem_vector')

        if refcheck == 2:
            ref_sch_vector = [SchDef.from_dict(schdict)
                              for schdict in data['ref_sch_vector']]
            ref_bem_vector = [BEMDef.from_dict(bemdict)
                              for bemdict in data['ref_bem_vector']]

        if refcheck == 2:
            model.ref_bem_vector, model.ref_sch_vector = \
                model._check_reference_data(ref_bem_vector, ref_sch_vector)

        return model

    @property
    def epw_path(self):
        """Get full path to rural .epw file to be morphed."""
        return self._epw_path

    @epw_path.setter
    def epw_path(self, value):
        if value is not None:
            assert os.path.exists(
                value), 'File: "{}" does not exist.'.format(value)
        self._epw_path = value

    @property
    def new_epw_path(self):
        """Get text string for new epw filepath."""
        if self._new_epw_path is None:
            if self._new_epw_dir is None or self._new_epw_name is None:
                epw_dir, epw_name = os.path.split(self.epw_path)
                if self._new_epw_dir is None:
                    self._new_epw_dir = epw_dir
                if self._new_epw_name is None:
                    self._new_epw_name = epw_name.strip('.epw') + '_UWG.epw'
            self._new_epw_path = os.path.join(
                self._new_epw_dir, self._new_epw_name)
        return self._new_epw_path

    @property
    def refBEM(self):
        """Get matrix of DOE reference BEMDefs defined by built type, era, and zone."""
        if self._refBEM is None:
            self._refBEM, self._refSchedule = UWG.load_refDOE()
        return self._refBEM

    @property
    def refSchedule(self):
        """Get matrix of DOE reference SchDefs defined by built type, era, and zone."""
        if self._refSchedule is None:
            self._refBEM, self._refSchedule = UWG.load_refDOE()
        return self._refSchedule

    @property
    def month(self):
        """Get or set number (1-12) as simulation start month."""
        return self._month

    @month.setter
    def month(self, value):
        self._month = int_in_range(value, 1, 12, 'month')

    @property
    def day(self):
        """Get or set number (1-31) as simulation start day."""
        return self._day

    @day.setter
    def day(self, value):
        self._day = int_in_range(value, 1, 31, 'day')

    @property
    def nday(self):
        """Get or set number of days to simulate."""
        return self._nday

    @nday.setter
    def nday(self, value):
        self._nday = int_positive(value, 'nday')

    @property
    def dtsim(self):
        """Get or set simulation timestep in seconds."""
        return self._dtsim

    @dtsim.setter
    def dtsim(self, value):
        self._dtsim = int_positive(value, 'dtsim')

    @property
    def dtweather(self):
        """Get or set weather data timestep in seconds."""
        return self._dtweather

    @dtweather.setter
    def dtweather(self, value):
        self._dtweather = float_positive(value, 'dtweather')

    @property
    def autosize(self):
        """Get or set boolean to autosize HVAC system."""
        return self._autosize

    @autosize.setter
    def autosize(self, value):
        assert isinstance(value, (bool, int, float)), 'Input autosize must be a ' \
            'boolean value or a number that can be cast as or boolean value. ' \
            'Got: {}.'.format(value)
        self._autosize = bool(value)

    @property
    def sensocc(self):
        """Get or set sensible heat from occupant [W]."""
        return self._sensocc

    @sensocc.setter
    def sensocc(self, value):
        self._sensocc = float_positive(value, 'sensocc')

    @property
    def latfocc(self):
        """Get or set latent heat fraction from occupant."""
        return self._latfocc

    @latfocc.setter
    def latfocc(self, value):
        self._latfocc = float_in_range(value, 0, 1, 'latfocc')

    @property
    def radfocc(self):
        """Get or set radiant heat fraction from occupant."""
        return self._radfocc

    @radfocc.setter
    def radfocc(self, value):
        self._radfocc = float_in_range(value, 0, 1, 'radfocc')

    @property
    def radfequip(self):
        """Get or set radiant heat fraction from equipment."""
        return self._radfequip

    @radfequip.setter
    def radfequip(self, value):
        self._radfequip = float_in_range(value, 0, 1, 'radfequip')

    @property
    def radflight(self):
        """Get or set radiant heat fraction from electric light."""
        return self._radflight

    @radflight.setter
    def radflight(self, value):
        self._radflight = float_in_range(value, 0, 1, 'radflight')

    @property
    def h_ubl1(self):
        """Get or set daytime urban boundary layer height in meters."""
        return self._h_ubl1

    @h_ubl1.setter
    def h_ubl1(self, value):
        self._h_ubl1 = float_positive(value, 'h_ubl1')

    @property
    def h_ubl2(self):
        """Get or set nighttime urban boundary layer height in meters."""
        return self._h_ubl2

    @h_ubl2.setter
    def h_ubl2(self, value):
        self._h_ubl2 = float_positive(value, 'h_ubl2')

    @property
    def h_ref(self):
        """Get or set inversion height in meters."""
        return self._h_ref

    @h_ref.setter
    def h_ref(self, value):
        self._h_ref = float_positive(value, 'h_ref')

    @property
    def h_temp(self):
        """Get or set temperature height in meters."""
        return self._h_temp

    @h_temp.setter
    def h_temp(self, value):
        self._h_temp = float_positive(value, 'h_temp')

    @property
    def h_wind(self):
        """Get or set wind height in meters."""
        return self._h_wind

    @h_wind.setter
    def h_wind(self, value):
        self._h_wind = float_positive(value, 'h_wind')

    @property
    def c_circ(self):
        """Get or set wind scaling coefficient."""
        return self._c_circ

    @c_circ.setter
    def c_circ(self, value):
        self._c_circ = float_positive(value, 'c_circ')

    @property
    def c_exch(self):
        """Get or set exchange velocity coefficient."""
        return self._c_exch

    @c_exch.setter
    def c_exch(self, value):
        self._c_exch = float_positive(value, 'c_exch')

    @property
    def maxday(self):
        """Get or set maximum heat flux threshold for daytime conditions [W/m]."""
        return self._maxday

    @maxday.setter
    def maxday(self, value):
        self._maxday = float_positive(value, 'maxday')

    @property
    def maxnight(self):
        """Get or set maximum heat flux threshold for nighttime conditions [W/m2]."""
        return self._maxnight

    @maxnight.setter
    def maxnight(self, value):
        self._maxnight = float_positive(value, 'maxnight')

    @property
    def windmin(self):
        """Get or set minimum wind speed in m/s."""
        return self._windmin

    @windmin.setter
    def windmin(self, value):
        self._windmin = float_positive(value, 'windmin')

    @property
    def h_obs(self):
        """Get or set rural average obstacle height in meters."""
        return self._h_obs

    @h_obs.setter
    def h_obs(self, value):
        self._h_obs = float_positive(value, 'h_obs')

    @property
    def bldheight(self):
        """Get or set average urban building height in meters."""
        return self._bldheight

    @bldheight.setter
    def bldheight(self, value):
        self._bldheight = float_positive(value, 'bldheight')

    @property
    def h_mix(self):
        """Get or set fraction of building HVAC waste heat released to street canyon.

        It is assumed the rest of building HVAC waste heat is released from the roof.
        """
        return self._h_mix

    @h_mix.setter
    def h_mix(self, value):
        self._h_mix = float_in_range(value, 0, 1, 'h_mix')

    @property
    def vertohor(self):
        """Get or set vertical-to-horizontal urban area ratio.

        The vertical-to-horizontal urban area ratio is calculated by dividing the
        urban facade area by total urban area.
        """
        return self._vertohor

    @vertohor.setter
    def vertohor(self, value):
        self._vertohor = float_positive(value, 'vertohor')

    @property
    def charlength(self):
        """Get or set the urban characteristic length in meters.

        The characteristic length is the dimension of a square that encompasses the
        whole neighborhood.
        """
        return self._charlength

    @charlength.setter
    def charlength(self, value):
        self._charlength = float_positive(value, 'charlength')

    @property
    def albroad(self):
        """Get or set urban road albedo."""
        return self._albroad

    @albroad.setter
    def albroad(self, value):
        self._albroad = float_in_range(value, 0, 1, 'albroad')

    @property
    def droad(self):
        """Get or set thickness of urban road pavement thickness in meters."""
        return self._droad

    @droad.setter
    def droad(self, value):
        self._droad = float_positive(value, 'droad')

    @property
    def sensanth(self):
        """Get or set street level anthropogenic sensible heat [W/m2].

        Street level anthropogenic heat is non-building heat like heat emitted from cars,
        pedestrians, and street cooking.
        """
        return self._sensanth

    @sensanth.setter
    def sensanth(self, value):
        self._sensanth = float_positive(value, 'sensanth')

    @property
    def bld(self):
        """Get or set list of building types, eras, and fractions of urban building stock.

        This property consists of a list of tuples, each containing a string for the the
        built era, and a number between 0 and 1, inclusive, defining built stock
        fraction, i.e ('LargeOffice', 'New', 0.4). The fractions should sum to one.

        768 predefined models are built referencing 16 building types for 3 built eras
        and 16 climate zones according to models from the Department of Energy (DOE).
        Choose from the following text identifiers to reference a DOE building type:

        * 'fullservicerestaurant'
        * 'hospital'
        * 'largehotel'
        * 'largeoffice'
        * 'medoffice'
        * 'midriseapartment'
        * 'outpatient'
        * 'primaryschool'
        * 'quickservicerestaurant'
        * 'secondaryschool'
        * 'smallhotel'
        * 'smalloffice'
        * 'standaloneretail'
        * 'stripmall'
        * 'supermarket'
        * 'warehouse'

        Choose from the following built eras:

        * 'pre80'
        * 'pst80'
        * 'new'

        Custom building types can also be referenced in this property. For example, a
        built stock consisting of 40% post-1980's large office, 30% new midrise
        apartment, and 30% of a pre-1980s custom building type (defined by the user)
        is referenced as follows:

        .. code-block:: python

            bld = [('largeoffice', 'pst80', 0.4),
                   ('midriseapartment', 'new', 0.3),
                   ('custombuilding', 'pre80', 0.3)]
        """
        return self._bld

    @bld.setter
    def bld(self, value):

        assert isinstance(value, (list, tuple)), 'The bld property must be a list. ' \
            'Got {}.'.format(value)

        # check values and fraction sum.
        total_frac = 0.0
        for bld_row in value:
            assert len(bld_row) == 3, 'Each bld tuple must contain three ' \
                'items defining building type, era and fraction. Got ' \
                '{} values.'.format(len(bld_row))

            bldtype, builtera, frac = \
                bld_row[0], bld_row[1], bld_row[2]
            assert isinstance(bldtype, str), 'The first item in the ' \
                'bld tuple must be text defining the reference building ' \
                'type. Got: {}.'.format(bldtype)
            assert isinstance(builtera, str) and builtera.lower() in REF_BUILTERA_SET, \
                'The second item in the bld tuple must be text defining the built ' \
                'era as one of {}. Got: {}.'.format(
                    REF_BUILTERA, builtera.lower())
            assert 0.0 <= frac <= 1.0, 'The third item in the bld tuple ' \
                'must be a value between 0 and 1, inclusive, defining the ' \
                'fraction of total built stock. Got: {}.'.format(frac)
            total_frac += frac

        assert abs(total_frac - 1.0) < 1e-10, 'The sum of reference building ' \
            'fractions defined in bld must equal one. Got: {}.'.format(
                total_frac)

        self._bld = value

    @property
    def lattree(self):
        """Get or set fraction of latent heat absorbed by urban trees."""
        return self._lattree

    @lattree.setter
    def lattree(self, value):
        self._lattree = float_in_range(value, 0, 1, 'lattree')

    @property
    def latgrss(self):
        """Get or set fraction of latent heat absorbed by urban grass."""
        return self._latgrss

    @latgrss.setter
    def latgrss(self, value):
        self._latgrss = float_in_range(value, 0, 1, 'latgrss')

    @property
    def zone(self):
        """Get or set text representing an ASHRAE climate zone.

        This value is used to specify climate zone-specific construction, and
        HVAC parameters for the DOE reference building types. This will not effect
        the simulation if only custom reference buildings are used.

        Choose from the following:

        * '1A' - (i.e Miami)
        * '2A' - (i.e Houston)
        * '2B' - (i.e Phoenix)
        * '3A' - (i.e Atlanta)
        * '3B-CA' - (i.e Los Angeles)
        * '3B' - (i.e Las Vegas)
        * '3C' - (i.e San Francisco)
        * '4A' - (i.e Baltimore)
        * '4B' - (i.e Albuquerque)
        * '4C' - (i.e Seattle)
        * '5A' - (i.e Chicago)
        * '5B' - (i.e Boulder)
        * '6A' - (i.e Minneapolis)
        * '6B' - (i.e Helena)
        * '7' - (i.e Duluth)
        * '8' - (i.e Fairbanks)
        """
        return self._zone

    @zone.setter
    def zone(self, value):
        assert isinstance(value, str), \
            'zone must be a string. Got: {}.'.format(value)
        value = value.upper()
        assert value in REF_ZONETYPE_SET, 'zone must be on of {}. Got: {}.'.format(
            REF_ZONETYPE, value)
        self._zone = value

    @property
    def vegstart(self):
        """Get or set value from 1 to 12 for month at which vegetation starts to evapotranspire.

        This month corresponds to when the leaves of vegetation are assumed to be out.
        """
        return self._vegstart

    @vegstart.setter
    def vegstart(self, value):
        self._vegstart = int_in_range(value, 1, 12, 'vegstart')

    @property
    def vegend(self):
        """Get or set value from 1 to 12 for month at which vegetation stops evapotranspiration.

        This month corresponds to when the leaves of vegetation are assumed to fall.
        """
        return self._vegend

    @vegend.setter
    def vegend(self, value):
        self._vegend = int_in_range(value, 1, 12, 'vegend')

    @property
    def blddensity(self):
        """Get or set building footprint density as fraction of urban area.

        The sum of blddensity, grasscover and treecover must be less than or equal to 1.
        """
        return self._blddensity

    @blddensity.setter
    def blddensity(self, value):
        try:
            assert self.vegcover + value <= 1, 'The sum of the blddensity, treecover '\
                ' and grasscover ratios must be less than one. Got: {}, {} and {}, ' \
                'respectively.'.format(value, self.treecover, self.grasscover)
        except AttributeError:
            pass  # attributes haven't been set
        self._blddensity = float_in_range(value, 0, 1, 'blddensity')

    @property
    def treecover(self):
        """Get or set fraction of urban area covered in trees.

        The sum of blddensity, grasscover and treecover must be less than or equal to 1.
        """
        return self._treecover

    @treecover.setter
    def treecover(self, value):
        try:
            assert self.grasscover + self.blddensity + value <= 1, 'The sum of the ' \
                'blddensity, treecover and grasscover ratios must be less than one. ' \
                'Got: {}, {} and {}, respectively.'.format(
                    self.blddensity, value, self.grasscover)
        except AttributeError:
            pass  # attributes haven't been set
        self._treecover = float_in_range(value, 0, 1, 'treecover')

    @property
    def grasscover(self):
        """Get or set fraction of urban area covered exclusively in grass.

        This value does not including grass under trees. The sum of blddensity,
        grasscover and treecover must be less than or equal to 1.
        """
        return self._grasscover

    @grasscover.setter
    def grasscover(self, value):
        try:
            assert self.treecover + self.blddensity + value <= 1, 'The sum of the ' \
                'blddensity, treecover and grasscover ratios must be less than one. ' \
                'Got: {}, {} and {}, respectively.'.format(
                    self.blddensity, self.treecover, value)
        except AttributeError:
            pass  # attributes haven't been set
        self._grasscover = float_in_range(value, 0, 1, 'grasscover')

    @property
    def vegcover(self):
        """Get fraction of urban ground covered by trees and grass."""
        return self.treecover + self.grasscover

    @property
    def albveg(self):
        """Get or set vegetation albedo."""
        return self._albveg

    @albveg.setter
    def albveg(self, value):
        self._albveg = float_in_range(value, 0, 1, 'albveg')

    @property
    def rurvegcover(self):
        """Get or set fraction of rural ground covered by vegetation."""
        return self._rurvegcover

    @rurvegcover.setter
    def rurvegcover(self, value):
        self._rurvegcover = float_in_range(value, 0, 1, 'rurvegcover')

    @property
    def kroad(self):
        """Get or set road pavement conductivity [W/mK]."""
        return self._kroad

    @kroad.setter
    def kroad(self, value):
        self._kroad = float_positive(value, 'kroad')

    @property
    def croad(self):
        """Get or set road pavement volumetric heat capacity [J/m^3K]."""
        return self._croad

    @croad.setter
    def croad(self, value):
        self._croad = float_positive(value, 'croad')

    @property
    def schtraffic(self):
        """Get or set matrix for schedule of fractional anthropogenic heat load.

        This property consists of a 3 x 24 matrix. Each row corresponding to a schedule
        for a weekday, Saturday, and Sunday, and each column corresponds to an hour in
        the day, for example:

        .. code-block:: python

            # Weekday schedule
            wkday = [0.2, 0.2, 0.2, 0.2, 0.2, 0.4, 0.7, 0.9, 0.9, 0.6, 0.6, 0.6, 0.6,
                     0.6, 0.7, 0.8, 0.9, 0.9, 0.8, 0.8, 0.7, 0.3, 0.2, 0.2]
            # Saturday schedule
            satday = [0.2, 0.2, 0.2, 0.2, 0.2, 0.3, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                     0.5, 0.6, 0.7, 0.7, 0.7, 0.7, 0.5, 0.4, 0.3, 0.2, 0.2]
            # Sunday schedule
            sunday = [0.2, 0.2, 0.2, 0.2, 0.2, 0.3, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
                     0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.3, 0.3, 0.2, 0.2]

            schtraffic = [wkday, satday, sunday]
        """
        return self._schtraffic

    @schtraffic.setter
    def schtraffic(self, value):

        value = SchDef.check_week_validity(value, 'schtraffic')
        self._schtraffic = [[0 for c in range(24)] for r in range(3)]

        # Check column number and add value
        for i in range(3):
            for j in range(24):
                self._schtraffic[i][j] = float_in_range(value[i][j])

    @property
    def shgc(self):
        """Get or set average building glazing Solar Heat Gain Coefficient.

        If value is None, a unique shgc is set for each building from the refBEM.
        (Default: None).
        """
        return self._shgc

    @shgc.setter
    def shgc(self, value):
        if value is None:
            self._shgc = value
        else:
            self._shgc = float_in_range(value, 0, 1, 'shgc')

    @property
    def albroof(self):
        """Get or set average building roof albedo.

        If value is None, a unique albroof is set for each building from the refBEM.
        (Default: None).
        """
        return self._albroof

    @albroof.setter
    def albroof(self, value):
        if value is None:
            self._albroof = value
        else:
            self._albroof = float_in_range(value, 0, 1, 'albroof')

    @property
    def glzr(self):
        """Get or set average building glazing ratio.

        If value is None, a unique glzr is set for each building from the refBEM.
        (Default: None).
        """
        return self._glzr

    @glzr.setter
    def glzr(self, value):
        if value is None:
            self._glzr = value
        else:
            self._glzr = float_in_range(value, 0, 1, 'glzr')

    @property
    def vegroof(self):
        """Get or set fraction of roofs covered in grass/shrubs.

        If value is None, a unique vegroof is set for each building from the refBEM.
        (Default: None).
        """
        return self._vegroof

    @vegroof.setter
    def vegroof(self, value):
        if value is None:
            self._vegroof = value
        else:
            self._vegroof = float_in_range(value, 0, 1, 'vegroof')

    @property
    def albwall(self):
        """Get or set average building albedo.

        If value is None, a unique albwall is set for each building from the refBEM.
        (Default: None).
        """
        return self._albwall

    @albwall.setter
    def albwall(self, value):
        if value is None:
            self._albwall = value
        else:
            self._albwall = float_in_range(value, 0, 1, 'albwall')

    @property
    def flr_h(self):
        """Get or set average building floor height in meters.

        If value is None, a unique flr_h is set for each building from the refBEM.
        (Default: None).
        """
        return self._flr_h

    @flr_h.setter
    def flr_h(self, value):
        if value is None:
            self._flr_h = value
        else:
            self._flr_h = float_positive(value, 'flr_h')

    @property
    def ref_bem_vector(self):
        """Get list of custom BEMDef objects to add to refBEM.

        If value is None, all BEMDef objects are referenced from the DOE typologies
        defined by default in the refBEM matrix. (Default: None).
        """
        return self._ref_bem_vector

    @ref_bem_vector.setter
    def ref_bem_vector(self, value):
        if self._ref_bem_vector is None:
            self._ref_bem_vector = value
        else:
            raise Exception('Cannot reset ref_bem_vector value.')

    @property
    def ref_sch_vector(self):
        """Get list of custom SchDef objects to add to refSchedule.

        If value is None, all SchDef objects are referenced from the DOE typologies
        defined by default in the refSch matrix. (Default: None).
        """
        return self._ref_sch_vector

    @ref_sch_vector.setter
    def ref_sch_vector(self, value):
        if self._ref_sch_vector is None:
            self._ref_sch_vector = value
        else:
            raise Exception('Cannot reset ref_sch_vector value.')

    def to_dict(self, include_refDOE=False):
        """UWG dictionary representation.

        Args:
            add_refDOE: Optional boolean to include custom reference BEMDef and SchDef
                objects from the ref_bem_vector and ref_sch_vector attributes.
                (Default: False).
        """

        base = {'type': 'UWG'}

        # Add UWG parameters
        for attr in self.PARAMETER_LIST:
            base[attr] = getattr(self, attr)

        # Add reference data
        if include_refDOE and self.ref_bem_vector and self.ref_sch_vector:
            base['ref_sch_vector'] = [s.to_dict() for s in self.ref_sch_vector]
            base['ref_bem_vector'] = [b.to_dict() for b in self.ref_bem_vector]

        return base

    def generate(self):
        """Generate all UWG objects after input parameters are set."""

        if not self.epw_path:
            raise Exception('Cannot generate the UWG while epw_path is None.')
        if self.ref_bem_vector:
            self._customize_reference_data()
        self._read_epw()
        self._compute_BEM()
        self._compute_input()
        self._hvac_autosize()

    def simulate(self):
        """Simulate UWG object and produce urban canyon weather timeseries data.

        This function will set the following attributes in the UWG object:

        * N - Total number of hours in simulation
        * ph - Number of simulation time steps per hour
        * dayType - Number representing day type: Sunday, Saturday or Weekday
        * ceil_time_step - sim timestep fitted to weather file timestep
        * solar - SolarCalcs object for current timestep solar calculation
        * WeatherData - N x 1 output vector of forc object instance
        * UCMData - N x 1 output vector of UCM object instance
        * UBLData - N x 1 output vector of UBL object instance
        * RSMData - N x 1 output vector of RSM object instace
        * USMData - N x 1 output vector of USM object instance
        """

        # total number of hours in simulation
        self.N = int(self.simTime.days * 24)
        n = 0  # weather time step counter
        self.ph = self.simTime.dt / 3600.  # dt (simulation time step) in hours

        # Data dump variables
        self.WeatherData = [None for x in range(self.N)]
        self.UCMData = [None for x in range(self.N)]
        self.UBLData = [None for x in range(self.N)]
        self.RSMData = [None for x in range(self.N)]
        self.USMData = [None for x in range(self.N)]

        print('Simulating new temperature and humidity values for '
              '{} days from {}/{}.'.format(self.nday, self.month, self.day))

        # iterate through every simulation time-step (i.e 5 min) defined by UWG
        for it in range(1, self.simTime.nt, 1):
            # Update water temperature (estimated)
            if self.nSoil < 3:
                # for BUBBLE/CAPITOUL/Singapore only
                self.forc.deepTemp = sum(
                    self.forcIP.temp) / float(len(self.forcIP.temp))
                self.forc.waterTemp = \
                    sum(self.forcIP.temp) / float(len(self.forcIP.temp)) - 10.0
            else:
                # soil temperature by depth, by month
                self.forc.deepTemp = self.Tsoil[self._soilindex1][self.simTime.month - 1]
                self.forc.waterTemp = self.Tsoil[2][self.simTime.month - 1]

            # There's probably a better way to update the weather...
            self.simTime.update_date()

            # simulation time increment raised to weather time step
            self.ceil_time_step = int(math.ceil(it * self.ph)) - 1
            # minus one to be consistent with forcIP list index
            # Updating forcing instance
            # horizontal Infrared Radiation Intensity (W m-2)
            self.forc.infra = self.forcIP.infra[self.ceil_time_step]
            # wind speed (m s-1)
            self.forc.wind = max(self.forcIP.wind[self.ceil_time_step],
                                 self.geoParam.windMin)
            # wind direction
            self.forc.uDir = self.forcIP.uDir[self.ceil_time_step]
            # specific humidty (kg kg-1)
            self.forc.hum = self.forcIP.hum[self.ceil_time_step]
            # Pressure (Pa)
            self.forc.pres = self.forcIP.pres[self.ceil_time_step]
            # air temp (C)
            self.forc.temp = self.forcIP.temp[self.ceil_time_step]
            self.forc.rHum = self.forcIP.rHum[self.ceil_time_step]  # RH (%)
            # Precip (mm h-1)
            self.forc.prec = self.forcIP.prec[self.ceil_time_step]
            # horizontal solar diffuse radiation (W m-2)
            self.forc.dif = self.forcIP.dif[self.ceil_time_step]
            # normal solar direct radiation (W m-2)
            self.forc.dir = self.forcIP.dir[self.ceil_time_step]
            # Canyon humidity (absolute) same as rural
            self.UCM.canHum = copy.copy(self.forc.hum)

            # Update solar flux
            self.solar = SolarCalcs(
                self.UCM, self.BEM, self.simTime, self.RSM, self.forc, self.geoParam,
                self.rural)
            self.rural, self.UCM, self.BEM = self.solar.solarcalcs()

            # Update building & traffic schedule
            # Assign day type (1 = weekday, 2 = sat, 3 = sun/other)
            if utilities.is_near_zero(self.simTime.julian % 7, 1e-10):
                self.dayType = 3  # Sunday
            elif utilities.is_near_zero(self.simTime.julian % 7 - 6., 1e-10):
                self.dayType = 2  # Saturday
            else:
                self.dayType = 1  # Weekday

            # Update anthropogenic heat load for each hour (building & UCM)
            self.UCM.sensAnthrop = \
                self.sensanth * \
                (self.schtraffic[self.dayType - 1][self.simTime.hourDay])

            # Update the energy components for building types defined in initialize.UWG
            for i in range(len(self.BEM)):

                di = self.dayType - 1
                hi = self.simTime.hourDay

                # Set temperature

                # add from temperature schedule for cooling
                self.BEM[i].building.cool_setpoint_day = \
                    self.Sch[i].cool[di][hi] + 273.15
                self.BEM[i].building.cool_setpoint_night = \
                    self.BEM[i].building.cool_setpoint_day
                # add from temperature schedule for heating
                self.BEM[i].building.heat_setpoint_day = \
                    self.Sch[i].heat[di][hi] + 273.15
                self.BEM[i].building.heat_setpoint_night = \
                    self.BEM[i].building.heat_setpoint_day

                # Internal Heat Load Schedule (W/m^2 of floor area for Q)

                # Qelec x elec fraction for day
                self.BEM[i].elec = self.Sch[i].q_elec * \
                    self.Sch[i].elec[di][hi]
                # Qlight x light fraction for day
                self.BEM[i].light = self.Sch[i].q_light * \
                    self.Sch[i].light[di][hi]
                # Number of occupants x occ fraction for day
                self.BEM[i].Nocc = self.Sch[i].n_occ * self.Sch[i].occ[di][hi]
                # Sensible Q occ * fraction occ sensible Q * number of occ
                self.BEM[i].Qocc = self.sensocc * \
                    (1 - self.latfocc) * self.BEM[i].Nocc

                # SWH and ventilation schedule

                # L/hr/m2 x SWH fraction for day
                self.BEM[i].swh = self.Sch[i].v_swh * self.Sch[i].swh[di][hi]
                # m^3/s/m^2 x Vent fraction for day
                self.BEM[i].building.vent = self.Sch[i].vent
                # Gas Equip Schedule, per m^2 of floor
                self.BEM[i].gas = self.Sch[i].q_gas * self.Sch[i].gas[di][hi]

                # This is quite messy, should update
                # Update internal heat and corresponding fractional loads
                intHeat = self.BEM[i].light + \
                    self.BEM[i].elec + self.BEM[i].Qocc
                # W/m2 from light, electricity, occupants
                self.BEM[i].building.int_heat_day = intHeat
                self.BEM[i].building.int_heat_night = intHeat
                # fraction of radiant heat from light/equipment of whole internal heat
                self.BEM[i].building.int_heat_f_rad = \
                    (self.radflight * self.BEM[i].light + self.radfequip *
                     self.BEM[i].elec) / intHeat

                # fraction of latent heat (from occupants) of whole internal heat
                self.BEM[i].building.int_heat_flat = \
                    self.latfocc * self.sensocc * self.BEM[i].Nocc / intHeat

                # Update envelope temperature layers
                self.BEM[i].T_wallex = self.BEM[i].wall.layerTemp[0]
                self.BEM[i].T_wallin = self.BEM[i].wall.layerTemp[-1]
                self.BEM[i].T_roofex = self.BEM[i].roof.layerTemp[0]
                self.BEM[i].T_roofin = self.BEM[i].roof.layerTemp[-1]

            # Update rural heat fluxes & update vertical diffusion model (VDM)
            self.rural.infra = self.forc.infra - self.rural.emissivity * self.SIGMA * \
                self.rural.layerTemp[0] ** 4.    # Infrared radiation from rural road
            self.rural.SurfFlux(self.forc, self.geoParam, self.simTime,
                                self.forc.hum, self.forc.temp, self.forc.wind, 2., 0.)
            self.RSM.vdm(self.forc, self.rural, self.geoParam, self.simTime)

            # Calculate urban heat fluxes, update UCM & UBL
            self.UCM, self.UBL, self.BEM = urbflux(
                self.UCM, self.UBL, self.BEM, self.forc, self.geoParam, self.simTime,
                self.RSM)
            self.UCM.UCModel(self.BEM, self.UBL.ublTemp,
                             self.forc, self.geoParam)
            self.UBL.ublmodel(
                self.UCM, self.RSM, self.rural, self.forc, self.geoParam, self.simTime)

            # Experimental code to run diffusion model in the urban area
            # N.B Commented out in python UWG because computed wind speed in urban
            # VDM: y = =0.84 * ln((2 - x / 20) / 0.51) results in negative log for
            # building heights >= 40m.

            # Uroad = copy.copy(self.UCM.road)
            # Uroad.sens = copy.copy(self.UCM.sensHeat)
            # Uforc = copy.copy(self.forc)
            # Uforc.wind = copy.copy(self.UCM.canWind)
            # Uforc.temp = copy.copy(self.UCM.canTemp)
            # self.USM.VDM(Uforc,Uroad,self.geoParam,self.simTime)

            istimestep = utilities.is_near_zero(
                self.simTime.secDay % self.simTime.timePrint, 1e-10)
            if istimestep and n < self.N:

                self.WeatherData[n] = copy.copy(self.forc)
                _Tdb, _w, self.UCM.canRHum, _h, self.UCM.Tdp, _v = psychrometrics(
                    self.UCM.canTemp, self.UCM.canHum, self.forc.pres)

                self.UBLData[n] = copy.copy(self.UBL)
                self.UCMData[n] = copy.copy(self.UCM)
                self.RSMData[n] = copy.copy(self.RSM)

                n += 1

    def write_epw(self):
        """Write new EPW file to new_epw_path property."""

        epw_prec = self.epw_precision  # precision of epw file input

        for iJ in range(len(self.UCMData)):
            # [iJ+self.simTime.timeInitial-8] - increments along weather timestep in epw
            # [6 to 21]                       - column data of epw

            # dry bulb temperature  [C]
            self.epwinput[iJ + self.simTime.timeInitial - 8][6] = \
                '{0:.{1}f}'.format(self.UCMData[iJ].canTemp - 273.15, epw_prec)
            # dew point temperature [?C]
            self.epwinput[iJ + self.simTime.timeInitial - 8][7] = \
                '{0:.{1}f}'.format(self.UCMData[iJ].Tdp, epw_prec)
            # relative humidity     [%]
            self.epwinput[iJ + self.simTime.timeInitial - 8][8] = \
                '{0:.{1}f}'.format(self.UCMData[iJ].canRHum, epw_prec)
            # wind speed [m/s]
            self.epwinput[iJ + self.simTime.timeInitial - 8][21] = \
                '{0:.{1}f}'.format(self.WeatherData[iJ].wind, epw_prec)

        # Writing new EPW file
        epw_new_id = open(self.new_epw_path, 'w')

        for i in range(8):
            new_epw_line = \
                '{}\n'.format(
                    reduce(lambda x, y: x + ',' + y, self._header[i]))
            epw_new_id.write(new_epw_line)

        for i in range(len(self.epwinput)):
            printme = ''
            for ei in range(34):
                printme += "{}".format(self.epwinput[i][ei]) + ','
            printme = printme + "{}".format(self.epwinput[i][ei])
            new_epw_line = '{0}\n'.format(printme)
            epw_new_id.write(new_epw_line)

        epw_new_id.close()

        print('New climate file is generated at {}.'.format(
              self.new_epw_path))

    def _read_input(self, param_path):
        """Read the parameter input file (.uwg file) and set as UWG attributes."""

        # Open .UWG file and feed csv data to initializeDataFile
        param_data = utilities.read_csv(param_path)

        # The initialize.UWG is read with a dictionary so that users changing
        # line endings or line numbers doesn't make reading input incorrect
        count = 0
        self._init_param_dict = {}
        while count < len(param_data):
            row = param_data[count]
            row = [c.replace(' ', '').lower() for c in row]

            # optional parameters might be empty so handle separately
            is_optional_parameter = \
                row[0] in self.OPTIONAL_PARAMETER_SET if len(
                    row) > 0 else False

            try:
                if row == [] or '#' in row[0]:
                    count += 1
                elif row[0] == 'schtraffic':
                    # SchTraffic: 3 x 24 matrix
                    trafficrows = param_data[count + 1:count + 4]
                    self._init_param_dict[row[0]] = \
                        [utilities.str2fl(r[:24]) for r in trafficrows]
                    count += 4
                elif row[0] == 'bld':
                    # bld: list of (bldtype, builtera, frac)
                    count += 1
                    bld = []
                    nextrow = [c.replace(' ', '').lower()
                               for c in param_data[count]]
                    while len(nextrow) > 0 and nextrow[0] in REF_BLDTYPE_SET:
                        bldtype, builtera, frac = nextrow[0], nextrow[1], nextrow[2]
                        bld.append((bldtype, builtera, float(frac)))
                        count += 1
                        nextrow = [c.replace(' ', '').lower()
                                   for c in param_data[count]]
                    self._init_param_dict[row[0]] = bld
                elif row[0] == 'zone':
                    self._init_param_dict[row[0]] = row[1]
                    count += 1
                elif is_optional_parameter:
                    self._init_param_dict[row[0]] = None if row[1] == '' else float(
                        row[1])
                    count += 1
                else:
                    self._init_param_dict[row[0]] = float(row[1])
                    count += 1
            except (ValueError, IndexError) as e:
                print(e)
                print('Error while reading parameter at row {}. Got: {}.'.format(
                    count, row))

        # Set UWG parameters
        for attr in self.PARAMETER_LIST:
            assert attr in self._init_param_dict, 'The {} attribute is not defined in ' \
                'the {} parameter file.'.format(
                    attr, os.path.split(param_path)[-1])
            setattr(self, attr, self._init_param_dict[attr])

    def _read_epw(self):
        """Read EPW file and sets corresponding UWG weather and site properties.

        This function will set the following attributes in the UWG object:

        * epwinput - EPW weather data as list
        * lat - latitude from EPW header
        * lon - longitude from EPW header
        * gmt - Greenwich Mean Time from EPW header
        * nSoil - Number of soil depths from EPW header
        * Tsoil - nSoil x 12 matrix for soil temperture from EPW header[1]
        * depth_soil - nSoil x 1 matrix for soil depth from EPW header[1]
        * new_epw_path - file path to modified EPW file

        Note:
            [1] EPW CSV Format (In/Out). Retrieved August 26, 2020,
            from https://bigladdersoftware.com/epx/docs/8-2/auxiliary-programs/epw-csv-format-inout.html
        """

        # Open epw file and feed csv data to climate_data
        climate_data = utilities.read_csv(self.epw_path)

        # Read header lines (1 to 8) from EPW and ensure TMY2 format.
        self._header = climate_data[0:8]

        # Read weather data from EPW for each time step in weather file. (lines 8 - end)
        self.epwinput = climate_data[8:]

        # Read Lat, Long (line 1 of EPW)
        self.lat = float(self._header[0][6])
        self.lon = float(self._header[0][7])
        self.gmt = float(self._header[0][8])

        # Read in soil temperature data (assumes this is always there)
        soilData = self._header[3]
        # Number of ground temperature depths
        self.nSoil = int(soilData[1])
        # nSoil x 12 matrix for soil temperture (K)
        self.Tsoil = [[0 for c in range(12)] for r in range(self.nSoil)]
        # nSoil x 1 matrix for soil depth (m)
        self.depth_soil = [[0] for r in range(self.nSoil)]

        # Read monthly data for each layer of soil from EPW file
        for i in range(self.nSoil):
            # get soil depth for each nSoil
            self.depth_soil[i][0] = float(soilData[2 + (i * 16)])
            # monthly data
            for j in range(12):
                # 12 months of soil T for specific depth
                self.Tsoil[i][j] = float(soilData[6 + (i * 16) + j]) + 273.15

    def _compute_BEM(self):
        """Define BEMDef objects for each archetype defined in the bld matrix.

        This function will set the following attributes in the UWG object:

        * r_glaze_total - Area-weighted average of glazing ratio from urban building
            stock [Km/W].
        * SHGC_total - Area-weighted average of SHGC from urban building stock.
        * alb_wall_total - Area-weighted average of wall albedo from urban building
            stock.
        * BEM - list of BEMDef objects extracted from readDOE
        * Sch - list of Schedule objects extracted from readDOE
        """

        # Define building energy models
        k = 0
        self.r_glaze_total = 0.
        self.SHGC_total = 0.
        self.alb_wall_total = 0.
        h_floor = self.flr_h or 3.05  # average floor height

        total_urban_bld_area = math.pow(self.charlength, 2) * self.blddensity * \
            self.bldheight / h_floor  # total building floor area

        self.BEM = []  # list of BEMDef objects
        self.Sch = []  # list of Schedule objects

        # Modify zone to be used as python index
        zone_idx = REF_ZONETYPE.index(self.zone)

        # Build unique key based on bldtype and builtera strings
        bld_dict = {bldg[0] + bldg[1]: (bldg[0], REF_BUILTERA.index(bldg[1]), bldg[2])
                    for bldg in self.bld}

        for i in range(len(self.refBEM)):  # ~16 building types (more w/ custom refs)
            for j in range(3):  # ~ 3 built eras

                if not self.refBEM[i][j][0]:
                    # when add custom types some matrix elements are None
                    continue

                ref_key = \
                    self.refBEM[i][j][0].bldtype + \
                    self.refBEM[i][j][0].builtera
                if ref_key in bld_dict:
                    # Add to BEM list
                    bldtype, builtera_idx, frac = bld_dict[ref_key]
                    self.BEM.append(self.refBEM[i][builtera_idx][zone_idx])
                    self.BEM[k].frac = frac
                    self.BEM[k].fl_area = frac * total_urban_bld_area

                    # Overwrite with optional parameters if provided
                    if self.glzr:
                        self.BEM[k].building.glazing_ratio = self.glzr
                    if self.albroof:
                        self.BEM[k].roof.albedo = self.albroof
                    if self.vegroof:
                        self.BEM[k].roof.vegcoverage = self.vegroof
                    if self.shgc:
                        self.BEM[k].building.shgc = self.shgc
                    if self.albwall:
                        self.BEM[k].wall.albedo = self.albwall
                    if self.flr_h:
                        self.BEM[k].building.floor_height = self.flr_h

                    # Keep track of total urban r_glaze, SHGC, and alb_wall for UCM model
                    self.r_glaze_total += \
                        self.BEM[k].frac * self.BEM[k].building.glazing_ratio
                    self.SHGC_total += self.BEM[k].frac * \
                        self.BEM[k].building.shgc
                    self.alb_wall_total += self.BEM[k].frac * \
                        self.BEM[k].wall.albedo
                    # Add to schedule list
                    self.Sch.append(
                        self.refSchedule[i][builtera_idx][zone_idx])
                    k += 1

    def _compute_input(self):
        """Create input objects from user-defined parameters.

        This function will set the following input objects as attributes in the UWG
        object:

        * simTime - SimParam object for time parameters
        * weather - Weather object for weather data
        * forcIP - Forcing object to estimate deep ground/water temperature
        * forc - empty Forcing object to estimate deep ground/water temperature
        * RSM - RSMDef object for rural site and vertical diffusion model
        * USM - RSMDef object for urban site and vertical diffusion model
        * geoParam - Param object for geographic parameters
        * UBL - UBL object for urban boundary layer model
        * UCM - UCMDef object for urban canopy model
        * road - Element object for urban road
        * rural - Element object for rural road
        """
        self.simTime = SimParam(self.dtsim, self.dtweather, self.month,
                                self.day, self.nday)  # simulation time parameters

        # weather file data for simulation time period
        self.weather = Weather(
            self.epw_path, self.simTime.timeInitial, self.simTime.timeFinal)
        self.forcIP = Forcing(self.weather.staTemp,
                              self.weather)  # init Forcing obj
        self.forc = Forcing()  # empty Forcing obj

        # Initialize geographic Param and Urban Boundary Layer Objects
        nightStart = 18.  # arbitrary values for begin/end hour for night setpoint
        nightEnd = 8.
        maxdx = 250.  # max dx (m)

        self.geoParam = Param(
            self.h_ubl1, self.h_ubl2, self.h_ref, self.h_temp, self.h_wind, self.c_circ,
            self.maxday, self.maxnight, self.lattree, self.latgrss, self.albveg,
            self.vegstart, self.vegend, nightStart, nightEnd, self.windmin, self.WGMAX,
            self.c_exch, maxdx, self.G, self.CP, self.VK, self.R, self.RV, self.LV,
            math.pi, self.SIGMA, self.WATERDENS, self.LVTT, self.TT, self.ESTT, self.CL,
            self.CPV, self.B, self.CM, self.COLBURN)

        self.UBL = UBLDef(
            'C', self.charlength, self.weather.staTemp[0], maxdx,
            self.geoParam.dayBLHeight, self.geoParam.nightBLHeight)

        # Defining road
        emis = 0.93
        asphalt = Material(self.kroad, self.croad, 'asphalt')
        road_T_init = 293.
        road_horizontal = 1
        # fraction of vegetation (tree & grass) coverage on unbuilt surface
        road_veg_coverage = self.vegcover / (1 - self.blddensity)
        road_grass_coverage = self.treecover / (1 - self.blddensity)
        road_tree_coverage = self.grasscover / (1 - self.blddensity)

        # define road layers
        road_layer_num = int(math.ceil(self.droad / 0.05))
        # 0.5/0.05 ~ 10 x 1 matrix of 0.05 thickness
        thickness_vector = [0.05 for r in range(road_layer_num)]
        material_vector = [asphalt for r in range(road_layer_num)]

        self.road = Element(
            self.albroad, emis, thickness_vector, material_vector, road_veg_coverage,
            road_T_init, road_horizontal, name='urban_road')
        self.road.treecoverage = road_tree_coverage
        self.road.grasscoverage = road_grass_coverage

        self.rural = Element(
            self.albroad, emis, thickness_vector, material_vector, self.rurvegcover,
            road_T_init, road_horizontal, name='rural_road')

        # Reference site class (also include VDM)
        self.RSM = RSMDef(
            self.lat, self.lon, self.gmt, self.h_obs, self.weather.staTemp[0],
            self.weather.staPres[0], self.geoParam, self.Z_MESO_PATH)
        self.USM = RSMDef(
            self.lat, self.lon, self.gmt, self.bldheight /
            10., self.weather.staTemp[0],
            self.weather.staPres[0], self.geoParam, self.Z_MESO_PATH)

        T_init = self.weather.staTemp[0]
        H_init = self.weather.staHum[0]
        wind_init = self.weather.staUmod[0]

        self.UCM = UCMDef(
            self.bldheight, self.blddensity, self.vertohor, self.treecover,
            self.sensanth, self.latanth, T_init, H_init, wind_init,
            self.geoParam, self.r_glaze_total, self.SHGC_total, self.alb_wall_total,
            self.road)

        self.UCM.h_mix = self.h_mix

        # Define Road Element & buffer to match ground temperature depth
        roadMat, newthickness = UWG._procmat(
            self.road, self.MAXTHICKNESS, self.MINTHICKNESS)

        for i in range(self.nSoil):
            # if soil depth is greater then the thickness of the road
            # we add new slices of soil at max thickness until road is greater or equal

            is_soildepth_equal = utilities.is_near_zero(
                self.depth_soil[i][0] - sum(newthickness), 1e-15)

            if is_soildepth_equal or (self.depth_soil[i][0] > sum(newthickness)):
                while self.depth_soil[i][0] > sum(newthickness):
                    newthickness.append(self.MAXTHICKNESS)
                    roadMat.append(self.SOIL)
                self._soilindex1 = i
                break

        self.road = Element(
            self.road.albedo, self.road.emissivity, newthickness, roadMat,
            self.road.vegcoverage, self.road.layerTemp[0], self.road.horizontal,
            self.road.name)

        # Define Rural Element
        ruralMat, newthickness = self._procmat(
            self.rural, self.MAXTHICKNESS, self.MINTHICKNESS)

        for i in range(self.nSoil):
            # if soil depth is greater then the thickness of the road
            # we add new slices of soil at max thickness until road is greater or equal

            is_soildepth_equal = utilities.is_near_zero(
                self.depth_soil[i][0] - sum(newthickness), 1e-15)

            if is_soildepth_equal or (self.depth_soil[i][0] > sum(newthickness)):
                while self.depth_soil[i][0] > sum(newthickness):
                    newthickness.append(self.MAXTHICKNESS)
                    ruralMat.append(self.SOIL)

                self._soilindex2 = i
                break

        self.rural = Element(
            self.rural.albedo, self.rural.emissivity, newthickness, ruralMat,
            self.rural.vegcoverage, self.rural.layerTemp[0], self.rural.horizontal,
            self.rural.name)

    def _hvac_autosize(self):
        """HVAC autosizing (unlimited cooling & heating) in BEM objects."""

        for i in range(len(self.BEM)):
            if self.autosize:
                self.BEM[i].building.coolcap = 9999.
                self.BEM[i].building.heat_cap = 9999.

    def _check_reference_data(self, ref_bem_vector, ref_sch_vector):
        """Check validity of reference BEMDef and SchDef lists.

        Args:
            ref_bem_vector: List of custom BEMDef objects to override or add to
                the refBEM matrix according to the BEMDef bldtype and builtera values.
            ref_sch_vector: List of custom SchDef objects to override or add to
                the refSchedule matrix according to the SchDef bldtype and builtera
                values.

        Returns:
            Tuple consisting of validated ref_bem_vector and ref_sch_vector.
        """
        # check for null vector
        if ref_bem_vector is None:
            return None, None

        assert len(ref_sch_vector) == len(ref_bem_vector), 'The ' \
            'ref_sch_vector and ref_bem_vector properties must be lists of equal ' \
            'length. Got lengths {} and {}, respectively.'.format(
                len(ref_sch_vector), len(ref_bem_vector))

        assert all(isinstance(v, SchDef) for v in ref_sch_vector), 'All items in the ' \
            'ref_sch_vector must be a SchDef object.'
        assert all(isinstance(v, BEMDef) for v in ref_bem_vector), 'All items in the ' \
            'ref_bem_vector must be a BEMDef object.'

        for ref_bem, ref_sch in zip(ref_bem_vector, ref_sch_vector):
            assert ref_bem.bldtype == ref_sch.bldtype, 'The bldtype for ' \
                'corresponding items in ref_bem_vector and ref_sch_vector must be the ' \
                'same. Got {} and {}.'.format(ref_bem.bldtype, ref_sch.bldtype)
            assert ref_bem.builtera == ref_sch.builtera, 'The builtera for ' \
                'corresponding items in ref_bem_vector and ref_sch_vector must be the ' \
                'same. Got {} and {}.'.format(
                    ref_bem.builtera, ref_sch.builtera)

        return ref_bem_vector, ref_sch_vector

    def _customize_reference_data(self):
        """Customize refBEM and refSchedule data by extending or overriding DOE reference data.

        The custom BEMDef and SchDef objects must contain bldtype and builtera
        identifiers referencing a nonzero fraction of urban area in the UWG bld property
        to be used in the UWG model. Also note that this method should be used before
        calling the generate method, in order to ensure the reference data gets
        transferred over to the UWG object BEM and Schedule properties.

        Args:
            ref_bem_vector: List of custom BEMDef objects to override or add to
                the refBEM matrix according to the BEMDef bldtype and builtera values.
            ref_sch_vector: List of custom SchDef objects to override or add to
                the refSchedule matrix according to the SchDef bldtype and builtera
                values.
        """
        zi = REF_ZONETYPE.index(self.zone)
        # Insert or extend refSchedule matrix
        for sch in self.ref_sch_vector:
            ei = REF_BUILTERA.index(sch.builtera)
            try:
                ti = REF_BLDTYPE.index(sch.bldtype)  # will fail if custom type
                print('Overwrite DOE reference schedule "{} {}" '
                      'with custom schedule.'.format(sch.builtera, sch.bldtype))
            except ValueError:
                # Add new rows based on type index in object
                ti = len(self.refSchedule)
                self.refSchedule.append([[None for c in range(16)]
                                         for r in range(3)])
                print('Add custom schedule for "{} {}".'.format(
                    sch.builtera, sch.bldtype))
            self.refSchedule[ti][ei][zi] = sch

        # Insert or extend refBEM matrix
        for bem in self.ref_bem_vector:
            ei = REF_BUILTERA.index(bem.builtera)
            try:
                ti = REF_BLDTYPE.index(bem.bldtype)  # will fail if custom type
                print('Overwrite DOE reference BEM "{} {}" '
                      'with custom schedule.'.format(bem.builtera, bem.bldtype))
            except ValueError:
                # Add new rows based on type index in object
                ti = len(self.refBEM)
                self.refBEM.append([[None for c in range(16)]
                                    for r in range(3)])
                print('Add custom bem for "{} {}".'.format(
                    bem.builtera, bem.bldtype))
            self.refBEM[ti][ei][zi] = bem

    @ staticmethod
    def load_refDOE(refDOE_path=REFDOE_PATH):
        """Static method to deserialize DOE reference data.

        Args:
            readDOE_path: Text string for full path to the refDOE pickle.
            (Default: the filepath specified in the UWG.REFDOE_PATH constant).

        Returns:
            Two 16 x 3 x 16 matrices of reference BEMDef and SchDef objects,
            respectively.
        """
        assert os.path.exists(refDOE_path), \
            'File: {} does not exist.'.format(refDOE_path)

        # open pickle file in binary form
        with open(refDOE_path, 'rb') as refDOE_file:
            refBEM = pickle.load(refDOE_file)
            refSchedule = pickle.load(refDOE_file)
        return refBEM, refSchedule

    @ staticmethod
    def _procmat(materials, max_thickness, min_thickness):
        """Processes material layer slices in Element objects based on layer number and depth.

        This function will divide a material with single layer thickness or that is too
        thick into subdivisions.
        """
        newmat = []
        newthickness = []
        k = materials.layerThermalCond
        Vhc = materials.layerVolHeat

        if len(materials.layer_thickness_lst) > 1:

            for j in range(len(materials.layer_thickness_lst)):
                # Break up each layer that's more than max thickness (0.05m)
                if materials.layer_thickness_lst[j] > max_thickness:
                    nlayers = math.ceil(
                        materials.layer_thickness_lst[j] / float(max_thickness))
                    for i in range(int(nlayers)):
                        newmat.append(
                            Material(k[j], Vhc[j], name=materials.name))
                        newthickness.append(
                            materials.layer_thickness_lst[j] / float(nlayers))
                # Material that's less then min_thickness is not added.
                elif materials.layer_thickness_lst[j] < min_thickness:
                    print('WARNING: Material layer too thin (less then 2 cm) to process.'
                          'Material {} is {:.2f} cm.'.format(
                              materials.name, min_thickness * 100))
                else:
                    newmat.append(Material(k[j], Vhc[j], name=materials.name))
                    newthickness.append(materials.layer_thickness_lst[j])

        else:

            # Divide single layer into two (UWG assumes at least 2 layers)
            if materials.layer_thickness_lst[0] > max_thickness:
                nlayers = math.ceil(
                    materials.layer_thickness_lst[0] / float(max_thickness))
                for i in range(int(nlayers)):
                    newmat.append(Material(k[0], Vhc[0], name=materials.name))
                    newthickness.append(
                        materials.layer_thickness_lst[0] / float(nlayers))
            # Material should be at least 1cm thick, so if we're here,
            # should give warning and stop. Only warning given for now.
            elif materials.layer_thickness_lst[0] < min_thickness * 2:
                newthickness = [min_thickness / 2., min_thickness / 2.]
                newmat = [Material(k[0], Vhc[0], name=materials.name),
                          Material(k[0], Vhc[0], name=materials.name)]
                print('WARNING: Material layer less then 2 cm is found.'
                      'Material {} is {:.2f} cm. May cause error.'.format(
                          materials.name, min_thickness * 100))
            else:
                newthickness = [materials.layer_thickness_lst[0] / 2.,
                                materials.layer_thickness_lst[0] / 2.]
                newmat = [Material(k[0], Vhc[0], name=materials.name),
                          Material(k[0], Vhc[0], name=materials.name)]

        return newmat, newthickness

    def ToString(self):
        """Overwrite .NET ToString method."""
        return self.__repr__()

    def __repr__(self):
        def _split_string(s):
            return s[0] + ':\n  ' + s[1].replace(',', '\n  ')

        def _tabbed(s):
            return _split_string(s.__repr__().split(':'))

        def _list_2_tabbed(b):
            return reduce(lambda a, b: a + '\n' + b, [_tabbed(_b) for _b in b])

        simtime_str = _tabbed(self.simTime) + '\n' \
            if hasattr(self, 'simTime') else 'No simTime attr.\n'
        weather_str = _tabbed(self.weather) + '\n' \
            if hasattr(self, 'weather') else 'No weather attr.\n'
        param_str = _tabbed(self.geoParam) + '\n' \
            if hasattr(self, 'geoParam') else 'No geoParam attr.\n'
        ubl_str = _tabbed(self.UBL) + '\n' \
            if hasattr(self, 'UBL') else 'No UBL attr.\n'
        rsm_str = 'Rural ' + _tabbed(self.RSM) + '\n' \
            if hasattr(self, 'RSM') else 'No Rural RSM attr.\n'
        usm_str = 'Urban ' + _tabbed(self.USM) + '\n' \
            if hasattr(self, 'USM') else 'No Urban RSM attr.\n'
        ucm_str = _tabbed(self.UCM) + '\n' \
            if hasattr(self, 'UCM') else 'No UCM attr.\n'
        bem_str = reduce(lambda a, b: a.__repr__() + '\n' + b.__repr__(), self.BEM) \
            if hasattr(self, 'BEM') else 'No BEM attr.'

        return 'UWG for {}:\n\n{}{}{}{}{}{}{}{}'.format(
            self.epw_path, simtime_str, weather_str, param_str, ubl_str, rsm_str,
            usm_str, ucm_str, bem_str)
