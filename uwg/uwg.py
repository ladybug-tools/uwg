"""Urban Weather Generator (UWG) Version 4.2

Original Author: B. Bueno[1]
Edited by A. Nakano & Lingfu Zhang
Modified by Joseph Yang (joeyang@mit.edu) - May, 2016
Translated to Python by Saeran Vasanthakumar - February, 2018

Note:
    [1] Bueno, Bruno; Norford, Leslie; Hidalgo, Julia; Pigeon, Gregoire (2013).
    The urban weather generator, Journal of Building Performance Simulation. 6:4,269-281.
    doi: 10.1080/19401493.2012.718797
"""

from __future__ import division, print_function
from functools import reduce
from honeybee.typing import int_in_range, float_in_range, int_positive, float_positive

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


class UWG(object):
    """Morph a rural EPW file to urban conditions based on defined urban parameters.

    Args:
        epw_path: Text string for the name of the rural epw file that will be morphed.
        param_path: Optional text string for the UWG parameter file (.UWG) path.
            If None the UWG input parameters must be manually set in the UWG object.
            (Default: None).
        new_epw_dir: Optional text string destination directory for the morphed
            EPW file. If None the morphed file will be written into the same directory
            as the rural EPW file (the epwDir). (Default: None).
        new_epw_name: Optional destination file name for the morphed EPW file.
            If None the morphed file will append '_UWG' to the original file name.
            (Default: None).

    Properties:
        * epw_path
        * new_epw_path
        * refBEM
        * refSchedule
        * month
        * day
        * nday
        * dtsim
        * dtweather
        * autosize
        * sensocc
        * latfocc
        * radfocc
        * radfequip
        * radflight
        * h_ubl1
        * h_ubl2
        * h_ref
        * h_temp
        * h_wind
        * c_circ
        * c_exch
        * maxday
        * maxnight
        * windmin
        * h_obs
        * bldheight
        * h_mix
        * blddensity
        * vertohor
        * charlength
        * albroad
        * droad
        * sensanth
        * zone
        * vegcover
        * treecoverage
        * vegstart
        * vegend
        * albveg
        * rurvegcover
        * latgrss
        * lattree
        * schtraffic
        * kroad
        * croad
        * bld
        * albroof
        * vegroof
        * glzr
        * albwall
        * shgc
        * flr_h
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
    PARAMETER_LIST = ['month', 'day', 'nday', 'dtsim', 'dtweather', 'autosize',
                      'sensocc', 'latfocc', 'radfocc', 'radfequip', 'radflight',
                      'h_ubl1', 'h_ubl2', 'h_ref', 'h_temp', 'h_wind', 'c_circ',
                      'c_exch', 'maxday', 'maxnight', 'windmin', 'h_obs', 'bldheight',
                      'h_mix', 'blddensity', 'vertohor', 'charlength', 'albroad',
                      'droad', 'sensanth', 'zone', 'vegcover', 'treecoverage',
                      'vegstart', 'vegend', 'albveg', 'rurvegcover', 'latgrss',
                      'lattree', 'schtraffic', 'kroad', 'croad', 'bld', 'shgc',
                      'albroof', 'glzr', 'vegroof', 'albwall', 'flr_h']
    OPTIONAL_PARAMETER_SET = {'shgc', 'albroof', 'glzr', 'vegroof', 'albwall', 'flr_h'}
    DEFAULT_BLD = [
        [0, 0, 0],  # FullServiceRestaurant
        [0, 0, 0],  # Hospital
        [0, 0, 0],  # LargeHotel
        [0, 0.4, 0],  # LargeOffice
        [0, 0, 0],  # MediumOffice
        [0, 0.6, 0],  # MidRiseApartment
        [0, 0, 0],  # OutPatient
        [0, 0, 0],  # PrimarySchool
        [0, 0, 0],  # QuickServiceRestaurant
        [0, 0, 0],  # SecondarySchool
        [0, 0, 0],  # SmallHotel
        [0, 0, 0],  # SmallOffice
        [0, 0, 0],  # Stand-aloneRetail
        [0, 0, 0],  # StripMall
        [0, 0, 0],  # SuperMarket
        [0, 0, 0]]  # Warehouse
    DEFAULT_SCHTRAFFIC = [
        [0.2, 0.2, 0.2, 0.2, 0.2, 0.4, 0.7, 0.9, 0.9, 0.6, 0.6, 0.6, 0.6, 0.6, 0.7, 0.8,
         0.9, 0.9, 0.8, 0.8, 0.7, 0.3, 0.2, 0.2],  # Weekday
        [0.2, 0.2, 0.2, 0.2, 0.2, 0.3, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.6, 0.7,
         0.7, 0.7, 0.7, 0.5, 0.4, 0.3, 0.2, 0.2],  # Saturday
        [0.2, 0.2, 0.2, 0.2, 0.2, 0.3, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
         0.4, 0.4, 0.4, 0.4, 0.3, 0.3, 0.2, 0.2]]  # Sunday

    # Constant file paths
    CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
    Z_MESO_PATH = os.path.join(CURRENT_PATH, 'refdata', 'z_meso.txt')
    READDOE_PATH = os.path.join(CURRENT_PATH, 'refdata', 'readDOE.pkl')

    def __init__(self, epw_path, new_epw_dir=None, new_epw_name=None):

        # Logger will be disabled by default unless explicitly called in tests
        self.logger = logging.getLogger(__name__)

        # set filepath data
        self.epw_path = epw_path
        self._new_epw_dir, self._new_epw_name = new_epw_dir, new_epw_name
        self._new_epw_path = None

        # set defaults for reference data
        self._refBEM = None
        self._refSchedule = None

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
    def from_param_args(cls, epw_path, bldheight, blddensity, vertohor, zone, month=1,
                        day=1, nday=31, dtsim=300, dtweather=3600, autosize=False,
                        h_mix=1, sensOcc=100, latfocc=0.3, radfocc=0.2, radfequip=0.5,
                        radflight=0.7, bld=DEFAULT_BLD, charlength=1000, albroad=0.1,
                        droad=0.5, sensanth=20, kroad=1, croad=1600000, vegcover=0.2,
                        treecoverage=0.1, vegstart=4, vegend=10, albveg=0.25,
                        rurvegcover=0.9, latgrss=0.4, lattree=0.6,
                        schtraffic=DEFAULT_SCHTRAFFIC, h_ubl1=1000, h_ubl2=80, h_ref=150,
                        h_temp=2, h_wind=10, c_circ=1.2, c_exch=1, maxday=150,
                        maxnight=20, windmin=1, h_obs=0.1, new_epw_dir=None,
                        new_epw_name=None):

        uwg_model = UWG(epw_path, new_epw_dir, new_epw_name)

        # Set defaults of parameters
        uwg_model.bldheight = bldheight
        uwg_model.blddensity = blddensity
        uwg_model.vertohor = vertohor
        uwg_model.zone = zone
        uwg_model.month = month
        uwg_model.day = day
        uwg_model.nday = nday
        uwg_model.dtsim = dtsim
        uwg_model.dtweather = dtweather
        uwg_model.autosize = autosize
        uwg_model.h_mix = h_mix
        uwg_model.sensOcc = sensOcc
        uwg_model.latfocc = latfocc
        uwg_model.radfocc = radfocc
        uwg_model.radfequip = radfequip
        uwg_model.radflight = radflight
        uwg_model.bld = bld
        uwg_model.charlength = charlength
        uwg_model.albroad = albroad
        uwg_model.droad = droad
        uwg_model.sensanth = sensanth
        uwg_model.kroad = kroad
        uwg_model.croad = croad
        uwg_model.vegcover = vegcover
        uwg_model.treecoverage = treecoverage
        uwg_model.vegstart = vegstart
        uwg_model.vegend = vegend
        uwg_model.albveg = albveg
        uwg_model.rurvegcover = rurvegcover
        uwg_model.latgrss = latgrss
        uwg_model.lattree = lattree
        uwg_model.schtraffic = schtraffic
        uwg_model.h_ubl1 = h_ubl1
        uwg_model.h_ubl2 = h_ubl2
        uwg_model.h_ref = h_ref
        uwg_model.h_temp = h_temp
        uwg_model.h_wind = h_wind
        uwg_model.c_circ = c_circ
        uwg_model.c_exch = c_exch
        uwg_model.maxday = maxday
        uwg_model.maxnight = maxnight
        uwg_model.windmin = windmin
        uwg_model.h_obs = h_obs

        return uwg_model

    @classmethod
    def from_param_file(cls, epw_path, param_path, new_epw_dir=None, new_epw_name=None):
        """Morph a rural EPW file to urban conditions based on .uwg file."""

        uwg_model = UWG(epw_path, new_epw_dir, new_epw_name)
        uwg_model._read_input(param_path)
        return uwg_model

    @classmethod
    def from_dict(cls, data):
        """Create an UWG object from a dictionary.

        Args:
            data: An UWG dictionary following the format below. Note that
                this example has been truncated for the sake of brevity. For
                the full list of required properties in the UWG, see the
                initialization docstrings.

        .. code-block:: python

            {
            "type": "UWG",
            "epw_path": "/path/to/epw/SGP_Singapore.486980_IWEC.epw",
            "new_epw_dir": null,
            "new_epw_name": null,
            "bldheight": 10,
            "blddensity": 0.5,
            "vertohor": 0.8,
            ...
            "h_obs": 0.1,
            "flr_h": 3.5,
            "shgc": None,
            "ref_sch_vector": [sch.to_dict()]  # Optional vector of SchDef dictionary.
            "ref_bem_vector": [bem.to_dict()]  # Optional vector of BEMDef dictionary.
            }
        """
        assert data['type'] == 'UWG', \
            'Expected UWG dictionary. Got {}.'.format(data['type'])
        assert 'epw_path' in data, \
            'The epw_path must be defined to create an UWG object.'

        uwg_model = UWG(data['epw_path'], data['new_epw_dir'], data['new_epw_name'])

        # set UWG parameters
        for attr in cls.PARAMETER_LIST:
            setattr(uwg_model, attr, data[attr])

        # check and add reference data
        refcheck = int('ref_sch_vector' in data) + int('ref_bem_vector' in data)
        assert refcheck == 2 or refcheck == 0, 'The ref_sch_vector and ref_bem_vector ' \
            'properties must both be defined in order to modify the UWG reference ' \
            'data. Only {} is defined.'.format(
                'ref_sch_vector' if 'ref_sch_vector' in data else 'ref_bem_vector')

        if refcheck == 2:
            ref_sch_vector = [SchDef.from_dict(schdict)
                              for schdict in data['ref_sch_vector']]
            ref_bem_vector = [BEMDef.from_dict(bemdict)
                              for bemdict in data['ref_bem_vector']]
            uwg_model.customize_reference_data(ref_bem_vector, ref_sch_vector)

        return uwg_model

    def customize_reference_data(self, ref_bem_vector, ref_sch_vector):
        """Customize refBEM and refSchedule data by extending or overriding DOE reference data.

        The custom BEMDef and SchDef objects must contain bldtype and builtera values
        referencing a nonzero fraction of urban area in the UWG bld matrix to be used in
        the UWG model. Also note that this method should be used before calling the
        generate method, in order to ensure the reference data gets transferred over to
        the UWG object BEM and Schedule properties.

        Args:
            ref_bem_vector: List of custom SchDef objects to add to the refSch matrix
                property.
            ref_sch_vector: List of custom BEMDef objects to add to the refBEM matrix
                property.
        """
        assert len(ref_sch_vector) == len(ref_bem_vector), 'The ref_sch_vector ' \
            'and ref_bem_vector properties must be lists of equal length. Got ' \
            'lengths {} and {}, respectively.'.format(
                len(ref_sch_vector), len(ref_bem_vector))

        zi = self.zone - 1

        # Insert or extend refSchedule matrix
        for sch in ref_sch_vector:
            ti, ei = sch.bldtype, sch.builtera
            try:
                self.refSchedule[ti][ei][zi] = sch
            except IndexError:
                # Add new rows based on type index in object
                new_rows_num = ti + 1 - len(self.refSchedule)
                for i in range(new_rows_num):
                    self.refSchedule.append(
                        [[None for c in range(16)] for r in range(3)])
                self.refSchedule[ti][ei][zi] = sch

        # Insert or extend refBEM matrix
        for bem in ref_bem_vector:
            ti, ei = bem.bldtype, bem.builtera
            try:
                self.refBEM[ti][ei][zi] = bem
            except IndexError:
                # Add new rows based on type index in object
                new_rows_num = ti + 1 - len(self.refBEM)
                for i in range(new_rows_num):
                    self.refBEM.append(
                        [[None for c in range(16)] for r in range(3)])
                self.refBEM[ti][ei][zi] = bem

    def to_dict(self, include_refDOE=False):
        """UWG dictionary representation.

        Args:
            add_refDOE: Optional boolean to include reference BEMDef and SchDef objects
                from the refBEM and refSch matrices. Only BEMDef and SchDef objects
                with a bldtype and builtera value referenced in the bld matrix will be
                included. Set this value to True if custom reference data has been
                added to the UWG object. (Default: False).
        """

        base = {'type': 'UWG'}
        base['epw_path'] = self.epw_path
        base['new_epw_dir'] = self._new_epw_dir
        base['new_epw_name'] = self._new_epw_path

        # Add UWG parameters
        for attr in self.PARAMETER_LIST:
            base[attr] = getattr(self, attr)

        # Add reference data
        if include_refDOE:

            zi = self.zone - 1
            type_num = len(self.bld)
            ref_sch_vec, ref_bem_vec = [], []
            for ti in range(type_num):
                for ei in range(3):
                    if utilities.is_near_zero(self.bld[ti][ei], 1e-5):
                        continue
                    ref_sch_vec.append(
                        self.refSchedule[ti][ei][zi].to_dict())
                    ref_bem_vec.append(
                        self.refBEM[ti][ei][zi].to_dict())

            base['ref_sch_vector'] = ref_sch_vec
            base['ref_bem_vector'] = ref_bem_vec

        return base

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
            self._new_epw_path = os.path.join(self._new_epw_dir, self._new_epw_name)
        return self._new_epw_path

    @property
    def refBEM(self):
        if self._refBEM is None:
            self._load_readDOE(self.READDOE_PATH)
        return self._refBEM

    @property
    def refSchedule(self):
        if self._refSchedule is None:
            self._load_readDOE(self.READDOE_PATH)
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
        """Get or set sensible heat in Watts from occupant."""
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
        """Get or set microclimate inversion height in meters."""
        return self._h_ref

    @h_ref.setter
    def h_ref(self, value):
        self._h_ref = float_positive(value, 'h_ref')

    @property
    def h_temp(self):
        """Get or set microclimate temperature height in meters."""
        return self._h_temp

    @h_temp.setter
    def h_temp(self, value):
        self._h_temp = float_positive(value, 'h_temp')

    @property
    def h_wind(self):
        """Get or set microclimate wind height in meters."""
        return self._h_wind

    @h_wind.setter
    def h_wind(self, value):
        self._h_wind = float_positive(value, 'h_wind')

    @property
    def c_circ(self):
        """Get or set microclimate circulation coefficient."""
        return self._c_circ

    @c_circ.setter
    def c_circ(self, value):
        self._c_circ = float_positive(value, 'c_circ')

    @property
    def c_exch(self):
        """Get or set microclimate exchange coefficient."""
        return self._c_exch

    @c_exch.setter
    def c_exch(self, value):
        self._c_exch = float_positive(value, 'c_exch')

    @property
    def maxday(self):
        """Get or set microclimate maximum day threshold."""
        return self._maxday

    @maxday.setter
    def maxday(self, value):
        try:
            self._maxday = float(value)
        except TypeError:
            raise TypeError('Input maxday must be an integer or float. Got: '
                            '{}.'.format(value))

    @property
    def maxnight(self):
        """Get or set microclimate maximum night threshold."""
        return self._maxnight

    @maxnight.setter
    def maxnight(self, value):
        try:
            self._maxnight = float(value)
        except TypeError:
            raise TypeError('Input maxnight must be an integer or float. Got: '
                            '{}.'.format(value))

    @property
    def windmin(self):
        """Get or set microclimate minimum wind speed in m/s."""
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
    def blddensity(self):
        """Get or set building footprint density relative to urban area."""
        return self._blddensity

    @blddensity.setter
    def blddensity(self, value):
        self._blddensity = float_in_range(value, 0, 1, 'blddensity')

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
        """Get or set street level anthropogenic sensible heat [W/m^2].

        Street level anthropogenic heat is non-building heat like heat emitted from cars,
        pedestrians, and street cooking.
        """
        return self._sensanth

    @sensanth.setter
    def sensanth(self, value):
        self._sensanth = float_positive(value, 'sensanth')

    @property
    def bld(self):
        """Get or set matrix representing fraction of urban building stock.

        This property consists of a 16 x 3 matrix referencing the fraction of the urban
        building stock from 16 building types and 3 built eras representing, in
        combination with 16 climate zones, 768 building archetypes generated from the
        Commercial Building Energy Survey. The sum of the fractional values in the bld
        matrix must sum to one. Each column represent a pre-1980's, post-1980's, or new
        construction era, and rows represent building types, for example:

        .. code-block:: python

            # Represent 40% post-1980's LargeOffice, and 60% new construction
            # MidRiseApartment.

            bld = [[0, 0, 0],  # FullServiceRestaurant
                   [0, 0, 0],  # Hospital
                   [0, 0, 0],  # LargeHotel
                   [0, 0. 4,0],  # LargeOffice
                   [0, 0, 0],  # MediumOffice
                   [0, 0, 0.6],  # MidRiseApartment
                   [0, 0, 0],  # OutPatient
                   [0, 0, 0],  # PrimarySchool
                   [0, 0, 0],  # QuickServiceRestaurant
                   [0, 0, 0],  # SecondarySchool
                   [0, 0, 0],  # SmallHotel
                   [0, 0, 0],  # SmallOffice
                   [0, 0, 0],  # Stand-aloneRetail
                   [0, 0, 0],  # StripMall
                   [0, 0, 0],  # SuperMarket
                   [0, 0, 0]]  # Warehouse
        """
        return self._bld

    @bld.setter
    def bld(self, value):

        assert isinstance(value, (list, tuple)), 'The bld property must be a list ' \
            'or tuple. Got {}.'.format(value)
        type_num = len(value)
        assert type_num >= 16, 'The bld property must have greater than or equal to ' \
            '16 rows. Got {} rows.'.format(type_num)

        self._bld = [[0 for c in range(3)] for r in range(type_num)]

        # Check column number and add value
        for i in range(type_num):
            assert len(value[i]) == 3, 'The bld property must be a 16 (or greater) ' \
                'x 3 matrix. Got {} columns for the row {}.'.format(len(value[i]), i)
            for j in range(3):
                self._bld[i][j] = float_in_range(value[i][j])

    @property
    def lattree(self):
        """Get or set fraction of latent heat absorbed by tree."""
        return self._lattree

    @lattree.setter
    def lattree(self, value):
        self._lattree = float_in_range(value, 0, 1, 'lattree')

    @property
    def latgrss(self):
        """Get or set fraction of latent heat absorbed by grass."""
        return self._latgrss

    @latgrss.setter
    def latgrss(self, value):
        self._latgrss = float_in_range(value, 0, 1, 'latgrss')

    @property
    def zone(self):
        """Get or set number representing an ASHRAE climate zone.

        Choose from the following:

            1 - 1A (Miami)
            2 - 2A (Houston)
            3 - 2B (Phoenix)
            4 - 3A (Atlanta)
            5 - 3B-CA (Los Angeles)
            6 - 3B (Las Vegas)
            7 - 3C (San Francisco)
            8 - 4A (Baltimore)
            9 - 4B (Albuquerque)
            10 - 4C (Seattle)
            11 - 5A (Chicago)
            12 - 5B (Boulder)
            13 - 6A (Minneapolis)
            14 - 6B (Helena)
            15 - 7 (Duluth)
            16 - 8 (Fairbanks)
        """
        return self._zone

    @zone.setter
    def zone(self, value):
        self._zone = int_in_range(value, 1, 16, 'zone')

    @property
    def vegstart(self):
        """Get or set number for the month in which vegetation starts to evapotranspire.

        This month corresponds to when the leaves of vegetation are assumed to be out.
        """
        return self._vegstart

    @vegstart.setter
    def vegstart(self, value):
        self._vegstart = int_in_range(value, 1, 12, 'vegstart')

    @property
    def vegend(self):
        """Get or set number for the month in which vegetation stops evapotranspiration.

        This month corresponds to when the leaves of vegetation are assumed to fall.
        """
        return self._vegend

    @vegend.setter
    def vegend(self, value):
        self._vegend = int_in_range(value, 1, 12, 'vegend')

    @property
    def vegcover(self):
        """Get or set fraction of urban ground covered in grass only."""
        return self._vegcover

    @vegcover.setter
    def vegcover(self, value):
        self._vegcover = float_in_range(value, 0, 1, 'vegcover')

    @property
    def treecoverage(self):
        """Get or set fraction of urban ground covered in trees."""
        return self._treecoverage

    @treecoverage.setter
    def treecoverage(self, value):
        self._treecoverage = float_in_range(value, 0, 1, 'treecoverage')

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
        """Get or set average building glazing Solar Heat Gain Coefficient."""
        return self._shgc

    @shgc.setter
    def shgc(self, value):
        if value is None:
            self._shgc = value
        else:
            self._shgc = float_in_range(value, 0, 1, 'shgc')

    @property
    def albroof(self):
        """Get or set average building roof albedo."""
        return self._albroof

    @albroof.setter
    def albroof(self, value):
        if value is None:
            self._albroof = value
        else:
            self._albroof = float_in_range(value, 0, 1, 'albroof')

    @property
    def glzr(self):
        """Get or set average building glazing ratio."""
        return self._glzr

    @glzr.setter
    def glzr(self, value):
        if value is None:
            self._glzr = value
        else:
            self._glzr = float_in_range(value, 0, 1, 'glzr')

    @property
    def vegroof(self):
        """Get or set fraction of roofs covered in grass/shrubs."""
        return self._vegroof

    @vegroof.setter
    def vegroof(self, value):
        if value is None:
            self._vegroof = value
        else:
            self._vegroof = float_in_range(value, 0, 1, 'vegroof')

    @property
    def albwall(self):
        """Get or set average building albedo."""
        return self._albwall

    @albwall.setter
    def albwall(self, value):
        if value is None:
            self._albwall = value
        else:
            self._albwall = float_in_range(value, 0, 1, 'albwall')

    @property
    def flr_h(self):
        """Get or set average building floor height in meters."""
        return self._flr_h

    @flr_h.setter
    def flr_h(self, value):
        if value is None:
            self._flr_h = value
        else:
            self._flr_h = float_positive(value, 'flr_h')

    def generate(self):
        """Generate all UWG objects after input parameters are set."""

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

        self.N = int(self.simTime.days * 24)  # total number of hours in simulation
        n = 0  # weather time step counter
        self.ph = self.simTime.dt / 3600.  # dt (simulation time step) in hours

        # Data dump variables
        self.WeatherData = [None for x in range(self.N)]
        self.UCMData = [None for x in range(self.N)]
        self.UBLData = [None for x in range(self.N)]
        self.RSMData = [None for x in range(self.N)]
        self.USMData = [None for x in range(self.N)]

        print('\nSimulating new temperature and humidity values for '
              '{} days from {}/{}.\n'.format(self.nday, self.month, self.day))
        self.logger.info('Start simulation')

        # iterate through every simulation time-step (i.e 5 min) defined by UWG
        for it in range(1, self.simTime.nt, 1):
            # Update water temperature (estimated)
            if self.nSoil < 3:
                # for BUBBLE/CAPITOUL/Singapore only
                self.forc.deepTemp = sum(self.forcIP.temp) / float(len(self.forcIP.temp))
                self.forc.waterTemp = \
                    sum(self.forcIP.temp) / float(len(self.forcIP.temp)) - 10.0
            else:
                # soil temperature by depth, by month
                self.forc.deepTemp = self.Tsoil[self._soilindex1][self.simTime.month - 1]
                self.forc.waterTemp = self.Tsoil[2][self.simTime.month-1]

            # There's probably a better way to update the weather...
            self.simTime.update_date()

            self.logger.info('\n{0} m={1}, d={2}, h={3}, s={4}'.format(
                __name__, self.simTime.month, self.simTime.day,
                self.simTime.secDay / 3600., self.simTime.secDay))

            # simulation time increment raised to weather time step
            self.ceil_time_step = int(math.ceil(it * self.ph))-1
            # minus one to be consistent with forcIP list index
            # Updating forcing instance
            # horizontal Infrared Radiation Intensity (W m-2)
            self.forc.infra = self.forcIP.infra[self.ceil_time_step]
            # wind speed (m s-1)
            self.forc.wind = max(self.forcIP.wind[self.ceil_time_step],
                                 self.geoParam.windMin)
            self.forc.uDir = self.forcIP.uDir[self.ceil_time_step]  # wind direction
            # specific humidty (kg kg-1)
            self.forc.hum = self.forcIP.hum[self.ceil_time_step]
            self.forc.pres = self.forcIP.pres[self.ceil_time_step]  # Pressure (Pa)
            self.forc.temp = self.forcIP.temp[self.ceil_time_step]  # air temp (C)
            self.forc.rHum = self.forcIP.rHum[self.ceil_time_step]  # RH (%)
            self.forc.prec = self.forcIP.prec[self.ceil_time_step]  # Precip (mm h-1)
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
            if utilities.is_near_zero(self.simTime.julian % 7, 1e-10):
                self.dayType = 3                                        # Sunday
            elif utilities.is_near_zero(self.simTime.julian % 7 - 6., 1e-10):
                self.dayType = 2                                        # Saturday
            else:
                self.dayType = 1                                        # Weekday

            # Update anthropogenic heat load for each hour (building & UCM)
            self.UCM.sensAnthrop = \
                self.sensanth * (self.schtraffic[self.dayType - 1][self.simTime.hourDay])

            # Update the energy components for building types defined in initialize.UWG
            for i in range(len(self.BEM)):

                di = self.dayType - 1
                hi = self.simTime.hourDay

                # Set temperature

                # add from temperature schedule for cooling
                self.BEM[i].building.coolSetpointDay = self.Sch[i].cool[di][hi] + 273.15
                self.BEM[i].building.coolSetpointNight = \
                    self.BEM[i].building.coolSetpointDay
                # add from temperature schedule for heating
                self.BEM[i].building.heatSetpointDay = self.Sch[i].heat[di][hi] + 273.15
                self.BEM[i].building.heatSetpointNight = \
                    self.BEM[i].building.heatSetpointDay

                # Internal Heat Load Schedule (W/m^2 of floor area for Q)

                # Qelec x elec fraction for day
                self.BEM[i].elec = self.Sch[i].Qelec * self.Sch[i].elec[di][hi]
                # Qlight x light fraction for day
                self.BEM[i].light = self.Sch[i].Qlight * self.Sch[i].light[di][hi]
                # Number of occupants x occ fraction for day
                self.BEM[i].Nocc = self.Sch[i].Nocc * self.Sch[i].occ[di][hi]
                # Sensible Q occ * fraction occ sensible Q * number of occ
                self.BEM[i].Qocc = self.sensocc * (1 - self.latfocc) * self.BEM[i].Nocc

                # SWH and ventilation schedule

                # litres per hour x SWH fraction for day
                self.BEM[i].swh = self.Sch[i].Vswh * self.Sch[i].swh[di][hi]
                # m^3/s/m^2 of floor
                self.BEM[i].building.vent = self.Sch[i].Vent
                # Gas Equip Schedule, per m^2 of floor
                self.BEM[i].gas = self.Sch[i].Qgas * self.Sch[i].gas[di][hi]

                # This is quite messy, should update
                # Update internal heat and corresponding fractional loads
                intHeat = self.BEM[i].light + self.BEM[i].elec + self.BEM[i].Qocc
                # W/m2 from light, electricity, occupants
                self.BEM[i].building.intHeatDay = intHeat
                self.BEM[i].building.intHeatNight = intHeat
                # fraction of radiant heat from light/equipment of whole internal heat
                self.BEM[i].building.intHeatFRad = \
                    (self.radflight * self.BEM[i].light + self.radfequip *
                     self.BEM[i].elec) / intHeat

                # fraction of latent heat (from occupants) of whole internal heat
                self.BEM[i].building.intHeatFLat = \
                    self.latfocc * self.sensocc * self.BEM[i].Nocc / intHeat

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
                self.UCM, self.UBL, self.BEM, self.forc, self.geoParam, self.simTime,
                self.RSM)
            self.UCM.UCModel(self.BEM, self.UBL.ublTemp, self.forc, self.geoParam)
            self.UBL.UBLModel(
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

            self.logger.info('dbT = {}'.format(self.UCM.canTemp-273.15))
            if n > 0:
                logging.info('dpT = {}'.format(self.UCM.Tdp))
                logging.info('RH  = {}'.format(self.UCM.canRHum))

            istimestep = utilities.is_near_zero(
                self.simTime.secDay % self.simTime.timePrint, 1e-10)
            if istimestep and n < self.N:

                self.logger.info(
                    '{0} ----sim time step = {1}----\n\n'.format(__name__, n))

                self.WeatherData[n] = copy.copy(self.forc)
                _Tdb, _w, self.UCM.canRHum, _h, self.UCM.Tdp, _v = psychrometrics(
                    self.UCM.canTemp, self.UCM.canHum, self.forc.pres)

                self.UBLData[n] = copy.copy(self.UBL)
                self.UCMData[n] = copy.copy(self.UCM)
                self.RSMData[n] = copy.copy(self.RSM)

                self.logger.info('dbT = {}'.format(self.UCMData[n].canTemp-273.15))
                self.logger.info('dpT = {}'.format(self.UCMData[n].Tdp))
                self.logger.info('RH  = {}'.format(self.UCMData[n].canRHum))

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
                '{}\n'.format(reduce(lambda x, y: x + ',' + y, self._header[i]))
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

        assert os.path.exists(param_path), 'Parameter file "{}" does not ' \
            'exist.'.format(param_path)

        # Open .UWG file and feed csv data to initializeDataFile
        param_data = utilities.read_csv(param_path)

        # The initialize.UWG is read with a dictionary so that users changing
        # line endings or line numbers doesn't make reading input incorrect
        count = 0
        self._init_param_dict = {}
        while count < len(param_data):
            row = param_data[count]
            row = [row[i].replace(' ', '').lower() for i in range(len(row))]

            # optional parameters might be empty so handle separately
            is_optional_parameter = \
                row[0] in self.OPTIONAL_PARAMETER_SET if len(row) > 0 else False

            try:
                if row == [] or '#' in row[0]:
                    count += 1
                    continue
                elif row[0] == 'schtraffic':
                    # SchTraffic: 3 x 24 matrix
                    trafficrows = param_data[count+1:count+4]
                    self._init_param_dict[row[0]] = \
                        [utilities.str2fl(r[:24]) for r in trafficrows]
                    count += 4
                elif row[0] == 'bld':
                    # bld: 17 x 3 matrix
                    bldrows = param_data[count+1:count+17]
                    self._init_param_dict[row[0]] = \
                        [utilities.str2fl(r[:3]) for r in bldrows]
                    count += 17
                elif is_optional_parameter:
                    self._init_param_dict[row[0]] = \
                        None if row[1] == '' else float(row[1])
                    count += 1
                else:
                    self._init_param_dict[row[0]] = float(row[1])
                    count += 1
            except ValueError:
                print('Error while reading parameter at row {}. Got: {}.'.format(
                    count, row))

        # Set UWG parameters
        for attr in self.PARAMETER_LIST:
            assert attr in self._init_param_dict, 'The {} attribute is not defined in ' \
                'the .UWG parameter file.'.format(attr)
            setattr(self, attr, self._init_param_dict[attr])

    def _read_epw(self):
        """Read EPW file and sets corresponding UWG weather and site properties.

        This function will set the following attributes in the UWG object:

        * epwinput - EPW weather data as list
        * lat - latitude from EPW header
        * lon - longitude from EPW header
        * GMT - Greenwich Mean Time from EPW header
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
        self.GMT = float(self._header[0][8])

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

        * r_glaze_total - Glazing ratio summation for total building stock
        * SHGC_total - SHGC summation for total building stock
        * alb_wall_total -  Albedo wall summation for total building stock
        * BEM - list of BEMDef objects extracted from readDOE
        * Sch - list of Schedule objects extracted from readDOE
        """

        # TODO add check in _compute_BEM to check lenghts of refSchedle/reFBEM w/ bld matrix


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
        zone_idx = self.zone - 1

        for i in range(len(self.refBEM)):  # ~16 building types
            for j in range(3):  # 3 built eras
                if self.bld[i][j] > 0.:
                    # Add to BEM list
                    self.BEM.append(self.refBEM[i][j][zone_idx])
                    self.BEM[k].frac = self.bld[i][j]
                    self.BEM[k].fl_area = self.bld[i][j] * total_urban_bld_area

                    # Overwrite with optional parameters if provided
                    if self.glzr:
                        self.BEM[k].building.glazingRatio = self.glzr
                    if self.albroof:
                        self.BEM[k].roof.albedo = self.albroof
                    if self.vegroof:
                        self.BEM[k].roof.vegcoverage = self.vegroof
                    if self.shgc:
                        self.BEM[k].building.shgc = self.shgc
                    if self.albwall:
                        self.BEM[k].wall.albedo = self.albwall
                    if self.flr_h:
                        self.BEM[k].building.floorHeight = self.flr_h

                    # Keep track of total urban r_glaze, SHGC, and alb_wall for UCM model
                    self.r_glaze_total += \
                        self.BEM[k].frac * self.BEM[k].building.glazingRatio
                    self.SHGC_total += self.BEM[k].frac * self.BEM[k].building.shgc
                    self.alb_wall_total += self.BEM[k].frac * self.BEM[k].wall.albedo
                    # Add to schedule list
                    self.Sch.append(self.refSchedule[i][j][zone_idx])
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
        self.forcIP = Forcing(self.weather.staTemp, self.weather)  # init Forcing obj
        self.forc = Forcing()  # empty Forcing obj

        # Initialize geographic Param and Urban Boundary Layer Objects
        nightStart = 18.        # arbitrary values for begin/end hour for night setpoint
        nightEnd = 8.
        maxdx = 250.            # max dx (m)

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
        # fraction of surface vegetation coverage
        road_veg_coverage = min(self.vegcover / (1 - self.blddensity), 1.)

        # define road layers
        road_layer_num = int(math.ceil(self.droad/0.05))
        # 0.5/0.05 ~ 10 x 1 matrix of 0.05 thickness
        thickness_vector = [0.05 for r in range(road_layer_num)]
        material_vector = [asphalt for r in range(road_layer_num)]

        self.road = Element(
            self.albroad, emis, thickness_vector, material_vector, road_veg_coverage,
            road_T_init, road_horizontal, name='urban_road')

        self.rural = copy.deepcopy(self.road)
        self.rural.vegcoverage = self.rurvegcover
        self.rural.name = 'rural_road'

        # Reference site class (also include VDM)
        self.RSM = RSMDef(
            self.lat, self.lon, self.GMT, self.h_obs, self.weather.staTemp[0],
            self.weather.staPres[0], self.geoParam, self.Z_MESO_PATH)
        self.USM = RSMDef(
            self.lat, self.lon, self.GMT, self.bldheight / 10., self.weather.staTemp[0],
            self.weather.staPres[0], self.geoParam, self.Z_MESO_PATH)

        T_init = self.weather.staTemp[0]
        H_init = self.weather.staHum[0]

        self.UCM = UCMDef(
            self.bldheight, self.blddensity, self.vertohor, self.treecoverage,
            self.sensanth, self.latanth, T_init, H_init, self.weather.staUmod[0],
            self.geoParam, self.r_glaze_total, self.SHGC_total, self.alb_wall_total,
            self.road)

        self.UCM.h_mix = self.h_mix

        # Define Road Element & buffer to match ground temperature depth
        roadMat, newthickness = UWG._procmat(
            self.road, self.MAXTHICKNESS, self.MINTHICKNESS)

        for i in range(self.nSoil):
            # if soil depth is greater then the thickness of the road
            # we add new slices of soil at max thickness until road is greater or equal

            is_soildepth_equal = \
                utilities.is_near_zero(self.depth_soil[i][0] - sum(newthickness), 1e-15)

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
        ruralMat, newthickness = \
            self._procmat(self.rural, self.MAXTHICKNESS, self.MINTHICKNESS)

        for i in range(self.nSoil):
            # if soil depth is greater then the thickness of the road
            # we add new slices of soil at max thickness until road is greater or equal

            is_soildepth_equal = \
                utilities.is_near_zero(self.depth_soil[i][0] - sum(newthickness), 1e-15)

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
                self.BEM[i].building.coolCap = 9999.
                self.BEM[i].building.heatCap = 9999.

    def _load_readDOE(self, readDOE_path):
        # Serialized DOE reference data
        assert os.path.exists(readDOE_path), \
            'File: {} does not exist.'.format(readDOE_path)

        # open pickle file in binary form
        with open(readDOE_path, 'rb') as readDOE_file:
            self._refBEM = pickle.load(readDOE_file)
            self._refSchedule = pickle.load(readDOE_file)

    @staticmethod
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
                    nlayers = math.ceil(materials.layer_thickness_lst[j] / float(max_thickness))
                    for i in range(int(nlayers)):
                        newmat.append(Material(k[j], Vhc[j], name=materials.name))
                        newthickness.append(materials.layer_thickness_lst[j] / float(nlayers))
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
                nlayers = math.ceil(materials.layer_thickness_lst[0] / float(max_thickness))
                for i in range(int(nlayers)):
                    newmat.append(Material(k[0], Vhc[0], name=materials.name))
                    newthickness.append(materials.layer_thickness_lst[0] / float(nlayers))
            # Material should be at least 1cm thick, so if we're here,
            # should give warning and stop. Only warning given for now.
            elif materials.layer_thickness_lst[0] < min_thickness*2:
                newthickness = [min_thickness/2., min_thickness/2.]
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
