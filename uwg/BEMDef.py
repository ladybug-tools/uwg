"""Building Energy Model (BEM) class definition."""

from .building import Building
from .element import Element
from .utilities import float_in_range


class BEMDef(object):
    """Building Energy Model (BEM) definition.

    Args:
        building: Building object.
        mass: Element object representing building internal mass.
        wall: Element object representing building wall.
        roof: Element object representing building roof.
        frac: Fraction of the urban floor space of this building typology.
        bldtype: Number between 0 and 15 corresponding to the following building
            types: FullServiceRestaurant (0), Hospital (1), LargeHotel (2),
            LargeOffice (3), MediumOffice (4), MidRiseApartment (5), OutPatient (6),
            PrimarySchool (7), QuickServiceRestaurant (8), SecondarySchool (9),
            SmallHotel (10), SmallOffice (11), StandaloneRetail (12), StripMall (13),
            SuperMarket (14), Warehouse (15). Additional building types can be defined
            with a number greater then 15. This value is used to reference the fraction
            of urban area the BEMDef object defines in the UWG bld matrix.
        builtera: Number between 0 and 2 corresponding to the following built eras:
            Pre-1980s (0), Post1980s (1), New construction (2). This value is used to
            reference the fraction of urban area the BEMDef object defines in the UWG
            bld matrix.

    Properties:
        * building
        * mass
        * wall
        * roof
        * frac
        * bldtype
        * builtera
        * fl_area -- building typology urban floor area [m2]
        * elec -- actual electricity consumption [W/m2]
        * gas -- actual gas consumption [W/m2]
        * light -- actual light [W/m2]
        * Qocc -- actual heat load from occupant [W/m2]
        * swh -- actual hot water usage
        * Nocc -- number of occupants
        * ElecTotal -- actual total electricity [W/m2]
        * T_wallex -- wall surface temp (ext) [K]
        * T_wallin -- wall surface temp (int) [K]
        * T_roofex -- roof surface temp (ext) [K]
        * T_roofin -- roof surface temp (int) [K]
    """

    def __init__(self, building, mass, wall, roof, frac, bldtype, builtera):

        # Initialization
        self.building = building
        self.mass = mass
        self.wall = wall
        self.roof = roof
        self.frac = frac

        # Properties to be set in readDOE
        self.bldtype = bldtype  # DOE reference building type
        self.builtera = builtera  # pre80, pst80, new
        self.zonetype = None  # climate zone number (only used in testing).

        # Properties to be computed during UWG simulation
        self.fl_area = 0
        self.elec = 0
        self.gas = 0
        self.light = 0
        self.Qocc = 0
        self.swh = 0
        self.Nocc = 0
        self.ElecTotal = 0
        self.T_wallex = 293
        self.T_wallin = 293
        self.T_roofex = 293
        self.T_roofin = 293

    @classmethod
    def from_dict(cls, data):
        """Create a BEMDef object from a dictionary.

        Args:
            data: A BEMDef dictionary following the format below.

        .. code-block:: python
            {
            "building": building.to_dict()  # dictionary representation of a Building.
            "mass": mass.to_dict()  # dictionary representation of mass Element.
            "wall": wall.to_dict()  # dictionary representation of wall Element.
            "roof": roof.to_dict()  # dictionary representation of roof Element.
            "frac": 0.4,  # fraction of urban floor space of this type
            "bldtype": 0,  # building type index
            "builtera": 1,  # built era index
            }
        """
        assert data['type'] == 'BEMDef', 'Expected ' \
            'BEMDef dictionary. Got {}.'.format(data['type'])

        building = Building.from_dict(data['building'])
        mass = Element.from_dict(data['mass'])
        wall = Element.from_dict(data['wall'])
        roof = Element.from_dict(data['roof'])
        frac = data['frac']
        bldtype, builtera = data['bldtype'], data['builtera']

        return cls(building, mass, wall, roof, frac, bldtype, builtera)

    @property
    def building(self):
        """Get or set Building object."""
        return self._building

    @building.setter
    def building(self, value):
        assert isinstance(value, Building), 'building must be a Building ' \
            'object. Got: {}.'.format(value)
        self._building = value

    @property
    def mass(self):
        """Get or set Element as representative building internal mass."""
        return self._mass

    @mass.setter
    def mass(self, value):
        assert isinstance(value, Element), 'mass must be an Element object. ' \
            'Got: {}.'.format(value)
        self._mass = value

    @property
    def roof(self):
        """Get or set Element as representative building roof."""
        return self._roof

    @roof.setter
    def roof(self, value):
        assert isinstance(value, Element), 'roof must be an Element object. ' \
            'Got: {}.'.format(value)
        self._roof = value

    @property
    def wall(self):
        """Get or set Element as representative building wall."""
        return self._wall

    @wall.setter
    def wall(self, value):
        assert isinstance(value, Element), 'wall must be an Element object. ' \
            'Got: {}.'.format(value)
        self._wall = value

    @property
    def frac(self):
        """Get or set fraction of the urban floor space of this building typology."""
        return self._frac

    @frac.setter
    def frac(self, value):
        self._frac = float_in_range(value, 0, 1, 'frac')

    def to_dict(self):
        """BEMDef dictionary representation."""
        base = {'type': 'BEMDef'}
        base['building'] = self.building.to_dict()
        base['mass'] = self.mass.to_dict()
        base['wall'] = self.wall.to_dict()
        base['roof'] = self.roof.to_dict()
        base['frac'] = self.frac
        base['bldtype'] = self.bldtype
        base['builtera'] = self.builtera
        return base

    def __repr__(self):
        return 'BEMDef,\n bldtype: {}\n builtera: {}\n mass: {}\n ' \
            'wall: {}\n roof: {}\n frac: {}'.format(
                self.bldtype, self.builtera, self.mass.name,
                self.wall.name, self.roof.name, self.frac)
