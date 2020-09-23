"""Building Energy Model (BEM) class definition."""

from .building import Building
from .element import Element


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
        * fl_area
        * elec
        * gas
        * light
        * Qocc
        * swh
        * Nocc
        * ElecTotal
        * T_wallex
        * T_wallin
        * T_roofex
        * T_roofin
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
        self.fl_area = 0  # building typology urban floor area [m2]
        self.elec = 0  # actual electricity consumption [W/m2]
        self.gas = 0  # actual gas consumption [W/m2]
        self.light = 0  # actual light [W/m2]
        self.Qocc = 0  # actual heat load from occupant [W/m2]
        self.swh = 0  # actual hot water usage
        self.Nocc = 0  # number of occupants
        self.ElecTotal = 0  # actual total electricity [W/m2]
        self.T_wallex = 293  # wall surface temp (ext) [K]
        self.T_wallin = 293  # wall surface temp (int) [K]
        self.T_roofex = 293  # roof surface temp (ext) [K]
        self.T_roofin = 293  # roof surface temp (int) [K]

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
