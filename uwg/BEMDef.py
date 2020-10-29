"""Building Energy Model (BEM) class definition."""

from .building import Building
from .element import Element
from .utilities import REF_BUILTERA_SET, REF_BUILTERA

try:
    str = basestring
except NameError:
    pass


class BEMDef(object):
    """Building Energy Model (BEM) definition.

    Args:
        building: Building object.
        mass: Element object representing building internal mass.
        wall: Element object representing building wall.
        roof: Element object representing building roof.
        bldtype: Text referring to a building type. By default, 16 building types are
            defined in the UWG according to models from the Department of Energy (DOE).
            Custom building types can also be defined with a new name. Note that this
            value along with the BEMDef builtera must exactly match the identifiers in
            the UWG bld list in order to specify the fraction of total built stock the
            building occupies in the UWG simulation. Choose from the following to
            reference or overwrite a BEM associated with a DOE reference building type:
            'fullservicerestaurant', 'hospital', 'largehotel', 'largeoffice',
            'mediumoffice', 'midriseapartment', 'outpatient', 'primaryschool',
            'quickservicerestaurant', 'secondaryschool', 'smallhotel', 'smalloffice',
            'standaloneretail', 'stripmall', 'supermarket', or 'warehouse'.
        builtera: Text defining building built era. Must be one of the following:
            "pre80" (pre-1980s), "pst80" (post-1980s), or "new" (new construction).
            This value along with the bldtype must exactly match the identifiers in
            the bld array in order to specify the fraction of total built stock the
            building occupies in the UWG simulation.

    Properties:
        * building -- Building object
        * mass -- Element object representing building internal mass.
        * wall -- Element object representing building wall.
        * roof -- Element object representing building roof.
        * frac -- fraction of the urban floor space of this building typology.
        * bldtype -- number for building type
        * builtera -- number for built era
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

    def __init__(self, building, mass, wall, roof, bldtype, builtera):

        # Initialization
        self.building = building
        self.mass = mass
        self.wall = wall
        self.roof = roof

        # Properties to be set in readDOE
        self.bldtype = bldtype  # DOE reference building type
        self.builtera = builtera  # pre80, pst80, new
        self.zonetype = None  # climate zone number (only used in testing).

        # Properties to be computed during UWG simulation
        self.frac = 0.0
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
        bldtype, builtera = data['bldtype'], data['builtera']

        return cls(building, mass, wall, roof, bldtype, builtera)

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
    def bldtype(self):
        """Get or set text for building type.

        By default, 16 building types are defined in the UWG according to models from
        the Department of Energy (DOE). Choose from the following to reference or
        overwrite a BEM associated with a DOE reference building type:

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

        Custom building types can also be defined with a new name. If a custom BEMDef is
        defined with the same name as a reference DOE building type from the list above,
        the reference BEMDef will be overwritten by the custom BEMDef. Note that this
        value along with the BEMDef builtera must exactly match the identifiers in the
        UWG bld list in order to specify the fraction of total built stock the building
        occupies in the UWG simulation.

        """
        return self._bldtype

    @bldtype.setter
    def bldtype(self, value):
        assert isinstance(value, str), 'The bldtype must be a string. ' \
            'Got: {}.'.format(value.lower())
        self._bldtype = value

    @property
    def builtera(self):
        """Get or set text for built era.

        Must be one of the following:

        * 'pre80' - pre-1980s
        * 'pst80' - post-1980s
        * 'new' - new construction

        This value along with the bldtype must exactly match the identifiers in
        the bld array in order to specify the fraction of total built stock the
        building occupies in the UWG simulation.
        """
        return self._builtera

    @builtera.setter
    def builtera(self, value):
        assert isinstance(value, str) and value in REF_BUILTERA_SET, \
            'The builtera must be one of {}.Got: {}.'.format(
                REF_BUILTERA, value.lower())
        self._builtera = value

    def to_dict(self):
        """BEMDef dictionary representation."""
        base = {'type': 'BEMDef'}
        base['building'] = self.building.to_dict()
        base['mass'] = self.mass.to_dict()
        base['wall'] = self.wall.to_dict()
        base['roof'] = self.roof.to_dict()
        base['bldtype'] = self.bldtype
        base['builtera'] = self.builtera
        return base

    def __repr__(self):
        return 'BEMDef,\n bldtype: {}\n builtera: {}\n mass: {}\n ' \
            'wall: {}\n roof: {}\n frac: {}'.format(
                self.bldtype, self.builtera, self.mass.name,
                self.wall.name, self.roof.name, self.frac)
