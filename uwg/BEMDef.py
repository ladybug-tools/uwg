"""Building Energy Model (BEM) class definition."""


class BEMDef(object):
    """Building Energy Model (BEM) definition.

    Args:
        building: Building object.
        mass: Element object representing building internal mass.
        wall: Element object representing building wall.
        roof: Element object representing building roof.
        frac: Fraction of the urban floor space of this building typology.

    Properties:
        building
        mass
        wall
        roof
        frac
        Type
        Era
        Zone
        fl_area
        elec
        gas
        light
        Qocc
        swh
        Nocc
        ElecTotal
        T_wallex
        T_wallin
        T_roofex
        T_roofin
    """

    def __init__(self, building, mass, wall, roof, frac):

        # Initialization
        self.building = building
        self.mass = mass
        self.wall = wall
        self.roof = roof
        self.frac = frac

        # Properties to be set in readDOE
        self.Type = ''  # DOE reference building type
        self.Era = ''  # pre80, pst80, new
        self.Zone = ''  # climate zone number

        # Properties to be computed during UWG simulation
        self.fl_area = 0  # building typology urban floor area [m2]
        self.elec = 0  # actual electricity consumption [W/m2]
        self.gas = 0  # actual gas consumption [W/m2]
        self.light = 0  # actual light [W/m2]
        self.Qocc = 0  # actual heat load from occupant [W/m2]
        self.swh = 0  # actual hot water usage
        self.Nocc = 0  # number of occupants
        self.ElecTotal = 0  #  actual total electricity [W/m2]
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
            "frac": 0.4  # fraction of urban floor space of this type
            }
        """
        pass

    def to_dict(self):
        """BEMDef dictionary representation."""
        base = {'type': 'BEMDef'}
        base['building'] = self.building.to_dict()
        base['mass'] = self.mass.to_dict()
        base['wall'] = self.wall.to_dict()
        base['roof'] = self.roof.to_dict()
        base['frac'] = self.frac
        return base

    def __repr__(self):
        return "BEMDef: Type = {}, Zone = {}, Era = {}, Construction = {}".format(
            self.building.Type, self.building.Zone, self.building.Era, self.wall._name)
