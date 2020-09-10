"""Material class"""


class Material(object):
    """Material

    Args:
        thermalCond: Number for material thermal conductivity (W/m K)
        volHeat: NUmber for volumetric heat capacity (J/m^3 K)
        name: Text string for name of the material.

    Properties:
        * thermalCond
        * volHeat
    """
    def __init__(self, thermalCond, volHeat, name=''):
        self.name = name
        self.thermalCond = thermalCond
        self.volHeat = volHeat

    @classmethod
    def from_dict(cls, data):
        """Create a Material object from a dictionary.

        Args:
            data: A Material dictionary following the format below.

        .. code-block:: python

            {
            "type": "Material",
            "name": "Concrete",
            "thermalCond": 1.311,  # thermal conductivity [W m-1 K-1]
            "volHeat": 1874432.0  # volumetric heat capacity [J m-3 K-1]
            }
        """
        assert data['type'] == 'Material', 'Expected ' \
            'Material dictionary. Got {}.'.format(data['type'])

        return cls(data['thermalCond'], data['volHeat'], data['name'])

    def to_dict(self):
        """Material dictionary representation."""
        base = {'type': 'Material'}
        base['name'] = self.name
        base['thermalCond'] = self.thermalCond
        base['volHeat'] = self.volHeat
        return base

    def __repr__(self):
        return "Material, name: {}\n thermalCond: {}\n volHeat: {}".format(
            self.name, self.thermalCond, self.volHeat)
