"""Material class"""

from .utilities import float_in_range_excl_incl

try:
    str = basestring
except NameError:
    pass


class Material(object):
    """Material class.

    Args:
        thermalcond: Number for thermal conductivity [W m-1 K-1].
        volheat: Number for volumetric heat capacity [J m-3 K-1].
        name: Text string for name of the Material.

    Properties:
        * thermalcond
        * volheat
    """

    def __init__(self, thermalcond, volheat, name):
        self._name = name
        self.thermalcond = thermalcond
        self.volheat = volheat

    @property
    def name(self):
        """Get or set text string for name of Material."""
        return self._name

    @property
    def thermalcond(self):
        """Get or set number for thermal conductivity [W/(m-K)]."""
        return self._thermalcond

    @thermalcond.setter
    def thermalcond(self, value):
        self._thermalcond = \
            float_in_range_excl_incl(value, mi=0, input_name='thermalcond')

    @property
    def volheat(self):
        """Get or set number for volumetric capacity [J/(m3-K)]."""
        return self._volheat

    @volheat.setter
    def volheat(self, value):
        self._volheat = float_in_range_excl_incl(
            value, mi=0, input_name='volheat')

    @classmethod
    def from_dict(cls, data):
        """Create a Material object from a dictionary.

        Args:
            data: A Material dictionary following the format below.

        .. code-block:: python

            {
            "type": "Material",
            "name": "Concrete",
            "thermalcond": 1.311,  # thermal conductivity [W m-1 K-1]
            "volheat": 1874432.0  # volumetric heat capacity [J m-3 K-1]
            }
        """
        assert data['type'] == 'Material', 'Expected ' \
            'Material dictionary. Got {}.'.format(data['type'])

        return cls(data['thermalcond'], data['volheat'], data['name'])

    def to_dict(self):
        """Material dictionary representation."""
        base = {'type': 'Material'}
        base['name'] = self.name
        base['thermalcond'] = self.thermalcond
        base['volheat'] = self.volheat
        return base

    def __repr__(self):
        return "Material, name: {}\n thermalcond: {}\n volheat: {}".format(
            self.name, self.thermalcond, self.volheat)
