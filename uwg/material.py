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
    def __init__(self, thermalCond, volHeat, name='noname'):
        self._name = name
        self.thermalCond = thermalCond
        self.volHeat = volHeat

    def __repr__(self):
        return "Material: {}, k={}, spec vol={}".format(
            self._name, self.thermalCond, self.volHeat)
