"""Material class"""


class Material(object):
    """uwg Material

    Attributes:
        thermalCond: Thermal conductivity (W/m K)
        volHeat: Volumetric heat capacity (J/m^3 K)
        name: Name of the material.
    """
    def __init__(self, thermalCond, volHeat, name='noname'):
        self._name = name  # purely for internal purpose
        self.thermalCond = thermalCond
        self.volHeat = volHeat

    def __repr__(self):
        return "Material: {a}, k={b}, spec vol={c}".format(
            a=self._name,
            b=self.thermalCond,
            c=self.volHeat
            )
