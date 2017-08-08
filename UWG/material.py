"""Material class"""

class Material(object):
    """UWG Material

    Attributes:
        thermalCond: Thermal conductivity (W/m K)
        volHeat: Volumetric heat capacity (J/m^3 K)
        name: Name of the material.
    """
    def __init__(self, thermalCond, volHeat, name=None):
        self.name = name
        self.thermalCond = thermalCond
        self.volHeat = volHeat

    def __repr__(self):
        return "Name: {a}, k={b}, spec vol={c}".format(a=self.name,b=self.thermalCond,c=self.volHeat)
