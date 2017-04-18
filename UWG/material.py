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
