"""Geographic Parameters class."""


class Param(object):
    """Geographic Parameters class.

    Args:
        dayBLHeight: Daytime mixing height [orig: 700] [m].
        nightBLHeight: Nighttime boundary-layer height
            (orig: 80, Sing: 80, Bub-Cap: 50) [m]
        refHeight: Reference height at which the vertical profile of potential
            temperature is vertical
        tempHeight: Temperature measuremnt height at the weather station [m]
        windHeight: Air velocity measuremnt height at the weather station [m]
        circCoeff: Wind scaling coefficient
        dayThreshold: Heat flux threshold for daytime conditions [W m-2]
        nightThreshold: Heat flux threshold for nighttime conditions [W m-2]
        treeFLat: Value between 0 and 1 for latent fraction of trees
        grassFLat: Value between 0 and 1 for latent fraction of grass
        vegAlbedo: Value between 0 and 1 for albedo of vegetation
        vegStart: Value between 1 and 12 for begin month for vegetation participation
        vegEnd: Value between 1 and 12 for end month for vegetation participation
        nightSetStart: Begin hour for night thermal set point schedule
        nightSetEnd: End hour for night thermal set point schedule
        windMin: Minimum wind speed [m s-1]
        wgmax: Maximum film water depth on horizontal surfaces [m]
        exCoeff: Exchange velocity coefficient
        maxdx: Maximum discretization length for the UBL model [m]
        g: Gravity [m s-2]
        cp: Heat capacity for air [J/kg K]
        vk: Von karman constant [dimensionless]
        r: Gas constant for dry air [J/kg K]
        rv: Gas constant for water vapor [J/kg K]
        lv: Latent heat of evaporation [J/kg]
        pi: Pi
        sigma: Stefan Boltzmann constant [W K-4 m-2]
        waterDens: Water density [kg m-3]
        lvtt: Undefined.
        tt: Undefined.
        estt: Undefined.
        cl: Undefined.
        cpv: Undefined.
        b  : Coefficients derived by Louis [1979]
        cm : Undefined.
        colburn: [Pr/Sc]^[2/3] for Colburn analogy in water evaporation

    Properties:
        * dayBLHeight
        * nightBLHeight
        * refHeight
        * tempHeight
        * windHeight
        * circCoeff
        * dayThreshold
        * nightThreshold
        * treeFLat
        * grassFLat
        * vegAlbedo
        * vegStart
        * vegEnd
        * nightSetStart
        * nightSetEnd
        * windMin
        * wgmax
        * exCoeff
        * maxdx
        * g
        * cp
        * vk
        * r
        * rv
        * lv
        * pi
        * sigma
        * waterDens
        * lvtt
        * tt
        * estt
        * cl
        * cpv
        * b
        * cm
        * colburn
    """

    def __init__(self, dayBLHeight, nightBLHeight, refHeight, tempHeight, windHeight,
                 circCoeff, dayThreshold, nightThreshold, treeFLat, grassFLat, vegAlbedo,
                 vegStart, vegEnd, nightSetStart, nightSetEnd, windMin, wgmax, exCoeff,
                 maxdx, g, cp, vk, r, rv, lv, pi, sigma, waterDens, lvtt, tt, estt, cl,
                 cpv, b, cm, colburn):

        self.dayBLHeight = dayBLHeight
        self.nightBLHeight = nightBLHeight
        self.refHeight = refHeight
        self.tempHeight = tempHeight
        self.windHeight = windHeight
        self.circCoeff = circCoeff
        self.dayThreshold = dayThreshold
        self.nightThreshold = nightThreshold
        self.treeFLat = treeFLat
        self.grassFLat = grassFLat
        self.vegAlbedo = vegAlbedo
        self.vegStart = vegStart
        self.vegEnd = vegEnd
        self.nightSetStart = nightSetStart
        self.nightSetEnd = nightSetEnd
        self.windMin = windMin
        self.wgmax = wgmax
        self.exCoeff = exCoeff
        self.maxdx = maxdx
        self.g = g
        self.cp = cp
        self.vk = vk
        self.r = r
        self.rv = rv
        self.lv = lv
        self.pi = pi
        self.sigma = sigma
        self.waterDens = waterDens
        self.lvtt = lvtt
        self.tt = tt
        self.estt = estt
        self.cl = cl
        self.cpv = cpv
        self.b = b
        self.cm = cm
        self.colburn = colburn

    def __repr__(self):
        return 'Parameter,\n dayBLHeight: {}\n nightBLHeight: {}\n refHeight: {}\n ' \
            'tempHeight: {}\n windHeight: {}\n circCoeff: {}\n dayThreshold: {}\n ' \
            'nightThreshold: {}\n treeFLat: {}\n grassFLat: {}\n vegAlbedo: {}\n ' \
            'vegStart: {}\n vegEnd: {}\n nightSetStart: {}\n nightSetEnd: {}\n ' \
            'windMin: {}\n wgmax: {}\n exCoeff: {}\n maxdx: {}\n'.format(
                self.dayBLHeight, self.nightBLHeight, self.refHeight, self.tempHeight,
                self.windHeight, self.circCoeff, self.dayThreshold, self.nightThreshold,
                self.treeFLat, self.grassFLat, self.vegAlbedo, self.vegStart,
                self.vegEnd, self.nightSetStart, self.nightSetEnd, self.windMin,
                self.wgmax, self.exCoeff, self.maxdx)
