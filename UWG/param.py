class Param(object):
    """
    PARAMETERS
    """

    def __init__(self, dayBLHeight, nightBLHeight, refHeight, tempHeight, windHeight,
                circCoeff, dayThreshold, nightThreshold, treeFLat, grassFLat, vegAlbedo, vegStart,
                vegEnd, nightSetStart, nightSetEnd, windMin, wgmax, exCoeff, maxdx, g,
                cp, vk, r, rv, lv, pi, sigma, waterDens, lvtt, tt, estt, cl, cpv, b, cm, colburn):
        self.dayBLHeight = dayBLHeight           # daytime mixing height, orig = 700
        self.nightBLHeight = nightBLHeight       # Sing: 80, Bub-Cap: 50, nighttime boundary-layer height (m); orig 80
        self.refHeight = refHeight               # Reference height at which the vertical profile of potential temperature is vertical
        self.tempHeight = tempHeight             # Temperature measuremnt height at the weather station (m)
        self.windHeight = windHeight             # Air velocity measuremnt height at the weather station (m)
        self.circCoeff =  circCoeff              # Wind scaling coefficient
        self.dayThreshold = dayThreshold         # heat flux threshold for daytime conditions (W m-2)
        self.nightThreshold = nightThreshold     # heat flux threshold for nighttime conditions (W m-2)
        self.treeFLat = treeFLat                 # latent fraction of trees
        self.grassFLat = grassFLat               # latent fraction of grass
        self.vegAlbedo = vegAlbedo               # albedo of vegetation
        self.vegStart = vegStart                 # begin month for vegetation participation
        self.vegEnd = vegEnd                     # end month for vegetation participation
        self.nightSetStart = nightSetStart       # begin hour for night thermal set point schedule
        self.nightSetEnd = nightSetEnd           # end hour for night thermal set point schedule
        self.windMin = windMin                   # minimum wind speed (m s-1)
        self.wgmax = wgmax                       # maximum film water depth on horizontal surfaces (m)
        self.exCoeff = exCoeff                   # exchange velocity coefficient
        self.maxdx = maxdx                       # maximum discretization length for the UBL model (m)
        self.g = g                               # gravity
        self.cp = cp                             # heat capacity for air ((J/kg.K))
        self.vk = vk                             # von karman constant
        self.r = r                               # gas constant (for dry air?)
        self.rv = rv                             # gas constant (for water vapor?)
        self.lv = lv                             # latent heat of evaporation
        self.pi = pi                             # pi
        self.sigma = sigma                       # Stefan Boltzmann constant
        self.waterDens = waterDens               # water density
        self.lvtt = lvtt                         # ?
        self.tt = tt                             # ?
        self.estt = estt                         # ?
        self.cl = cl                             # ?
        self.cpv = cpv                           # ?
        self.b   = b                             # coefficients derived by Louis (1979)
        self.cm  = cm                            # ?
        self.colburn = colburn                   # (Pr/Sc)^(2/3) for Colburn analogy in water evaporation

    def __repr__(self):
        return "Param w/ dayBLht {a}, treeFlat {b} ".format(a=self.dayBLHeight,b=self.treeFLat)
