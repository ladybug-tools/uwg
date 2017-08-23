class Param(object):
    """
    PARAMETERS
    """

    def __init__(self, dayBLHeight, nightBLHeight, refHeight, tempHeight, windHeight,
                circCoeff, dayThreshold, nightThreshold, treeFLat, grassFLat, vegAlbedo, vegStart,
                vegEnd, nightSetStart, nightSetEnd, windMin, wgmax, exCoeff, maxdx, g,
                cp, vk, r, rv, lv, pi, sigma, waterDens, lvtt, tt, estt, cl, cpv, b, cm, colburn):
        obj.dayBLHeight = dayBLHeight           # daytime mixing height, orig = 700
        obj.nightBLHeight = nightBLHeight       # Sing: 80, Bub-Cap: 50, nighttime boundary-layer height (m); orig 80
        obj.refHeight = refHeight               # Reference height at which the vertical profile of potential temperature is vertical
        obj.tempHeight = tempHeight             # Temperature measuremnt height at the weather station (m)
        obj.windHeight = windHeight             # Air velocity measuremnt height at the weather station (m)
        obj.circCoeff =  circCoeff              # Wind scaling coefficient
        obj.dayThreshold = dayThreshold         # heat flux threshold for daytime conditions (W m-2)
        obj.nightThreshold = nightThreshold     # heat flux threshold for nighttime conditions (W m-2)
        obj.treeFLat = treeFLat                 # latent fraction of trees
        obj.grassFLat = grassFLat               # latent fraction of grass
        obj.vegAlbedo = vegAlbedo               # albedo of vegetation
        obj.vegStart = vegStart                 # begin month for vegetation participation
        obj.vegEnd = vegEnd                     # end month for vegetation participation
        obj.nightSetStart = nightSetStart       # begin hour for night thermal set point schedule
        obj.nightSetEnd = nightSetEnd           # end hour for night thermal set point schedule
        obj.windMin = windMin                   # minimum wind speed (m s-1)
        obj.wgmax = wgmax                       # maximum film water depth on horizontal surfaces (m)
        obj.exCoeff = exCoeff                   # exchange velocity coefficient
        obj.maxdx = maxdx                       # maximum discretization length for the UBL model (m)
        obj.g = g                               # gravity
        obj.cp = cp                             # heat capacity for air ((J/kg.K))
        obj.vk = vk                             # von karman constant
        obj.r = r                               # gas constant (for dry air?)
        obj.rv = rv                             # gas constant (for water vapor?)
        obj.lv = lv                             # latent heat of evaporation
        obj.pi = pi                             # pi
        obj.sigma = sigma                       # Stefan Boltzmann constant
        obj.waterDens = waterDens               # water density
        obj.lvtt = lvtt                         # ?
        obj.tt = tt                             # ?
        obj.estt = estt                         # ?
        obj.cl = cl                             # ?
        obj.cpv = cpv                           # ?
        obj.b   = b                             # coefficients derived by Louis (1979)
        obj.cm  = cm                            # ?
        obj.colburn = colburn                   # (Pr/Sc)^(2/3) for Colburn analogy in water evaporation

    def __repr__(self):
        return "Param w/ dayBLht {:a}, nightBLht {:b} ".format(a=self.dayBLHeight,b=self.nightBLHeight)
