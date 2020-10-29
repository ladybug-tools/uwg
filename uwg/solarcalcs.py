"""Class for solar calculations."""
from __future__ import division

try:
    range = xrange
except NameError:
    pass

import math


class SolarCalcs(object):
    """Solar calculations.

    Args:
        UCM: Urban Canopy - Building Energy Model object.
        BEM: Building Energy Model object.
        simTime: Simulation time object
        RSM: Rural Site & Vertical Diffusion Model Object.
        forc: Forcing object.
        parameter: Param Object.
        rural: Rural road Element object.

    Properties:
        * # Properties added from solarcalcs
        * dir -- direct sunlight (perpendicular to the sun's ray)
        * dif -- diffuse sunlight
        * tanzen -- tangente of solar zenithal angle
        * critOrient -- critical canyon angle for which solar radiation reaches the road
        * horSol -- direct horizontal radiation
        * Kw_term -- fractional terms for wall
        * Kr_term -- fractional terms for road
        * mr -- fractional terms for wall
        * mw -- fractional terms for road
        * ut -- elapsed hours on current day
        * ad -- fractional year in radians
        * eqtime -- time variable
        * decsol -- solar declination angle
        * zenith -- Angle between normal to earth's surface and sun position
    """

    def __init__(self, UCM, BEM, simTime, RSM, forc, parameter, rural):
        self.UCM = UCM
        self.BEM = BEM
        self.simTime = simTime
        self.RSM = RSM
        self.forc = forc
        self.parameter = parameter
        self.rural = rural

    def solarcalcs(self):
        """ Solar Calculation.

        Mutates RSM, BEM, and UCM objects based on following parameters:

        * UCM -- Urban Canopy - Building Energy Model object
        * BEM -- Building Energy Model object
        * simTime -- Simulation time object
        * RSM -- Rural Site & Vertical Diffusion Model Object
        * forc -- Forcing object
        * parameter -- Geo Param Object
        * rural -- Rural road Element object
        """

        # Direct sunlight (perpendicular to the sun's ray)
        self.dir = self.forc.dir
        self.dif = self.forc.dif  # Diffuse sunlight

        if (self.dir + self.dif) > 0.:

            # calculate zenith tangent, and critOrient solar angles
            self.solarangles()

            # Direct horizontal radiation
            self.horSol = max(math.cos(self.zenith) * self.dir, 0.0)
            # Fractional terms for wall & road
            self.Kw_term = (
                min(abs(1. / self.UCM.canAspect*(0.5 - self.critOrient / math.pi) + 1 /
                        math.pi * self.tanzen * (1 - math.cos(self.critOrient))), 1.))
            self.Kr_term = (
                min(abs(2. * self.critOrient / math.pi -
                        (2 / math.pi * self.UCM.canAspect * self.tanzen) *
                        (1 - math.cos(self.critOrient))),
                    1 - 2 * self.UCM.canAspect * self.Kw_term))

            # Direct and diffuse solar radiation
            # Assume trees are shorter than buildings
            self.bldSol = self.horSol * self.Kw_term + self.UCM.wallConf * self.dif
            self.roadSol = self.horSol * self.Kr_term + self.UCM.roadConf * self.dif

            # Solar reflections. Add diffuse radiation from vegetation to alb_road if
            # in season
            if self.simTime.month < self.parameter.vegStart or \
               self.simTime.month > self.parameter.vegEnd:
                alb_road = self.UCM.road.albedo
            else:
                alb_road = self.UCM.road.albedo * (1. - self.UCM.road.vegcoverage) + \
                    self.parameter.vegAlbedo * self.UCM.road.vegcoverage

            # First set of reflections
            rr = alb_road * self.roadSol
            rw = self.UCM.alb_wall * self.bldSol

            # bounces
            fr = (
                1. - (1. - 2. * self.UCM.wallConf) * self.UCM.alb_wall +
                (1. - self.UCM.roadConf) * self.UCM.wallConf * alb_road *
                self.UCM.alb_wall)

            # (1.0-self.UCM.roadConf) road to wall view
            self.mr = (rr + (1.0 - self.UCM.roadConf) * alb_road *
                       (rw + self.UCM.wallConf * self.UCM.alb_wall * rr)) / fr
            self.mw = (rw + self.UCM.wallConf * self.UCM.alb_wall * rr) / fr

            # Receiving solar, including bounces (W m-2)
            # Includes road covered by vegetation.
            self.UCM.road.solRec = self.roadSol + \
                (1 - self.UCM.roadConf) * self.mw

            for j in range(len(self.BEM)):
                self.BEM[j].roof.solRec = self.horSol + self.dif
                self.BEM[j].wall.solRec = (
                    self.bldSol + (1 - 2*self.UCM.wallConf)*self.mw +
                    self.UCM.wallConf * self.mr)

            self.rural.solRec = self.horSol + self.dif  # Solar received by rural
            self.UCM.SolRecRoof = self.horSol + self.dif  # Solar received by roof
            self.UCM.SolRecRoad = self.UCM.road.solRec  # Solar received by road
            # Solar received by wall
            self.UCM.SolRecWall = (self.bldSol + (1 - 2 * self.UCM.wallConf) *
                                   self.UCM.road.albedo * self.roadSol)

            # Modification from UWG_Matlab on 09/20:
            # Consolidate and incorporate the treeCoverage fraction into the tree
            # sensible and latent heat (W-m2) calculation. Previously this was
            # factored (inconsistently) at the urbflux and UCMDef.UCModel functions.

            # Vegetation heat (tree) (per m^2 of veg)
            self.UCM.treeSensHeat = (
                (1 - self.parameter.vegAlbedo) * (1 - self.parameter.treeFLat) *
                self.UCM.SolRecRoad * self.UCM.treeCoverage)
            self.UCM.treeLatHeat = (
                (1 - self.parameter.vegAlbedo) * self.parameter.treeFLat *
                self.UCM.SolRecRoad * self.UCM.treeCoverage)

            # Modification from UWG_Matlab on 09/20:
            # Add the sensible and latent heat fraction of grass (vegetation not
            # accounted for by tree fraction). Note that the self.UCM.road.sens property
            # also contains a calculation of sensible heat from vegetation and surface
            # convection, but was not used here because the UCModel function in the
            # UCMDef module already accounts for the entire road convective heat transfer
            # so using the UCM.road.sens property would double-count convection. This
            # sensible heat therefore just accounts for the absorbed solar radiation
            # split into it's sensible heat fraction.

            # Vegetation heat (grass) (per m^2 of veg)
            grasscover = self.UCM.vegcover - self.UCM.treeCoverage
            self.UCM.treeSensHeat += (
                (1 - self.parameter.vegAlbedo) * (1 - self.parameter.grassFLat) *
                self.UCM.SolRecRoad * grasscover)
            self.UCM.treeLatHeat += (
                (1 - self.parameter.vegAlbedo) * self.parameter.grassFLat *
                self.UCM.SolRecRoad * grasscover)

        else:  # No Sun

            self.UCM.road.solRec = 0.
            self.rural.solRec = 0.

            for j in range(len(self.BEM)):
                self.BEM[j].roof.solRec = 0.
                self.BEM[j].wall.solRec = 0.

            self.UCM.SolRecRoad = 0.  # Solar received by road
            self.UCM.SolRecRoof = 0.  # Solar received by roof
            self.UCM.SolRecWall = 0.  # Solar received by wall
            self.UCM.treeSensHeat = 0.
            self.UCM.treeLatHeat = 0.

        return self.rural, self.UCM, self.BEM

    def solarangles(self):
        """Solves for zenith angle, tangent of zenithal angle, and critical canyon angle.

        Calculation based on NOAA from:

        * canAspect -- aspect Ratio of canyon
        * simTime -- simulation time parameters
        * RSM.lon -- longitude (deg)
        * RSM.lat -- latitude (deg)
        * RSM.GMT -- GMT hour correction
        """

        month = self.simTime.month
        day = self.simTime.day
        secDay = self.simTime.secDay  # Total elapsed seconds in simulation
        # inobis: total days for first of month
        # i.e [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]
        inobis = self.simTime.inobis

        canAspect = self.UCM.canAspect
        lon = self.RSM.lon
        lat = self.RSM.lat
        GMT = self.RSM.gmt

        # Get elapsed hours on current day
        self.ut = (24. + (int(secDay) / 3600. % 24.)) % 24.
        ibis = list(range(len(inobis)))

        for JI in range(1, 12):
            ibis[JI] = inobis[JI] + 1

        date = day + inobis[month-1] - 1  # Julian day of the year
        # divide circle by 365 days, multiply by elapsed days + hours
        # ad: fractional year (radians)
        self.ad = 2.0 * math.pi / 365. * (date - 1 + (self.ut - (12 / 24.)))

        self.eqtime = (
            229.18 * (0.000075 + 0.001868 * math.cos(self.ad) - 0.032077 *
                      math.sin(self.ad) - 0.01461 * math.cos(2 * self.ad) - 0.040849 *
                      math.sin(2 * self.ad)))

        # Declination angle (angle of sun with equatorial plane)
        self.decsol = (
            0.006918 - 0.399912*math.cos(self.ad) + 0.070257*math.sin(self.ad) -
            0.006758 * math.cos(2. * self.ad) + 0.000907 * math.sin(2. * self.ad) -
            0.002697 * math.cos(3. * self.ad) + 0.00148 * math.sin(3. * self.ad))

        time_offset = self.eqtime - 4. * lon + 60 * GMT
        tst = secDay + time_offset * 60

        ha = (tst / 4. / 60. - 180.) * math.pi / 180.
        zlat = lat * (math.pi / 180.)  # change angle units to radians

        # Calculate zenith solar angle
        self.zenith = (
            math.acos(math.sin(zlat) * math.sin(self.decsol) + math.cos(zlat) *
                      math.cos(self.decsol) * math.cos(ha)))

        # tangente of solar zenithal angle
        if abs(0.5 * math.pi - self.zenith) < 1e-6:
            if 0.5 * math.pi - self.zenith > 0.:
                self.tanzen = math.tan(0.5 * math.pi - 1e-6)

            elif 0.5 * math.pi - self.zenith <= 0.:
                self.tanzen = math.tan(0.5 * math.pi + 1e-6)

        elif abs(self.zenith) < 1e-6:
            # lim x->0 tan(x) -> 0 which results in division by zero error
            # when calculating the critical canyon angle
            # so set tanzen to 1e-6 which will result in critical canyon angle = 90
            self.tanzen = 1e-6

        else:
            self.tanzen = math.tan(self.zenith)

        # critical canyon angle for which solar radiation reaches the road
        # has to do with street canyon orientation for given solar angle
        self.critOrient = math.asin(min(abs(1. / self.tanzen) / canAspect, 1.))
