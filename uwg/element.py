"""Element class."""
from __future__ import division

try:
    range = xrange
except NameError:
    pass

import math
from .utilities import is_near_zero
from .material import Material


class Element(object):
    """Element object defines wall, and roof constructions.

    Args:
        alb: Number for outer surface albedo.
        emis: Number for outer surface emissivity.
        thicknessLst: List of layer thickness in meters.
        materialLst: List of Material objects in Element.
        vegCoverage: Number for fraction of surface vegetation coverage.
        T_init: Element initial temperature [K].
        horizontal: Boolean indicating if Element is horizontal or not (vertical).
        name: Text string for name of Element.

    Properties:
        * albedo
        * emissivity
        * layerThickness
        * layerThermalCond
        * layerVolHeat
        * vegCoverage
        * layerTemp
        * waterStorage
        * horizontal
        * solRec
        * infra
        * lat
        * sens
        * solAbs
        * aeroCond
        * T_ext
        * T_int
        * flux
    """

    def __init__(self, alb, emis, thicknessLst, materialLst, vegCoverage, T_init,
                 horizontal, name):

        assert len(thicknessLst) == len(materialLst), 'The number of layer thickness ' \
            'must match the number of layer materials. Got {} and {}, ' \
            'respectively.'.format(len(thicknessLst), len(materialLst))

        self.albedo = alb  # outer surface albedo
        self.emissivity = emis  # outer surface emissivity.
        self.layerThickness = thicknessLst  # list of layer thickness in meters
        self.materialLst = materialLst  # material objects in Element.
        self.vegCoverage = vegCoverage  # surface vegetation coverage
        self.T_init = T_init  # element initial temperature [K].
        self.horizontal = int(horizontal)  # 1-horizontal, 0-vertical
        self.name = name

        # layerThermaCond: vector of layer thermal conductivities [W m-1 K-1]
        self.layerThermalCond = [material.thermalCond for material in materialLst]
        # layerVolHeat: vector of layer volumetric heat [J m-3 K-1]
        self.layerVolHeat = [material.volHeat for material in materialLst]
        # layerTemp: # vector of layer temperatures [K]
        self.layerTemp = [T_init] * len(thicknessLst)
        # waterStorage: # thickness of water film [m] for horizontal surfaces only
        self.waterStorage = 0
        self.infra = 0  # net longwave radiation [W m-2]
        self.sens = 0  # surface sensible heat flux [W m-2]
        self.solRec = 0  # solar radiation received [W m-2]
        self.lat = 0  # surface latent heat flux [W m-2]
        self.solAbs = 0  # solar radiation absorbed [W m-2]
        self.aeroCond = 0  # convective heat transfer
        self.T_ext = 293  # external surface temperature
        self.T_int = 293  # internal surface temperature
        self.flux = 0  # external surface heat flux

    @classmethod
    def from_dict(cls, data):
        """Create a Element object from a dictionary.

        Args:
            data: A Element dictionary following the format below.

        .. code-block:: python

            {
            "albedo": 0.08,
            "emissivity": 0.92,
            "layerThickness": [0.0254, 0.0508, 0.0254],
            "materialLst": material_dict_lst,  # List of material dictionaries
            "vegCoverage": 0.5,
            "T_init": 293,
            "horizontal": False,
            "name": MassWall
            }
        """
        assert data['type'] == 'Element', 'Expected ' \
            'Element dictionary. Got {}.'.format(data['type'])

        materials = [Material.from_dict(m) for m in data['materialLst']]
        return cls(data['albedo'], data['emissivity'], data['layerThickness'],
                   materials, data['vegCoverage'], data['T_init'],
                   data['horizontal'], data['name'])

    def to_dict(self):
        """Element dictionary representation."""
        base = {'type': 'Element'}
        base['albedo'] = self.albedo
        base['emissivity'] = self.emissivity
        base['layerThickness'] = self.layerThickness
        base['materialLst'] = [m.to_dict() for m in self.materialLst]
        base['vegCoverage'] = self.vegCoverage
        base['T_init'] = self.T_init
        base['horizontal'] = self.horizontal
        base['name'] = self.name
        return base

    def SurfFlux(self, forc, parameter, simTime, humRef, tempRef, windRef, boundCond,
                 intFlux):
        """ Calculate net heat flux, and update element layer temperatures."""

        # Calculated per unit area [m^2]

        # dens: air density (kgd m-3)
        dens = forc.pres / (1000 * 0.287042 * tempRef * (1. + 1.607858 * humRef))
        self.aeroCond = 5.8 + 3.7 * windRef  # convection coef (ref: uwg, eq. 12))

        if self.horizontal:
            # For roof, mass, road
            if not is_near_zero(self.waterStorage) and self.waterStorage > 0.0:
                # Evaporation [m s-1], Film water & soil latent heat
                # Note: in the current uwg code, latent heat from evapotranspiration,
                # stagnant water, or anthropogenic sources is not modelled due to the
                # difficulty of validation, and lack of reliability of precipitation
                # data from EPW files.Therefore this condition is never run because all
                # elements have had their waterStorage hardcoded to 0.
                qtsat = self.qsat([self.layerTemp[0]], [forc.pres], parameter)[0]
                eg = self.aeroCond * parameter.colburn * dens * \
                    (qtsat - humRef) / parameter.waterDens / parameter.cp
                self.waterStorage = \
                    min(self.waterStorage + simTime.dt * (forc.prec - eg),
                        parameter.wgmax)
                self.waterStorage = max(self.waterStorage, 0.)  # [m]
            else:
                eg = 0.
            soilLat = eg * parameter.waterDens * parameter.lv

            if simTime.month < parameter.vegStart and simTime.month > parameter.vegEnd:
                # Winter, no veg
                self.solAbs = (1.0 - self.albedo) * self.solRec  # (W m-2)
                vegLat = 0.
                vegSens = 0.
            else:
                # Summer, veg
                self.solAbs = ((1.0 - self.vegCoverage) * (1. - self.albedo) +
                               self.vegCoverage * (1.0 - parameter.vegAlbedo)) * \
                                   self.solRec
                vegLat = self.vegCoverage * parameter.grassFLat * \
                    (1. - parameter.vegAlbedo) * self.solRec
                vegSens = self.vegCoverage * (1. - parameter.grassFLat) * \
                    (1. - parameter.vegAlbedo) * self.solRec
            self.lat = soilLat + vegLat

            # Sensible & net heat flux
            self.sens = vegSens + self.aeroCond * (self.layerTemp[0] - tempRef)
            self.flux = -self.sens + self.solAbs + self.infra - self.lat  # [W m-2]

        else:
            # For vertical surfaces (wall)
            self.solAbs = (1.0 - self.albedo) * self.solRec
            self.lat = 0.0

            # Sensible & net heat flux
            self.sens = self.aeroCond * (self.layerTemp[0] - tempRef)
            self.flux = -self.sens + self.solAbs + self.infra - self.lat  # (W m-2)

        self.layerTemp = \
            self.Conduction(simTime.dt, self.flux, boundCond, forc.deepTemp, intFlux)
        self.T_ext = self.layerTemp[0]
        self.T_int = self.layerTemp[-1]

    def Conduction(self, dt, flx1, bc, temp2, flx2):
        """Solve the conductance of heat based on of the element layers.

        Args:
            flx1: net heat flux on surface [W m-2]
            bc: boundary condition parameter (1 or 2)
            temp2: deep soil temperature (ave of air temperature) [K]
            flx2: surface flux (sum of absorbed, emitted, etc.) [W m-2]

        Returns:
            A 1d vector of element layer temperatures.
        """
        t = self.layerTemp          # vector of layer temperatures (K)
        hc = self.layerVolHeat      # vector of layer volumetric heat (J m-3 K-1)
        tc = self.layerThermalCond  # vector of layer thermal conductivities (W m-1 K-1)
        d = self.layerThickness     # vector of layer thicknesses (m)

        fimp = 0.5                  # implicit coefficient
        fexp = 0.5                  # explicit coefficient
        num = len(t)                # number of layers

        # Mean thermal conductivity over distance between 2 layers (W/mK)
        tcp = [0 for x in range(num)]
        # Thermal capacity times layer depth (J/m2K)
        hcp = [0 for x in range(num)]
        # lower, main, and upper diagonals
        za = [[0 for y in range(3)] for x in range(num)]
        # RHS
        zy = [0 for x in range(num)]

        # --------------------------------------------------------------------------
        # Define the column vectors for heat capactiy and conductivity
        hcp[0] = hc[0] * d[0]
        for j in range(1, num):
            tcp[j] = 2. / (d[j-1] / tc[j-1] + d[j] / tc[j])
            hcp[j] = hc[j] * d[j]

        # --------------------------------------------------------------------------
        # Define the first row of za matrix, and RHS column vector
        za[0][0] = 0.
        za[0][1] = hcp[0] / dt + fimp * tcp[1]
        za[0][2] = -fimp * tcp[1]
        zy[0] = hcp[0] / dt*t[0] - fexp * tcp[1] * (t[0] - t[1]) + flx1

        # --------------------------------------------------------------------------
        # Define other rows
        for j in range(1, num - 1):
            za[j][0] = fimp * (-tcp[j])
            za[j][1] = hcp[j] / dt + fimp * (tcp[j] + tcp[j+1])
            za[j][2] = fimp * (-tcp[j+1])
            zy[j] = hcp[j] / dt * t[j] + fexp * \
                (tcp[j] * t[j - 1] - tcp[j] * t[j] - tcp[j+1] * t[j] + tcp[j+1] * t[j+1])

        # --------------------------------------------------------------------------
        # Boundary conditions
        if is_near_zero(bc - 1.0):
            # heat flux
            za[num-1][0] = fimp * (-tcp[num-1])
            za[num-1][1] = hcp[num-1] / dt + fimp * tcp[num-1]
            za[num-1][2] = 0.
            zy[num-1] = hcp[num-1] / dt * t[num-1] + fexp * tcp[num-1] * \
                (t[num-2] - t[num-1]) + flx2
        elif is_near_zero(bc - 2.0):
            # deep-temperature
            za[num-1][0] = 0.
            za[num-1][1] = 1.
            za[num-1][2] = 0.
            zy[num-1] = temp2
        else:
            raise Exception('Error during conduction calculation. Check input '
                            'parameters in the Conduction routine.')

        # --------------------------------------------------------------------------
        zx = Element.invert(num, za, zy)
        # t(:) = zx(:)
        return zx  # return zx as 1d vector of templayers

    def qsat(self, temp, pres, parameter):
        """Calculate vector of saturation humidity from air pressure and layer temperatures.

        Args:
            temp: List of layer temperatures [K].
            pres: Pressure (at current timestep) [Pa].
            parameter: Parameter object with geographic parameters.

        Returns:
            List of saturation humidity values.
        """
        gamw = (parameter.cl - parameter.cpv) / parameter.rv
        betaw = (parameter.lvtt / parameter.rv) + (gamw * parameter.tt)
        alpw = math.log(parameter.estt) + (betaw / parameter.tt) + \
            (gamw * math.log(parameter.tt))
        work2 = parameter.r / parameter.rv
        foes_lst = [0 for i in range(len(temp))]
        work1_lst = [0 for i in range(len(temp))]
        qsat_lst = [0 for i in range(len(temp))]

        for i in range(len(temp)):
            # saturation vapor pressure
            foes_lst[i] = math.exp(alpw - (betaw / temp[i]) - (gamw * math.log(temp[i])))
            work1_lst[i] = foes_lst[i] / pres[i]
            # saturation humidity
            qsat_lst[i] = work2 * work1_lst[i] / (1. + ((work2 - 1.) * work1_lst[i]))

        return qsat_lst

    @staticmethod
    def invert(nz, A, C):
        """Invert and solve tridiagonal matrix.

        Given A * X = C, solves for X, where
            a(*, 1) lower diagonal (Ai, i-1)
            a(*, 2) principal diagonal (Ai, i)
            a(*, 3) upper diagonal (Ai, i+1)

        Args:
            nz: number of layers.
            A: A matrix to invert.
            C: Result from A * X

        Returns:
            Solved matrix X from matrix inversion.
        """

        X = [0 for i in range(nz)]

        for i in reversed(range(nz-1)):
            C[i] = C[i] - A[i][2] * C[i+1] / A[i+1][1]
            A[i][1] = A[i][1] - A[i][2] * A[i+1][0] / A[i+1][1]

        for i in range(1, nz, 1):
            C[i] = C[i] - A[i][0] * C[i-1] / A[i-1][1]

        for i in range(nz):
            X[i] = C[i] / A[i][1]

        return X

    def __repr__(self):
        rval = round(sum([1.0 / t for t in self.layerThermalCond]), 2)
        return 'Element,\n name: {}\n emissivity: {}\n albedo: {}\n R value: ' \
            '{}\n vegCoverage: {}\n horizontal: {}\n layerThickness: {}\n ' \
            'layerTemp: {}'.format(
                self.name, self.emissivity, self.albedo, rval, self.vegCoverage,
                bool(self.horizontal), self.layerThickness, self.layerTemp)
