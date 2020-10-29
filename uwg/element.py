"""Element class."""
from __future__ import division

try:
    range = xrange
except NameError:
    pass

import math
from .utilities import is_near_zero, float_in_range, float_positive
from .material import Material

try:
    str = basestring
except NameError:
    pass


class Element(object):
    """Element object defines wall, and roof constructions.

    Args:
        albedo: Number for outer surface albedo.
        emissivity: Number for outer surface emissivity.
        layer_thickness_lst: List of layer thickness in meters.
        material_lst: List of Material objects in Element.
        vegcoverage: Number for fraction of surface vegetation coverage.
        t_init: Element initial temperature [K].
        horizontal: Boolean indicating if Element is horizontal or not (vertical).
        name: Text string for name of Element.

    Properties:
        * albedo -- outer surface albedo
        * emissivity -- outer surface emissivity.
        * layer_thickness_lst -- list of layer thickness [m]
        * material_lst -- material objects in Element.
        * layerThermalCond -- vector of layer thermal conductivities [W m-1 K-1]
        * layerVolHeat -- vector of layer volumetric heat [J m-3 K-1]
        * vegcoverage -- surface grass coverage
        * t_init -- element initial temperature [K].
        * layerTemp -- vector of layer temperatures [K]
        * waterStorage -- thickness of water film [m] for horizontal surfaces only
        * horizontal -- 1-horizontal, 0-vertical
        * solRec -- solar radiation received [W m-2]
        * infra -- net longwave radiation [W m-2]
        * lat -- surface latent heat flux [W m-2]
        * sens -- surface sensible heat flux [W m-2]
        * solAbs -- solar radiation absorbed [W m-2]
        * aeroCond -- convective heat transfer coefficient [-]
        * T_ext -- external surface temperature [K]
        * T_int -- internal surface temperature [K]
        * flux -- external surface heat flux [W m-2]
    """

    def __init__(self, albedo, emissivity, layer_thickness_lst, material_lst,
                 vegcoverage, t_init, horizontal, name):

        assert len(layer_thickness_lst) == len(material_lst), 'The number of layer ' \
            'thickness must match the number of layer materials. Got {} and {}, ' \
            'respectively.'.format(len(layer_thickness_lst), len(material_lst))

        self.albedo = albedo  #
        self.emissivity = emissivity  # outer surface emissivity.
        # list of layer thickness [m]
        self.layer_thickness_lst = layer_thickness_lst
        self.material_lst = material_lst  #
        self.vegcoverage = vegcoverage  #
        self.t_init = t_init  #
        self.horizontal = int(horizontal)  #
        self._name = name

        # layerThermaCond:
        self.layerThermalCond = [
            material.thermalcond for material in material_lst]
        # layerVolHeat:
        self.layerVolHeat = [material.volheat for material in material_lst]
        # layerTemp: #
        self.layerTemp = [t_init] * len(layer_thickness_lst)
        # waterStorage: #
        self.waterStorage = 0
        self.infra = 0  #
        self.sens = 0  #
        self.solRec = 0  #
        self.lat = 0  #
        self.solAbs = 0  #
        self.aeroCond = 0  #
        self.T_ext = 293  #
        self.T_int = 293  #
        self.flux = 0  #

    @property
    def albedo(self):
        """Get or set a value between 0 and 1 for outer surface albedo."""
        return self._albedo

    @albedo.setter
    def albedo(self, value):
        self._albedo = float_positive(value, 'albedo')

    @property
    def emissivity(self):
        """Get or set a value between 0 and 1 for outer surface emissivity."""
        return self._emissivity

    @emissivity.setter
    def emissivity(self, value):
        self._emissivity = float_positive(value, 'emissivity')

    @property
    def layer_thickness_lst(self):
        """Get or set list of thickness in meters of each Material in Element.

        The order of thickness should correspond to the order of the Material objects in
        material_lst.
        """
        return self._layer_thickness_lst

    @layer_thickness_lst.setter
    def layer_thickness_lst(self, value):
        assert all(v > 0 for v in value), 'Every value in layer_thickness_lst '
        'must be greater than 0.'
        self._layer_thickness_lst = value

    @property
    def material_lst(self):
        """Get or set list of Material objects in the element.

        The order of Material objects should correspond to the order of the thickness in
        layer_thickness_lst.
        """
        return self._material_lst

    @material_lst.setter
    def material_lst(self, value):
        assert all(isinstance(v, Material) for v in value), 'Every item in '
        'in material_lst must be a Material object.'
        self._material_lst = value

    @property
    def vegcoverage(self):
        """Get or set fraction of vegetation coverage on Element."""
        return self._vegcoverage

    @vegcoverage.setter
    def vegcoverage(self, value):
        self._vegcoverage = float_in_range(value, 0, 1, 'vegcoverage')

    @property
    def t_init(self):
        """Get or set initial temperature of Element [K]."""
        return self._t_init

    @t_init.setter
    def t_init(self, value):
        self._t_init = float_in_range(value, mi=0, input_name='t_init')

    @property
    def horizontal(self):
        """Get or set boolean value indicating if Element is horizontal or not."""
        return self._horizontal

    @horizontal.setter
    def horizontal(self, value):
        self._horizontal = bool(value)

    @property
    def name(self):
        """Get or set text string for name of Element."""
        return self._name

    @classmethod
    def from_dict(cls, data):
        """Create a Element object from a dictionary.

        Args:
            data: A Element dictionary following the format below.

        .. code-block:: python

            {
            "albedo": 0.08,
            "emissivity": 0.92,
            "layer_thickness_lst": [0.0254, 0.0508, 0.0254],
            "material_lst": material_dict_lst,  # List of material dictionaries
            "vegcoverage": 0.5,
            "t_init": 293,
            "horizontal": False,
            "name": MassWall
            }
        """
        assert data['type'] == 'Element', 'Expected ' \
            'Element dictionary. Got {}.'.format(data['type'])

        materials = [Material.from_dict(m) for m in data['material_lst']]
        return cls(data['albedo'], data['emissivity'], data['layer_thickness_lst'],
                   materials, data['vegcoverage'], data['t_init'],
                   data['horizontal'], data['name'])

    def to_dict(self):
        """Element dictionary representation."""
        base = {'type': 'Element'}
        base['albedo'] = self.albedo
        base['emissivity'] = self.emissivity
        base['layer_thickness_lst'] = self.layer_thickness_lst
        base['material_lst'] = [m.to_dict() for m in self.material_lst]
        base['vegcoverage'] = self.vegcoverage
        base['t_init'] = self.t_init
        base['horizontal'] = self.horizontal
        base['name'] = self.name
        return base

    def SurfFlux(self, forc, parameter, simTime, humRef, tempRef, windRef, boundCond,
                 intFlux):
        """Calculate net heat flux, and update element layer temperatures."""

        # Calculated per unit area [m^2]

        # dens: air density (kgd m-3)
        dens = forc.pres / (1000 * 0.287042 * tempRef *
                            (1. + 1.607858 * humRef))
        # convection coef (ref: uwg, eq. 12))
        self.aeroCond = 5.8 + 3.7 * windRef

        if self.horizontal:
            # For roof, mass, road
            if not is_near_zero(self.waterStorage) and self.waterStorage > 0.0:
                # Evaporation [m s-1], Film water & soil latent heat
                # Note: in the current uwg code, latent heat from evapotranspiration,
                # stagnant water, or anthropogenic sources is not modelled due to the
                # difficulty of validation, and lack of reliability of precipitation
                # data from EPW files.Therefore this condition is never run because all
                # elements have had their waterStorage hardcoded to 0.
                qtsat = self.qsat([self.layerTemp[0]], [
                                  forc.pres], parameter)[0]
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
                vegSen = 0.
            else:
                # Summer, veg
                self.solAbs = (
                    ((1.0 - self.vegcoverage) * (1. - self.albedo) + self.vegcoverage *
                     (1.0 - parameter.vegAlbedo)) * self.solRec)
                try:
                    # if road compute grass/tree fractions seperately
                    vegLat = (
                        self.grasscoverage * (1.0 - parameter.vegAlbedo) *
                        parameter.grassFLat * self.solRec)
                    vegLat += (
                        self.treecoverage * (1.0 - parameter.vegAlbedo) *
                        parameter.treeFLat * self.solRec)
                    vegSen = (
                        self.grasscoverage * (1.0 - parameter.vegAlbedo) *
                        (1.0 - parameter.grassFLat) * self.solRec)
                    vegSen += (
                        self.treecoverage * (1.0 - parameter.vegAlbedo) *
                        (1.0 - parameter.treeFLat) * self.solRec)
                except AttributeError:
                    # for all other Elements use veg fraction w/ grassFLat
                    vegLat = (
                        self.vegcoverage * (1.0 - parameter.vegAlbedo) *
                        parameter.grassFLat * self.solRec)
                    vegSen = (
                        self.vegcoverage * (1.0 - parameter.vegAlbedo) *
                        (1.0 - parameter.grassFLat) * self.solRec)

            self.lat = soilLat + vegLat

            # Sensible & net heat flux
            self.sens = vegSen + self.aeroCond * (self.layerTemp[0] - tempRef)
            self.flux = self.solAbs + self.infra - \
                self.lat - self.sens  # [W m-2]

        else:
            # For vertical surfaces (wall)
            self.solAbs = (1.0 - self.albedo) * self.solRec
            self.lat = 0.0

            # Sensible & net heat flux
            self.sens = self.aeroCond * (self.layerTemp[0] - tempRef)
            self.flux = self.solAbs + self.infra - \
                self.lat - self.sens  # [W m-2]

        self.layerTemp = \
            self.Conduction(simTime.dt, self.flux, boundCond,
                            forc.deepTemp, intFlux)
        self.T_ext = self.layerTemp[0]
        self.T_int = self.layerTemp[-1]

    def Conduction(self, dt, flx1, bc, temp2, flx2):
        """Solve the conductance of heat based on of the element layers.

        Args:
            dt: Simulation time step in seconds.
            flx1: Net heat flux on surface [W m-2]
            bc: Boundary condition parameter (1 or 2)
            temp2: Deep soil temperature (ave of air temperature) [K]
            flx2: Surface flux (sum of absorbed, emitted, etc.) [W m-2]

        Returns:
            A 1d vector of element layer temperatures.
        """
        t = self.layerTemp          # vector of layer temperatures (K)
        # vector of layer volumetric heat (J m-3 K-1)
        hc = self.layerVolHeat
        # vector of layer thermal conductivities (W m-1 K-1)
        tc = self.layerThermalCond
        d = self.layer_thickness_lst     # vector of layer thicknesses (m)

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
                (tcp[j] * t[j - 1] - tcp[j] * t[j] -
                 tcp[j+1] * t[j] + tcp[j+1] * t[j+1])

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
            foes_lst[i] = math.exp(
                alpw - (betaw / temp[i]) - (gamw * math.log(temp[i])))
            work1_lst[i] = foes_lst[i] / pres[i]
            # saturation humidity
            qsat_lst[i] = work2 * work1_lst[i] / \
                (1. + ((work2 - 1.) * work1_lst[i]))

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
            '{}\n vegcoverage: {}\n horizontal: {}\n layer_thickness_lst: {}\n ' \
            'layerTemp: {}'.format(
                self.name, self.emissivity, self.albedo, rval, self.vegcoverage,
                bool(self.horizontal), self.layer_thickness_lst, self.layerTemp)
