"""Class for rural site and Vertical Diffusion Model."""
from __future__ import division, print_function

try:
    range = xrange
except NameError:
    pass

import os
from math import sqrt, pow, log
from .utilities import is_near_zero


class RSMDef(object):
    """Rural Site & Vertical Diffusion Model (VDM).

    This class calculates the vertical profiles of air temperature above the weather
    station[1].

    Note:
        [1] 'The uwg' (2012) Eq. 4, 5, 6.

    Args:
        lat: Number for latitude in degrees.
        lon: Number for longitude in degrees.
        GMT: Number for GMT hour correction.
        height: Number for rural average obstacle height in meters.
        T_init: Number for initial dry bulb temperature.
        P_init: Number for initial pressure.
        parameter: Param object for geographic parameters.
        z_meso_path: Text string for height reference data filepath for mesoscale
            component.

    Properties
        * lat
        * lon
        * GMT
        * height
        * z0r
        * disp
        * z
        * dz
        * nz0
        * nzref
        * nzfor
        * nz10
        * nzi
        * tempProf
        * presProf
        * tempRealProf
        * densityProfC
        * densityProfS
        * windProf
        * ublPres
    """

    def __init__(self, lat, lon, GMT, height, T_init, P_init, parameter, z_meso_path):

        # z_meso: list of mesoscale heights
        self.z_meso = RSMDef.load_z_meso(z_meso_path)
        self.lat = lat  # latitude [deg]
        self.lon = lon  # longitude [deg]
        self.GMT = GMT  # GMT hour correction
        self.height = height  # average obstacle height [m]
        self.z0r = 0.1 * height  # rural roughness length [m]
        self.disp = 0.5 * height  # rural displacement lenght [m]

        # vertical grid at the rural site
        # z: vertical height [m]
        self.z = [0 for x in range(len(self.z_meso)-1)]
        # dz: vertical discretization [m]
        self.dz = [0 for x in range(len(self.z_meso)-1)]

        # Initialize empty parameters to compute for RSMDef
        self.nz0 = None  # layer number at zmt [m]
        self.nzref = None  # layer number at zref [m]
        self.nzfor = None  # layer number at zfor [m]
        self.nz10 = None  # layer number at zmu [m]
        self.nzi = None  # layer number at zi_d [m]
        # tempProf: potential temperature profile at the rural site [K]
        self.tempProf = None
        # presProf: pressure profile at the rural site [Pa]
        self.presProf = None
        # tempRealProf: real temperature profile at the rural site [K]
        self.tempRealProf = None
        # densityProfC: density profile at the center of layers [kg m-3]
        self.densityProfC = None
        # densityProfS: density profile at the sides of layers [kg m-3]
        self.densityProfS = None
        # windProf: wind profile at the rural site [m s-1]
        self.windProf = None
        self.ublPres = None  # Average pressure at UBL [Pa]

        for zi in range(len(self.z_meso)-1):
            self.z[zi] = 0.5 * (self.z_meso[zi] + self.z_meso[zi+1])
            self.dz[zi] = self.z_meso[zi+1] - self.z_meso[zi]

        # Define initial booleans
        ll = True
        mm = True
        nn = True
        oo = True
        pp = True

        # Define self.nz0, self.nzref, self.nzfor, self.nz10, self.nzi
        for iz in range(len(self.z_meso)-1):
            # self.nz0: self.z index >= reference height for weather station
            eq_th = is_near_zero(self.z[iz] - parameter.tempHeight)
            if (eq_th is True or self.z[iz] > parameter.tempHeight) and ll is True:
                self.nz0 = iz + 1  # layer number at zmt (m)
                ll = False

            # self.nzref: self.z index >= reference inversion height
            eq_rh = is_near_zero(self.z[iz] - parameter.refHeight)
            if (eq_rh is True or self.z[iz] > parameter.refHeight) and mm is True:
                self.nzref = iz + 1  # layer number at zref (m)
                mm = False

            # self.nzfor: self.z index >= nighttime boundary layer height
            eq_nh = is_near_zero(self.z[iz] - parameter.nightBLHeight)
            if (eq_nh is True or self.z[iz] > parameter.nightBLHeight) and nn is True:
                self.nzfor = iz + 1   # layer number at zfor (m)
                nn = False

            # self.nz10: self.z index >= wind height
            eq_wh = is_near_zero(self.z[iz] - parameter.windHeight)
            if (eq_wh is True or self.z[iz] > parameter.windHeight) and oo is True:
                self.nz10 = iz + 1  # layer number at zmu (m)
                oo = False

            eq_dh = is_near_zero(self.z[iz] - parameter.dayBLHeight)
            if (eq_dh is True or self.z[iz] > parameter.dayBLHeight) and pp is True:
                self.nzi = iz + 1  # layer number at zi_d (m)
                pp = False

        # Define temperature, pressure and density vertical profiles
        self.tempProf = [T_init for x in range(self.nzref)]
        self.presProf = [P_init for x in range(self.nzref)]

        for iz in range(1, self.nzref):
            self.presProf[iz] = \
                (self.presProf[iz-1] ** (parameter.r / parameter.cp) - parameter.g /
                 parameter.cp * (P_init ** (parameter.r / parameter.cp)) *
                 (1. / self.tempProf[iz] + 1. / self.tempProf[iz-1]) * 0.5 *
                 self.dz[iz]) ** (1. / (parameter.r / parameter.cp))

        self.tempRealProf = [T_init for x in range(self.nzref)]
        for iz in range(self.nzref):
            self.tempRealProf[iz] = \
                self.tempProf[iz] * (self.presProf[iz] / P_init) ** \
                (parameter.r / parameter.cp)

        self.densityProfC = [None for x in range(self.nzref)]
        for iz in range(self.nzref):
            self.densityProfC[iz] = \
                self.presProf[iz] / parameter.r / self.tempRealProf[iz]

        self.densityProfS = [self.densityProfC[0] for x in range(self.nzref + 1)]

        for iz in range(1, self.nzref):
            self.densityProfS[iz] = \
                (self.densityProfC[iz] * self.dz[iz-1] + self.densityProfC[iz-1] *
                 self.dz[iz]) / (self.dz[iz-1]+self.dz[iz])

        self.densityProfS[self.nzref] = self.densityProfC[self.nzref-1]
        self.windProf = [1 for x in range(self.nzref)]

    def VDM(self, forc, rural, parameter, simTime):
        """# Ref: The uwg (2012), Eq. (4)"""

        self.tempProf[0] = forc.temp    # Lower boundary condition

        # compute pressure profile
        for iz in reversed(list(range(self.nzref))[1:]):
            self.presProf[iz-1] = \
                (pow(self.presProf[iz], parameter.r / parameter.cp) +
                 parameter.g / parameter.cp *
                 (pow(forc.pres, parameter.r / parameter.cp)) *
                 (1./self.tempProf[iz] + 1./self.tempProf[iz-1]) *
                 0.5 * self.dz[iz]) ** (1. / (parameter.r / parameter.cp))

        # compute the real temperature profile
        for iz in range(self.nzref):
            self.tempRealProf[iz] = self.tempProf[iz] * \
                (self.presProf[iz] / forc.pres) ** (parameter.r / parameter.cp)

        # compute the density profile
        for iz in range(self.nzref):
            self.densityProfC[iz] = self.presProf[iz] / parameter.r / \
                self.tempRealProf[iz]
        self.densityProfS[0] = self.densityProfC[0]

        for iz in range(1, self.nzref):
            self.densityProfS[iz] = \
                (self.densityProfC[iz] * self.dz[iz-1] + self.densityProfC[iz-1] *
                 self.dz[iz]) / (self.dz[iz-1] + self.dz[iz])

        self.densityProfS[self.nzref] = self.densityProfC[self.nzref-1]

        # Ref: The uwg (2012), Eq. (5)
        # compute diffusion coefficient
        cd, ustarRur = self.DiffusionCoefficient(
            self.densityProfC[0], self.z, self.dz, self.z0r, self.disp, self.tempProf[0],
            rural.sens, self.nzref, forc.wind, self.tempProf, parameter)

        # solve diffusion equation
        self.tempProf = RSMDef.DiffusionEquation(
            self.nzref, simTime.dt, self.tempProf, self.densityProfC, self.densityProfS,
            cd, self.dz)

        # compute wind profile
        # N.B In Matlab, negative values are converted to complex values.
        # log(-x) = log(x) + log(-1) = log(x) + i*pi
        # Python will throw an exception. Negative value occurs here if
        # VDM is run for average obstacle height ~ 4m.
        for iz in range(self.nzref):
            self.windProf[iz] = \
                ustarRur / parameter.vk * log((self.z[iz] - self.disp) / self.z0r)

        # Average pressure
        self.ublPres = 0.
        for iz in range(self.nzfor):
            self.ublPres = \
                self.ublPres + self.presProf[iz] * self.dz[iz] / \
                (self.z[self.nzref-1] + self.dz[self.nzref-1] / 2.)

    def DiffusionCoefficient(self, rho, z, dz, z0, disp, tempRur, heatRur, nz, uref, th,
                             parameter):
        # Initialization
        Kt = [0 for x in range(nz+1)]
        ws = [0 for x in range(nz)]
        te = [0 for x in range(nz)]
        # Friction velocity (Louis 1979)
        ustar = parameter.vk * uref/log((10.-disp)/z0)

        # Monin-Obukhov length
        lengthRur = max(-rho * parameter.cp * ustar ** 3 * tempRur / parameter.vk /
                        parameter.g / heatRur, -50.)

        if heatRur > 1e-2:
            # Unstable conditions

            # Convective velocity scale
            wstar = (parameter.g * heatRur * parameter.dayBLHeight / rho /
                     parameter.cp / tempRur) ** (1 / 3.)
            # Wind profile function
            phi_m = (1-8. * 0.1 * parameter.dayBLHeight / lengthRur) ** (-1. / 3.)

            for iz in range(nz):
                # Mixed-layer velocity scale
                ws[iz] = (ustar ** 3 + phi_m * parameter.vk * wstar ** 3 * z[iz] /
                          parameter.dayBLHeight) ** (1 / 3.)
                # TKE approximation
                te[iz] = max(ws[iz] ** 2., 0.01)

        else:
            # Stable and neutral conditions
            for iz in range(nz):
                # TKE approximation
                te[iz] = max(ustar ** 2., 0.01)

        # lenght scales (l_up, l_down, l_k, l_eps)
        self.dlu, self.dld = RSMDef.DissipationBougeault(parameter.g, nz, z, dz, te, th)
        self.dld, dls, dlk = RSMDef.LengthBougeault(nz, self.dld, self.dlu, z)

        # Boundary-layer diffusion coefficient
        for iz in range(nz):
            Kt[iz] = 0.4 * dlk[iz] * sqrt(te[iz])

        Kt[nz] = Kt[nz-1]

        return Kt, ustar

    @staticmethod
    def DiffusionEquation(nz, dt, co, da, daz, cd, dz):

        cddz = [0 for i in range(nz + 2)]
        a = [[0 for j in range(3)] for i in range(nz)]
        c = [0 for i in range(nz)]

        cddz[0] = daz[0] * cd[0] / dz[0]
        for iz in range(1, nz):
            cddz[iz] = 2. * daz[iz] * cd[iz] / (dz[iz] + dz[iz-1])
        cddz[nz] = daz[nz] * cd[nz] / dz[nz]

        a[0][0] = 0.
        a[0][1] = 1.
        a[0][2] = 0.
        c[0] = co[0]

        for iz in range(1, nz - 1):
            dzv = dz[iz]
            a[iz][0] = -cddz[iz] * dt / dzv / da[iz]
            a[iz][1] = 1 + dt * (cddz[iz] + cddz[iz+1]) / dzv / da[iz]
            a[iz][2] = -cddz[iz+1] * dt / dzv / da[iz]
            c[iz] = co[iz]

        a[nz-1][0] = -1.
        a[nz-1][1] = 1.
        a[nz-1][2] = 0.
        c[nz-1] = 0.

        return RSMDef.invert(nz, a, c)

    @staticmethod
    def DissipationBougeault(g, nz, z, dz, te, pt):
        # N.B on translation from UWG_Matlab
        # list length (i.e nz) != list indexing (i.e dlu[0] in python
        # wherease in matlab it is

        dlu = [0 for x in range(nz)]
        dld = [0 for x in range(nz)]

        for iz in range(nz):
            zup = 0.
            dlu[iz] = z[nz] - z[iz] - (dz[iz] / 2.)
            zzz = 0.
            zup_inf = 0.
            beta = g / pt[iz]

            for izz in range(iz, nz - 1):
                dzt = (dz[izz+1] + dz[izz]) / 2.
                zup = zup - beta * pt[iz] * dzt
                zup = zup + beta * (pt[izz+1] + pt[izz]) * dzt / 2.
                zzz = zzz + dzt

                _chkzero = is_near_zero(te[iz] - zup_inf)
                if (te[iz] < zup) and ((te[iz] > zup_inf) or _chkzero):
                    bbb = (pt[izz+1] - pt[izz]) / dzt

                    if not is_near_zero(bbb - 0.):
                        _t1 = sqrt(max(0., (beta * (pt[izz] - pt[iz])) ** 2. +
                                   2. * bbb * beta * (te[iz] - zup_inf)))
                        tl = (-beta * (pt[izz] - pt[iz]) + _t1) / bbb / beta
                    else:
                        tl = (te[iz] - zup_inf) / (beta * (pt[izz] - pt[iz]))
                    dlu[iz] = max(1., zzz - dzt + tl)
                zup_inf = zup

            zdo = 0.
            zdo_sup = 0.
            dld[iz] = z[iz] + dz[iz] / 2.
            zzz = 0.

            for izz in range(iz, 0, -1):
                dzt = (dz[izz-1] + dz[izz]) / 2.
                zdo = zdo + beta * pt[iz] * dzt
                zdo = zdo - beta * (pt[izz - 1] + pt[izz]) * dzt / 2.
                zzz = zzz + dzt

                _chkzero = is_near_zero(te[iz] - zdo_sup)
                if (te[iz] < zdo) and ((te[iz] > zdo_sup) or _chkzero):

                    bbb = (pt[izz] - pt[izz-1]) / dzt

                    if not is_near_zero(bbb - 0.):
                        tl = (beta * (pt[izz] - pt[iz]) +
                              sqrt(max(0., (beta * (pt[izz] - pt[iz])) ** 2. + 2. *
                                   bbb * beta * (te[iz] - zdo_sup)))) / bbb / beta
                    else:
                        tl = (te[iz] - zdo_sup) / (beta * (pt[izz] - pt[iz]))
                    dld[iz] = max(1., zzz - dzt + tl)
                zdo_sup = zdo

        return dlu, dld

    @staticmethod
    def LengthBougeault(nz, dld, dlu, z):

        dlg = [0 for x in range(nz)]
        dls = [0 for x in range(nz)]
        dlk = [0 for x in range(nz)]

        for iz in range(nz):
            dlg[iz] = (z[iz] + z[iz+1]) / 2.

        for iz in range(nz):
            dld[iz] = min(dld[iz], dlg[iz])
            dls[iz] = sqrt(dlu[iz] * dld[iz])
            dlk[iz] = min(dlu[iz], dld[iz])

        return dld, dls, dlk

    @staticmethod
    def invert(nz, A, C):
        """Inversion and resolution of a tridiagonal matrix A X = C.

        nz number of layers
        a(*,1) lower diagonal (Ai, i-1)
        a(*,2) principal diagonal (Ai, i)
        a(*,3) upper diagonal (Ai, i+1)
        """

        X = [0 for i in range(nz)]

        for i in reversed(range(nz - 1)):
            C[i] = C[i] - A[i][2] * C[i+1] / A[i+1][1]
            A[i][1] = A[i][1] - A[i][2] * A[i+1][0] / A[i+1][1]

        for i in range(1, nz, 1):
            C[i] = C[i] - A[i][0] * C[i-1] / A[i-1][1]

        for i in range(nz):
            X[i] = C[i] / A[i][1]

        return X

    @staticmethod
    def load_z_meso(z_meso_path):
        """Get list of mesoscale height reference data."""

        assert os.path.exists(z_meso_path), \
            "File: '{}' does not exist.".format(z_meso_path)

        z_meso = []
        with open(z_meso_path, 'r') as f:
            for txtline in f:
                # Strip white spaces, change to float
                z = float(''.join(txtline.split()))
                z_meso.append(z)

        return z_meso

    def __repr__(self):
        return 'RSM: obstacle ht = {}m, surface roughness length = {}m, ' \
            'displacement length = {}m'.format(self.height, self.z0r, self.disp)
