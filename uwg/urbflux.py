"""Function for surface heat flux calculations."""
from __future__ import division

try:
    range = xrange
except NameError:
    pass

from .infracalcs import infracalcs
from math import log


def urbflux(UCM, UBL, BEM, forc, parameter, simTime, RSM):
    """Calculate the surface heat fluxes.

    Args:
        UCM: UCMDef object.
        UBL: UBLDef object.
        BEM: List of BEMDef objects.
        forc: Forcing object.
        parameter: Param object.
        simTime: SimParam object.
        RSM: RSMDef object.

    Returns:
        Tuple with updated UCMDef, UBLDef, and BEM objects.
    """

    T_can = UCM.canTemp
    Cp = parameter.cp
    UCM.Q_roof = 0.
    sigma = 5.67e-8  # Stephan-Boltzman constant
    UCM.roofTemp = 0.  # Average urban roof temperature
    UCM.wallTemp = 0.  # Average urban wall temperature

    for j in range(len(BEM)):
        # Building energy model
        BEM[j].building.BEMCalc(UCM, BEM[j], forc, parameter, simTime)
        BEM[j].ElecTotal = BEM[j].building.ElecTotal * BEM[j].fl_area  # W m-2

        # Update roof infra calc
        e_roof = BEM[j].roof.emissivity
        T_roof = BEM[j].roof.layerTemp[0]
        BEM[j].roof.infra = e_roof * (forc.infra - sigma * T_roof ** 4.)

        # update wall infra calc (road done later)
        e_wall = BEM[j].wall.emissivity
        T_wall = BEM[j].wall.layerTemp[0]
        # calculates the infrared radiation for wall, taking into account radiation
        # exchange from road
        _infra_road_, BEM[j].wall.infra = \
            infracalcs(UCM, forc, UCM.road.emissivity, e_wall, UCM.roadTemp, T_wall)

        # Update element temperatures
        BEM[j].mass.layerTemp = \
            BEM[j].mass.Conduction(
                simTime.dt, BEM[j].building.fluxMass, 1., 0., BEM[j].building.fluxMass)
        BEM[j].roof.SurfFlux(
            forc, parameter, simTime, UCM.canHum, T_can, max(forc.wind, UCM.canWind), 1.,
            BEM[j].building.fluxRoof)
        BEM[j].wall.SurfFlux(
            forc, parameter, simTime, UCM.canHum, T_can, UCM.canWind, 1.,
            BEM[j].building.fluxWall)

        # Note the average wall & roof temperature
        UCM.wallTemp = UCM.wallTemp + BEM[j].frac * BEM[j].wall.layerTemp[0]
        UCM.roofTemp = UCM.roofTemp + BEM[j].frac * BEM[j].roof.layerTemp[0]

    # Update road infra calc (assume walls have similar emissivity, so use the last one)
    UCM.road.infra, _wall_infra = \
        infracalcs(UCM, forc, UCM.road.emissivity, e_wall, UCM.roadTemp, UCM.wallTemp)
    UCM.road.SurfFlux(forc, parameter, simTime, UCM.canHum, T_can, UCM.canWind, 2., 0.)
    UCM.roadTemp = UCM.road.layerTemp[0]

    # Sensible & latent heat flux (total)
    try:
        UCM.latHeat += \
            UCM.latAnthrop + UCM.treeLatHeat + UCM.road.lat * (1. - UCM.bldDensity)
    except TypeError:
        pass  # latheat is None

    # ---------------------------------------------------------------------
    # Advective heat flux to UBL from VDM
    #
    # Note: UWG_Matlab code here is modified to compensate for rounding errors
    # that occur when recursively adding forDens, intAdv1, and intAdv2.
    # This causes issues in the UBL.advHeat calculation when large (1e5)
    # numbers are subtracted to produce small numbers (1e-10) that can
    # differ from equivalent matlab calculations by a factor of 2.
    # Values this small are ~ 0, but for consistency's sake Kahan Summation
    # algorithm is applied to keep margin of difference from UWG_Matlab low.
    # ---------------------------------------------------------------------

    forDens = 0.0
    intAdv1 = 0.0
    intAdv2 = 0.0

    # c1 & c2 stores values truncated by floating point rounding for values < 10^-16
    c1 = 0.0
    c2 = 0.0
    c3 = 0.0

    for iz in range(RSM.nzfor):
        # At c loss of precision at at low order of magnitude, that we need in
        # UBL.advHeat calc. Algebraically t is 0, but with floating pt numbers
        # c will accumulate truncated values
        y = (RSM.densityProfC[iz] * RSM.dz[iz] / (RSM.z[RSM.nzfor - 1] +
             RSM.dz[RSM.nzfor-1] / 2.))
        t = forDens + y
        c1 += (t - forDens) - y
        forDens = t

        y = RSM.windProf[iz] * RSM.tempProf[iz] * RSM.dz[iz]
        t = intAdv1 + y
        c2 += (t - intAdv1) - y
        intAdv1 = t

        y = RSM.windProf[iz] * RSM.dz[iz]
        t = intAdv2 + y
        c3 += (t - intAdv2) - y
        intAdv2 = t

    # Add the truncated values back
    forDens -= c1
    intAdv1 -= c2
    intAdv2 -= c3
    UBL.advHeat = (
        UBL.paralLength * Cp * forDens * (intAdv1 - (UBL.ublTemp * intAdv2)) /
        UBL.urbArea)

    # ---------------------------------------------------------------------
    # Convective heat flux to UBL from UCM (see Appendix - Bueno (2014))
    # ---------------------------------------------------------------------
    zrUrb = 2 * UCM.bldHeight
    zref = RSM.z[RSM.nzref-1]  # Reference height

    # Reference wind speed & canyon air density
    windUrb = (
        forc.wind*log(zref / RSM.z0r) / log(parameter.windHeight / RSM.z0r) *
        log(zrUrb / UCM.z0u) / log(zref / UCM.z0u))
    dens = forc.pres / (1000 * 0.287042 * T_can * (1. + 1.607858 * UCM.canHum))

    # Friction velocity
    UCM.ustar = parameter.vk * windUrb / log((zrUrb - UCM.l_disp) / UCM.z0u)

    # Convective scaling velocity
    wstar = (parameter.g * max(UCM.sensHeat, 0.0) * zref / dens / Cp / T_can) ** (1 / 3.)
    UCM.ustarMod = max(UCM.ustar, wstar)  # Modified friction velocity
    UCM.uExch = parameter.exCoeff * UCM.ustarMod  # Exchange velocity

    # Canyon wind speed, Eq. 27 Chp. 3 Hanna and Britter, 2002
    # assuming CD = 1 and lambda_f = verToHor/4
    UCM.canWind = UCM.ustarMod * (UCM.verToHor / 8.) ** (-1 / 2.)

    # Canyon turbulent velocities
    UCM.turbU = 2.4 * UCM.ustarMod
    UCM.turbV = 1.9 * UCM.ustarMod
    UCM.turbW = 1.3 * UCM.ustarMod

    # Urban wind profile
    for iz in range(RSM.nzref):
        UCM.windProf.append(
            UCM.ustar / parameter.vk * log((RSM.z[iz] + UCM.bldHeight - UCM.l_disp) /
                                           UCM.z0u))

    return UCM, UBL, BEM
