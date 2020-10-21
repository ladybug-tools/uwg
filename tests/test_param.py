"""Tests for Param object."""

import os
from uwg import Param, UWG


def test_param_init():
    """Param init test."""

    epw_path = \
        os.path.join(os.path.dirname(__file__), 'epw',
                     'SGP_Singapore.486980_IWEC.epw')
    testuwg = UWG.from_param_args(
        10, 0.5, 0.5, 0.1, 0.1, '1A', epw_path=epw_path)

    nightStart = 18.
    nightEnd = 8.
    maxdx = 250.

    geoparam = Param(
        testuwg.h_ubl1, testuwg.h_ubl2, testuwg.h_ref, testuwg.h_temp, testuwg.h_wind,
        testuwg.c_circ, testuwg.maxday, testuwg.maxnight, testuwg.lattree,
        testuwg.latgrss, testuwg.albveg, testuwg.vegstart, testuwg.vegend, nightStart,
        nightEnd, testuwg.windmin, testuwg.WGMAX, testuwg.c_exch, maxdx, testuwg.G,
        testuwg.CP, testuwg.VK, testuwg.R, testuwg.RV, testuwg.LV, 3.14, testuwg.SIGMA,
        testuwg.WATERDENS, testuwg.LVTT, testuwg.TT, testuwg.ESTT, testuwg.CL,
        testuwg.CPV, testuwg.B, testuwg.CM, testuwg.COLBURN)

    geoparam.__repr__()
