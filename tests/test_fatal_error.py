"""Test for uwg from df"""
from __future__ import division, print_function

try:
    import cPickle as pickle
except ImportError:
    import pickle

import os
import pytest
from .test_base import auto_setup_uwg


TEST_DIR = os.path.abspath(os.path.dirname(__file__))
EPW_PATH = os.path.join(TEST_DIR, 'epw', 'SGP_Singapore.486980_IWEC.epw')
PARAM_DIR = os.path.join(TEST_DIR, 'parameters')
PARAM_PATH = os.path.join(PARAM_DIR, 'initialize_fatal_error.uwg')
BUILDINGPARAMS = [
    'floor_height', 'int_heat', 'int_heat_night', 'int_heat_day', 'int_heat_frad',
    'int_heat_flat', 'infil', 'vent', 'glazing_ratio', 'u_value', 'shgc',
    'condtype', 'cop', 'cool_setpoint_day', 'cool_setpoint_night', 'heat_setpoint_day',
    'heat_setpoint_night', 'coolcap', 'heateff', 'msys', 'indoor_temp', 'indoor_hum',
    'heat_cap', 'cop_adj', 'canyon_fraction']


def check_obj_attr(obj1, obj2, attr_lst):
    for attr in attr_lst:
        v1, v2 = [getattr(obj, attr) for obj in [obj1, obj2]]

        # if v1 != v2:
        #     print('{}: {} != {}'.format(attr, v1, v2))
        #     is_equal = False
        # else:
        #     print('{}: {} == {}'.format(attr, v1, v2))

        if isinstance(v1, (int, float)) and isinstance(v2, (int, float)):
            assert v1 == pytest.approx(v2, abs=1e-10)
        else:
            assert v1 == v2


def unpickle(cpickle_path):
    with open(cpickle_path, 'rb') as f:
        try:
            return pickle.load(f)
        except EOFError:
            return None


def uwg_from_df():

    # spelling mistake
    # lathAnth = in DF = latAnth
    testuwg = auto_setup_uwg(EPW_PATH, PARAM_PATH)
    # Error on Aug 22 if dt at 300
    testuwg.month = 8
    testuwg.day = 1
    testuwg.nday = 31

    bem = unpickle(os.path.join(PARAM_DIR, 'initialize_fatal_error_bem.pkl'))
    testuwg._compute_BEM()
    testuwg.BEM = bem

    testuwg._read_epw()
    testuwg._compute_input()
    testuwg._hvac_autosize()

    return testuwg


def uwg_manual():
    # spelling mistake
    # lathAnth = in DF = latAnth
    testuwg = auto_setup_uwg(EPW_PATH, PARAM_PATH)
    testuwg._compute_BEM()

    # # Add custom typology parameters from DF
    df_typology_params = [[4.0, 0.5, 0.65, 0.35, 0.08, 0.5, 0.5],
                          [3.05, 0.5, 0.1499, 0.54, 0.15, 0.2, 0.2],
                          [4.0, 0.5, 0.0058, 0.251, 0.08, 0.2, 0.2]]

    # lathAnth type BEMDef:
    # Type = LargeOffice, Zone = 1A (Miami), Era = Pst80, Construction = MassWall
    # type BEMDef:
    #   Type = MidRiseApartment, Zone = 1A (Miami), Era = Pre80, Construction=SteelFrame
    # type BEMDef:
    #   Type = WareHouse, Zone = 1A (Miami), Era = Pst80, Construction = MassWall

    for i in range(len(testuwg.BEM)):
        testuwg.BEM[i].building.floor_height = df_typology_params[i][0]
        testuwg.BEM[i].building.canyon_fraction = df_typology_params[i][1]
        testuwg.BEM[i].building.glazing_ratio = df_typology_params[i][2]
        testuwg.BEM[i].building.shgc = df_typology_params[i][3]
        testuwg.BEM[i].wall.albedo = df_typology_params[i][4]
        testuwg.BEM[i].roof.albedo = df_typology_params[i][5]
        testuwg.BEM[i].roof.vegcoverage = df_typology_params[i][6]

    # Error on Aug 22 if dt at 300
    testuwg.month = 8
    testuwg.day = 1
    testuwg.nday = 31

    testuwg._read_epw()
    testuwg._compute_input()
    testuwg._hvac_autosize()

    return testuwg


# TODO: Repickle once variable names have been finished PEP8 conversion
# def test_compare():
#     uwgpkl = uwg_from_df()
#     uwgfile = uwg_manual()

#     for i in range(3):
#         check_obj_attr(uwgpkl.BEM[i].building, uwgfile.BEM[i].building, BUILDINGPARAMS)


def test_fatal_error_manual():
    # These will fail if simTime.dt at 300 o Aug 22nd
    # Will raise the FATAL ERROR exception
    with pytest.raises(Exception):
        uwgpkl = uwg_manual()
        uwgpkl.simulate()


def test_fatal_error_df():
    # These will fail if simTime.dt at 300 o Aug 22nd
    # Will raise the FATAL ERROR exception
    with pytest.raises(Exception):
        uwgpkl = uwg_from_df()
        uwgpkl.simulate()
