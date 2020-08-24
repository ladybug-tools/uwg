from __future__ import division, print_function
import os
import pytest
import uwg
import math
from .test_base import TestBase

try:
    import cPickle as pickle
except ImportError:
    import pickle

try:
    range = xrange
except NameError:
    pass

from pprint import pprint
from decimal import Decimal
pp = pprint
dd = Decimal.from_float


class TestDragonfly(TestBase):
    """Test for uwg from df
    """
    BUILDINGPARAMS = ['floorHeight',
                     'intHeat',
                     'intHeatNight',
                     'intHeatDay',
                     'intHeatFRad',
                     'intHeatFLat',
                     'infil',
                     'vent',
                     'glazingRatio',
                     'uValue',
                     'shgc',
                     'condType',
                     'cop',
                     'coolSetpointDay',
                     'coolSetpointNight',
                     'heatSetpointDay',
                     'heatSetpointNight',
                     'coolCap',
                     'heatEff',
                     'mSys',
                     'indoorTemp',
                     'indoorHum',
                     'heatCap',
                     'copAdj',
                     'canyon_fraction',
                     'Type',
                     'Era',
                     'Zone']

    @staticmethod
    def check_obj_attr(obj1, obj2, attr_lst):
        for attr in attr_lst:
            v1, v2 = [getattr(obj, attr) for obj in [obj1, obj2]]

            # if v1 != v2:
            #     print("{}: {} != {}".format(attr, v1, v2))
            #     is_equal = False
            # else:
            #     print("{}: {} == {}".format(attr, v1, v2))

            if isinstance(v1, (int, float)) and isinstance(v2, (int, float)):
                assert v1 == pytest.approx(v2, abs=1e-10)
            else:
                assert v1 == v2

    @staticmethod
    def unpickle(cpickle_path):
        with open(cpickle_path, "rb") as f:
            try:
                return pickle.load(f)
            except EOFError:
                return None

    def uwg_from_df(self):
        # spelling mistake
        # lathAnth = in DF = latAnth
        self.setup_uwg_integration("SGP_Singapore.486980_IWEC.epw", "initialize_fatal_error.uwg")
        self.uwg.set_input()
        bem = self.unpickle(os.path.join(self.uwg.uwgParamDir, "initialize_fatal_error_bem.pkl"))
        self.uwg.init_BEM_obj()
        self.uwg.BEM = bem

        # Error on Aug 22 if dt at 300
        self.uwg.Month = 8
        self.uwg.Day = 1
        self.uwg.nDay = 31

        self.uwg.read_epw()
        self.uwg.init_input_obj()
        self.uwg.hvac_autosize()

        return self.uwg

    def uwg_manual(self):
        # spelling mistake
        # lathAnth = in DF = latAnth
        self.setup_uwg_integration("SGP_Singapore.486980_IWEC.epw", "initialize_fatal_error.uwg")
        self.uwg.set_input()
        self.uwg.init_BEM_obj()

        # # Add custom typology parameters from DF
        df_typology_params = [[4.0, 0.5, 0.65, 0.35, 0.08, 0.5, 0.5],
                              [3.05, 0.5, 0.1499, 0.54, 0.15, 0.2, 0.2],
                              [4.0, 0.5, 0.0058, 0.251, 0.08, 0.2, 0.2]]

        # lathAnth type BEMDef: Type = LargeOffice, Zone = 1A (Miami), Era = Pst80, Construction = MassWall
        # type BEMDef: Type = MidRiseApartment, Zone = 1A (Miami), Era = Pre80, Construction = SteelFrame
        # type BEMDef: Type = WareHouse, Zone = 1A (Miami), Era = Pst80, Construction = MassWall

        for i in range(len(self.uwg.BEM)):
            self.uwg.BEM[i].building.floorHeight = df_typology_params[i][0]
            self.uwg.BEM[i].building.canyon_fraction = df_typology_params[i][1]
            self.uwg.BEM[i].building.glazingRatio = df_typology_params[i][2]
            self.uwg.BEM[i].building.shgc = df_typology_params[i][3]
            self.uwg.BEM[i].wall.albedo = df_typology_params[i][4]
            self.uwg.BEM[i].roof.albedo = df_typology_params[i][5]
            self.uwg.BEM[i].roof.vegCoverage = df_typology_params[i][6]

        # Error on Aug 22 if dt at 300
        self.uwg.Month = 8
        self.uwg.Day = 1
        self.uwg.nDay = 31

        self.uwg.read_epw()
        self.uwg.init_input_obj()
        self.uwg.hvac_autosize()

        return self.uwg

    def test_compare(self):
        uwgpkl = self.uwg_from_df()
        uwgfile = self.uwg_manual()

        for i in range(3):
            self.check_obj_attr(uwgpkl.BEM[i].building, uwgfile.BEM[i].building, self.BUILDINGPARAMS)

    def test_fatal_error_manual(self):
        # These will fail if simTime.dt at 300 o Aug 22nd
        # Will raise the FATAL ERROR exception
        with pytest.raises(Exception):
            uwgpkl = self.uwg_manual()
            uwgpkl.simulate()

    def test_fatal_error_df(self):
        # These will fail if simTime.dt at 300 o Aug 22nd
        # Will raise the FATAL ERROR exception
        with pytest.raises(Exception):
            uwgpkl = self.uwg_from_df()
            uwgpkl.simulate()


if __name__ == "__main__":
    pass
    #test = TestDragonfly()
    #test.test_compare()
    #test.test_fatal_error_df()
    #test.test_fatal_error_manual()
