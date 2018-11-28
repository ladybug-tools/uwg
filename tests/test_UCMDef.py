try:
    range = xrange
except NameError:
    pass

import pytest
import uwg
import os
import math
import pprint
import decimal
from .test_base import TestBase

dd = decimal.Decimal.from_float

class TestUCMDef(TestBase):

    def test_ucm_init(self):
        """ test ucm constructor """

        self.setup_uwg_integration()
        self.uwg.read_epw()
        self.uwg.set_input()
        self.uwg.init_BEM_obj()
        self.uwg.init_input_obj()

        # Because we are not using latAnthrop, this value can be None in uwg.
        # So we will hardcode value for now.
        self.uwg.UCM.latAnthrop = 2.0

        # Get uwg_python values
        uwg_python_val = [
            self.uwg.UCM.h_mix,         # waste heat mix into canyon ratio
            self.uwg.UCM.bldHeight,     # average building height (m)
            self.uwg.UCM.verToHor,      # vertical-to-horizontal urban area ratio
            self.uwg.UCM.bldDensity,    # horizontal building density (footprint)
            self.uwg.UCM.treeCoverage,  # horizontal tree coverage (footprint)
            self.uwg.UCM.sensAnthrop,   # sensible heat anthropogenic (W m-2)
            self.uwg.UCM.latAnthrop,    # latent heat anthropogenic (W m-2)
            self.uwg.UCM.roadShad,      # shadow fraction of road
            self.uwg.UCM.bldWidth,      # building width (m)
            self.uwg.UCM.canWidth,      # canyon width (m)
            self.uwg.UCM.canAspect,     # canyon aspect ratio
            self.uwg.UCM.roadConf,      # road sky view factor
            self.uwg.UCM.wallConf,      # wall sky view factor
            self.uwg.UCM.facArea,       # facade area (m2)
            self.uwg.UCM.roadArea,      # road area (m2)
            self.uwg.UCM.roofArea,      # roof area (m2)
            self.uwg.UCM.canTemp,       # canyon air temperature (db) (K)
            self.uwg.UCM.roadTemp,      # avg road temperature
            self.uwg.UCM.canHum,        # canyon specific humidity (kgv kgda-1)
            self.uwg.UCM.ublWind,       # urban boundary layer wind velocity (m s-1)
            self.uwg.UCM.canWind,       # urban canyon wind velocity (m s-1)
            self.uwg.UCM.ustar,         # friction velocity (m s-1)
            self.uwg.UCM.ustarMod,      # modified friction velocity (m s-1)
            self.uwg.UCM.z0u,           # urban roughness length (m)
            self.uwg.UCM.l_disp,        # urban displacement length (m)
            self.uwg.UCM.alb_wall,      # albedo of wall
            self.uwg.UCM.facAbsor,      # absorptivity of wall
            self.uwg.UCM.roadAbsor,     # average road absorptivity
            self.uwg.UCM.sensHeat       # urban sensible heat [W m-2]
        ]

        uwg_matlab_val = self.setup_open_matlab_ref("matlab_ucm","matlab_ref_ucm_init.txt")


        # matlab ref checking
        assert len(uwg_matlab_val) == len(uwg_python_val)
        for i in range(len(uwg_matlab_val)):
            #print uwg_python_val[i], uwg_matlab_val[i]
            tol = self.CALCULATE_TOLERANCE(uwg_python_val[i],15.0)
            assert uwg_python_val[i] == pytest.approx(uwg_matlab_val[i], tol), "error at index={}".format(i)

    def test_ucm_ucmodel(self):
        """ test ucm ucmodel """

        self.setup_uwg_integration()
        self.uwg.read_epw()
        self.uwg.set_input()

        # Test Jan 1 (winter, no vegetation coverage)
        self.uwg.Month = 1
        self.uwg.Day = 1
        self.uwg.nDay = 1

        # set_input
        self.uwg.init_BEM_obj()
        self.uwg.init_input_obj()


        # In order to avoid integration effects. Test only first time step
        # Subtract timestep to stop at 300 sec
        self.uwg.simTime.nt -= (23*12 + 11)

        # Run simulation
        self.uwg.hvac_autosize()
        self.uwg.simulate()

        # check date
        #print self.uwg.simTime
        assert self.uwg.simTime.month == 1
        assert self.uwg.simTime.day == 1
        assert self.uwg.simTime.secDay == pytest.approx(300.0,abs=1e-15)

        # Get uwg values
        uwg_python_val = [
            # heat load building
            self.uwg.UCM.Q_wall,        # convective sensible heat flux from building roof
            self.uwg.UCM.Q_window,      # sensible heat flux from window (via U-factor)
            self.uwg.UCM.Q_hvac,        # sensible heat flux from HVAC waste
            self.uwg.UCM.ElecTotal,     # total electricity consumption of urban area
            self.uwg.UCM.GasTotal,      # total gas consumption of urban area
            # heat load road/canyon/ubl
            self.uwg.UCM.Q_road,        # convective sensible heat flux from road
            self.uwg.UCM.Q_traffic,     # net sensible heat flux from traffic
            self.uwg.UCM.Q_vent,        # convective heat exchange from ventilation/infiltration
            self.uwg.UCM.Q_ubl,         # convective heat exchange with UBL layer
            # Building temperature
            self.uwg.UCM.roofTemp,      # average roof temp (K)
            self.uwg.UCM.wallTemp,      # average wall temp (K)
            # Sensible heat
            self.uwg.UCM.treeSensHeat,  # sensible heat from trees
            self.uwg.UCM.sensHeat,      # urban sensible heat
            self.uwg.UCM.canTemp        # canyon air temp (K)
        ]

        # Get uwg_matlab values
        uwg_matlab_val = self.setup_open_matlab_ref("matlab_ucm","matlab_ref_ucm_ucmodel.txt")

        # matlab ref checking
        assert len(uwg_matlab_val) == len(uwg_python_val)

        for i in range(len(uwg_matlab_val)):
            tol = self.CALCULATE_TOLERANCE(abs(uwg_python_val[i]),15.0)
            tol = 10**(math.log10(tol)+2) # Lowering tolerance by two orders of magnitude
            #print '-'*int(-(-15-math.log10(tol)))+'.123456789012345', tol
            assert uwg_python_val[i] == pytest.approx(uwg_matlab_val[i], abs=tol), "error at index={}".format(i)


if __name__ == "__main__":
    tucm = TestUCMDef()
    #tucm.test_ucm_init()
    tucm.test_ucm_ucmodel()
