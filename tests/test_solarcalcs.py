try:
    range = xrange
except NameError:
    pass

import pytest
import uwg
import os
import math
import pprint
from .test_base import TestBase

class TestSolarCalcs(TestBase):

    def test_solarangles(self):
        """ test solar angles """

        self.setup_uwg_integration()
        self.uwg.read_epw()
        self.uwg.set_input()
        self.uwg.init_BEM_obj()
        self.uwg.init_input_obj()


        solar = uwg.SolarCalcs(self.uwg.UCM, self.uwg.BEM, self.uwg.simTime, self.uwg.RSM, \
            self.uwg.forc, self.uwg.geoParam, self.uwg.rural)

        #timestep every 5 minutes (300s)
        for i in range(int(12*24*1.5)): #1.5 days, Jan 2nd, 12 noon
            solar.simTime.UpdateDate()
        solar.solarangles()

        # test simtime for solar after 1.5 days
        assert solar.ut == pytest.approx(12.0,abs=1e-15)
        assert solar.simTime.day == pytest.approx(2.0,abs=1e-15)
        assert solar.ad == pytest.approx(0.197963373,abs=1e-8)

        # recalculate solar angles with new time simulation, add to previous time increment
        for i in range(12*24*20 + 12*1 + 6): # Increment time to 22 days, 13hrs, 30min = Jan 22 at 1330
            solar.simTime.UpdateDate()

        # Run simulation
        solar.solarangles()

        uwg_python_val = [
        solar.simTime.month,
        solar.simTime.secDay,
        solar.simTime.day,
        solar.simTime.hourDay,
        solar.ut,
        solar.ad,
        solar.eqtime,
        solar.decsol,
        solar.zenith,
        solar.tanzen,
        solar.critOrient
        ]

        # open matlab ref file
        uwg_matlab_val = self.setup_open_matlab_ref("matlab_solarcalcs","matlab_ref_solarangles.txt")

        # matlab ref checking
        assert len(uwg_matlab_val) == len(uwg_python_val)
        for i in range(len(uwg_matlab_val)):
            #print uwg_python_val[i], uwg_matlab_val[i]
            assert uwg_python_val[i] == pytest.approx(uwg_matlab_val[i], abs=1e-15), "error at index={}".format(i)

    def test_solarcalcs(self):
        """ test solar calculation """

        self.setup_uwg_integration()
        self.uwg.read_epw()
        self.uwg.set_input()
        self.uwg.init_BEM_obj()
        self.uwg.init_input_obj()

        # We subtract 11 hours from total timestep so
        # we can stop simulation while we still have sun!
        # New time: Jan 31, 1300
        self.uwg.simTime.nt -= 12*11

        self.uwg.hvac_autosize()
        self.uwg.simulate()

        # check date
        assert self.uwg.simTime.month == 1
        assert self.uwg.simTime.day == 31
        assert self.uwg.simTime.secDay/3600. == pytest.approx(13.0,abs=1e-15)

        # Manual checking

        # Total Radiation (W m-2) at 13:00 on Jan 31 w/o reflected sunlight
        assert self.uwg.solar.bldSol == pytest.approx(170.4189961148073, abs=1e-12)
        assert self.uwg.solar.roadSol == pytest.approx(307.4032644249752, abs=1e-12)

        # Matlab Checking
        uwg_matlab_val = self.setup_open_matlab_ref("matlab_solarcalcs","matlab_ref_solarcalcs.txt")

        uwg_python_val = [
        self.uwg.solar.horSol,
        self.uwg.solar.Kw_term,
        self.uwg.solar.Kr_term,
        self.uwg.solar.mr,
        self.uwg.solar.mw,
        self.uwg.solar.BEM[0].roof.solRec,
        self.uwg.solar.BEM[0].wall.solRec,
        self.uwg.solar.BEM[1].roof.solRec,
        self.uwg.solar.BEM[1].wall.solRec,
        self.uwg.solar.rural.solRec,
        self.uwg.solar.UCM.SolRecRoof,
        self.uwg.solar.UCM.SolRecRoad,
        self.uwg.solar.UCM.SolRecWall,
        self.uwg.solar.UCM.treeSensHeat,
        self.uwg.solar.UCM.treeLatHeat
        ]

        # matlab ref checking
        assert len(uwg_matlab_val) == len(uwg_python_val)
        for i in range(len(uwg_matlab_val)):
            #print uwg_python_val[i], uwg_matlab_val[i]
            assert uwg_python_val[i] == pytest.approx(uwg_matlab_val[i], abs=1e-15), "error at index={}".format(i)


if __name__ == "__main__":
    tsc = TestSolarCalcs()
    tsc.test_solarangles()
    tsc.test_solarcalcs()
