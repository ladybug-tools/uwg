
import os
import pytest
import uwg
import math
from .test_base import TestBase

from pprint import pprint
from decimal import Decimal
pp = pprint
dd = Decimal.from_float

class TestUWG(TestBase):
    """Test for uwg.py
    """

    def test_read_epw(self):

        self.setup_uwg_integration("SGP_Singapore.486980_IWEC.epw", "initialize_singapore.uwg")
        self.uwg.read_epw()

        # test header
        assert self.uwg._header[0][0] == "LOCATION"
        assert self.uwg._header[0][1] == "SINGAPORE"
        assert self.uwg.lat == pytest.approx(1.37, abs=1e-3)
        assert self.uwg.lon == pytest.approx(103.98, abs=1e-3)
        assert self.uwg.GMT == pytest.approx(8, abs=1e-3)
        # test soil data
        assert self.uwg.nSoil == pytest.approx(3, abs=1e-2)
        # test soil depths
        assert self.uwg.depth_soil[0][0] == pytest.approx(0.5, abs=1e-3)
        assert self.uwg.depth_soil[1][0] == pytest.approx(2., abs=1e-3)
        assert self.uwg.depth_soil[2][0] == pytest.approx(4., abs=1e-3)
        # test soil temps over 12 months
        assert self.uwg.Tsoil[0][0] == pytest.approx(27.55+273.15, abs=1e-3)
        assert self.uwg.Tsoil[1][2] == pytest.approx(28.01+273.15, abs=1e-3)
        assert self.uwg.Tsoil[2][11] == pytest.approx(27.07+273.15, abs=1e-3)
        # test time step in weather file
        assert self.uwg.epwinput[0][0] == "1989"
        assert float(self.uwg.epwinput[3][6]) == pytest.approx(24.1,abs=1e-3)

    def test_read_input(self):

        self.setup_uwg_integration("SGP_Singapore.486980_IWEC.epw", "initialize_singapore.uwg")
        self.uwg.read_epw()
        self.uwg.set_input()
        self.uwg.init_BEM_obj()
        self.uwg.init_input_obj()

        #test uwg param dictionary first and last
        assert 'bldHeight' in self.uwg._init_param_dict
        assert 'h_obs' in self.uwg._init_param_dict

        assert self.uwg._init_param_dict['bldHeight'] == pytest.approx(10., abs=1e-6)
        assert self.uwg._init_param_dict['vegEnd'] == pytest.approx(10, abs=1e-6)
        assert self.uwg._init_param_dict['albRoof'] == None
        assert self.uwg._init_param_dict['h_ubl1'] == pytest.approx(1000., abs=1e-6)
        assert self.uwg._init_param_dict['h_ref'] == pytest.approx(150., abs=1e-6)

        # test SchTraffic schedule
        assert self.uwg._init_param_dict['SchTraffic'][0][0] == pytest.approx(0.2, abs=1e-6) # first
        assert self.uwg._init_param_dict['SchTraffic'][2][23] == pytest.approx(0.2, abs=1e-6) # last
        assert self.uwg._init_param_dict['SchTraffic'][0][19] == pytest.approx(0.8, abs=1e-6)
        assert self.uwg._init_param_dict['SchTraffic'][1][21] == pytest.approx(0.3, abs=1e-6)
        assert self.uwg._init_param_dict['SchTraffic'][2][6] == pytest.approx(0.4, abs=1e-6)

        # test bld fraction list
        assert self.uwg._init_param_dict['bld'][0][0] == pytest.approx(0., abs=1e-6)
        assert self.uwg._init_param_dict['bld'][3][1] == pytest.approx(0.4, abs=1e-6)
        assert self.uwg._init_param_dict['bld'][5][1] == pytest.approx(0.6, abs=1e-6)
        assert self.uwg._init_param_dict['bld'][15][2] == pytest.approx(0.0, abs=1e-6)

        # test BEMs
        assert len(self.uwg.BEM) == pytest.approx(2.,abs=1e-6)
        # test BEM office (BLD4 in DOE)
        assert self.uwg.BEM[0].building.Type == "LargeOffice"
        assert self.uwg.BEM[0].building.Zone == "1A (Miami)"
        assert self.uwg.BEM[0].building.Era == "Pst80"
        assert self.uwg.BEM[0].frac == 0.4

        # test BEM apartment
        assert self.uwg.BEM[1].building.Type == "MidRiseApartment"
        assert self.uwg.BEM[1].building.Zone == "1A (Miami)"
        assert self.uwg.BEM[1].building.Era == "Pst80"
        assert self.uwg.BEM[1].frac == 0.6

        # Check that schedules are called correctly
        assert self.uwg.Sch[0].Light[0][8] == pytest.approx(0.9, abs=1e-6)   #9am on Weekday for Office
        assert self.uwg.Sch[0].Light[0][7] == pytest.approx(0.3, abs=1e-6)   #9am on Weekday for Office
        assert self.uwg.Sch[1].Occ[1][11] == pytest.approx(0.25, abs=1e-6)     #12 noon on Weekend for apt

        # Check that soil ground depth is set correctly
        assert self.uwg.depth_soil[self.uwg.soilindex1][0] == pytest.approx(0.5, abs=1e-6)
        assert self.uwg.depth_soil[self.uwg.soilindex2][0] == pytest.approx(0.5, abs=1e-6)

        #self.road
        # Check the road layer splitting
        assert len(self.uwg.road.layerThickness) == pytest.approx(11., abs=1e-15)
        assert self.uwg.road.layerThickness[0] == pytest.approx(0.05, abs=1e-15)

        # Check the road layer splitting for rural
        assert len(self.uwg.rural.layerThickness) == pytest.approx(11., abs=1e-15)
        assert self.uwg.rural.layerThickness[0] == pytest.approx(0.05, abs=1e-6)

    def test_optional_blank_parameters(self):

        self.setup_uwg_integration("SGP_Singapore.486980_IWEC.epw", "initialize_singapore.uwg")

        # From blank inputs will be from DOE
        self.uwg.read_epw()
        self.uwg.set_input()
        self.uwg.init_BEM_obj()
        self.uwg.init_input_obj()

        assert self.uwg.BEM[0].building.glazingRatio == pytest.approx(0.38, abs=1e-15)
        assert self.uwg.BEM[0].roof.albedo == pytest.approx(0.2, abs=1e-15)
        assert self.uwg.BEM[0].roof.vegCoverage == pytest.approx(0.0, abs=1e-15)
        assert self.uwg.BEM[1].roof.albedo == pytest.approx(0.2, abs=1e-15)
        assert self.uwg.BEM[1].building.glazingRatio == pytest.approx(0.1499, abs=1e-15)
        assert self.uwg.BEM[1].roof.vegCoverage == pytest.approx(0.0, abs=1e-15)

    def test_optional_inputted_parameters(self):

        self.setup_uwg_integration("SGP_Singapore.486980_IWEC.epw", None)
        self.set_input_manually()

        self.uwg.albRoof = .5
        self.uwg.vegRoof = .1
        self.uwg.glzR = .5
        self.uwg.albWall = 0.91
        self.uwg.SHGC = 0.65
        self.uwg.flr_h = 4.5

        # From blank inputs will be from DOE
        self.uwg.read_epw()
        self.uwg.set_input()
        self.uwg.init_BEM_obj()
        self.uwg.init_input_obj()

        assert self.uwg.BEM[0].building.glazingRatio == pytest.approx(0.5, abs=1e-15)
        assert self.uwg.BEM[0].roof.albedo == pytest.approx(0.5, abs=1e-15)
        assert self.uwg.BEM[0].roof.vegCoverage == pytest.approx(0.1, abs=1e-15)
        assert self.uwg.BEM[1].building.glazingRatio == pytest.approx(0.5, abs=1e-15)
        assert self.uwg.BEM[1].roof.albedo == pytest.approx(0.5, abs=1e-15)
        assert self.uwg.BEM[1].roof.vegCoverage == pytest.approx(0.1, abs=1e-15)
        assert self.uwg.BEM[0].wall.albedo == pytest.approx(0.91, abs=1e-15)
        assert self.uwg.BEM[1].building.shgc == pytest.approx(0.65, abs=1e-15)
        assert self.uwg.BEM[0].building.floorHeight == pytest.approx(4.5, abs=1e-15)

    def test_procMat(self):
        """
        Test different max/min layer depths that generate different diffrent road layer
        thicknesses (to account for too deep elements with inaccurate heat transfer).
        """

        self.setup_uwg_integration("SGP_Singapore.486980_IWEC.epw", "initialize_singapore.uwg")
        self.uwg.read_epw()
        self.uwg.set_input()
        self.uwg.init_BEM_obj()
        self.uwg.init_input_obj()

        #test a 0.5m road split into 10 slices of 0.05m
        # base case; min=0.01, max=0.05, stays the same
        roadMat, newthickness = uwg.procMat(self.uwg.road, 0.05, 0.01)
        assert len(roadMat) == pytest.approx(11, abs=1e-6)
        assert len(newthickness) == pytest.approx(11, abs=1e-6)
        assert sum(newthickness) == pytest.approx(0.05*11, abs=1e-6)

        # modify to one layer for tests
        self.uwg.road.layerThickness = [0.05]
        self.uwg.road.layerThermalCond = self.uwg.road.layerThermalCond[:1]
        self.uwg.road.layerVolHeat = self.uwg.road.layerVolHeat[:1]

        #0.05 layer, will split in two
        roadMat, newthickness = uwg.procMat(self.uwg.road, 0.05, 0.01)
        assert len(roadMat) == pytest.approx(2, abs=1e-6)
        assert len(newthickness) == pytest.approx(2, abs=1e-6)
        assert sum(newthickness) == pytest.approx(0.025*2, abs=1e-6)

        #0.015 layer, will split in min thickness in two
        self.uwg.road.layerThickness = [0.015]
        roadMat, newthickness = uwg.procMat(self.uwg.road, 0.05, 0.01)
        assert len(roadMat) == pytest.approx(2, abs=1e-6)
        assert len(newthickness) == pytest.approx(2, abs=1e-6)
        assert sum(newthickness) == pytest.approx(0.005*2, abs=1e-6)

        #0.12 layer, will split into 3 layers b/c > max_thickness
        self.uwg.road.layerThickness = [0.12]
        roadMat, newthickness = uwg.procMat(self.uwg.road, 0.05, 0.01)
        assert len(roadMat) == pytest.approx(3, abs=1e-6)
        assert len(newthickness) == pytest.approx(3, abs=1e-6)
        assert sum(newthickness) == pytest.approx(0.04*3, abs=1e-6)

        # Old Tests
        # min=0.01, max=0.04, 0.05 cut into two = 0.025
        #roadMat, newthickness = uwg.procMat(self.uwg.road, 0.04, 0.01)
        #assert len(roadMat) == pytest.approx(20, abs=1e-6)
        #assert len(newthickness) == pytest.approx(20, abs=1e-6)
        #assert sum(newthickness) == pytest.approx(0.025*20, abs=1e-6)

        # min=0.0001, max=0.1, should stay the same
        #roadMat, newthickness = uwg.procMat(self.uwg.road,0.1,0.0001)
        #assert len(roadMat) == pytest.approx(11, abs=1e-6)
        #assert len(newthickness) == pytest.approx(11, abs=1e-6)
        #assert sum(newthickness) == pytest.approx(0.5, abs=1e-6)

        # min=0.06, max=0.1, should make new material at 0.06 thickness
        #roadMat, newthickness = uwg.procMat(self.uwg.road, 0.1, 0.06)
        #assert len(roadMat) == pytest.approx(11, abs=1e-6)
        #assert len(newthickness) == pytest.approx(11, abs=1e-6)
        #assert sum(newthickness) == pytest.approx(0.06*11, abs=1e-6)


    def test_hvac_autosize(self):

        self.setup_uwg_integration("SGP_Singapore.486980_IWEC.epw", "initialize_singapore.uwg")
        self.uwg.read_epw()
        self.uwg.set_input()
        self.uwg.init_BEM_obj()
        self.uwg.init_input_obj()
        self.uwg.hvac_autosize()

        assert self.uwg.autosize == pytest.approx(0.0, abs=1e-3)
        assert len(self.uwg.BEM) == pytest.approx(2, abs=1e-6)

        # coolCap and heatCap don't retain high accuracy when extracted from the
        # DOE reference csv, so we will reduce the tolerance here
        assert self.uwg.BEM[0].building.coolCap == pytest.approx((3525.66904*1000.0)/46320.0, abs=1e-3)
        assert self.uwg.BEM[0].building.heatCap == pytest.approx((2875.97378*1000.0)/46320.0, abs=1e-3)
        assert self.uwg.BEM[1].building.coolCap == pytest.approx((252.20895*1000.0)/3135., abs=1e-2)
        assert self.uwg.BEM[1].building.heatCap == pytest.approx((132.396*1000.0)/3135., abs=1e-2)

    def test_simulate(self):

        self.setup_uwg_integration("SGP_Singapore.486980_IWEC.epw", "initialize_singapore.uwg")
        self.uwg.read_epw()
        self.uwg.set_input()
        self.uwg.init_BEM_obj()
        self.uwg.init_input_obj()
        self.uwg.hvac_autosize()
        self.uwg.simulate()

        # Parameters from initialize.uwg
        # Month = 1;              % starting month (1-12)
        # Day = 1;                % starting day (1-31)
        # nDay = 31;              % number of days
        # dtSim = 300;            % simulation time step (s)
        # dtWeather = 3600;       % weather time step (s)

        assert self.uwg.N == pytest.approx(744., abs=1e-6)       # total hours in simulation
        assert self.uwg.ph == pytest.approx(0.083333, abs=1e-6)  # dt (simulation time step) in hours

        #test the weather data time series is equal to time step
        assert len(self.uwg.forcIP.infra) == pytest.approx((self.uwg.simTime.nt-1)/12., abs=1e-3)
        # check that simulation time is happening every 5 minutes 8928
        assert self.uwg.simTime.nt-1 == pytest.approx(31*24*3600/300., abs=1e-3)
        # check that weather step time is happening every 1 hour = 744
        assert len(self.uwg.forcIP.dif) ==  pytest.approx(31 * 24, abs=1e-3)

        # check that final day of timestep is at correct dayType
        assert self.uwg.dayType == pytest.approx(1., abs=1e-3)
        assert self.uwg.SchTraffic[self.uwg.dayType-1][self.uwg.simTime.hourDay] == pytest.approx(0.2, abs=1e-6)

if __name__ == "__main__":
    test = TestUWG()
    #test.test_read_epw()
    #test.test_read_input()
    #test.test_procMat()
    #test.test_hvac_autosize()
