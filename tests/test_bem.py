try:
    range = xrange
except NameError:
    pass

import pytest
import os
import math
import uwg
from .test_base import TestBase
import pprint

pp = pprint.pprint


class TestBEM(TestBase):

    def test_bem_building_init_largeoffice(self):
        """test for BEM.building init"""

        self.setup_uwg_integration()
        self.uwg.read_epw()
        self.uwg.set_input()
        self.uwg.init_BEM_obj()
        self.uwg.init_input_obj()

        uwg_python_val = [
            self.uwg.BEM[0].building.floorHeight,        # floor height (m)
            self.uwg.BEM[0].building.intHeat,            # timestep internal heat gains (W m-2 bld) (sensible only)
            self.uwg.BEM[0].building.intHeatNight,       # nighttime internal heat gains (W m-2 floor)
            self.uwg.BEM[0].building.intHeatDay,         # daytime internal heat gains (W m-2 floor)
            self.uwg.BEM[0].building.intHeatFRad,        # radiant fraction of internal gains
            self.uwg.BEM[0].building.intHeatFLat,        # latent fraction of internal gains
            self.uwg.BEM[0].building.infil,              # Infiltration (ACH)
            self.uwg.BEM[0].building.vent,               # Ventilation (ACH)
            self.uwg.BEM[0].building.glazingRatio,       # glazing ratio
            self.uwg.BEM[0].building.uValue,             # window U-value (W m-2 K-1) (including film coeff)
            self.uwg.BEM[0].building.shgc,               # window SHGC
            self.uwg.BEM[0].building.condType,           # cooling condensation system type {'AIR', 'WATER'}
            self.uwg.BEM[0].building.cop,                # COP of the cooling system (nominal)
            self.uwg.BEM[0].building.coolSetpointDay,    # daytime indoor cooling set-point (K)
            self.uwg.BEM[0].building.coolSetpointNight,  # nighttime indoor cooling set-point (K)
            self.uwg.BEM[0].building.heatSetpointDay,    # daytime indoor heating set-point (K)
            self.uwg.BEM[0].building.heatSetpointNight,  # nighttime indoor heating set-point (K)
            self.uwg.BEM[0].building.coolCap,            # rated cooling system capacity (W m-2)
            self.uwg.BEM[0].building.heatCap,            # default heating system capacity (W m-2)
            self.uwg.BEM[0].building.heatEff,            # heating system efficiency (-)
            self.uwg.BEM[0].building.mSys,               # HVAC supply mass flowrate (kg s-1 m-2)
            self.uwg.BEM[0].building.indoorTemp,         # indoor air temperature (K)
            self.uwg.BEM[0].building.indoorHum,          # indoor specific humidity (kg / kg)
            self.uwg.BEM[0].building.FanMax,             # max fan flow rate (m^3/s) per DOE
            self.uwg.BEM[0].building.Type,               # DOE reference building type
            self.uwg.BEM[0].building.Era,                # PRE80, PST80, NEW
            self.uwg.BEM[0].building.Zone                # Climate zone number
        ]
        #TODO: track down where canyon fraciton is comding from
        #self.uwg.BEM[0].building.canyon_fraction,    # fraction of waste heat released to canyon, default = 1

        uwg_matlab_val = self.setup_open_matlab_ref("matlab_bem","matlab_ref_bem_building_init_largeoffice.txt")

        # matlab ref checking
        assert len(uwg_matlab_val) == len(uwg_python_val)

        for i in range(len(uwg_matlab_val)):
            #print uwg_python_val[i], uwg_matlab_val[i]
            if type(uwg_python_val[i]) == type("a"):
                assert "".join(uwg_python_val[i].split()) == "".join(uwg_matlab_val[i].split()), "error at index={}".format(i)
            else:
                tol = self.CALCULATE_TOLERANCE(uwg_python_val[i],14.0)
                assert uwg_python_val[i] == pytest.approx(uwg_matlab_val[i], abs=tol), "error at index={}".format(i)

    def test_bem_building_bemcalc_largeoffice_cooling(self):
        """
        test for bem.building bemcalc during cooling period
        """

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


        # Calculated values
        # commneted out properties are instantiated building variables that are never actually set in UWG_Matlab
        uwg_python_val = [
            # changes from constructor
            self.uwg.BEM[0].building.intHeat,            # timestep internal heat gains (W m-2 bld) (sensible only)
            self.uwg.BEM[0].building.intHeatNight,       # nighttime internal heat gains (W m-2 floor)
            self.uwg.BEM[0].building.intHeatDay,         # daytime internal heat gains (W m-2 floor)
            self.uwg.BEM[0].building.indoorTemp,         # indoor air temperature (K)
            self.uwg.BEM[0].building.indoorHum,          # indoor specific humidity (kg / kg)
            # new values
            #self.uwg.BEM[0].building.Tdp,               # dew point
            self.uwg.BEM[0].building.indoorRhum,         # indoor relative humidity
            self.uwg.BEM[0].building.nFloor,             # number of floors
            #self.uwg.BEM[0].building.RadFOcc,            # Radiant fraction of occupant
            #self.uwg.BEM[0].building.LatFOcc,            # Latent fraction of occupant
            #self.uwg.BEM[0].building.RadFEquip,          # Radiant fraction of equipment
            #self.uwg.BEM[0].building.RadFLight,          # Radiant fraction of light
            self.uwg.BEM[0].building.sensCoolDemand,     # building sensible cooling demand (W m-2)
            self.uwg.BEM[0].building.sensHeatDemand,     # building sensible heating demand (W m-2)
            self.uwg.BEM[0].building.copAdj,             # adjusted COP per temperature
            self.uwg.BEM[0].building.dehumDemand,        # dehumidification energy (W m-2)
            self.uwg.BEM[0].building.coolConsump,        # cooling energy consumption (W m-2)
            self.uwg.BEM[0].building.heatConsump,        # heating energy consumption (W m-2)
            self.uwg.BEM[0].building.sensWaste,          # sensible waste heat (W m-2)
            #self.uwg.BEM[0].building.latWaste,           # lat waste heat (W m-2)
            self.uwg.BEM[0].building.fluxMass,           # mass surface heat flux (W m-2) (mass to indoor air)
            self.uwg.BEM[0].building.fluxWall,           # wall surface heat flux (W m-2) (wall to inside)
            self.uwg.BEM[0].building.fluxRoof,           # roof surface heat flux (W m-2) (roof to inside)
            self.uwg.BEM[0].building.fluxSolar,          # solar heat gain (W m-2) through window (SHGC)
            self.uwg.BEM[0].building.fluxWindow,         # heat gain/loss from window (U-value)
            self.uwg.BEM[0].building.fluxInterior,       # internal heat gain adjusted for latent/LW heat (W m-2)
            self.uwg.BEM[0].building.fluxInfil,          # heat flux from infiltration (W m-2)
            self.uwg.BEM[0].building.fluxVent,           # heat flux from ventilation (W m-2)
            self.uwg.BEM[0].building.ElecTotal,          # total electricity consumption - (W/m^2) of floor
            self.uwg.BEM[0].building.GasTotal,   ##        # total gas consumption - (W/m^2) of floor
            self.uwg.BEM[0].building.Qhvac,              # total heat removed (sensible + latent)
            self.uwg.BEM[0].building.Qheat               # total heat added (sensible only)
        ]

        uwg_matlab_val = self.setup_open_matlab_ref("matlab_bem","matlab_ref_bem_building_bemcalc_largeoffice_cooling.txt")

        #matlab ref checking
        assert len(uwg_matlab_val) == len(uwg_python_val)

        for i in range(len(uwg_matlab_val)):
            #print i, uwg_python_val[i], ' == ', uwg_matlab_val[i]
            tol = self.CALCULATE_TOLERANCE(uwg_python_val[i],15.0)
            assert uwg_python_val[i] == pytest.approx(uwg_matlab_val[i], abs=tol), "error at index={}".format(i)

if __name__ == "__main__":
    b = TestBEM()
    #b.test_bem_building_init_largeoffice()
    b.test_bem_building_bemcalc_largeoffice_cooling()
