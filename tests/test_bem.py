import pytest
import os
import math
import UWG

import pprint
pp = lambda x: pprint.pprint(x)

class TestBEM(object):

    DIR_UP_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    DIR_EPW_PATH = os.path.join(DIR_UP_PATH,"resources/epw")
    DIR_MATLAB_PATH = os.path.join(os.path.abspath(os.path.dirname(__file__)), "matlab_ref","matlab_bem")
    CALCULATE_TOLERANCE = lambda s,x: 1*10**-(15.0 - (int(math.log10(x)) + 1)) if (x > 1. or int(x)==1) else 1e-15


    def setup_uwg_integration(self):
        """ set up initialize.uwg """

        epw_dir = self.DIR_EPW_PATH
        epw_file_name = "SGP_Singapore.486980_IWEC.epw"
        uwg_param_dir = os.path.join(self.DIR_UP_PATH,"resources")
        uwg_param_file_name = "initialize.uwg"

        self.uwg = UWG.UWG(epw_dir, epw_file_name, uwg_param_dir, uwg_param_file_name)

    def setup_open_matlab_ref(self,matlab_ref_file_path):
        """ open the matlab reference file """
        def convert_type(x):
            try:
                return float(x)
            except ValueError:
                return x
        matlab_path = os.path.join(self.DIR_MATLAB_PATH,matlab_ref_file_path)
        if not os.path.exists(matlab_path):
            raise Exception("Failed to open {}!".format(matlab_path))
        matlab_file = open(matlab_path,'r')
        uwg_matlab_val = [convert_type(x) for x in matlab_file.readlines()]
        matlab_file.close()

        return uwg_matlab_val


    def test_bem_building_init_largeoffice(self):
        """test for BEM.building init"""

        self.setup_uwg_integration()
        self.uwg.read_epw()
        self.uwg.read_input()
        self.uwg.set_input()

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

        uwg_matlab_val = self.setup_open_matlab_ref("matlab_ref_bem_building_init_largeoffice.txt")

        # matlab ref checking
        assert len(uwg_matlab_val) == len(uwg_python_val)

        for i in xrange(len(uwg_matlab_val)):
            #print uwg_python_val[i], uwg_matlab_val[i]
            if type(uwg_python_val[i]) == type("a"):
                assert "".join(uwg_python_val[i].split()) == "".join(uwg_matlab_val[i].split()), "error at index={}".format(i)
            else:
                tol = self.CALCULATE_TOLERANCE(uwg_python_val[i])
                assert uwg_python_val[i] == pytest.approx(uwg_matlab_val[i], abs=tol), "error at index={}".format(i)

    def test_bem_building_bemcalc_largeoffice(self):
        """
        test for bem.building bemcalc
        """

        self.setup_uwg_integration()
        self.uwg.read_epw()
        self.uwg.read_input()

        # Test Jan 1 (winter, no vegetation coverage)
        self.uwg.Month = 1
        self.uwg.Day = 1
        self.uwg.nDay = 1

        # set_input
        self.uwg.set_input()

        # In order to avoid integration effects. Test only first time step
        # Subtract timestep to stop at 300 sec
        self.uwg.simTime.nt -= (23*12 + 11)

        # Run simulation
        self.uwg.hvac_autosize()
        self.uwg.uwg_main()

        # check date
        #print self.uwg.simTime
        assert self.uwg.simTime.month == 1
        assert self.uwg.simTime.day == 1
        assert self.uwg.simTime.secDay == pytest.approx(300.0,abs=1e-15)


        # Calculated values
        uwg_python_val = [
            # changes from constructor
            self.uwg.BEM[0].building.intHeat,            # timestep internal heat gains (W m-2 bld) (sensible only)
            self.uwg.BEM[0].building.intHeatNight,       # nighttime internal heat gains (W m-2 floor)
            self.uwg.BEM[0].building.intHeatDay,         # daytime internal heat gains (W m-2 floor)
            self.uwg.BEM[0].building.indoorTemp,         # indoor air temperature (K)
            self.uwg.BEM[0].building.indoorHum,          # indoor specific humidity (kg / kg)
            # new values
            """
            self.uwg.BEM[0].building.Twb,                # wetbulb temperature
            self.uwg.BEM[0].building.Tdp,                # dew point
            self.uwg.BEM[0].building.indoorRhum,         # indoor relative humidity
            self.uwg.BEM[0].building.area_floor,         # total floor space of the BEM
            self.uwg.BEM[0].building.nFloor,             # number of floors
            self.uwg.BEM[0].building.RadFOcc,            # Radiant fraction of occupant
            self.uwg.BEM[0].building.LatFOcc,            # Latent fraction of occupant
            self.uwg.BEM[0].building.RadFEquip,          # Radiant fraction of equipment
            self.uwg.BEM[0].building.RadFLight,          # Radiant fraction of light
            self.uwg.BEM[0].building.sensCoolDemand,     # building sensible cooling demand (W m-2)
            self.uwg.BEM[0].building.sensHeatDemand,     # building sensible heating demand (W m-2)
            self.uwg.BEM[0].building.copAdj,             # adjusted COP per temperature
            self.uwg.BEM[0].building.dehumDemand,        # dehumidification energy (W m-2)
            self.uwg.BEM[0].building.coolConsump,        # cooling energy consumption (W m-2)
            self.uwg.BEM[0].building.heatConsump,        # heating energy consumption (W m-2)
            self.uwg.BEM[0].building.sensWaste,          # sensible waste heat (W m-2)
            self.uwg.BEM[0].building.latWaste,           # lat waste heat (W m-2)
            self.uwg.BEM[0].building.fluxMass,           # mass surface heat flux (W m-2) (mass to indoor air)
            self.uwg.BEM[0].building.fluxWall,           # wall surface heat flux (W m-2) (wall to inside)
            self.uwg.BEM[0].building.fluxRoof,           # roof surface heat flux (W m-2) (roof to inside)
            self.uwg.BEM[0].building.fluxSolar,          # solar heat gain (W m-2) through window (SHGC)
            self.uwg.BEM[0].building.fluxWindow,         # heat gain/loss from window (U-value)
            self.uwg.BEM[0].building.fluxInterior,       # internal heat gain adjusted for latent/LW heat (W m-2)
            self.uwg.BEM[0].building.fluxInfil,          # heat flux from infiltration (W m-2)
            self.uwg.BEM[0].building.fluxVent,           # heat flux from ventilation (W m-2)
            self.uwg.BEM[0].building.ElecTotal,          # total electricity consumption - (W/m^2) of floor
            self.uwg.BEM[0].building.GasTotal,           # total gas consumption - (W/m^2) of floor
            self.uwg.BEM[0].building.Qhvac,              # total heat removed (sensible + latent)
            self.uwg.BEM[0].building.Qheat               # total heat added (sensible only)
            """
        ]

        uwg_matlab_val = self.setup_open_matlab_ref("matlab_ref_bem_building_bemcalc_largeoffice.txt")

        # matlab ref checking
        #assert len(uwg_matlab_val) == len(uwg_python_val)

        for i in xrange(len(uwg_matlab_val)):
            if i<len(uwg_python_val):
                print uwg_python_val[i], uwg_matlab_val[i]
            #assert uwg_python_val[i] == pytest.approx(uwg_matlab_val[i], abs=tol), "error at index={}".format(i)

if __name__ == "__main__":
    b = TestBEM()
    b.test_bem_building_init_largeoffice()
    b.test_bem_building_bemcalc_largeoffice()
