import pytest
import UWG
import os
import math
import pprint

class TestElement(object):

    DIR_UP_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    DIR_EPW_PATH = os.path.join(DIR_UP_PATH,"resources/epw")
    DIR_MATLAB_PATH = os.path.join(os.path.abspath(os.path.dirname(__file__)), "matlab_ref","matlab_element")
    CALCULATE_TOLERANCE = lambda s,x: 1*10**-(15.0 - (int(math.log10(x)) + 1)) if (x > 1. or int(x)==1) else 1e-15

    def setup_uwg_integration(self):
        """ set up uwg object from initialize.uwg """

        epw_dir = self.DIR_EPW_PATH
        epw_file_name = "SGP_Singapore.486980_IWEC.epw"
        uwg_param_dir = os.path.join(self.DIR_UP_PATH,"resources")
        uwg_param_file_name = "initialize.uwg"

        self.uwg = UWG.UWG(epw_dir, epw_file_name, uwg_param_dir, uwg_param_file_name)

    def setup_open_matlab_ref(self,matlab_ref_file_path):
        """ open the matlab reference file """

        matlab_path = os.path.join(self.DIR_MATLAB_PATH,matlab_ref_file_path)
        print matlab_path
        if not os.path.exists(matlab_path):
            raise Exception("Failed to open {}!".format(matlab_path))
        matlab_file = open(matlab_path,'r')
        uwg_matlab_val_ = [float(x) for x in matlab_file.readlines()]
        matlab_file.close()
        return uwg_matlab_val_

    def test_SurfFlux_with_waterStorage(self):
        """ Edge case: test element SurfFlux against matlab reference
        when waterStorage > 0.0.
        This has to be hardcoded b/c doesn't get used otherwise

        """
        self.setup_uwg_integration()
        self.uwg.read_epw()
        self.uwg.read_input()
        self.uwg.set_input()

        # We subtract 11 hours from total timestep so
        # we can stop simulation while we still have sun!
        # New time: Jan 16, 1300
        self.uwg.simTime.nt -= (12*24*15 + 12*13)
        self.uwg.simTime.nt += 1
        # turn rural road waterStorage to 1.
        self.uwg.rural.waterStorage = 0.005 # .5cm thick film (from wgmax constant)

        # run simulation
        self.uwg.hvac_autosize()
        self.uwg.uwg_main()

        # check date
        assert self.uwg.simTime.month == 1
        assert self.uwg.simTime.day == 16
        assert self.uwg.simTime.secDay == pytest.approx(39900.,abs=1e-15)

        # Check waterStorage
        assert 0.005 == pytest.approx(self.uwg.rural.waterStorage, 1e-15)


    def test_SurfFlux(self):
        """ test element SurfFlux against matlab references

        TODO: Description of method

        Example of application:
        for it in xrange(1,simTime.nt,1):
            ...
            rural.SurfFlux(forc,parameter,simTime,humRef,tempRef,windRef,boundCond,intFlux)

        albedo;          % outer surface albedo
        emissivity;      % outer surface emissivity
        layerThickness;  % vector of layer thicknesses (m)
        layerThermalCond;% vector of layer thermal conductivities (W m-1 K-1)
        layerVolHeat;    % vector of layer volumetric heat (J m-3 K-1)
        vegCoverage;     % surface vegetation coverage
        layerTemp;       % vector of layer temperatures (K)
        waterStorage;    % thickness of water film (m) (only for horizontal surfaces)
        horizontal;      % 1-horizontal, 0-vertical
        solRec;          % solar radiation received (W m-2)
        infra;           % net longwave radiation (W m-2)
        lat;             % surface latent heat flux (W m-2)
        sens;            % surface sensible heat flux (W m-2)
        solAbs;          % solar radiation absorbed (W m-2)
        aeroCond;        % convective heat transfer
        T_ext;           % external surface temperature
        T_int;           % internal surface temperature
        flux;            % external surface heat flux

        """
        #X: add all the properties from surflux
        #TODO: Schematizie surfflux from thesis
        #TODO: output matlab properties
        #TODO: test/debug final tdd
        #TODO: describe method
        #TODO: start translation
        #TODO: write manual tests to aid understanding

        self.setup_uwg_integration()

        self.uwg.read_epw()
        self.uwg.read_input()
        self.uwg.set_input()
        self.uwg.hvac_autosize()

        # Change vegStart so that vegetation cover can be included in heat calculation
        self.uwg.geoParam.vegStart = 2
        """
        # Reset simulation time so we can capture sun and vegetation covered
        # February 15th, 13:00
        startmonth = 1
        startday = 1
        days = 31 + 15 # February 15th
        self.uwg.simTime = UWG.SimParam(self.uwg.simTime.dt,self.uwg.simTime.timeForcing,startmonth,startday,days)
        # We subtract 11 hours from total timestep so still have sun. New time: 1300
        self.uwg.simTime.nt -= 12*11

        # check date
        assert self.uwg.simTime.month == 2
        assert self.uwg.simTime.day == 15
        assert self.uwg.simTime.secDay/3600. == pytest.approx(13.0,abs=1e-15)

        # run simulation
        self.uwg.uwg_main()


        uwg_python_val = [
        self.uwg.rural.aeroCond,        # Convection coef (refL UWG, eq.12)
        self.uwg.rural.waterStorage,    # thickness of water film (m) (only for horizontal surfaces)
        self.uwg.rural.solAbs,          # solar radiation absorbed (W m-2)
        self.uwg.rural.lat,             # surface latent heat flux (W m-2)
        self.uwg.rural.sens,            # surface sensible heat flux (W m-2)
        self.uwg.rural.flux             # external surface heat flux (W m-2)
        ]

        # Matlab Checking for rural road
        uwg_matlab_val = self.setup_open_matlab_ref("matlab_ref_element_surfflux.txt")

        # matlab ref checking
        assert len(uwg_matlab_val) == len(uwg_python_val)

        for i in xrange(len(uwg_matlab_val)):
            print uwg_python_val[i], uwg_matlab_val[i]
            #tol = self.CALCULATE_TOLERANCE(uwg_matlab_val[i])
            #assert uwg_python_val[i] == pytest.approx(uwg_matlab_val[i], abs=tol), "error at index={}".format(i)
        """

if __name__ == "__main__":
    te = TestElement()
    te.test_SurfFlux()
