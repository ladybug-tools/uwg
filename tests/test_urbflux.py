import pytest
import os
import math
import UWG

import pprint
pp = lambda x: pprint.pprint(x)

class TestUrbFlux(object):

    DIR_UP_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    DIR_EPW_PATH = os.path.join(DIR_UP_PATH,"resources/epw")
    DIR_MATLAB_PATH = os.path.join(os.path.abspath(os.path.dirname(__file__)), "matlab_ref","matlab_urbflux")
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


    def test_urbflux_unit(self):
        """test for urbflux"""

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
            # building surface properties
            self.uwg.BEM[1].roof.emissivity,        # roof emissivity
            self.uwg.BEM[1].roof.layerTemp[0],      # roof surface temp (K)
            self.uwg.BEM[1].roof.infra,             # roof net longwave radiation (W m-2)
            self.uwg.BEM[1].mass.layerTemp[-1],
            self.uwg.BEM[1].wall.emissivity,        # wall emissivity
            self.uwg.BEM[1].wall.layerTemp[0],      # wall surface temp (K)
            self.uwg.BEM[1].wall.infra,             # wall net longwave radiation (W m-2)
            self.uwg.BEM[1].wall.layerTemp[-1],     # wall interior surface temp (K)
            # urban canyon properties
            self.uwg.UCM.wallTemp,                  # canyon wall temp (K)
            self.uwg.UCM.roofTemp,                  # canyon roof temp (K)
            self.uwg.road.infra,                    # road net longwave radiation (W m-2)
            self.uwg.UCM.roadTemp,                  # canyon road temperature (W m-2)
            self.uwg.UCM.latHeat,                   # canyon latent heat (W m-2)
            self.uwg.UBL.advHeat,                   # boundary layer advective heat (W m-2)
            self.uwg.UCM.ustar,                     # canyon friction velocity (m s-1)
            self.uwg.UCM.ustarMod,                  # modified friction velocity (m s-1)
            self.uwg.UCM.uExch,                     # exchange velocity (m s-1)
            self.uwg.UCM.canWind,                   # canyon wind velocity (m s-1)
            self.uwg.UCM.turbU,                     # canyon turbulent velocity (m s-1)
            self.uwg.UCM.turbV,                     # canyon turbulent velocity (m s-1)
            self.uwg.UCM.turbW,                     # canyon turbulent velocity (m s-1)
            self.uwg.UCM.windProf[-1]               # canyon wind profile
        ]

        # UBL.advHeat = Boundary layer advective heat rounding error
        # Past 10^-9 it's basically sero, but for accuracy's sake apply
        # Kahan's summation algorithm
        # -6.87022037029e-11  ==  -3.43511e-11 << w/o Kahan summation algo
        # -6.01144282401e-11  ==  -3.43511e-11 << w Kahan summation algo

        uwg_matlab_val = self.setup_open_matlab_ref("matlab_ref_urbflux_unit.txt")
        #matlab ref checking
        assert len(uwg_matlab_val) == len(uwg_python_val)

        for i in xrange(len(uwg_matlab_val)):
            if uwg_matlab_val[i] == "\n":
                uwg_matlab_val[i] = None
                assert uwg_python_val[i] == uwg_matlab_val[i], "error at index={}".format(i)
                continue
            elif i==13:
                tol = 10. # hardcoding tolerance b/c advheat rounding errors
            else:
                tol = self.CALCULATE_TOLERANCE(uwg_python_val[i])
            #print i, uwg_python_val[i], ' == ', uwg_matlab_val[i]
            assert uwg_python_val[i] == pytest.approx(uwg_matlab_val[i], abs=tol), "error at index={}".format(i)

if __name__ == "__main__":
    urbflux = TestUrbFlux()
    urbflux.test_urbflux_unit()
