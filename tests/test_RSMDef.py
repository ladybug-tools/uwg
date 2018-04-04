import os
import pytest
import UWG
import math

class TestRSMDef(object):
    """Test for RSMDef.py - Rural Site & Vertical Diffusion Model (RSM & VDM)

    Naming: Test prefixed test classes (without an __init__ method)
    for test autodetection by pytest
    """
    #RESOURCE_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', "resources"))
    #DIR_EPW_PATH = os.path.join(RESOURCE_PATH,"epw")


    """
    """
    DIR_UP_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    DIR_EPW_PATH = os.path.join(DIR_UP_PATH,"resources/epw")
    DIR_MATLAB_PATH = os.path.join(os.path.abspath(os.path.dirname(__file__)), "matlab_ref","matlab_rsmdef")
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

    def test_rsm_init(self):
        """
        Initialize RSM instance for rural and urban site parameters

        Reference construction:

        lat = 1.37
        lon = 103.98
        GMT = 8.0
        rural_height = 0.1            # average obstacle height from initialize.uwg
        urban_height = 1.0            # average building height / 10.
        T_init = 297.85               # initial dry bulb
        P_init = 100900.0             # initial pressure

        self.uwg.geoParam # geographic parameters

        RSM = RSMDef(lat,lon,GMT,rural_height,T_init,P_init,geo_param,self.RESOURCE_PATH)
        USM = RSMDef(lat,lon,GMT,urban_height,T_init,P_init,geo_param,self.RESOURCE_PATH)

        """
        self.setup_uwg_integration()

        self.uwg.read_epw()
        self.uwg.read_input()
        self.uwg.set_input()

        # check date
        assert self.uwg.simTime.month == 1
        assert self.uwg.simTime.day == 1
        assert self.uwg.simTime.secDay/3600. == pytest.approx(0.0,abs=1e-15)

        # Test z_meso list lenght
        assert len(self.uwg.RSM.z_meso) == pytest.approx(56., abs=1e-6)
        assert self.uwg.RSM.dz[1] == pytest.approx(4.4, abs=1e-6)
        assert self.uwg.RSM.z[1] == pytest.approx(6.2, abs=1e-6)

        # Test self.nz0, self.nzref, self.nzfor, self.nz10, self.nzi
        assert self.uwg.RSM.nz0 == pytest.approx(1, abs=1e-6)
        assert self.uwg.RSM.nzref == pytest.approx(17, abs=1e-6)
        assert self.uwg.RSM.nzfor == pytest.approx(13., abs=1e-6)
        assert self.uwg.RSM.nz10 == pytest.approx(3, abs=1e-6)
        assert self.uwg.RSM.nzi == pytest.approx(35, abs=1e-6)

        # Test the tempProfile
        assert len(self.uwg.RSM.tempProf) == pytest.approx(17, abs=1e-6)
        assert self.uwg.RSM.tempProf[16] == pytest.approx(297.8500, abs=1e-6)

        # Test the presProfile with values from UWG_Matlab to 15 digits of precision
        matlab_presProf = [1.009000000000000,1.008490604570617,1.007930481749694,
        1.007314603288614,1.006637447436641,1.005892951541777,1.005074460319214,
        1.004174669438980,1.003185564067109,1.002098351979773,1.000903390854001,
        0.999590109338784,0.998146921496037,0.996561134212841,0.994818847169896,
        0.992904845102566,0.990802481868362] # * 1.01e+05

        # Test varying precision manually
        assert self.uwg.RSM.presProf[0] == pytest.approx(1.009e5, abs=1e-15)
        assert self.uwg.RSM.presProf[-1] == pytest.approx(0.990802481868362e5, abs=1e-10)
        # Test all
        for i in xrange(1,len(matlab_presProf)-1):
            tol = self.CALCULATE_TOLERANCE(matlab_presProf[i]*1e5)
            assert self.uwg.RSM.presProf[i] == pytest.approx(matlab_presProf[i]*1e5, abs=tol)

        # Test the tempRealProf
        # Add dynamic tolerance
        assert self.uwg.RSM.tempRealProf[0] == pytest.approx(2.978499999999999*1e2, abs=1e-10)
        assert self.uwg.RSM.tempRealProf[6] == pytest.approx(2.975182902489641*1e2, abs=1e-10)
        assert self.uwg.RSM.tempRealProf[-1] == pytest.approx(2.963044480678944*1e2, abs=1e-10)

        # Test the densityProfS
        matlab_densityProfS = [1.180352339267655,1.180149676936383,1.179703870430437,
        1.179213600207828,1.178674444441198,1.178081544270841,1.177429561181760,
        1.176712630344852,1.175924309566055,1.175057523461647,1.174104502451500,
        1.173056716134468,1.171904800590439,1.170638479118457,1.169246475901031,
        1.167716422060640,1.166034753626033,1.165110236615915]

        for i in xrange(len(matlab_densityProfS)):
            tol = self.CALCULATE_TOLERANCE(matlab_densityProfS[i]*1e5)
            assert self.uwg.RSM.densityProfS[i] == pytest.approx(matlab_densityProfS[i], abs=tol)

        assert len(self.uwg.RSM.densityProfC) == pytest.approx(self.uwg.RSM.nzref, abs=1e-15)
        assert len(self.uwg.RSM.densityProfS) == pytest.approx(self.uwg.RSM.nzref+1, abs=1e-15)

    def test_rsm_vdm(self):
        """ test RSM VDM against matlab references
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
        self.uwg.simulate()

        # check date
        #print self.uwg.simTime
        assert self.uwg.simTime.month == 1
        assert self.uwg.simTime.day == 1
        assert self.uwg.simTime.secDay == pytest.approx(300.0,abs=1e-15)


        # Matlab Checking for RSM.VDM
        # 2d matrix = 7 x 16
        uwg_python_val = [
            self.uwg.RSM.presProf,
            self.uwg.RSM.tempRealProf,
            self.uwg.RSM.densityProfC,
            self.uwg.RSM.densityProfS,
            self.uwg.RSM.tempProf,
            self.uwg.RSM.windProf,
            [self.uwg.RSM.ublPres]
        ]

        # Flatten 2d matrix into 1d vector
        uwg_python_val = reduce(lambda x,y: x+y, uwg_python_val)

        uwg_matlab_val = self.setup_open_matlab_ref("matlab_ref_rsmdef_vdm.txt")

        # matlab ref checking
        assert len(uwg_matlab_val) == len(uwg_python_val)

        for i in xrange(len(uwg_matlab_val)):
            #print uwg_python_val[i], uwg_matlab_val[i]
            tol = self.CALCULATE_TOLERANCE(uwg_python_val[i])
            assert uwg_python_val[i] == pytest.approx(uwg_matlab_val[i], abs=tol), "error at index={}".format(i)



if __name__ == "__main__":
    test_rsm = TestRSMDef()
    test_rsm.test_rsm_init()
    test_rsm.test_rsm_vdm()
