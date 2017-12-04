import os
import pytest
import UWG

class TestRSMDef(object):
    """Test for RSMDef.py - Rural Site & Vertical Diffusion Model (RSM & VDM)

    Naming: Test prefixed test classes (without an __init__ method)
    for test autodetection by pytest
    """
    RESOURCE_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', "resources"))
    DIR_EPW_PATH = os.path.join(RESOURCE_PATH,"epw")


    def setup_init(self):
        """ Initialize rsm class """

        # Make UWG instance
        epw_dir = self.DIR_EPW_PATH
        epw_file_name = "SGP_Singapore.486980_IWEC.epw"
        uwg_param_dir = self.RESOURCE_PATH
        uwg_param_file_name = "initialize.uwg"
        self.uwg = UWG.UWG(epw_dir, epw_file_name, uwg_param_dir, uwg_param_file_name)
        self.uwg.read_epw()
        self.uwg.read_input()
        print '-----'
        # Make RSM instance for rural site parameters
        lat = 1.37
        lon = 103.98
        GMT = 8.0
        rural_height = 0.1            # average obstacle height from initialize.uwg
        urban_height = 1.0            # average building height / 10.
        T_init = 297.85               # initial dry bulb
        P_init = 100900.0             # initial pressure

        geo_param = self.uwg.geoParam # geographic parameters
        self.rural_rsm = UWG.RSMDef(lat,lon,GMT,rural_height,T_init,P_init,geo_param,self.RESOURCE_PATH)
        #self.urban_rsm = UWG.RSMDef(lat,lon,GMT,urban_height,T_init,P_init,geo_param,self.RESOURCE_PATH)

    def test_rsm(self):
        self.setup_init()

        # Test z_meso list lenght
        assert len(self.rural_rsm.z_meso) == pytest.approx(56., abs=1e-6)
        assert self.rural_rsm.dz[1] == pytest.approx(4.4, abs=1e-6)
        assert self.rural_rsm.z[1] == pytest.approx(6.2, abs=1e-6)

        # Test self.nz0, self.nzref, self.nzfor, self.nz10, self.nzi
        assert self.rural_rsm.nz0 == pytest.approx(0, abs=1e-6)
        assert self.rural_rsm.nzref == pytest.approx(16, abs=1e-6)
        assert self.rural_rsm.nzfor == pytest.approx(12., abs=1e-6)
        assert self.rural_rsm.nz10 == pytest.approx(2, abs=1e-6)
        assert self.rural_rsm.nzi == pytest.approx(34, abs=1e-6)

        # Test the tempProfile
        assert len(self.rural_rsm.tempProf) == pytest.approx(17, abs=1e-6)
        assert self.rural_rsm.tempProf[16] == pytest.approx(297.8500, abs=1e-6)

        # Test the presProfile with values from UWG_Matlab to 15 digits of precision
        # Note: 16 digits of precision == max for python 64bit float (2^-53 = 1.11e-16)
        # So for numbers > 0. we have to subtract the 16 significant digits (significand)
        # from the integer part when testing the numbers

        matlab_presProf = [
            1.009000000000000,
            1.008490604570617,
            1.007930481749694,
            1.007314603288614,
            1.006637447436641,
            1.005892951541777,
            1.005074460319214,
            1.004174669438980,
            1.003185564067109,
            1.002098351979773,
            1.000903390854001,
            0.999590109338784,
            0.998146921496037,
            0.996561134212841,
            0.994818847169896,
            0.992904845102566,
            0.990802481868362
            ]

        #assert self.rural_rsm.presProf[0] == pytest.approx(1.009e5, abs=1e-15)
        #assert self.rural_rsm.presProf[-1] == pytest.approx(0.990802481868362e5, abs=1e-10) #
        #for i in xrange(1,len(matlab_presProf)-1):
        #    assert self.rural_rsm.presProf[i] == pytest.approx(matlab_presProf[i]*1e5, abs=1e-6)


if __name__ == "__main__":
    test_rsm = TestRSMDef()
    test_rsm.test_rsm()
