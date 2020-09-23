"""Test for RSMDef.py - Rural Site & Vertical Diffusion Model (RSM & VDM)"""

import os
import pytest
from .test_base import auto_setup_uwg, calculate_tolerance, setup_open_matlab_ref
from functools import reduce


TEST_DIR = os.path.abspath(os.path.dirname(__file__))


def test_rsm_init():
    """Test Initialize RSM instance for rural site parameters

    Reference init:

    lat = 1.37
    lon = 103.98
    GMT = 8.0
    rural_height = 0.1            # average obstacle height from initialize.uwg
    urban_height = 1.0            # average building height / 10.
    T_init = 297.85               # initial dry bulb
    P_init = 100900.0             # initial pressure

    testuwg.geoParam # geographic parameters

    RSM = RSMDef(
        lat, lon, GMT, rural_height, T_init, P_init, geo_param, self.RESOURCE_PATH)
    """

    testuwg = auto_setup_uwg()
    testuwg.generate()
    testuwg.RSM.__repr__()

    # check date
    assert testuwg.simTime.month == 1
    assert testuwg.simTime.day == 1
    assert testuwg.simTime.secDay / 3600. == pytest.approx(0.0, abs=1e-15)

    # Test z_meso list lenght
    assert len(testuwg.RSM.z_meso) == pytest.approx(56., abs=1e-6)
    assert testuwg.RSM.dz[1] == pytest.approx(4.4, abs=1e-6)
    assert testuwg.RSM.z[1] == pytest.approx(6.2, abs=1e-6)

    # Test self.nz0, self.nzref, self.nzfor, self.nz10, self.nzi
    assert testuwg.RSM.nz0 == pytest.approx(1, abs=1e-6)
    assert testuwg.RSM.nzref == pytest.approx(17, abs=1e-6)
    assert testuwg.RSM.nzfor == pytest.approx(13., abs=1e-6)
    assert testuwg.RSM.nz10 == pytest.approx(3, abs=1e-6)
    assert testuwg.RSM.nzi == pytest.approx(35, abs=1e-6)

    # Test the tempProfile
    assert len(testuwg.RSM.tempProf) == pytest.approx(17, abs=1e-6)
    assert testuwg.RSM.tempProf[16] == pytest.approx(297.8500, abs=1e-6)

    # Test the presProfile with values from UWG_Matlab to 15 digits of precision
    matlab_presProf = [
        1.009000000000000, 1.008490604570617, 1.007930481749694, 1.007314603288614,
        1.006637447436641, 1.005892951541777, 1.005074460319214, 1.004174669438980,
        1.003185564067109, 1.002098351979773, 1.000903390854001, 0.999590109338784,
        0.998146921496037, 0.996561134212841, 0.994818847169896, 0.992904845102566,
        0.990802481868362]  # * 1.01e+05

    # Test varying precision manually
    assert testuwg.RSM.presProf[0] == pytest.approx(1.009e5, abs=1e-15)
    assert testuwg.RSM.presProf[-1] == pytest.approx(0.990802481868362e5, abs=1e-10)
    # Test all
    for i in range(1, len(matlab_presProf) - 1):
        tol = calculate_tolerance(matlab_presProf[i] * 1e5, 15.0)
        assert testuwg.RSM.presProf[i] == pytest.approx(matlab_presProf[i]*1e5, abs=tol)

    # Test the tempRealProf
    # Add dynamic tolerance
    assert testuwg.RSM.tempRealProf[0] == pytest.approx(2.978499999999999*1e2, abs=1e-10)
    assert testuwg.RSM.tempRealProf[6] == pytest.approx(2.975182902489641*1e2, abs=1e-10)
    assert testuwg.RSM.tempRealProf[-1] == \
        pytest.approx(2.963044480678944*1e2, abs=1e-10)

    # Test the densityProfS
    matlab_densityProfS = [
        1.180352339267655, 1.180149676936383, 1.179703870430437, 1.179213600207828,
        1.178674444441198, 1.178081544270841, 1.177429561181760, 1.176712630344852,
        1.175924309566055, 1.175057523461647, 1.174104502451500, 1.173056716134468,
        1.171904800590439, 1.170638479118457, 1.169246475901031, 1.167716422060640,
        1.166034753626033, 1.165110236615915]

    for i in range(len(matlab_densityProfS)):
        tol = calculate_tolerance(matlab_densityProfS[i] * 1e5, 15.0)
        assert testuwg.RSM.densityProfS[i] == \
            pytest.approx(matlab_densityProfS[i], abs=tol)

    assert len(testuwg.RSM.densityProfC) == pytest.approx(testuwg.RSM.nzref, abs=1e-15)
    assert len(testuwg.RSM.densityProfS) == pytest.approx(testuwg.RSM.nzref+1, abs=1e-15)


def test_rsm_vdm():
    """ test RSM VDM against matlab references."""

    testuwg = auto_setup_uwg()

    # Test Jan 1 (winter, no vegetation coverage)
    testuwg.month = 1
    testuwg.day = 1
    testuwg.nday = 1

    testuwg.generate()

    # In order to avoid integration effects. Test only first time step
    # Subtract timestep to stop at 300 sec
    testuwg.simTime.nt -= (23 * 12 + 11)

    # Run simulation
    testuwg.simulate()

    # check date
    assert testuwg.simTime.month == 1
    assert testuwg.simTime.day == 1
    assert testuwg.simTime.secDay == pytest.approx(300.0, abs=1e-15)

    # Matlab Checking for RSM.VDM
    # 2d matrix = 7 x 16
    uwg_python_val = [
        testuwg.RSM.presProf,
        testuwg.RSM.tempRealProf,
        testuwg.RSM.densityProfC,
        testuwg.RSM.densityProfS,
        testuwg.RSM.tempProf,
        testuwg.RSM.windProf,
        [testuwg.RSM.ublPres]]

    # Flatten 2d matrix into 1d vector
    uwg_python_val = reduce(lambda x, y: x + y, uwg_python_val)

    # Matlab checking
    uwg_matlab_val = setup_open_matlab_ref("matlab_rsmdef", "matlab_ref_rsmdef_vdm.txt")

    # Matlab ref checking
    assert len(uwg_matlab_val) == len(uwg_python_val)

    for i in range(len(uwg_matlab_val)):
        # print uwg_python_val[i], uwg_matlab_val[i]
        tol = calculate_tolerance(uwg_python_val[i], 15.0)
        assert uwg_python_val[i] == pytest.approx(uwg_matlab_val[i], abs=tol), \
            "error at index={}".format(i)


def test_rsm_dissipation_bougeault():
    """Test RSM VDM against matlab references."""

    epw_path = os.path.join(TEST_DIR, 'epw', 'CAN_ON_Toronto.716240_CWEC.epw')
    param_path = os.path.join(TEST_DIR, 'parameters', 'initialize_toronto.uwg')
    testuwg = auto_setup_uwg(epw_path=epw_path, param_path=param_path)
    testuwg.month = 6
    testuwg.day = 1
    testuwg.nday = 30

    # run simulation
    testuwg.generate()
    testuwg.simulate()

    # check date
    # print testuwg.simTime
    assert testuwg.simTime.month == 7
    assert testuwg.simTime.day == 1
    assert testuwg.simTime.secDay == pytest.approx(0.0, abs=1e-15)

    # dld, dlu values after 1 day
    uwg_python_val = [
        testuwg.RSM.dlu,
        testuwg.RSM.dld]

    # Flatten 2d matrix into 1d vector
    uwg_python_val = reduce(lambda x, y: x + y, uwg_python_val)

    # Matlab checking
    uwg_matlab_val = setup_open_matlab_ref(
        "matlab_rsmdef", "matlab_rsmdef_dissipation_bougeault.txt")

    # Matlab ref checking
    assert len(uwg_matlab_val) == len(uwg_python_val)

    for i in range(len(uwg_matlab_val)):
        tol = calculate_tolerance(uwg_python_val[i], 11.0)
        assert uwg_python_val[i] == pytest.approx(uwg_matlab_val[i], abs=tol), \
            "error at index={}".format(i)
