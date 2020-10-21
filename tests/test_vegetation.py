"""Test tree and vegetation coverage."""
from __future__ import division

import os
import pytest
from uwg import UWG, solarcalcs


EPW_PATH = \
    os.path.join(os.path.dirname(__file__), 'epw',
                 'SGP_Singapore.486980_IWEC.epw')


def test_tree_vegetation_ratio():
    """Test road vegetation fractions."""

    model = UWG.from_param_args(
        bldheight=10, blddensity=0.5, vertohor=0.8, zone='1A', grasscover=0,
        treecover=0, epw_path=EPW_PATH)

    with pytest.raises(AssertionError):
        model.blddensity = 0.5
        model.grasscover = 0.5
        model.treecover = 0.1  # > grasscover

    with pytest.raises(AssertionError):
        model.blddensity = 0.5
        model.treecover = 0.3
        model.grasscover = 0.3  # > treecover

    with pytest.raises(AssertionError):
        model.treecover = 0.51
        model.grasscover = 0.49
        model.blddensity = 0.001  # > veg + grass

    model.grasscover = 0.125
    model.treecover = 0.125

    # generate road
    model.generate()

    assert model.vegcover == pytest.approx(0.25, 1e-10)
    # test road coverage (vegcover / (1 - density))
    assert model.blddensity == pytest.approx(0.5, 1e-10)
    # (0.125 + 0.125)/0.5 = 0.5
    assert model.road.vegcoverage == pytest.approx(0.5, 1e-10)
    # 0.125/0.5 = 1/8 * 2 = 1/4
    assert model.UCM.roadShad == pytest.approx(0.25, 1e-10)


def test_tree_drybulb():
    """Test if drybulb decreases when tree cover increases."""

    # Double tree fraction, and increase veg fraction.
    nday = 31
    model1 = UWG.from_param_args(
        bldheight=10, blddensity=0.5, vertohor=0.8, zone='1A', nday=nday,
        vegstart=1, vegend=12, latgrss=0.4, lattree=0.6, albveg=0.5,
        treecover=0, grasscover=0, epw_path=EPW_PATH)

    model2 = UWG.from_param_args(
        bldheight=10, blddensity=0, vertohor=0.8, zone='1A', nday=nday,
        vegstart=1, vegend=12, latgrss=0.4, lattree=0.6, albveg=0.5,
        treecover=0, grasscover=0, epw_path=EPW_PATH)

    # override model1
    model1.blddensity = 0.5
    model1.treecover = 0.125  # 50% tree
    model1.grasscover = 0.125
    # model1.vegcover = 0.25

    # override model2
    model2.blddensity = 0.5
    model2.treecover = 0.25  # 200% tree
    model2.grasscover = 0.125
    # model1.vegcover = 0.375

    model1.generate()
    model2.generate()

    # check vegcover
    assert model1.vegcover == pytest.approx(0.25, 1e-10)
    assert model1.UCM.vegcover == pytest.approx(0.25, 1e-10)
    assert model2.vegcover == pytest.approx(0.375, 1e-10)
    assert model2.UCM.vegcover == pytest.approx(0.375, 1e-10)

    assert model1.road.vegcoverage == pytest.approx(0.5, 1e-10)
    assert model2.road.vegcoverage == pytest.approx(0.75, 1e-10)

    # check tree coverage
    assert model1.UCM.roadShad == pytest.approx(0.25, 1e-10)
    assert model2.UCM.roadShad == pytest.approx(0.5, 1e-10)

    # run models for 7 days starting from Jan 1
    model1.simulate()
    model2.simulate()

    # calc drybulb mean
    temps1 = [ucm.canTemp - 273.15 for ucm in model1.UCMData]
    temps2 = [ucm.canTemp - 273.15 for ucm in model2.UCMData]

    mean_temps1, mean_temps2 = sum(
        temps1) / len(temps1), sum(temps2) / len(temps2)

    # confirm if increasing trees result in lower temps
    assert mean_temps1 > mean_temps2


def test_veg_drybulb():
    """Test if drybulb decreases when grasscover increases."""

    # Double grass fraction, and increase veg fraction.
    model1 = UWG.from_param_args(
        bldheight=10, blddensity=0.5, vertohor=0.8, zone='1A', nday=7,
        vegstart=1, vegend=12, latgrss=0.4, lattree=0.6, albveg=0.5,
        treecover=0, grasscover=0, epw_path=EPW_PATH)

    model2 = UWG.from_param_args(
        bldheight=10, blddensity=0, vertohor=0.8, zone='1A', nday=7,
        vegstart=1, vegend=12, latgrss=0.4, lattree=0.6, albveg=0.5,
        treecover=0, grasscover=0, epw_path=EPW_PATH)

    # override model1
    model1.blddensity = 0.5
    model1.treecover = 0.125  # 50% tree
    model1.grasscover = 0.125  # TODO: make into grasscover
    # model1.vegcover = 0.25

    # override model2
    model2.blddensity = 0.5
    model2.treecover = 0.125  # 200% tree
    model2.grasscover = 0.25  # TODO: make the grasscover
    # model1.vegcover = 0.375

    model1.generate()
    model2.generate()

    # check vegcover
    assert model1.vegcover == pytest.approx(0.25, 1e-10)
    assert model1.UCM.vegcover == pytest.approx(0.25, 1e-10)
    assert model2.vegcover == pytest.approx(0.375, 1e-10)
    assert model2.UCM.vegcover == pytest.approx(0.375, 1e-10)

    assert model1.road.vegcoverage == pytest.approx(0.5, 1e-10)
    assert model2.road.vegcoverage == pytest.approx(0.75, 1e-10)

    # check tree coverage
    assert model1.UCM.roadShad == pytest.approx(0.25, 1e-10)
    assert model2.UCM.roadShad == pytest.approx(0.25, 1e-10)

    # run models for 7 days starting from Jan 1
    model1.simulate()
    model2.simulate()

    # calc drybulb mean
    temps1 = [ucm.canTemp - 273.15 for ucm in model1.UCMData]
    temps2 = [ucm.canTemp - 273.15 for ucm in model2.UCMData]

    mean_temps1, mean_temps2 = sum(
        temps1) / len(temps1), sum(temps2) / len(temps2)

    # confirm if increasing grass result in lower temps
    assert mean_temps1 > mean_temps2


def test_tree_only_drybulb():
    """Test drybulb impact of treecover change but vegcover is constant.

    Note: while counterintuitive that increasing tree coverage increases dry bulb,
    this is primarily due to the tree canopy blocking long-wave radiation.
    """

    # Double tree fraction, and subtract grass fraction.
    model1 = UWG.from_param_args(
        bldheight=10, blddensity=0, vertohor=0.8, zone='1A', nday=7,
        vegstart=1, vegend=12, albveg=0.5, treecover=0, grasscover=0,
        epw_path=EPW_PATH)

    model2 = UWG.from_param_args(
        bldheight=10, blddensity=0, vertohor=0.8, zone='1A', nday=7,
        vegstart=1, vegend=12, albveg=0.5, treecover=0, grasscover=0,
        epw_path=EPW_PATH)

    # override model1
    model1.blddensity = 0.5
    model1.treecover = 0.125  # 50% tree
    model1.grasscover = 0.125
    model1.latgrss = 0.6
    model1.lattree = 0.6
    # model1.vegcover = 0.25

    # override model2
    model2.blddensity = 0.5
    model2.treecover = 0.25  # 200% tree
    model2.grasscover = 0.0
    model2.latgrss = 0.6
    model2.lattree = 0.6
    # model2.vegcover = 0.25

    model1.generate()
    model2.generate()

    # check vegcover
    assert model1.vegcover == pytest.approx(0.25, 1e-10)
    assert model1.UCM.vegcover == pytest.approx(0.25, 1e-10)
    assert model2.vegcover == pytest.approx(0.25, 1e-10)
    assert model2.UCM.vegcover == pytest.approx(0.25, 1e-10)

    assert model1.road.vegcoverage == pytest.approx(0.5, 1e-10)
    assert model2.road.vegcoverage == pytest.approx(0.5, 1e-10)

    # check tree coverage
    assert model1.UCM.roadShad == pytest.approx(0.25, 1e-10)
    assert model2.UCM.roadShad == pytest.approx(0.5, 1e-10)

    # run models for 7 days starting from Jan 1
    model1.simulate()
    model2.simulate()

    # calc drybulb mean
    temps1 = [ucm.canTemp - 273.15 for ucm in model1.UCMData]
    temps2 = [ucm.canTemp - 273.15 for ucm in model2.UCMData]

    mean_temps1, mean_temps2 = sum(
        temps1) / len(temps1), sum(temps2) / len(temps2)

    # increasing trees result in higher temps due to long-wave radiation
    # print(mean_temps1, mean_temps2)
    assert mean_temps1 < mean_temps2


def test_veg_only_drybulb():
    """Test if drybulb decreases when grasscover increases but vegcover is constant."""

    # Double grass fraction, and maintain veg fraction.
    model1 = UWG.from_param_args(
        bldheight=10, blddensity=0.5, vertohor=0.8, zone='1A', nday=7,
        vegstart=1, vegend=12, latgrss=0.4, lattree=0.6, albveg=0.5,
        treecover=0, grasscover=0, epw_path=EPW_PATH)

    model2 = UWG.from_param_args(
        bldheight=10, blddensity=0, vertohor=0.8, zone='1A', nday=7,
        vegstart=1, vegend=12, albveg=0.5, treecover=0, grasscover=0,
        epw_path=EPW_PATH)

    # override model1
    model1.blddensity = 0.5
    model1.treecover = 0.125  # 50% tree
    model1.grasscover = 0.125
    model1.latgrss = 0.4
    model1.lattree = 0.6
    # model1.vegcover = 0.25

    # override model2
    model2.blddensity = 0.5
    model2.treecover = 0.0  # 200% tree
    model2.grasscover = 0.25
    model2.latgrss = 0.4
    model2.lattree = 0.6
    # model1.vegcover = 0.25

    model1.generate()
    model2.generate()

    # check vegcover
    assert model1.vegcover == pytest.approx(0.25, 1e-10)
    assert model1.UCM.vegcover == pytest.approx(0.25, 1e-10)
    assert model2.vegcover == pytest.approx(0.25, 1e-10)
    assert model2.UCM.vegcover == pytest.approx(0.25, 1e-10)

    assert model1.road.vegcoverage == pytest.approx(0.5, 1e-10)
    assert model2.road.vegcoverage == pytest.approx(0.5, 1e-10)

    # check tree coverage
    assert model1.UCM.roadShad == pytest.approx(0.25, 1e-10)
    assert model2.UCM.roadShad == pytest.approx(0.00, 1e-10)

    # run models for 7 days starting from Jan 1
    model1.simulate()
    model2.simulate()

    # calc drybulb mean
    temps1 = [ucm.canTemp - 273.15 for ucm in model1.UCMData]
    temps2 = [ucm.canTemp - 273.15 for ucm in model2.UCMData]

    mean_temps1, mean_temps2 = sum(
        temps1) / len(temps1), sum(temps2) / len(temps2)

    # confirm if increasing grass result in lower temps
    assert mean_temps1 > mean_temps2


def test_tree_sens_heat():
    """Test treeSensHeat."""

    model = UWG.from_param_args(
        bldheight=10, blddensity=0, vertohor=0.8, zone='1A', nday=3,
        vegstart=1, vegend=12, latgrss=0.4, lattree=0.6, albveg=0.5,
        treecover=0, grasscover=0, epw_path=EPW_PATH)

    # override model1
    model.blddensity = 0.5
    model.treecover = 0.1  # 20% tree
    model.grasscover = 0.4
    # model1.vegcover = 0.5

    # run models for 3 days starting from Jan 1
    model.generate()
    model.simulate()

    # check tree coverage
    assert model.UCM.roadShad == pytest.approx(0.2, 1e-10)

    hr = (24 * 2) + 12  # afternoon on last day

    # treeSensHeat is just Solar on road, even when fraction is 0
    sol = model.UCMData[hr].SolRecRoad
    sens_chk = sol * 0.4 * 0.5 * 0.1  # sensFrac*alb*treeFrac = tree sens
    sens_chk += sol * 0.6 * 0.5 * 0.4  # sensFrac*alb*grassFrac = grass sens
    # Note that sensible heat accounts for fraction of tree/grass relative to
    # entire urban area, not relative to road, since SolRecRoad doesn't account
    # for road fraction (just road view factor)
    assert model.UCMData[hr].treeSensHeat == pytest.approx(sens_chk, abs=1e-10)

    # with zero coerage
    model = UWG.from_param_args(
        bldheight=10, blddensity=0, vertohor=0.8, zone='1A', nday=3,
        vegstart=1, vegend=12, latgrss=0.4, lattree=0.6, albveg=0.5,
        treecover=0, grasscover=0, epw_path=EPW_PATH)

    # override model1
    model.blddensity = 0.5
    model.treecover = 0
    model.grasscover = 0  # TODO: make into grasscover
    # model1.vegcover = 0

    # run models for 3 days starting from Jan 1
    model.generate()
    model.simulate()

    # check tree coverage
    assert model.UCM.roadShad == pytest.approx(0.0, 1e-10)

    hr = (24 * 2) + 12  # afternoon on last day

    # treeSensHeat is just Solar on road, even when fraction is 0
    sol = model.UCMData[hr].SolRecRoad
    assert model.UCMData[hr].treeSensHeat == pytest.approx(0, abs=1e-10)


def test_veg_element():
    """Test veg Element calcs"""

    model = UWG.from_param_args(
        bldheight=10, blddensity=0.5, vertohor=0.8, zone='1A', nday=1,
        treecover=0.125, grasscover=0.125, epw_path=EPW_PATH)

    model.vegstart = 1
    model.vegend = 12
    model.lattree = 0.5
    model.albveg = 0.5
    model.latgrss = 0.5
    model.lattree = 0.5
    model.albroad = 0.25

    # simulate
    model.generate()
    model.simulate()

    assert model.vegcover == pytest.approx(0.25, 1e-10)
    assert model.road.vegcoverage == pytest.approx(0.5, 1e-10)
    assert model.UCM.vegcover == pytest.approx(model.vegcover, 1e-10)

    # Summer, veg
    solRec = 10
    vegcoverage = 0.5
    albedo = 0.25
    grassFLat = 0.5
    vegAlbedo = 0.5

    solAbs = (1.0 - vegcoverage) * (1.0 - albedo) * solRec
    solAbs += vegcoverage * (1.0 - vegAlbedo) * solRec
    vegSens = vegcoverage * (1.0 - vegAlbedo) * solRec * (1.0 - grassFLat)

    solref = (1.0 - vegcoverage) * (albedo) * solRec
    solref += vegcoverage * (vegAlbedo) * solRec

    # check sensible & net heat flux
    sens = vegSens  # + aeroCond * (layerTemp[0] - tempRef)
    flux = -sens + solAbs  # + infra - lat  # [W m-2]

    # print to understand interactions
    # print('solrefl:', solref)
    # print('% solAbs:', (0.5 * 0.75) + (0.5 * 0.75))
    # print('solAbs: ', solAbs)
    # print('% vegSens:', 0.5 * 0.5 * 0.5)
    # print('vegSens: ', vegSens)
    # print('sens: ', sens)
    # print('flux: ', flux)

    UCM = model.UCMData[12]  # afternoon
    sol = model.UCMData[12].SolRecRoad
    treeSensHeat = sol * 0.5 * 0.5 * model.vegcover
    assert treeSensHeat == pytest.approx(UCM.treeSensHeat, 1e-10)


def test_veg_roof_sens():
    """Test veg roof calcs.
    """

    model1 = UWG.from_param_args(
        bldheight=10, blddensity=0.5, vertohor=0.8, zone='1A', month=6,
        nday=31, treecover=0.0, grasscover=0.0, vegstart=1, vegend=12, dtsim=300,
        epw_path=EPW_PATH)
    model1.vegroof = 0.1
    model1.latgrss = 0.4

    model2 = UWG.from_param_args(
        bldheight=10, blddensity=0.5, vertohor=0.8, zone='1A', month=6,
        nday=31, treecover=0.0, grasscover=0.0, vegstart=1, vegend=12, dtsim=300,
        epw_path=EPW_PATH)
    model2.vegroof = 1
    model2.latgrss = 0.4

    # generate
    model1.generate()
    model2.generate()

    print(model1.BEM)

    # check fractions
    for i in range(2):
        assert model1.BEM[i].roof.vegcoverage == pytest.approx(0.1, 1e-10)
        assert model2.BEM[i].roof.vegcoverage == pytest.approx(1, 1e-10)

    # simulate
    model1.simulate()
    model2.simulate()

    # Roof vegetation sensible heat is greater with more vegetation, but convection
    # here mitigates that. This only works b/c latgrss is very low so that veg
    # sens heat is > then convection. For checking only, no test.
    # qroof1 = [ucm.Q_roof for ucm in model1.UCMData]
    # qroof2 = [ucm.Q_roof for ucm in model2.UCMData]
    # mean_qroof1, mean_qroof2 = sum(qroof1) / len(qroof1), sum(qroof2) / len(qroof2)
    # print(mean_qroof1, mean_qroof2)

    # calc drybulb mean
    temps1 = [ucm.canTemp - 273.15 for ucm in model1.UCMData]
    temps2 = [ucm.canTemp - 273.15 for ucm in model2.UCMData]

    mean_temps1, mean_temps2 = sum(
        temps1) / len(temps1), sum(temps2) / len(temps2)

    # Increasing roof vegetation results in slightly lower temps, as long as
    # latgrss is very low, which means sens heat is > then convection.
    # print(mean_temps1, mean_temps2)
    assert mean_temps1 > mean_temps2


def test_grass_sens():
    """Test changing grass sens makes a difference.

    This tests if changing grass sensible heat fraction on tree sensible heat,
    urban sensible heat, and drybulb.

    While rural grass sensible heat is accounted for in UWG, the code from
    the UWG_Matlab doesn't account for grass sensible heat on the urban
    road.

    To verify that the UHI doesn't change with changes to the sensible heat
    fraction of grass, comment out the following code to eliminate the
    impact of the rural road grass, and revert treeSensHeat to original
    UWG_Matlab code:
        At UBLDef.ublmodel  // set heatDif = max(self.sensHeat - rural.sens, 0) to 0
        At RSMDef.vdm  //set rural.sens to 1.0
    """

    model1 = UWG.from_param_args(
        bldheight=10, blddensity=0.0, vertohor=0.8, grasscover=0.1, treecover=0.1,
        zone='1A', nday=31, vegstart=1, vegend=12, dtsim=600, epw_path=EPW_PATH)
    # 100% of open space is grass with 0.9 sensible heat
    model1.vegroof = 0.0
    model1.blddensity = 0.5
    model1.treecover = 0.0
    model1.grasscover = 0.4
    model1.latgrss = 0.3

    model2 = UWG.from_param_args(
        bldheight=10, blddensity=0.0, vertohor=0.8, grasscover=0.1, treecover=0.1,
        zone='1A', nday=31, vegstart=1, vegend=12, dtsim=600, epw_path=EPW_PATH)
    # 100% of open space is grass with 0.1 sensible heat
    model2.vegroof = 0.0
    model2.blddensity = 0.5
    model2.treecover = 0.0
    model2.grasscover = 0.4
    model2.latgrss = 0.7

    # model1 should have produce more vegetation sensible heat.

    # generate & simulate
    model1.generate()
    model1.simulate()
    model2.generate()
    model2.simulate()

    # confirm tree/veg coverage has more sensible veg contribution with more sens factor.
    # The solar blocked by veg will be more then sensible heat fraction, so net effect
    # is less heat to canopy.
    qtree1 = [ucm.treeSensHeat for ucm in model1.UCMData]
    qtree2 = [ucm.treeSensHeat for ucm in model2.UCMData]
    mean_qtree1, mean_qtree2 = sum(
        qtree1) / len(qtree1), sum(qtree2) / len(qtree2)
    # print(mean_qtree1, mean_qtree2)
    assert mean_qtree1 > mean_qtree2

    # confirm canopy sens heat
    qsens1 = [ucm.sensHeat for ucm in model1.UCMData]
    qsens2 = [ucm.sensHeat for ucm in model2.UCMData]
    mean_qsens1, mean_qsens2 = sum(
        qsens1) / len(qsens1), sum(qsens2) / len(qsens2)
    # print(mean_qsens1, mean_qsens2)
    assert mean_qsens1 > mean_qsens2

    # calc drybulb mean
    temps1 = [ucm.canTemp - 273.15 for ucm in model1.UCMData]
    temps2 = [ucm.canTemp - 273.15 for ucm in model2.UCMData]
    mean_temps1, mean_temps2 = sum(
        temps1) / len(temps1), sum(temps2) / len(temps2)
    # print(mean_temps1, mean_temps2)
    assert mean_temps1 > mean_temps2  # higher latent fraction makes cooler temps
