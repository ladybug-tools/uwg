"""Test the readDOE module."""

import os
import math
import copy
from functools import reduce
import pytest

from .test_base import calculate_tolerance, MATLAB_DIR
from uwg.readDOE import readDOE


MATLAB_DIR = os.path.join(MATLAB_DIR, 'matlab_readDOE')

# Test matrix of ref data as nested nested lists [16, 3, 16]:

# refDOE = Building objs
# Schedule = SchDef objs
# refBEM = BEMDef

# where:
#     [16,3,16] is Type = 1-16, Era = 1-3, climate zone = 1-16
#     i.e.
#     Type: FullServiceRestaurant, Era: Pre80, Zone: 6A Minneapolis
# Nested tree:
# [TYPE_1:
#     ERA_1:
#         CLIMATE_ZONE_1
#         ...
#         CLIMATE_ZONE_16
#     ERA_2:
#         CLIMATE_ZONE_1
#         ...
#         CLIMATE_ZONE_16
#     ...
#     ERA_3:
#         CLIMATE_ZONE_1
#         ...
#         CLIMATE_ZONE_16]


def test_refDOE():
    """Tests for refDOE (Building class)"""

    # Run readDOE.py
    refBEM, Schedule = readDOE(serialize_output=False)

    bldlst = ['floorHeight', 'intHeat', 'intHeatNight', 'intHeatDay', 'intHeatFRad',
              'intHeatFLat', 'infil', 'vent', 'glazingRatio', 'uValue', 'shgc',
              'condType', 'cop', 'coolSetpointDay', 'coolSetpointNight',
              'heatSetpointDay', 'heatSetpointNight', 'coolCap', 'heatEff', 'mSys',
              'indoorTemp', 'indoorHum', 'heatCap', 'copAdj', 'fanMax']

    for bi in range(len(bldlst)):
        bldid = bldlst[bi]
        matlab_path = os.path.join(MATLAB_DIR, 'matlab_ref_{}_.txt'.format(bldid))

        # check file exists
        if not os.path.exists(matlab_path):
            raise Exception('Failed to open {}!'.format(matlab_path))

        matlab_file = open(matlab_path, 'r')

        refDOE = refBEM

        assert len(refDOE) == pytest.approx(16.0, abs=1e-3)

        for bldType in range(16):         # bldType
            assert len(refDOE[bldType]) == pytest.approx(3.0, abs=1e-3)

            for bldEra in range(3):        # bltEra
                assert len(refDOE[bldType][bldEra]) == pytest.approx(16.0, abs=1e-3)

                for climateZone in range(16):  # ZoneType

                    matlab_ref_value = matlab_file.readline()

                    if bldid != 'condType':
                        matlab_ref_value = float(matlab_ref_value)
                        tol = calculate_tolerance(matlab_ref_value, 14.0)
                    else:
                        matlab_ref_value = 'AIR'

                    bldg = refDOE[bldType][bldEra][climateZone].building

                    # run tests
                    if bldid == 'floorHeight':
                        assert bldg.floor_height == pytest.approx(matlab_ref_value, abs=tol), \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                    elif bldid == 'intHeat':
                        assert bldg.int_heat == pytest.approx(matlab_ref_value, abs=tol), \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                    elif bldid == 'intHeatNight':
                        assert bldg.int_heat_night == pytest.approx(matlab_ref_value, abs=tol), \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                    elif bldid == 'intHeatDay':
                        assert bldg.int_heat_day == pytest.approx(matlab_ref_value, abs=tol), \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                    elif bldid == 'intHeatFRad':
                        assert bldg.int_heat_frad == pytest.approx(matlab_ref_value, abs=tol), \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                    elif bldid == 'intHeatFLat':
                        assert bldg.int_heat_flat == pytest.approx(matlab_ref_value, abs=tol), \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                    elif bldid == 'infil':
                        assert bldg.infil == pytest.approx(matlab_ref_value,  abs=tol),  \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                    elif bldid == 'vent':
                        assert bldg.vent == pytest.approx(matlab_ref_value,  abs=tol), \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                    elif bldid == 'glazingRatio':
                        assert bldg.glazing_ratio == pytest.approx(matlab_ref_value, abs=tol), \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                    elif bldid == 'uValue':
                        assert bldg.u_value == pytest.approx(matlab_ref_value, abs=tol), \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                    elif bldid == 'shgc':
                        assert bldg.shgc == pytest.approx(matlab_ref_value,  abs=tol), \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                    elif bldid == 'condType':
                        assert bldg.condtype == matlab_ref_value,  \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                    elif bldid == 'cop':
                        assert bldg.cop == pytest.approx(matlab_ref_value,  abs=tol), \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                    elif bldid == 'coolSetpointDay':
                        assert bldg.cool_setpoint_day == pytest.approx(matlab_ref_value, abs=tol), \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                    elif bldid == 'coolSetpointNight':
                        assert bldg.cool_setpoint_night == pytest.approx(matlab_ref_value, abs=tol), \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                    elif bldid == 'heatSetpointDay':
                        assert bldg.heat_setpoint_day == pytest.approx(matlab_ref_value, abs=tol), \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                    elif bldid == 'heatSetpointNight':
                        assert bldg.heat_setpoint_night == pytest.approx(matlab_ref_value, abs=tol), \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                    elif bldid == 'coolCap': # knocking the testing tolerance one order of magnitude down due to excel extraction problem
                        assert bldg.coolcap == pytest.approx(matlab_ref_value, abs=10 ** (math.log10(tol) + 1)), \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                    elif bldid == 'heatEff':
                        assert bldg.heateff == pytest.approx(matlab_ref_value, abs=tol), \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                    elif bldid == 'mSys':
                        assert bldg.msys == pytest.approx(matlab_ref_value, abs=tol), \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                    elif bldid == 'indoorTemp':
                        assert bldg.indoor_temp == pytest.approx(matlab_ref_value, abs=tol), \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                    elif bldid == 'indoorHum':
                        assert bldg.indoor_hum == pytest.approx(matlab_ref_value, abs=tol), \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                    elif bldid == 'heatCap': # knocking the testing tolerance one order of magnitude down due to excel extraction problem
                        assert bldg.heat_cap == pytest.approx(matlab_ref_value, abs=10 ** (math.log10(tol) + 1)), \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                    elif bldid == 'copAdj':
                        assert bldg.cop_adj == pytest.approx(matlab_ref_value, abs=tol), \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                    elif bldid == 'canyon_fraction':
                        assert bldg.canyon_fraction == pytest.approx(matlab_ref_value,  abs=tol), \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                    elif bldid == 'fanMax':
                        assert refBEM[bldType][bldEra][climateZone].building.FanMax == pytest.approx(matlab_ref_value,  abs=tol), \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

        matlab_file.close()


def test_refBEM():
    """Tests for refBEM (BEMDef class)."""

    refBEM, Schedule = readDOE(serialize_output=False)

    elementlst = [
        'albedo', 'emissivity', 'layerThickness', 'layerThermalCond', 'layerVolHeat',
        'vegCoverage', 'layerTemp', 'horizontal']

    bemlst = [
        ('mass', copy.deepcopy(elementlst)),
        ('wall', copy.deepcopy(elementlst)),
        ('roof', copy.deepcopy(elementlst))]

    for bemi in range(len(bemlst)):
        for ei in range(len(elementlst)):
            bemid = bemlst[bemi][0] + '_' + bemlst[bemi][1][ei]

            matlab_path = \
                os.path.join(MATLAB_DIR, 'matlab_ref_bemdef_{}.txt'.format(bemid))

            # check file exists
            if not os.path.exists(matlab_path):
                raise Exception('Failed to open {}!'.format(matlab_path))

            matlab_file = open(matlab_path, 'r')

            assert len(refBEM) == pytest.approx(16.0, abs=1e-3)

            for bldType in range(16):         # bldType
                assert len(refBEM[bldType]) == pytest.approx(3.0, abs=1e-3)

                for bldEra in range(3):        # bltEra
                    assert len(refBEM[bldType][bldEra]) == pytest.approx(16.0, abs=1e-3)

                    for climateZone in range(16):  # ZoneType
                        # next line
                        matlab_ref_value = float(matlab_file.readline())
                        tol = calculate_tolerance(matlab_ref_value, 15.0)

                        # run tests
                        if bemid == 'mass_albedo':
                            assert refBEM[bldType][bldEra][climateZone].mass.albedo == pytest.approx(matlab_ref_value,  abs=tol), \
                                'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bemid == 'mass_emissivity':
                            assert refBEM[bldType][bldEra][climateZone].mass.emissivity == pytest.approx(matlab_ref_value,  abs=tol), \
                                'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bemid == 'mass_layerThickness':
                            assert refBEM[bldType][bldEra][climateZone].mass.layer_thickness_lst[0] == pytest.approx(matlab_ref_value, abs=tol), \
                                'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bemid == 'mass_layerThermalCond':
                            assert refBEM[bldType][bldEra][climateZone].mass.layerThermalCond[0] == pytest.approx(matlab_ref_value,  abs=tol), \
                                'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bemid == 'mass_layerVolHeat':
                            assert refBEM[bldType][bldEra][climateZone].mass.layerVolHeat[0] == pytest.approx(matlab_ref_value,  abs=tol), \
                                'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bemid == 'mass_vegCoverage':
                            assert refBEM[bldType][bldEra][climateZone].mass.vegcoverage == pytest.approx(matlab_ref_value, abs=tol), \
                                'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bemid == 'mass_layerTemp':
                            assert refBEM[bldType][bldEra][climateZone].mass.layerTemp[0] == pytest.approx(matlab_ref_value,  abs=tol), \
                                'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bemid == 'mass_horizontal':
                            assert refBEM[bldType][bldEra][climateZone].mass.horizontal == pytest.approx(matlab_ref_value,  abs=tol), \
                                'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bemid == 'wall_albedo':
                            assert refBEM[bldType][bldEra][climateZone].wall.albedo == pytest.approx(matlab_ref_value,  abs=tol), \
                                'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bemid == 'wall_emissivity':
                            assert refBEM[bldType][bldEra][climateZone].wall.emissivity == pytest.approx(matlab_ref_value,  abs=tol), \
                                'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bemid == 'wall_layerThickness':
                            assert refBEM[bldType][bldEra][climateZone].wall.layer_thickness_lst[0] == pytest.approx(matlab_ref_value, abs=tol), \
                                'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bemid == 'wall_layerThermalCond':
                            assert refBEM[bldType][bldEra][climateZone].wall.layerThermalCond[0] == pytest.approx(matlab_ref_value,  abs=tol), \
                                'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bemid == 'wall_layerVolHeat':
                            assert refBEM[bldType][bldEra][climateZone].wall.layerVolHeat[0] == pytest.approx(matlab_ref_value,  abs=tol), \
                                'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bemid == 'wall_vegCoverage':
                            assert refBEM[bldType][bldEra][climateZone].wall.vegcoverage == pytest.approx(matlab_ref_value, abs=tol), \
                                'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bemid == 'wall_layerTemp':
                            assert refBEM[bldType][bldEra][climateZone].wall.layerTemp[0] == pytest.approx(matlab_ref_value,  abs=tol), \
                                'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bemid == 'wall_horizontal':
                            assert refBEM[bldType][bldEra][climateZone].wall.horizontal == pytest.approx(matlab_ref_value,  abs=tol), \
                                'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bemid == 'roof_albedo':
                            assert refBEM[bldType][bldEra][climateZone].roof.albedo == pytest.approx(matlab_ref_value,  abs=tol), \
                                'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bemid == 'roof_emissivity':
                            assert refBEM[bldType][bldEra][climateZone].roof.emissivity == pytest.approx(matlab_ref_value,  abs=tol), \
                                'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bemid == 'roof_layerThickness':
                            assert refBEM[bldType][bldEra][climateZone].roof.layer_thickness_lst[0] == pytest.approx(matlab_ref_value, abs=tol), \
                                'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bemid == 'roof_layerThermalCond':
                            assert refBEM[bldType][bldEra][climateZone].roof.layerThermalCond[0] == pytest.approx(matlab_ref_value,  abs=tol), \
                                'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bemid == 'roof_layerVolHeat':
                            assert refBEM[bldType][bldEra][climateZone].roof.layerVolHeat[0] == pytest.approx(matlab_ref_value,  abs=tol), \
                                'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bemid == 'roof_vegCoverage':
                            assert refBEM[bldType][bldEra][climateZone].roof.vegcoverage == pytest.approx(matlab_ref_value, abs=tol), \
                                'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bemid == 'roof_layerTemp':
                            assert refBEM[bldType][bldEra][climateZone].roof.layerTemp[0] == pytest.approx(matlab_ref_value,  abs=tol), \
                                'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bemid == 'roof_horizontal':
                            assert refBEM[bldType][bldEra][climateZone].roof.horizontal == pytest.approx(1, abs=tol), \
                                'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        else:
                            assert refBEM[bldType][bldEra][climateZone].building.FanMax == pytest.approx(matlab_ref_value,  abs=tol), \
                                'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

    matlab_file.close()


def test_Schedule():
    """Tests for Schedule (Schdef class)."""

    refBEM, Schedule = readDOE(serialize_output=False)

    schlst = ['Elec', 'Light', 'Gas', 'Occ', 'Cool', 'Heat', 'SWH', 'Qelec', 'Qlight',
              'Nocc', 'Qgas', 'Vent', 'Vswh']

    for si in range(len(schlst)):
        schid = schlst[si]

        # Define file path
        matlab_path = os.path.join(MATLAB_DIR,'matlab_ref_sch_{}_.txt'.format(schid))
        # check file exists
        if not os.path.exists(matlab_path):
            raise Exception('Failed to open {}!'.format(matlab_path))

        matlab_file = open(matlab_path,'r')

        assert len(Schedule) == pytest.approx(16.0, abs=1e-3)

        for bldType in range(1): # bldType
            assert len(Schedule[bldType]) == pytest.approx(3.0, abs=1e-3)

            for bldEra in range(1): # bltEra
                assert len(Schedule[bldType][bldEra]) == pytest.approx(16.0, abs=1e-3)

                for climateZone in range(1): # ZoneType

                    if schid in schlst[:7]: #3x24 matrix of schedule for SWH (WD,Sat,Sun)
                        matlab_ref_str = matlab_file.readline()
                        # replace [,] -> '' && split at ' '
                        matlab_ref_str = matlab_ref_str.replace('[','').replace(']','').replace(';','_').replace(' ', '_')
                        matlab_ref_str = ''.join(matlab_ref_str.split())
                        matlab_ref_str = matlab_ref_str.split('_')
                        matlab_ref_value = [float(x) for x in matlab_ref_str]
                        # can compare equality of flat lists [0.12, 4, 5,1.12, 5, 6] == pytest.approx([0.124, 4., 5, 1.125, 5, 6], abs=1e-2)
                        assert len(matlab_ref_value) == pytest.approx(24.0*3.,abs=1e-3)

                        # Calculate tolerances and set minimum one for all
                        tol = min([calculate_tolerance(x,15.0) for x in matlab_ref_value])
                    else:
                        matlab_ref_value = float(matlab_file.readline())
                        tol = calculate_tolerance(matlab_ref_value, 15.0)

                    if schid == 'Elec': #3x24 matrix of schedule for SWH (WD,Sat,Sun)
                        flatten_sch = reduce(lambda x,y: x+y, Schedule[bldType][bldEra][climateZone].elec)
                        assert flatten_sch == pytest.approx(matlab_ref_value,  abs=tol), \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                    elif schid == 'Light': #3x24 matrix of schedule for SWH (WD,Sat,Sun)
                        flatten_sch = reduce(lambda x,y: x+y, Schedule[bldType][bldEra][climateZone].light)
                        assert flatten_sch == pytest.approx(matlab_ref_value,  abs=tol), \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                    elif schid == 'Gas': #3x24 matrix of schedule for SWH (WD,Sat,Sun)
                        flatten_sch = reduce(lambda x,y: x+y, Schedule[bldType][bldEra][climateZone].gas)
                        assert flatten_sch == pytest.approx(matlab_ref_value,  abs=tol), \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                    elif schid == 'Occ': #3x24 matrix of schedule for SWH (WD,Sat,Sun)
                        flatten_sch = reduce(lambda x,y: x+y, Schedule[bldType][bldEra][climateZone].occ)
                        assert flatten_sch == pytest.approx(matlab_ref_value,  abs=tol), \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                    elif schid == 'Cool': #3x24 matrix of schedule for SWH (WD,Sat,Sun)
                        flatten_sch = reduce(lambda x,y: x+y, Schedule[bldType][bldEra][climateZone].cool)
                        assert flatten_sch == pytest.approx(matlab_ref_value,  abs=tol), \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                    elif schid == 'Heat': #3x24 matrix of schedule for SWH (WD,Sat,Sun)
                        flatten_sch = reduce(lambda x,y: x+y, Schedule[bldType][bldEra][climateZone].heat)
                        assert flatten_sch == pytest.approx(matlab_ref_value,  abs=tol), \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                    elif schid == 'SWH': #3x24 matrix of schedule for SWH (WD,Sat,Sun)
                        flatten_sch = reduce(lambda x,y: x+y, Schedule[bldType][bldEra][climateZone].swh)
                        assert flatten_sch == pytest.approx(matlab_ref_value,  abs=tol), \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                    elif schid == 'Qelec':
                        assert Schedule[bldType][bldEra][climateZone].q_elec == pytest.approx(matlab_ref_value,  abs=tol), \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                    elif schid == 'Qlight':
                        assert Schedule[bldType][bldEra][climateZone].q_light == pytest.approx(matlab_ref_value,  abs=tol), \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                    elif schid == 'Nocc':
                        assert Schedule[bldType][bldEra][climateZone].n_occ == pytest.approx(matlab_ref_value,  abs=tol), \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                    elif schid == 'Qgas':
                        assert Schedule[bldType][bldEra][climateZone].q_gas == pytest.approx(matlab_ref_value,  abs=tol), \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                    elif schid == 'Vent':
                        assert Schedule[bldType][bldEra][climateZone].vent == pytest.approx(matlab_ref_value,  abs=tol), \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                    elif schid == 'Vswh':
                        assert Schedule[bldType][bldEra][climateZone].v_swh == pytest.approx(matlab_ref_value,  abs=tol), \
                            'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

    matlab_file.close()

