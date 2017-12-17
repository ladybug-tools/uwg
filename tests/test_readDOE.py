import pytest
import os
import UWG
from decimal import Decimal
import copy
dd = lambda x: Decimal.from_float(x)

class TestReadDOE(object):
    """Test for readDOE.py

    Test matrix of ref data as nested nested lists [16, 3, 16]:

    refDOE = Building objs
    Schedule = SchDef objs
    refBEM = BEMDef

    where:
        [16,3,16] is Type = 1-16, Era = 1-3, climate zone = 1-16
        i.e.
        Type: FullServiceRestaurant, Era: Pre80, Zone: 6A Minneapolis
    Nested tree:
    [TYPE_1:
        ERA_1:
            CLIMATE_ZONE_1
            ...
            CLIMATE_ZONE_16
        ERA_2:
            CLIMATE_ZONE_1
            ...
            CLIMATE_ZONE_16
        ...
        ERA_3:
            CLIMATE_ZONE_1
            ...
            CLIMATE_ZONE_16]
    """

    DIR_CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
    DIR_MATLAB_PATH = os.path.join(DIR_CURRENT_PATH,"matlab_ref","matlab_readDOE")

    def test_refDOE(self):
        """ Tests for refDOE (Building class)"""

        refDOE, refBEM, Schedule = UWG.readDOE(serialize_output=False)

        bldprop = [
            'floorHeight',
            'intHeat',
            'intHeatNight',
            'intHeatDay',
            'intHeatFRad',
            'intHeatFLat',
            'infil',
            'vent',
            'glazingRatio',
            'uValue',
            'shgc',
            'condType',
            'cop',
            'coolSetpointDay',
            'coolSetpointNight',
            'heatSetpointDay',
            'heatSetpointNight',
            'coolCap',
            'heatEff',
            'mSys',
            'indoorTemp',
            'indoorHum',
            'heatCap',
            'copAdj'
            ]

        for bpi in xrange(len(bldprop)):
            bpropid = bldprop[bpi]
            matlab_path = os.path.join(self.DIR_MATLAB_PATH,"matlab_ref_{}_.txt".format(bpropid))
            # check file exists
            if not os.path.exists(matlab_path):
                raise Exception("Failed to open {}!".format(matlab_path))

            matlab_file = open(matlab_path,'r')

            assert len(refDOE) == pytest.approx(16.0, abs=1e-3)
            for bldType in xrange(16):         # bldType
                assert len(refDOE[bldType]) == pytest.approx(3.0, abs=1e-3)

                for bldEra in xrange(3):        # bltEra
                    assert len(refDOE[bldType][bldEra]) == pytest.approx(16.0, abs=1e-3)

                    for climateZone in xrange(16):  # ZoneType
                        # next line
                        matlab_ref_value = float(matlab_file.next())
                        # run tests
                        if bpropid == 'floorHeight':
                        	assert refDOE[bldType][bldEra][climateZone].floorHeight == pytest.approx(matlab_ref_value, abs=1e-15),\
                        		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bpropid == 'intHeat':
                        	assert refDOE[bldType][bldEra][climateZone].intHeat == pytest.approx(matlab_ref_value, abs=1e-15),\
                        		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bpropid == 'intHeatNight':
                        	assert refDOE[bldType][bldEra][climateZone].intHeatNight == pytest.approx(matlab_ref_value, abs=1e-15),\
                        		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bpropid == 'intHeatDay':
                        	assert refDOE[bldType][bldEra][climateZone].intHeatDay == pytest.approx(matlab_ref_value, abs=1e-15),\
                        		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bpropid == 'intHeatFRad':
                        	assert refDOE[bldType][bldEra][climateZone].intHeatFRad == pytest.approx(matlab_ref_value, abs=1e-15),\
                        		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bpropid == 'intHeatFLat':
                        	assert refDOE[bldType][bldEra][climateZone].intHeatFLat == pytest.approx(matlab_ref_value, abs=1e-15),\
                        		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        #TODO elif bpropid == 'infil':
                        #	assert refDOE[bldType][bldEra][climateZone].infil == pytest.approx(matlab_ref_value, abs=1e-15),\
                        #		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bpropid == 'vent':
                        	assert refDOE[bldType][bldEra][climateZone].vent == pytest.approx(matlab_ref_value, abs=1e-15),\
                        		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bpropid == 'glazingRatio':
                        	assert refDOE[bldType][bldEra][climateZone].glazingRatio == pytest.approx(matlab_ref_value, abs=1e-15),\
                        		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bpropid == 'uValue':
                        	assert refDOE[bldType][bldEra][climateZone].uValue == pytest.approx(matlab_ref_value, abs=1e-15),\
                        		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bpropid == 'shgc':
                        	assert refDOE[bldType][bldEra][climateZone].shgc == pytest.approx(matlab_ref_value, abs=1e-15),\
                        		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        #TODO elif bpropid == 'condType':
                        #	assert refDOE[bldType][bldEra][climateZone].condType == pytest.approx(matlab_ref_value, abs=1e-15),\
                        #		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        #TODO elif bpropid == 'cop':
                        #	assert refDOE[bldType][bldEra][climateZone].cop == pytest.approx(matlab_ref_value, abs=1e-15),\
                        #		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bpropid == 'coolSetpointDay':
                        	assert refDOE[bldType][bldEra][climateZone].coolSetpointDay == pytest.approx(matlab_ref_value, abs=1e-15),\
                        		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bpropid == 'coolSetpointNight':
                        	assert refDOE[bldType][bldEra][climateZone].coolSetpointNight == pytest.approx(matlab_ref_value, abs=1e-15),\
                        		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bpropid == 'heatSetpointDay':
                        	assert refDOE[bldType][bldEra][climateZone].heatSetpointDay == pytest.approx(matlab_ref_value, abs=1e-15),\
                        		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bpropid == 'heatSetpointNight':
                        	assert refDOE[bldType][bldEra][climateZone].heatSetpointNight == pytest.approx(matlab_ref_value, abs=1e-15),\
                        		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        #TODO elif bpropid == 'coolCap':
                        #	assert refDOE[bldType][bldEra][climateZone].coolCap == pytest.approx(matlab_ref_value, abs=1e-15),\
                        #		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bpropid == 'heatEff':
                        	assert refDOE[bldType][bldEra][climateZone].heatEff == pytest.approx(matlab_ref_value, abs=1e-15),\
                        		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bpropid == 'mSys':
                        	assert refDOE[bldType][bldEra][climateZone].mSys == pytest.approx(matlab_ref_value, abs=1e-15),\
                        		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bpropid == 'indoorTemp':
                        	assert refDOE[bldType][bldEra][climateZone].indoorTemp == pytest.approx(matlab_ref_value, abs=1e-15),\
                        		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bpropid == 'indoorHum':
                        	assert refDOE[bldType][bldEra][climateZone].indoorHum == pytest.approx(matlab_ref_value, abs=1e-15),\
                        		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        #TODO elif bpropid == 'heatCap':
                        #	assert refDOE[bldType][bldEra][climateZone].heatCap == pytest.approx(matlab_ref_value, abs=1e-15),\
                        #		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        #TODO elif bpropid == 'copAdj':
                        #	assert refDOE[bldType][bldEra][climateZone].copAdj == pytest.approx(matlab_ref_value, abs=1e-15),\
                        #		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                        elif bpropid == 'canyon_fraction':
                        	assert refDOE[bldType][bldEra][climateZone].canyon_fraction == pytest.approx(matlab_ref_value, abs=1e-15),\
                        		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

            matlab_file.close()

    def test_refBEM(self):
        """ Tests for refBEM (BEMDef class) """

        refDOE, refBEM, Schedule = UWG.readDOE(serialize_output=False)

        elementlst=[
        'albedo',
        'emissivity',
        'layerThickness',
        'layerThermalCond',
        'layerVolHeat',
        'vegCoverage',
        'layerTemp',
        'horizontal'
        ]

        bemlst=[
        ("mass", copy.deepcopy(elementlst)),
        ("wall", copy.deepcopy(elementlst)),
        ("roof", copy.deepcopy(elementlst))
        ]

        for bemi in xrange(len(bemlst)):
            for ei in xrange(len(elementlst)):
                bemid = bemlst[bemi][0] + "_" + bemlst[bemi][1][ei]
                matlab_path = os.path.join(self.DIR_MATLAB_PATH,"matlab_ref_bemdef_{}.txt".format(bemid))
                # check file exists
                if not os.path.exists(matlab_path):
                    raise Exception("Failed to open {}!".format(matlab_path))

                matlab_file = open(matlab_path,'r')

                for bldType in xrange(16):         # bldType

                    for bldEra in xrange(3):        # bltEra

                        for climateZone in xrange(16):  # ZoneType
                            # next line
                            matlab_ref_value = float(matlab_file.next())
                            # run tests
                            if bemid == 'mass_albedo':
                            	assert refBEM[bldType][bldEra][climateZone].mass.albedo == pytest.approx(matlab_ref_value, abs=1e-15),\
                            		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                            elif bemid == 'mass_emissivity':
                            	assert refBEM[bldType][bldEra][climateZone].mass.emissivity == pytest.approx(matlab_ref_value, abs=1e-15),\
                            		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                            elif bemid == 'mass_layerThickness':
                            	assert refBEM[bldType][bldEra][climateZone].mass.layerThickness[0] == pytest.approx(matlab_ref_value, abs=1e-15),\
                            		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                            elif bemid == 'mass_layerThermalCond':
                            	assert refBEM[bldType][bldEra][climateZone].mass.layerThermalCond[0] == pytest.approx(matlab_ref_value, abs=1e-15),\
                            		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                            elif bemid == 'mass_layerVolHeat':
                            	assert refBEM[bldType][bldEra][climateZone].mass.layerVolHeat[0] == pytest.approx(matlab_ref_value, abs=1e-15),\
                            		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                            elif bemid == 'mass_vegCoverage':
                            	assert refBEM[bldType][bldEra][climateZone].mass.vegCoverage == pytest.approx(matlab_ref_value, abs=1e-15),\
                            		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                            elif bemid == 'mass_layerTemp':
                            	assert refBEM[bldType][bldEra][climateZone].mass.layerTemp[0] == pytest.approx(matlab_ref_value, abs=1e-15),\
                            		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                            elif bemid == 'mass_horizontal':
                            	assert refBEM[bldType][bldEra][climateZone].mass.horizontal == pytest.approx(matlab_ref_value, abs=1e-15),\
                            		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                            elif bemid == 'wall_albedo':
                            	assert refBEM[bldType][bldEra][climateZone].wall.albedo == pytest.approx(matlab_ref_value, abs=1e-15),\
                            		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                            elif bemid == 'wall_emissivity':
                            	assert refBEM[bldType][bldEra][climateZone].wall.emissivity == pytest.approx(matlab_ref_value, abs=1e-15),\
                            		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                            elif bemid == 'wall_layerThickness':
                            	assert refBEM[bldType][bldEra][climateZone].wall.layerThickness[0] == pytest.approx(matlab_ref_value, abs=1e-15),\
                            		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                            elif bemid == 'wall_layerThermalCond':
                            	assert refBEM[bldType][bldEra][climateZone].wall.layerThermalCond[0] == pytest.approx(matlab_ref_value, abs=1e-15),\
                            		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                            #TODO elif bemid == 'wall_layerVolHeat':
                            #	assert refBEM[bldType][bldEra][climateZone].wall.layerVolHeat[0] == pytest.approx(matlab_ref_value, abs=1e-15),\
                            #		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                            elif bemid == 'wall_vegCoverage':
                            	assert refBEM[bldType][bldEra][climateZone].wall.vegCoverage == pytest.approx(matlab_ref_value, abs=1e-15),\
                            		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                            elif bemid == 'wall_layerTemp':
                            	assert refBEM[bldType][bldEra][climateZone].wall.layerTemp[0] == pytest.approx(matlab_ref_value, abs=1e-15),\
                            		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                            elif bemid == 'wall_horizontal':
                            	assert refBEM[bldType][bldEra][climateZone].wall.horizontal == pytest.approx(matlab_ref_value, abs=1e-15),\
                            		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                            elif bemid == 'roof_albedo':
                            	assert refBEM[bldType][bldEra][climateZone].roof.albedo == pytest.approx(matlab_ref_value, abs=1e-15),\
                            		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                            elif bemid == 'roof_emissivity':
                            	assert refBEM[bldType][bldEra][climateZone].roof.emissivity == pytest.approx(matlab_ref_value, abs=1e-15),\
                            		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                            elif bemid == 'roof_layerThickness':
                            	assert refBEM[bldType][bldEra][climateZone].roof.layerThickness[0] == pytest.approx(matlab_ref_value, abs=1e-15),\
                            		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                            elif bemid == 'roof_layerThermalCond':
                            	assert refBEM[bldType][bldEra][climateZone].roof.layerThermalCond[0] == pytest.approx(matlab_ref_value, abs=1e-15),\
                            		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                            elif bemid == 'roof_layerVolHeat':
                            	assert refBEM[bldType][bldEra][climateZone].roof.layerVolHeat[0] == pytest.approx(matlab_ref_value, abs=1e-15),\
                            		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                            elif bemid == 'roof_vegCoverage':
                            	assert refBEM[bldType][bldEra][climateZone].roof.vegCoverage == pytest.approx(matlab_ref_value, abs=1e-15),\
                            		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                            elif bemid == 'roof_layerTemp':
                            	assert refBEM[bldType][bldEra][climateZone].roof.layerTemp[0] == pytest.approx(matlab_ref_value, abs=1e-15),\
                            		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                            elif bemid == 'roof_horizontal':
                            	assert refBEM[bldType][bldEra][climateZone].roof.horizontal == pytest.approx(matlab_ref_value, abs=1e-15),\
                            		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

                            elif bemid == 'roof_horizontal':
                            	assert refBEM[bldType][bldEra][climateZone].building.FanMax == pytest.approx(matlab_ref_value, abs=1e-15),\
                            		'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)

        matlab_file.close()

        # Check fanmax
        matlab_path = os.path.join(self.DIR_MATLAB_PATH,"matlab_ref_bemdef_fanmax.txt")
        # check file exists
        if not os.path.exists(matlab_path):
            raise Exception("Failed to open {}!".format(matlab_path))

        matlab_file = open(matlab_path,'r')

        for bldType in xrange(16):         # bldType
            for bldEra in xrange(3):        # bltEra
                for climateZone in xrange(16):  # ZoneType
                    #TODO
                    #assert refBEM[bldType][bldEra][climateZone].building.FanMax == pytest.approx(matlab_ref_value, abs=1e-15),\
                    #    'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)
                    pass

        matlab_file.close()

    def test_Schedule(self):
        """ Tests for Schedule (Schdef class) """

        refDOE, refBEM, Schedule = UWG.readDOE(serialize_output=False)

        elementlst=[
        'albedo',
        'emissivity',
        'layerThickness',
        'layerThermalCond',
        'layerVolHeat',
        'vegCoverage',
        'layerTemp',
        'horizontal'
        ]

if __name__ == "__main__":
    test_read_doe = TestReadDOE()
    test_read_doe.test_refDOE()
    test_read_doe.test_refBEM()
    test_read_doe.test_Schedule()
