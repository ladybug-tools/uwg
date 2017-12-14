import pytest
import os
import UWG

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

    def test_refDOE(self):
        """ Tests for refDOE (Building class)

        B = Building(
            hCeiling[j],                        # floorHeight by era
            1,                                  # intHeatNight
            1,                                  # intHeatDay
            0.1,                                # intHeatFRad
            0.1,                                # intHeatFLat
            Infil[j],                           # infil (ACH) by era
            Vent[j]/1000,                       # vent (m^3/s/m^2) by era, converted from liters
            glazing[j],                         # glazing ratio by era
            Uwindow[j][k],                      # uValue by era, by climate type
            SHGC[j][k],                         # SHGC, by era, by climate type
            'AIR',                              # a/c type
            COP[j][k],                          # cop by era, climate type
            297,                                # coolSetpointDay = 24 C
            297,                                # coolSetpointNight
            293,                                # heatSetpointDay = 20 C
            293,                                # heatSetpointNight
            (HVAC[j][k]*1000.0)/AreaFloor[j],   # coolCap converted to W/m2 by era, climate type
            EffHeat[j][k],                      # heatEff by era, climate type
            293)                                # initialTemp at 20 C

        # heatCap = (HEAT[j][k]*1000.0)/AreaFloor[j]  # heating Capacity converted to W/m2 by era, climate type
        # Type = bldType[i]
        # Era = builtEra[j]
        # Zone = zoneType[k]
        # B.FanMax = FanFlow[j][k]

        """

        refDOE, refBEM, Schedule = UWG.readDOE(serialize_output=False)



        """
        #Test sheet 3
        matlab_rvalue_bld_1 = ((1.757469244288225,2.375296912114014,2.793296089385475),
            (1.757469244288225,2.666666666666667,2.793296089385475))
        #if i==1:
        #    pprint.pprint(RvalRoof)
        
        assert refDOE[0][0][0].vent == pytest.approx(0.00533624898002, abs=1e-15)
        assert refDOE[0][1][15].uValue == pytest.approx(2.956, abs=1e-6)
        assert refDOE[0][2][2].heatEff == pytest.approx(0.784213988722061, abs=1e-15)
        #if i==0 and j==0: test_treeDOE.test_equality_tol(B.vent,5.34/1000.0,toggle=True)
        #if i==1 and j==2 and k==15: test_treeDOE.test_equality_tol(refDOE[i][j][k].vent,1.8/1000.0,toggle=True)
        """
        # Location data
        """
        test_readDOE.test_equality(list_doe3[0][2],"Location Summary")
        test_readDOE.test_equality(16,len(TypeWall[0]),toggle=True)
        test_readDOE.test_equality(16,len(TypeWall[1]),toggle=True)
        test_readDOE.test_equality(16,len(TypeWall[2]),toggle=True)
        test_readDOE.test_equality(16,len(RvalWall[0]),toggle=True)
        if i==1: test_readDOE.test_in_string('MassWall', TypeWall[1][2])
        if i==1:
            for i_ in xrange(2):
                for i__ in xrange(3):
                    test_readDOE.test_equality_tol(RvalRoof[i_][i__], matlab_rvalue_bld_1[i_][i__],toggle=True)
        if i==0: test_readDOE.test_in_string('SteelFrame',TypeWall[0][0],toggle=True)
        if i==0: test_readDOE.test_equality_tol(RvalWall[0][0],0.765696784074,toggle=True)
        if i==0: test_readDOE.test_equality_tol(Uwindow[0][0],5.835,toggle=True)
        if i==0: test_readDOE.test_equality_tol(SHGC[0][11],0.407,toggle=True)
        if i==0: test_readDOE.test_equality_tol(HEAT[0][0],174.49709,toggle=True)
        if i==0: test_readDOE.test_equality_tol(FanFlow[2][1],5.67,toggle=True)
        """

    def test_refBEM(self):
        """
        """
        # test_readDOE.test_equality(len(Elec),3,toggle=True)
        pass
        """
        if i==1 and j==1 and k==15: test_treeDOE.test_equality_tol(refBEM[i][j][k].building.FanMax,101.52)

        """

    def test_Schedule(self):
        """
        """
        pass
        """
        #Test sheet 4
        test_readDOE.test_equality(list_doe4[0][2],"Schedule")
        test_readDOE.test_equality(list_doe4[0][2],"Schedule")
        test_readDOE.test_equality_tol(len(SchEquip[0]),24)
        if i==0: test_readDOE.test_equality_tol(SchEquip[1][0],0.1)
        if i==0: test_readDOE.test_equality_tol(SchSWH[2][23],0.2)
        """

if __name__ == "__main__":
    test_read_doe = TestReadDOE()
    test_read_doe.test_refDOE()
    test_read_doe.test_refBEM()
    test_read_doe.test_Schedule()
