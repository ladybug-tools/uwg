# -------------------------------------------------------------------------
# Element/Material Definitions
# -------------------------------------------------------------------------

from uwg_test import UWG_Test
from building import Building
from material import Material
from element import Element
from BEMDef import BEMDef
from schdef import SchDef
from simparam import SimParam
from weather import Weather

def sim_singapore():
    test_uwg_param = UWG_Test("singapore_test", True)

    # Material: [conductivity (W m-1 K-1), Vol heat capacity (J m-3 K-1)]
    bldMat = Material(0.67,1.2e6)      # material (concrete? reference?)
    roadMat = Material(1.0,1.6e6)      # material (asphalt? reference?)


    # Define & build base elements
    # Element: [albedo, emissivity, thicknesses (m)(outer layer first),
    #  materials, vegetation coverage, initial temperature (K),
    #  inclination (horizontal - 1, vertical - 0) ]
    wall = Element(0.2,0.9,[0.01,0.05,0.1,0.05,0.01],\
        [bldMat,bldMat,bldMat,bldMat,bldMat],0.,300.,0)
    roof = Element(0.2,0.9,[0.01,0.05,0.1,0.05,0.01],\
        [bldMat,bldMat,bldMat,bldMat,bldMat],0.,300.,1)
    road = Element(0.5,0.95,[0.05,0.1,0.1,0.5,0.5],\
        [roadMat,roadMat,roadMat,roadMat,roadMat],0.2,300.,1)
    rural = Element(0.1,0.95,[0.05,0.1,0.1,0.5,0.5],\
        [roadMat,roadMat,roadMat,roadMat,roadMat],0.73,300.,1)
    mass = Element(0.7,0.9,[0.05,0.05],[bldMat,bldMat],0.,300.,0)

    # -------------------------------------------------------------------------
    # Simulation Parameters
    # -------------------------------------------------------------------------
    cityName = 'Singapore'   # For plot/printing
    LAT = 1.37
    LON = 103.98
    ELEV = 0.1
    dtSim = 300              # Sim time step
    dtWeather = 3600         # Weather data time-step
    monthName = 'July'       # For plot/printing
    MONTH = 7                # Begin month
    DAY = 30                 # Begin day of the month
    NUM_DAYS = 7             # Number of days of simulation
    autosize = 0             # Autosize HVAC
    ##CityBlock (8,3) = Block(NUM_DAYS * 24,wall,roof,mass,road)

    #climate_file = "rural_weather_data_changi.epw
    climate_file = "SGP_Singapore.486980_IWEC.epw"

    # Create simulation class (SimParam.m)
    simTime = SimParam(dtSim,dtWeather,MONTH,DAY,NUM_DAYS)

    # Simulation Parameters tests
    test_uwg_param.test_equality_tol(simTime.timeSim,168,True)
    test_uwg_param.test_equality_tol(simTime.timeMax,604800,True)
    test_uwg_param.test_equality_tol(simTime.nt,2017,True)

    # Read Rural weather data (EPW file - http://apps1.eere.energy.gov/)
    weather_ = Weather(climate_file,simTime.timeInitial,simTime.timeFinal)
    print weather_
    # Weather Tests
    test_uwg_param.test_equality_tol(len(weather_.staDif),simTime.timeFinal - simTime.timeInitial + 1,False)
    test_uwg_param.test_equality_tol(len(weather_.staHum),simTime.timeFinal - simTime.timeInitial + 1,False)
    test_uwg_param.test_equality_tol(len(weather_.staTemp),simTime.timeFinal - simTime.timeInitial + 1,False)
    if climate_file == "SGP_Singapore.486980_IWEC.epw":
        test_uwg_param.test_equality_tol(weather_.staTemp[3],24.+273.15,False)
        test_uwg_param.test_equality_tol(weather_.staTemp[-1],27.+273.15,False)
        test_uwg_param.test_equality_tol(weather_.staPres[10],100600.,False)
        test_uwg_param.test_equality_tol(weather_.staInfra[13],428.,False)
        test_uwg_param.test_equality_tol(weather_.staDif[6],0.,False)
        test_uwg_param.test_equality_tol(weather_.staDif[8],95.,False)
        test_uwg_param.test_equality_tol(weather_.staUdir[2],270,False)  # 270 deg
        test_uwg_param.test_equality_tol(weather_.staUmod[4],.5,False)   # 0.5 m/s
        test_uwg_param.test_equality_tol(weather_.staRobs[8],0.0,False)  # 0. mm/hr


    # Building definitions
    # Residential building with AC
    res_wAC = Building(3.0,     # floorHeight
        4.0,                    # nighttime internal heat gains (W m-2 floor)
        4.0,                    # daytime internal heat gains (W m-2 floor)
        0.2,                    # radiant fraction of internal gains
        0.2,                    # latent fraction of internal gains
        0.5,                    # Infiltration (ACH)
        0.0,                    # Ventilation (ACH)
        0.3,                    # glazing ratio
        2.715,                  # window U-value (W m-2 K)
        0.75,                   # SHGC window solar heat gain coefficient
        'AIR',                  # cooling condensation system type {'AIR','WATER'}
        2.5,                    # COP of the cooling system
        297.,                   # 24 C daytime indoor cooling set-point (K)
        297.,                   # 24 C nighttime indoor cooling set-point (K)
        293.,                   # 20 C daytime indoor heating set-point (K)
        293.,                   # 20 C nighttime indoor heating set-point (K)
        225.,                   # rated cooling system capacity (W m-2 bld)
        0.9,                    # heating system efficiency (-)
        300.)                   # 26.85 C intial indoor temp (K)

    #Add this b/c doesn't appear in current building.py
    res_wAC.canyon_fraction = 1.0     # fraction of waste heat released into the canyon
    print res_wAC

    # -------------------------------------------------------------------------
    # Urban Area Definitions (re-do this?)
    # -------------------------------------------------------------------------

    # Define Reference (RSMDef(lat,lon,height,initialTemp,initialPres,Param))
    #RSM = RSMDef(LAT,LON,ELEV,weather.staTemp(1),weather.staPres(1),Param),

    T_init = weather.staTemp(1)     # start dry bulb
    Hum_init = weather.staHum(1),   # start relative humidity
    Wind_init = weather.staUmod(1), # wind speed

    #TODO: need to add parameter to __init__ here
    UCM = UCMDef(bldHeight,bldDensity,verToHor,treeCoverage,sensAnthrop,latAnthrop,
        T_init,Hum_init,Wind_init,r_glaze,SHGC,alb_wall,road,rural),
    #UBL = UBLDef('C',1000.,weather.staTemp(1),Param.maxdx),

    """
    #res_wAC.BEMCalc(UCM,res_wAC,forc,parameter,simTime)

    print test_uwg_param.test_results()
    """
if __name__ == "__main__":

    sim_singapore()
