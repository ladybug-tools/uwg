
from uwg_test import UWG_Test
from building import Building
from material import Material
from element import Element
from BEMDef import BEMDef
from schdef import SchDef
from simparam import SimParam
from weather import Weather
from param import Param
from UCMDef import UCMDef

from math import pow, pi


# Physical constants (moved here from UWG.m)
g = 9.81                # gravity
cp = 1004.              # heat capacity for air (J/kg.K)
vk = 0.40               # von karman constant
r = 287.                # gas constant for air
rv = 461.5              # gas constant for vapour
lv = 2.26e6             # latent heat of evaporation
sigma = 5.67e-08        # Stefan Boltzmann constant
waterDens = 1000.       # water density (kg/m^3)
lvtt = 2.5008e6         #
tt = 273.16             #
estt = 611.14           #
cl = 4.218e3            #
cpv = 1846.1            #
b = 9.4                 # Coefficients derived by Louis (1979)
cm = 7.4                #
colburn = pow((0.713/0.621),(2./3.)) # (Pr/Sc)^(2/3) for Colburn analogy in water evaporation

def sim_singapore():
    test_uwg_param = UWG_Test("singapore_test", True)

    # -------------------------------------------------------------------------
    # Material
    # -------------------------------------------------------------------------
    # Material: [conductivity (W m-1 K-1), Vol heat capacity (J m-3 K-1)]
    bldMat = Material(0.67,1.2e6)      # material (concrete? reference?)
    roadMat = Material(1.0,1.6e6)      # material (asphalt? reference?)

    if test_uwg_param.run_test==True:
        print "INIT MATERIALS"
        print '\t', bldMat
        print '\t', roadMat

    # -------------------------------------------------------------------------
    # Elements
    # -------------------------------------------------------------------------
    # Element: [albedo, emissivity, thicknesses (m)(outer layer first),
    # materials, vegetation coverage, initial temperature (K),
    # inclination (horizontal - 1, vertical - 0) ]
    wall = Element(0.2,0.9,[0.01,0.05,0.1,0.05,0.01],\
        [bldMat,bldMat,bldMat,bldMat,bldMat],0.,300.,0)
    roof = Element(0.2,0.9,[0.01,0.05,0.1,0.05,0.01],\
        [bldMat,bldMat,bldMat,bldMat,bldMat],0.,300.,1)
    road = Element(0.5,0.95,[0.05,0.1,0.1,0.5,0.5],\
        [roadMat,roadMat,roadMat,roadMat,roadMat],0.2,300.,1)
    #rural = Element(0.1,0.95,[0.05,0.1,0.1,0.5,0.5],\
    #    [roadMat,roadMat,roadMat,roadMat,roadMat],0.73,300.,1)
    mass = Element(0.7,0.9,[0.05,0.05],[bldMat,bldMat],0.,300.,0)

    if test_uwg_param.run_test==True:
        print "INIT ELEMENTS"
        print '\t', 'wall', wall
        print '\t', 'roof', roof
        print '\t', 'road', road
        #print '\t', 'rural', rural
        print '\t', 'mass', mass
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

    if test_uwg_param.run_test==True:
        print "INIT SIMPARAM"
        print '\t', simTime

    # Simulation Parameters tests
    test_uwg_param.test_equality_tol(simTime.timeSim,168,True)
    test_uwg_param.test_equality_tol(simTime.timeMax,604800,True)
    test_uwg_param.test_equality_tol(simTime.nt,2017,True)

    # -------------------------------------------------------------------------
    # Weather
    # -------------------------------------------------------------------------

    # Read Rural weather data (EPW file - http://apps1.eere.energy.gov/)
    weather_ = Weather(climate_file,simTime.timeInitial,simTime.timeFinal)

    if test_uwg_param.run_test==True:
        print "INIT WEATHER"
        print '\t', weather_

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

    # -------------------------------------------------------------------------
    # Building
    # -------------------------------------------------------------------------

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
    if test_uwg_param.run_test==True:
        print "INIT BUILDING"
        print '\t', res_wAC

    # -------------------------------------------------------------------------
    # Urban MicroClimate Parameters
    # -------------------------------------------------------------------------

    # Urban microclimate parameters from initialize.uwg
    h_ubl1 = 1000.          # ubl height - day (m)
    h_ubl2 = 80.            # ubl height - night (m)
    h_ref = 150.            # inversion height (m)
    h_temp = 2.             # temperature height (m)
    h_wind = 10.            # wind height (m)
    c_circ = 1.2            # circulation coefficient
    c_exch = 1.0            # exchange coefficient
    maxDay = 150.           # max day threshhold heat flux (W/m^2)
    maxNight = 20.          # max night threshhold (W/m^2)
    windMin = 1.0           # min wind speed (m/s)
    h_obs = 0.1             # rural average obstacle height (m)

     # Vegetatin parameters
    vegCover = 0.2          # urban area veg coverage ratio
    treeCoverage = 0.1      # urban area tree coverage ratio
    vegStart = 4.           # vegetation start month
    vegEnd = 10.            # vegetation end month
    albVeg = 0.25           # Vegetation albedo
    latGrss = 0.5           # latent fraction of grass
    latTree = 0.5           # latent fraction of tree
    rurVegCover = 0.9       # rural vegetation cover

    nightStart = 18         # begin hour for night thermal set point schedule
    nightEnd = 8            # end hour for night thermal set point schedule

    # Site-specific parameters
    wgmax = 0.005           # maximum film water depth on horizontal surfaces (m)
    maxdx = 250             # Max Dx (m)

    geoParam = Param(h_ubl1,h_ubl2,h_ref,h_temp,h_wind,c_circ,maxDay,maxNight,
        latTree,latGrss,albVeg,vegStart,vegEnd,nightStart,nightEnd,windMin,wgmax,c_exch,maxdx,
        g, cp, vk, r, rv, lv, pi, sigma, waterDens, lvtt, tt, estt, cl, cpv, b, cm, colburn)

    if test_uwg_param.run_test==True:
        print "INIT PARAM"
        print '\t', geoParam

    # -------------------------------------------------------------------------
    # RSM/UCM/UBL
    # -------------------------------------------------------------------------

    # Define Reference (RSMDef(lat,lon,height,initialTemp,initialPres,Param))
    #RSM = RSMDef(LAT,LON,ELEV,weather.staTemp(1),weather.staPres(1),Param)

    T_init = weather_.staTemp[0]     # start dry bulb
    Hum_init = weather_.staHum[0]    # start relative humidity
    Wind_init = weather_.staUmod[0]  # wind speed

    # Urban characteristics
    bldHeight = 10          # average building height (m)
    h_mix = 1               # fraction of waste heat to canyon
    bldDensity = 0.5        # urban area building plan density (0-1)
    verToHor = 0.8          # urban area vertical to horizontal ratio
    h_floor = 3.05          # average floor height
    charLength = 1000       # urban area characteristic length (m)
    alb_road = 0.2          # road albedo (0 - 1)
    d_road = 0.5            # road pavement thickness (m)
    sensAnth = 20           # non-building sens heat (W/m^2)
    latAnth = 2             # non-building latent heat (W/m^2) (currently not used)

    # In UWG.py this is done as average of all refDOE types
    r_glaze = 0.3     # glazing ratio from Building
    SHGC = 0.75       # from Building
    alb_wall = 0.2    # from wall Element

    UCM_ = UCMDef(bldHeight,bldDensity,verToHor,treeCoverage,sensAnth,latAnth,
        T_init,Hum_init,Wind_init,geoParam,r_glaze,SHGC,alb_wall,road)#,rural)

    if test_uwg_param.run_test==True:
        print "INIT UCM"
        print '\t', UCM_


    #UBL = UBLDef('C',1000.,weather.staTemp(1),Param.maxdx),

    # -------------------------------------------------------------------------
    # Method Testing
    # -------------------------------------------------------------------------
    """
    #res_wAC.BEMCalc(UCM,res_wAC,forc,parameter,simTime)

    print test_uwg_param.test_results()
    """
if __name__ == "__main__":

    sim_singapore()
