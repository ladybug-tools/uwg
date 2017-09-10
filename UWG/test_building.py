import pytest
import os

from building import Building
from material import Material
from element import Element
from BEMDef import BEMDef
from schdef import SchDef
from simparam import SimParam
from weather import Weather
from param import Param
from UCMDef import UCMDef
from forcing import Forcing

import math

REL_EPW_PATH = "data/epw/SGP_Singapore.486980_IWEC.epw"
debug_log = True

def setup_simparam():
    """setup simparam instance"""

    dtSim = 300             # Sim time step
    dtWeather = 3600        # Weather data time-step
    MONTH = 7               # Begin month
    DAY = 30                # Begin day of the month
    NUM_DAYS = 7            # Number of days of simulation

    simparam = SimParam(dtSim,dtWeather,MONTH,DAY,NUM_DAYS)

    if debug_log==True:
        print "INIT SIMPARAM"
        print '\t', simparam

    return simparam

def setup_weather(simTime,climate_file):
    """setup weather instance"""
    weather = Weather(climate_file,simTime.timeInitial,simTime.timeFinal)

    if debug_log==True:
        print "INIT WEATHER"
        print '\t', weather

    return weather

def setup_material():
    """ setup material.py """
    # Material: [conductivity (W m-1 K-1), Vol heat capacity (J m-3 K-1)]
    bldMat = Material(0.67,1.2e6,"Concrete")      # material (concrete? reference?)
    roadMat = Material(1.0,1.6e6, "Ashphalt")     # material (asphalt? reference?)

    if debug_log==True:
        print "INIT MATERIALS"
        print '\t', bldMat
        print '\t', roadMat

    return bldMat, roadMat

def setup_element(bldMat, roadMat):
    """ setup element.py """
    # Element: [albedo, emissivity, thicknesses (m)(outer layer first),
    # materials, vegetation coverage, initial temperature (K),
    # inclination (horizontal - 1, vertical - 0) ]
    wall = Element(0.2,0.9,[0.01,0.05,0.1,0.05,0.01],\
        [bldMat,bldMat,bldMat,bldMat,bldMat],0.,300.,0,"MassWall")
    roof = Element(0.2,0.9,[0.01,0.05,0.1,0.05,0.01],\
        [bldMat,bldMat,bldMat,bldMat,bldMat],0.,300.,1,"MassRoof")
    road = Element(0.5,0.95,[0.05,0.1,0.1,0.5,0.5],\
        [roadMat,roadMat,roadMat,roadMat,roadMat],0.2,300.,1,"MassRoad")
    #rural = Element(0.1,0.95,[0.05,0.1,0.1,0.5,0.5],\
    #    [roadMat,roadMat,roadMat,roadMat,roadMat],0.73,300.,1)
    mass = Element(0.7,0.9,[0.05,0.05],[bldMat,bldMat],0.,300.,0,"MassFloor")

    if debug_log==True:
        print "INIT ELEMENTS"
        print '\t', wall
        print '\t', roof
        print '\t', road
        #print '\t', 'rural', rural
        print '\t', mass

    return wall, roof, road, mass


def setup_param():
    """ setup geoparam parameters from initialize.uwg """
    # Physical constants (moved here from UWG.m)
    # TODO: Inquire about moving this to param class
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
    colburn = math.pow((0.713/0.621),(2./3.)) # (Pr/Sc)^(2/3) for Colburn analogy in water evaporation

    # variables
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

     # Vegetation parameters
    vegCover = 0.2          # urban area veg coverage ratio
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

    geo_param = Param(h_ubl1,h_ubl2,h_ref,h_temp,h_wind,c_circ,maxDay,maxNight,
        latTree,latGrss,albVeg,vegStart,vegEnd,nightStart,nightEnd,windMin,wgmax,c_exch,maxdx,
        g, cp, vk, r, rv, lv, math.pi, sigma, waterDens, lvtt, tt, estt, cl, cpv, b, cm, colburn)

    if debug_log==True:
        print "INIT PARAM"
        print '\t', geo_param

    return geo_param

def setup_UCM(weather_,geoParam,road):
    """ setup UCM.py """
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

    # Vegetation parameters
    treeCoverage = 0.1      # urban area tree coverage ratio

    #UCM needs to be tested
    UCM = UCMDef(bldHeight,bldDensity,verToHor,treeCoverage,sensAnth,latAnth,
        T_init,Hum_init,Wind_init,geoParam,r_glaze,SHGC,alb_wall,road)#,rural)

    if debug_log==True:
        print "INIT UCM"
        print '\t', UCM

    #UBL = UBLDef('C',1000.,weather.staTemp(1),Param.maxdx),
    return UCM
    weather, param,
def setup_forcing(weather_, simTime, geoParam, UCM):
    """ setup forcing.py"""
    forcIP = Forcing(weather_.staTemp, weather_)
    forc = Forcing()

    ph = simTime.dt/3600.0      # per hour
    it = range(simTime.nt)[0]   # simTime incrment

    #print 'ph', ph
    #print 'it', it

    #TODO: unittests to understand ph it
    sim_dt_index = int(math.ceil(it*ph))
    #print 'simdtindex', sim_dt_index

    # Update the weather per UWG

    forc.infra = forcIP.infra[sim_dt_index]
    forc.wind = max(forcIP.wind[sim_dt_index], geoParam.windMin)
    forc.uDir = forcIP.uDir[sim_dt_index]
    forc.hum = forcIP.hum[sim_dt_index]
    forc.pres = forcIP.pres[sim_dt_index]
    forc.temp = forcIP.temp[sim_dt_index]
    forc.rHum = forcIP.rHum[sim_dt_index]
    forc.prec = forcIP.prec[sim_dt_index]
    forc.dir = forcIP.dir[sim_dt_index]
    forc.dif = forcIP.dif[sim_dt_index]
    UCM.canHum = forc.hum      # Canyon humidity (absolute) same as rural

    if debug_log==True:
        print "INIT FORCING"
        print '\t', forc

    return forc

def setup_building():
    """ setup building.py """
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
    res_wAC.FanMax = 10.22

    res_wAC.Type = 'MidRiseApartment'
    res_wAC.Zone = "1A (Miami)"
    res_wAC.Era = "Pst80"

    if debug_log==True:
        print "INIT BUILDING"
        print '\t', res_wAC
    return res_wAC

def setup_BEMDef(building, mass, wall, roof):
    """ setup BEMDef.py """

    res_wAC_BEM = BEMDef(building, mass, wall, roof, 0.0)

    res_wAC_BEM.Elec = 0.        # Actual electricity consumption(W/m^2)
    res_wAC_BEM.Light = 0.       # Actual light (W/m^2)
    res_wAC_BEM.Nocc = 0.        # Actual gas consumption(W/m^2)
    res_wAC_BEM.Qocc = 0.        # Actual heat load from occupant (W/m^2)
    res_wAC_BEM.SWH = 0.         # Actual hot water usage
    res_wAC_BEM.Gas = 0.         # Actual gas consumption(W/m^2)

    res_wAC_BEM.T_wallex = res_wAC_BEM.wall.layerTemp[0]    # Wall surface temp (ext)
    res_wAC_BEM.T_wallin = res_wAC_BEM.wall.layerTemp[-1]   # Wall surface temp (int)
    res_wAC_BEM.T_roofex = res_wAC_BEM.roof.layerTemp[0]    # Roof surface temp (ext)
    res_wAC_BEM.T_roofin = res_wAC_BEM.roof.layerTemp[-1]   # Roof surface temp (int)

    if debug_log==True:
        print "INIT BEMDef"
        print '\t', res_wAC_BEM
    return res_wAC_BEM

def test_building():
    """test for Building Class"""

    simparam = setup_simparam()
    dir_path = os.path.dirname(__file__)
    cf = os.path.join(dir_path,REL_EPW_PATH)
    weather = setup_weather(simparam,cf)

    geo_param = setup_param()

    bldMat, roadMat = setup_material()
    wall, roof, road, mass = setup_element(bldMat, roadMat)
    ucm = setup_UCM(weather, geo_param, road)

    forc = setup_forcing(weather, simparam, geo_param, ucm)

    # bldg defs
    building = setup_building()
    bemdef = setup_BEMDef(building, mass, wall, roof)

    # TODO: Hold off on this until finished test_UWG.py
    # tests
    #print '\n-- TESTING BEMCALC --\n'
    #bemdef.building.BEMCalc(ucm,bemdef,forc,geo_param,simparam)
    #print '\n---------------------'


if __name__ == "__main__":

    test_building()


# -------------------------------------------------------------------------
# RSM
# -------------------------------------------------------------------------
# Define Reference (RSMDef(lat,lon,height,initialTemp,initialPres,Param))
# RSM = RSMDef(LAT,LON,ELEV,weather.staTemp(1),weather.staPres(1),Param)

# -------------------------------------------------------------------------
# SolarCalcs
# -------------------------------------------------------------------------
# Update solar flux
# rural,UCM,BEM = SolarCalcs(UCM,BEM,simTime,RSM,forc,geoParam,rural)
