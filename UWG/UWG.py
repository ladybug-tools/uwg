"""
=========================================================================
 THE URBAN WEATHER GENERATOR (UWG)
=========================================================================
Version 4.2

Original Author: B. Bueno
Edited by A. Nakano & Lingfu Zhang
Modified by Joseph Yang (joeyang@mit.edu) - May, 2016
Translated to Python by Chris Mackey (chris@mackeyarchitecture.com) - May, 2017

Original Pulbication on the UWG's Methods:
Bueno, Bruno; Norford, Leslie; Hidalgo, Julia; Pigeon, Gregoire (2013).
The urban weather generator, Journal of Building Performance Simulation. 6:4,269-281.
doi: 10.1080/19401493.2012.718797
=========================================================================
"""

import os
import math
import cPickle
import copy
import pprint

import utilities
from material import Material
from simparam import SimParam
from weather import Weather
from forcing import Forcing
from element import Element
from param import Param
from RSMDef import RSMDef
from UCMDef import UCMDef

class UWG(object):
    """Morph a rural EPW file to urban conditions using a file with a list of urban parameters.

    args:
        epwDir: The directory in which the rural EPW file sits.
        epwFileName: The name of the rural epw file that will be morphed.
        uwgParamDir: The directory in which the UWG Parameter File (.uwg) sits.
        uwgParamFileName: The name of the UWG Parameter File (.uwg).
        destinationDir: Optional destination directory for the morphed EPW file.
            If left blank, the morphed file will be written into the same directory
            as the rural EPW file (the epwDir).

    returns:
        newClimateFile: the path to a new EPW file that has been morphed to account
            for uban conditions.
    """

    """ Section 1 - Definitions for constants / other parameters """
    #TODO: capitalize for constant covnention?
    minThickness = 0.01    # Minimum layer thickness (to prevent crashing) (m)
    maxThickness = 0.05    # Maximum layer thickness (m)
    soilTcond = 1          # http://web.mit.edu/parmstr/Public/NRCan/nrcc29118.pdf (Figly & Snodgrass)
    soilvolHeat = 2e6      # http://www.europment.org/library/2013/venice/bypaper/MFHEEF/MFHEEF-21.pdf (average taken from Table 1)
    soil = Material(soilTcond, soilvolHeat)  # Soil material used for soil-depth padding

    # Physical constants
    g = 9.81               # gravity
    cp = 1004.             # heat capacity for air (J/kg.K)
    vk = 0.40              # von karman constant
    r = 287                # gas constant
    rv = 461.5             #
    lv = 2.26e6            # latent heat of evaporation
    sigma = 5.67e-08       # Stefan Boltzmann constant
    waterDens = 1000       # water density (kg/m^3)
    lvtt = 2.5008e6        #
    tt = 273.16            #
    estt = 611.14          #
    cl = 4.218e3           #
    cpv = 1846.1           #
    b = 9.4                # Coefficients derived by Louis (1979)
    cm = 7.4               #
    colburn = math.pow((0.713/0.621), (2/3)) # (Pr/Sc)^(2/3) for Colburn analogy in water evaporation

    # Site-specific parameters
    wgmax = 0.005          # maximum film water depth on horizontal surfaces (m)

    # File path parameter
    DIR_UP_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

    def __init__(self, epwDir, epwFileName, uwgParamDir, uwgParamFileName, destinationDir=None, destinationFile=None):
        self.epwDir = epwDir
        self.epwFileName = epwFileName
        self.uwgParamDir = uwgParamDir
        self.uwgParamFileName = uwgParamFileName
        self.destinationDir = destinationDir
        self.destinationFile = destinationFile
        self._init_param_dict = None

    def __repr__(self):
        return "UWG: {} ".format(self.epwFileName)

    def is_near_zero(self,num,eps=1e-10):
        return abs(float(num)) < eps

    def read_epw(self):
        """Section 2 - Read EPW file
        properties:
            self.newPathName
            self._header    # header data
            self.epwinput   # timestep data for weather
            self.lat        # latitude
            self.lon        # longitude
            self.GMT        # GMT
            self.nSoil      # Number of soil depths
            self.Tsoil      # nSoil x 12 matrix for soil temperture (K)
            self.depth_soil # nSoil x 1 matrix for soil depth (m)
        """

        # Revise epw file name if not end with epw
        if not self.epwFileName.lower().endswith('.epw'):
            self.epwFileName = self.epwFileName + '.epw'

        # Make dir path to epw file
        climateDataPath = os.path.join(self.epwDir, self.epwFileName)

        # Open epw file and feed csv data to climate_data
        try:
            climate_data = utilities.read_csv(climateDataPath)
        except Exception as e:
            raise Exception("Failed to read epw file! {}".format(e.message))

        # Read header lines (1 to 8) from EPW and ensure TMY2 format.
        self._header = climate_data[0:8]

        # Read weather data from EPW for each time step in weather file. (lines 8 - end)
        self.epwinput = climate_data[8:]

        # Read Lat, Long (line 1 of EPW)
        self.lat = float(self._header[0][6])
        self.lon = float(self._header[0][7])
        self.GMT = float(self._header[0][8])

        # Read in soil temperature data (assumes this is always there)
        # ref: http://bigladdersoftware.com/epx/docs/8-2/auxiliary-programs/epw-csv-format-inout.html
        soilData = self._header[3]
        self.nSoil = int(soilData[1])           # Number of ground temperature depths
        self.Tsoil = utilities.zeros(self.nSoil,12)  # nSoil x 12 matrix for soil temperture (K)
        self.depth_soil = utilities.zeros(self.nSoil,1)   # nSoil x 1 matrix for soil depth (m)

        # Read monthly data for each layer of soil from EPW file
        for i in xrange(self.nSoil):
            self.depth_soil[i][0] = float(soilData[2 + (i*16)]) # get soil depth for each nSoil
            # Monthly data
            for j in xrange(12):
                self.Tsoil[i][j] = float(soilData[6 + (i*16) + j]) + 273.15 # 12 months of soil T for specific depth

        # Set new directory path for the moprhed EPW file.
        if self.destinationDir is None:
            destinationDir = self.epwDir
        if self.destinationFile is None:
            destinationFile = self.epwFileName.lower().strip('.epw') + '_UWG.epw'
        self.newPathName = destinationDir + destinationFile

    def read_input(self):
        """Section 3 - Read Input File (.m, file)
        Note: UWG_Matlab input files are xlsm, XML, .m, file.
        properties:
            self._init_param_dict    # dictionary of simulation initialization parameters
            self.BEM                # list of BEMDef objects extracted from readDOE
            self.Sch                # list of Schedule objects extracted from readDOE

            self.simTime            # simulation time parameter obj
            self.weather            # weather obj for simulation time period
            self.forcIP             # Forcing obj
            self.forc               # Empty forcing obj
            self.geoParam           # geographic parameters obj
            self.RSM                # Rural site & vertical diffusion model obj
            self.USM                # Urban site & vertical diffusion model obj
            self.UCM                # Urban canopy model obj

            self.road               # urban road element
            self.rural              # rural road element

            self.soilindex1         # soil index for urban rsoad depth
            self.soilindex2         # soil index for rural road depth
        """

        uwg_param_file_path = os.path.join(self.uwgParamDir,self.uwgParamFileName)
        climate_file_path = os.path.join(self.epwDir, self.epwFileName)

        if not os.path.exists(uwg_param_file_path):
            raise Exception("Param file: '{}' does not exist.".format(uwg_param_file))

        # Open .uwg file and feed csv data to initializeDataFile
        try:
            uwg_param_data = utilities.read_csv(uwg_param_file_path)
        except Exception as e:
            raise Exception("Failed to read .uwg file! {}".format(e.message))

        # The initialize.uwg is read with a dictionary so that users changing
        # line endings or line numbers doesn't make reading input incorrect
        # It may make sense to change .uwg into json or something for more control over i/o
        self._init_param_dict = {}
        count = 0
        while  count < len(uwg_param_data):
            row = uwg_param_data[count]
            if row == [] or "#" in row[0]:
                count += 1
                continue
            elif row[0] == "SchTraffic":
                # SchTraffic: 3 x 24 matrix
                trafficrows = uwg_param_data[count+1:count+4]
                self._init_param_dict[row[0]] = map(lambda r: utilities.str2fl(r[:24]),trafficrows)
                count += 4
            elif row[0] == "bld":
                #bld: 17 x 3 matrix
                bldrows = uwg_param_data[count+1:count+17]
                self._init_param_dict[row[0]] = map(lambda r: utilities.str2fl(r[:3]),bldrows)
                count += 17
            else:
                count += 1
                self._init_param_dict[row[0]] = float(row[1])

        #rename initial parameters for local scope
        ipd = self._init_param_dict
        # Define Simulation and Weather parameters
        Month = ipd['Month']                # starting month (1-12)
        Day = ipd['Day']                    # starting day (1-31)
        nDay = ipd['nDay']                  # number of days
        dtSim = ipd['dtSim']                # simulation time step (s)

        dtWeather = ipd['dtWeather']        # seconds (s)
        autosize = ipd['autosize']          # autosize HVAC (1 or 0)
        sensOcc = ipd['sensOcc']            # Sensible heat from occupant
        LatFOcc = ipd['LatFOcc']            # Latent heat fraction from occupant (normally 0.3)
        RadFOcc = ipd['RadFOcc']            # Radiant heat fraction from occupant (normally 0.2)
        RadFEquip = ipd['RadFEquip']        # Radiant heat fraction from equipment (normally 0.5)
        RadFLight = ipd['RadFLight']        # Radiant heat fraction from light (normally 0.7)

        self.simTime = SimParam(dtSim,dtWeather,Month,Day,nDay)  # simulation time parametrs
        self.weather = Weather(climate_file_path,self.simTime.timeInitial,
        self.simTime.timeFinal) # weather file data for simulation time period
        self.forcIP = Forcing(self.weather.staTemp,self.weather) # initialized Forcing class
        self.forc = Forcing() # empty forcing class

        # Define Urban microclimate parameters
        h_ubl1 = ipd['h_ubl1']              # ubl height - day (m)
        h_ubl2 = ipd['h_ubl2']              # ubl height - night (m)
        h_ref = ipd['h_ref']                # inversion height
        h_temp = ipd['h_temp']              # temperature height
        h_wind = ipd['h_wind']              # wind height
        c_circ = ipd['c_circ']              # circulation coefficient
        c_exch = ipd['c_exch']              # exchange coefficient
        maxDay = ipd['maxDay']              # max day threshhold
        maxNight = ipd['maxNight']          # max night threshhold
        windMin = ipd['windMin']            # min wind speed (m/s)
        h_obs = ipd['h_obs']                # rural average obstacle height

        # Urban characteristics
        bldHeight = ipd['bldHeight']        # average building height (m)
        h_mix = ipd['h_mix']                 # mixing height (m)
        bldDensity = ipd['bldDensity']      # building density (0-1)
        verToHor = ipd['verToHor']          # building aspect ratio
        charLength = ipd['charLength']      # radius defining the urban area of study [aka. characteristic length] (m)
        alb_road = ipd['albRoad']           # road albedo
        d_road = ipd['dRoad']               # road pavement thickness
        sensAnth = ipd['sensAnth']          # non-building sensible heat (W/m^2)
        latAnth = ipd['latAnth']            # non-building latent heat heat (W/m^2)

        # climate Zone
        zone = int(ipd['zone'])-1

        # Vegetation parameters
        vegCover = ipd['vegCover']          # urban area veg coverage ratio
        treeCoverage = ipd['treeCoverage']  # urban area tree coverage ratio
        vegStart = ipd['vegStart']          # vegetation start month
        vegEnd = ipd['vegEnd']              # vegetation end month
        albVeg = ipd['albVeg']              # Vegetation albedo
        latGrss = ipd['latGrss']            # latent fraction of grass
        latTree = ipd['latTree']            # latent fraction of tree
        rurVegCover = ipd['rurVegCover']    # rural vegetation cover

        # Initialize geographic Param and Urban Boundary Layer Objects
        nightStart = 18.        # arbitrary values for begin/end hour for night setpoint
        nightEnd = 8.
        maxdx = 250;            # max dx (m)

        self.geoParam = Param(h_ubl1,h_ubl2,h_ref,h_temp,h_wind,c_circ,maxDay,maxNight,latTree,latGrss,albVeg,vegStart,vegEnd,\
            nightStart,nightEnd,windMin,self.wgmax,c_exch,maxdx,self.g,self.cp,self.vk,self.r,self.rv,self.lv,math.pi,\
            self.sigma,self.waterDens,self.lvtt,self.tt,self.estt,self.cl,self.cpv,self.b, self.cm,self.colburn)

        #TODO:  write UBLDef
        #UBL = UBLDef('C',charLength,weather.staTemp(1),maxdx,geoParam.dayBLHeight,geoParam.nightBLHeight);

        # Define Traffic schedule
        SchTraffic = ipd['SchTraffic']

        # Define Road (Assume 0.5m of asphalt)
        kRoad = ipd['kRoad']                 # road pavement conductivity (W/m K)
        cRoad = ipd['cRoad']                 # road volumetric heat capacity (J/m^3 K)

        emis = 0.93
        asphalt = Material(kRoad,cRoad,'asphalt')
        road_T_init = 293.
        road_horizontal = 1
        road_veg_coverage = min(vegCover/(1-bldDensity),1.) # fraction of surface vegetation coverage

        # define road layers
        road_layer_num = int(math.ceil(d_road/0.05))
        thickness_vector = map(lambda r: 0.05, range(road_layer_num)) # 0.5/0.05 ~ 10 x 1 matrix of 0.05 thickness
        material_vector = map(lambda n: asphalt, range(road_layer_num))

        self.road = Element(alb_road,emis,thickness_vector,material_vector,road_veg_coverage,\
            road_T_init,road_horizontal,name="urban_road")

        self.rural = copy.deepcopy(self.road)
        self.rural.vegCoverage = rurVegCover
        self.rural._name = "rural_road"


        # Define BEM for each DOE type (read the fraction)
        readDOE_file_path = os.path.join(self.DIR_UP_PATH,"resources","readDOE.pkl")
        if not os.path.exists(readDOE_file_path):
            raise Exception("readDOE.pkl file: '{}' does not exist.".format(readDOE_file_path))

        readDOE_file = open(readDOE_file_path, 'rb') # open pickle file in binary form
        refDOE = cPickle.load(readDOE_file)
        refBEM = cPickle.load(readDOE_file)
        refSchedule = cPickle.load(readDOE_file)
        readDOE_file.close()

        #TODO: Include optional parameters from intialize.uwg here after testing
        bld = ipd['bld']                    # fraction of building type/era
        albRoof = ipd['albRoof']            # roof albedo (0 - 1)
        vegRoof = ipd['vegRoof']            # Fraction of the roofs covered in grass/shrubs (0-1)
        glzR = ipd['glzR']               # Glazing Ratio. If not provided, all buildings are assumed to have 40% glazing ratio
        hvac = ipd['hvac']                  # HVAC TYPE; 0 = Fully Conditioned (21C-24C); 1 = Mixed Mode Natural Ventilation (19C-29C + windows open >22C); 2 = Unconditioned (windows open >22C)

        # Define building energy models
        k = 0
        r_glaze = 0             # Glazing ratio for total building stock
        SHGC = 0                # SHGC addition for total building stock
        alb_wall = 0            # albedo wall addition for total building stock
        h_floor = 3.05          # average floor height

        total_urban_bld_area = math.pow(charLength,2)*bldDensity*bldHeight/h_floor  # total building floor area
        area_matrix = utilities.zeros(16,3)

        self.BEM = []           # list of BEMDef objects
        self.Sch = []           # list of Schedule objects

        for i in xrange(16):    # 16 building types
            for j in xrange(3): # 3 built eras
                if bld[i][j] > 0.:
                    # Add to BEM list
                    self.BEM.append(refBEM[i][j][zone])
                    self.BEM[k].frac = bld[i][j]
                    self.BEM[k].fl_area = bld[i][j] * total_urban_bld_area

                    # Keep track of total urban r_glaze, SHGC, and alb_wall for UCM model
                    r_glaze = r_glaze + self.BEM[k].frac * self.BEM[k].building.glazingRatio
                    SHGC = SHGC + self.BEM[k].frac * self.BEM[k].building.shgc
                    alb_wall = alb_wall + self.BEM[k].frac * self.BEM[k].wall.albedo;

                    # Add to schedule list
                    self.Sch.append(refSchedule[i][j][zone])
                    k += 1

        #TODO: Finish RSM Class
        # Reference site class (also include VDM)
        self.RSM = RSMDef(self.lat,self.lon,self.GMT,h_obs,self.weather.staTemp[1],self.weather.staPres[1],self.geoParam)
        self.USM = RSMDef(self.lat,self.lon,self.GMT,bldHeight/10.,self.weather.staTemp[1],self.weather.staPres[1],self.geoParam)

        T_init = self.weather.staTemp[1]
        H_init = self.weather.staHum[1]

        self.UCM = UCMDef(bldHeight,bldDensity,verToHor,treeCoverage,sensAnth,latAnth,T_init,H_init,\
        self.weather.staUmod[1],self.geoParam,r_glaze,SHGC,alb_wall,self.road)
        self.UCM.h_mix = h_mix

        # Define Road Element & buffer to match ground temperature depth
        roadMat, newthickness = procMat(self.road,self.maxThickness,self.minThickness)

        for i in xrange(self.nSoil):
            # if soil depth is greater then the thickness of the road
            # we add new slices of soil at max thickness until road is greater or equal
            if self.depth_soil[i][0] >= sum(newthickness):
                #to avoid floating point precision problems compare equality
                is_greater = self.depth_soil[i][0] > sum(newthickness)
                is_equal = self.is_near_zero(self.depth_soil[i][0] - sum(newthickness),1e-10)
                while(not is_equal and is_greater):
                    newthickness.append(self.maxThickness)
                    roadMat.append(self.soil)
                self.soilindex1 = i
                break

        self.road = Element(self.road.albedo, self.road.emissivity, newthickness, roadMat,\
            self.road.vegCoverage, self.road.layerTemp[0], self.road.horizontal, self.road._name)


        # Define Rural Element
        ruralMat, newthickness = procMat(self.rural,self.maxThickness,self.minThickness)
        for i in xrange(self.nSoil):
            # if soil depth is greater then the thickness of the road
            # we add new slices of soil at max thickness until road is greater or equal
            if self.depth_soil[i][0] >= sum(newthickness):
                #to avoid floating point precision problems compare equality
                is_greater = self.depth_soil[i][0] > sum(newthickness)
                is_equal = self.is_near_zero(self.depth_soil[i][0] - sum(newthickness),1e-10)
                while(not is_equal and is_greater):
                    newthickness.append(self.maxThickness)
                    ruralMat.append(self.soil)
                self.soilindex2 = i
                break

        self.rural = Element(self.rural.albedo, self.rural.emissivity, newthickness,\
            ruralMat,self.rural.vegCoverage,self.rural.layerTemp[0],self.rural.horizontal, self.rural._name)


    def hvac_autosize(self):
        """ Section 6 - HVAC Autosizing (unlimited cooling & heating) """

        for i in xrange(len(self.BEM)):
            self.BEM[i].building.coolCap = 9999.
            self.BEM[i].building.heatCap = 9999.

    def uwg_main(self):
        """ Section 7 - UWG main section

            self.N         #
            self.ph        # per hour
        """


        self.N = int(self.simTime.days * 24)       # total number of hours in simulation
        n = 0
        self.ph = self.simTime.dt/3600.       # dt (simulation time step) in hours

        # Data dump variables
        time = range(self.N)
        #TODO: Figure out what this means
        #WeatherData(self.N, 1) = Forcing() #Empty object
        #UCMData(N,1) = UCMDef
        #UBLData (N,1) = UBLDef;
        #RSMData (N,1) = RSMDef;
        #USMData (N,1) = RSMDef;


        bTemp = utilities.zeros(self.N,len(self.BEM))
        bRHum = utilities.zeros(self.N,len(self.BEM))
        bPelec = utilities.zeros(self.N,len(self.BEM))
        bQgas = utilities.zeros(self.N,len(self.BEM))
        bPequip = utilities.zeros(self.N,len(self.BEM))
        bPlight = utilities.zeros(self.N,len(self.BEM))
        bQocc = utilities.zeros(self.N,len(self.BEM))
        bFluxMass = utilities.zeros(self.N,len(self.BEM))
        bFluxRoof = utilities.zeros(self.N,len(self.BEM))
        bFluxWall = utilities.zeros(self.N,len(self.BEM))
        bFluxSolar = utilities.zeros(self.N,len(self.BEM))
        bFluxWindow = utilities.zeros(self.N,len(self.BEM))
        bFluxInfil = utilities.zeros(self.N,len(self.BEM))
        bFluxVent = utilities.zeros(self.N,len(self.BEM))
        bCoolConsump = utilities.zeros(self.N,len(self.BEM))
        bHeatConsump = utilities.zeros(self.N,len(self.BEM))
        bCoolDemand = utilities.zeros(self.N,len(self.BEM))
        bHeatDemand = utilities.zeros(self.N,len(self.BEM))
        bTwallext = utilities.zeros(self.N,len(self.BEM))
        bTroofext = utilities.zeros(self.N,len(self.BEM))
        bTwallin = utilities.zeros(self.N,len(self.BEM))
        bTroofin = utilities.zeros(self.N,len(self.BEM))
        bTmassin = utilities.zeros(self.N,len(self.BEM))
        bCOP = utilities.zeros(self.N,len(self.BEM))
        bVent = utilities.zeros(self.N,len(self.BEM))


        #TODO: keep at :15 slice for now
        print 31 * 24 * 3600/300., self.simTime.nt-1 # =8928, simulation time (every 5 minutes)
        print 31 * 24, (self.simTime.nt-1)/12., len(self.forcIP.infra) #=744, epw weather file (every 60 minutes)

        #f = open(os.path.join(DIR_UP_PATH ,"check_simtime.txt"),'w')

        for it in range(1,self.simTime.nt,1)[:]: #for every simulation time-step (i.e 5 min) defined by uwg
            # Update water temperature (estimated)
            if self.is_near_zero(self.nSoil):
                self.forc.deepTemp = sum(self.forcIP.temp)/float(len(self.forcIP.temp))             # for BUBBLE/CAPITOUL/Singapore only
                self.forc.waterTemp = sum(self.forcIP.temp)/float(len(self.forcIP.temp)) - 10.      # for BUBBLE/CAPITOUL/Singapore only
            else:
                self.forc.deepTemp = self.Tsoil[self.soilindex1][self.simTime.month] #soil temperature by depth, by month
                self.forc.waterTemp = self.Tsoil[2][self.simTime.month]

            # There's probably a better way to update the weather...
            self.simTime.UpdateDate()
            # Update forcing parameters, by simulation timestep * per hour
            ceil_time_step = int(math.ceil(it * self.ph)) - 1  # simulation time increment raised to weather time step

            #print 'it', it
            #print 'c-it*ph', int(math.ceil(time_increment_in_hours))
            #print 'it*ph', time_increment_in_hours
            #print '---'
            #fchk = "it: {a}\ncitph: {b}\nitph: {c}\n----\n".format(
            #    a=it,
            #    b= ceil_time_step,
            #    c= it * self.ph
            #    )

            self.forc.infra = self.forcIP.infra[ceil_time_step]        # horizontal Infrared Radiation Intensity (W m-2)
            self.forc.wind = max(self.forcIP.wind[ceil_time_step], self.geoParam.windMin) # wind speed (m s-1)
            self.forc.uDir = self.forcIP.uDir[ceil_time_step]          # wind direction
            self.forc.hum = self.forcIP.hum[ceil_time_step]            # specific humidty (kg kg-1)
            self.forc.pres = self.forcIP.pres[ceil_time_step]          # Pressure (Pa)
            self.forc.temp = self.forcIP.temp[ceil_time_step]          # air temperature (C)
            self.forc.rHum = self.forcIP.rHum[ceil_time_step]          # Relative humidity (%)
            self.forc.prec = self.forcIP.prec[ceil_time_step]          # Precipitation (mm h-1)
            self.forc.dir = self.forcIP.dir[ceil_time_step]            # normal solar direct radiation (W m-2)
            self.forc.dif = self.forcIP.dif[ceil_time_step]            # horizontal solar diffuse radiation (W m-2)
            self.UCM.canHum = self.forc.hum                            # Canyon humidity (absolute) same as rural

        #f.close()
        """

            % Update solar flux
            [rural,UCM,BEM] = SolarCalcs(UCM,BEM,simTime,RSM,forc,geoParam,rural);

            % Update buildling & traffic schedule
            if strcmp(ext,'.xlsm') || strcmp(ext,'.m')

                % Assign day type (1 = weekday, 2 = sat, 3 = sun/other)
                if mod (simTime.julian,7) == 0      % Sunday
                    dayType = 3;
                elseif mod (simTime.julian,7) == 6  % Saturday
                    dayType = 2;
                else                                % Weekday
                    dayType = 1;
                end

                % Update anthropogenic heat load for each hour (building & UCM)
                UCM.sensAnthrop = sensAnth*(SchTraffic(dayType,simTime.hourDay+1));

                for i = 1:numel(BEM)

                    % Set temperature
                    BEM(i).building.coolSetpointDay = Sch(i).Cool(dayType,simTime.hourDay+1) + 273.15;
                    BEM(i).building.coolSetpointNight = BEM(i).building.coolSetpointDay;
                    BEM(i).building.heatSetpointDay = Sch(i).Heat(dayType,simTime.hourDay+1) + 273.15;
                    BEM(i).building.heatSetpointNight = BEM(i).building.heatSetpointDay;

                    % Internal Heat Load Schedule (W/m^2 of floor area for Q)
                    BEM(i).Elec = Sch(i).Qelec*Sch(i).Elec(dayType,simTime.hourDay+1);
                    BEM(i).Light = Sch(i).Qlight*Sch(i).Light(dayType,simTime.hourDay+1);
                    BEM(i).Nocc = Sch(i).Nocc*Sch(i).Occ(dayType,simTime.hourDay+1);
                    BEM(i).Qocc = sensOcc*(1-LatFOcc)*BEM(i).Nocc;

                    % SWH and ventilation schedule
                    BEM(i).SWH = Sch(i).Vswh*Sch(i).SWH(dayType,simTime.hourDay+1);     % litres per hour / m^2 of floor space
                    BEM(i).building.vent = Sch(i).Vent;                                 % m^3/s/m^2 of floor
                    BEM(i).Gas = Sch(i).Qgas * Sch(i).Gas(dayType,simTime.hourDay+1);   % Gas Equip Schedule, per m^2 of floor

                    % This is quite messy, should update
                    intHeat = BEM(i).Light+BEM(i).Elec+BEM(i).Qocc;
                    BEM(i).building.intHeatDay = intHeat;
                    BEM(i).building.intHeatNight = intHeat;
                    BEM(i).building.intHeatFRad = (RadFLight *BEM(i).Light + RadFEquip*BEM(i).Elec)/intHeat;
                    BEM(i).building.intHeatFLat = LatFOcc*sensOcc*BEM(i).Nocc/intHeat;

                    BEM(i).T_wallex = BEM(i).wall.layerTemp(1);
                    BEM(i).T_wallin = BEM(i).wall.layerTemp(end);
                    BEM(i).T_roofex = BEM(i).roof.layerTemp(1);
                    BEM(i).T_roofin = BEM(i).roof.layerTemp(end);
                end
            end

            % Update rural heat fluxes & update vertical diffusion model (VDM)
            rural.infra = forc.infra-rural.emissivity*sigma*rural.layerTemp(1)^4.;
            rural = SurfFlux(rural,forc,geoParam,simTime,forc.hum,forc.temp,forc.wind,2,0.);
            RSM = VDM(RSM,forc,rural,geoParam,simTime);

            % Calculate urban heat fluxes, update UCM & UBL
            [UCM,UBL,BEM] = UrbFlux(UCM,UBL,BEM,forc,geoParam,simTime,RSM);
            UCM = UCModel(UCM,BEM,UBL.ublTemp,forc,geoParam);
            UBL = UBLModel(UBL,UCM,RSM,rural,forc,geoParam,simTime);

            % Experimental code to run diffusion model in the urban area
            Uroad = UCM.road;
            Uroad.sens = UCM.sensHeat;
            Uforc = forc;
            Uforc.wind = UCM.canWind;
            Uforc.temp = UCM.canTemp;
            USM = VDM(USM,Uforc,Uroad,geoParam,simTime);

            % Update variables to output data dump
            if mod(simTime.secDay,simTime.timePrint) == 0 && n < N
                n = n + 1;
                WeatherData (n) = forc;
                [~,~,UCM.canRHum,~,UCM.Tdp,~] = Psychrometrics (UCM.canTemp, UCM.canHum, forc.pres);
                UBLData (n) = UBL;
                UCMData (n) = UCM;
                USMData (n) = USM;
                RSMData (n) = RSM;

                for i = 1:numel(BEM)
                    bTemp(n,i) = BEM(i).building.indoorTemp;
                    bVent(n,i) = BEM(i).building.vent;
                    bRHum(n,i) = BEM(i).building.indoorRhum;
                    bPelec(n,i) = BEM(i).building.ElecTotal;    % HVAC + Lighting + Elec Equip
                    bQgas(n,i) = BEM(i).building.GasTotal;
                    bPequip(n,i) = BEM(i).Elec;                 % Electric equipment only
                    bPlight(n,i) = BEM(i).Light;
                    bQocc(n,i) = BEM(i).Qocc;
                    bFluxMass(n,i) = -BEM(i).building.fluxMass*2;    % Assume floor & ceiling
                    bFluxWall(n,i) = -BEM(i).building.fluxWall*UCM.verToHor/UCM.bldDensity/BEM(i).building.nFloor;
                    bFluxRoof(n,i) = -BEM(i).building.fluxRoof/BEM(i).building.nFloor;
                    bFluxSolar(n,i) = BEM(i).building.fluxSolar;
                    bFluxWindow(n,i) = BEM(i).building.fluxWindow;
                    bFluxInfil(n,i) = BEM(i).building.fluxInfil;
                    bFluxVent(n,i) = BEM(i).building.fluxVent;
                    bCoolConsump(n,i) = BEM(i).building.coolConsump;
                    bHeatConsump(n,i) = BEM(i).building.sensHeatDemand/BEM(i).building.heatEff;
                    bCoolDemand(n,i) = BEM(i).building.sensCoolDemand;
                    bHeatDemand(n,i) = BEM(i).building.sensHeatDemand;
                    bTwallext(n,i) = BEM(i).T_wallex;
                    bTroofext(n,i) = BEM(i).T_roofex;
                    bTwallin(n,i) = BEM(i).T_wallin;
                    bTroofin(n,i) = BEM(i).T_roofin;
                    bTmassin(n,i) = BEM(i).mass.layerTemp(1);
                    bCOP(n,i) = BEM(i).building.copAdj;
                end
                progressbar(it/simTime.nt); % Print progress
            end

        end
        progressbar(1); % Close progress bar

    # =========================================================================
    # Section 8 - Writing new EPW file
    # =========================================================================
    if strcmp('Yes',writeEPW)
        disp('Calculating new Temperature and humidity values')
        for iJ = 1:numel(UCMData)
            epwinput.values{iJ+simTime.timeInitial-8,7}{1,1} = num2str(UCMData(iJ).canTemp- 273.15,'%0.1f'); % dry bulb temperature  [C]
            epwinput.values{iJ+simTime.timeInitial-8,8}{1,1} = num2str(UCMData(iJ).Tdp,'%0.1f'); % dew point temperature [C]
            epwinput.values{iJ+simTime.timeInitial-8,9}{1,1} = num2str(UCMData(iJ).canRHum,'%0.0f'); % relative humidity     [%]
            epwinput.values{iJ+simTime.timeInitial-8,22}{1,1} = num2str(WeatherData(iJ).wind,'%0.1f'); % wind speed [m/s]
        end
        disp('writing new EPW file');

        % Writing new EPW file
        new_climate_file = strcat(newPathName,'\',newFileName,'.epw');
        epwnewid = fopen(new_climate_file,'w');

        for i = 1:8
            fprintf(epwnewid,'%s\r\n',header{i});
        end

        for i = 1:size(epwinput.values,1)
            printme = [];
            for e = 1:34
                printme = [printme epwinput.values{i,e}{1,1} ','];
            end
            printme = [printme epwinput.values{i,e}{1,1}];
            fprintf(epwnewid,'%s\r\n',printme);
        end
        disp(['New climate file generated: ',new_climate_file]);
    end

    return None
    """

def procMat(materials,max_thickness,min_thickness):
    """ Processes material layer so that a material with single
    layer thickness is divided into two and material layer that is too
    thick is subdivided
    """
    newmat = []
    newthickness = []
    k = materials.layerThermalCond
    Vhc = materials.layerVolHeat

    if len(materials.layerThickness) > 1:

        for j in xrange(len(materials.layerThickness)):
            # Break up each layer that's more than max thickness (0.05m)
            if materials.layerThickness[j] > max_thickness:
                nlayers = math.ceil(materials.layerThickness[j]/float(max_thickness))
                for i in xrange(int(nlayers)):
                    newmat.append(Material(k[j],Vhc[j]))
                    newthickness.append(materials.layerThickness[j]/float(nlayers))
            # Material that's less then min_thickness is added at min_thickness.
            elif materials.layerThickness[j] < min_thickness:
                newmat.append(Material(k[j],Vhc[j]))
                newthickness.append(min_thickness)
                print 'WARNING: Material layer found too thin (<{:.2f}cm), added at new minimum thickness'.format(min_thickness*100)
            else:
                newmat.append(Material(k[j],Vhc[j]))
                newthickness.append(materials.layerThickness[j])

    else:

        # Divide single layer into two (UWG assumes at least 2 layers)
        if materials.layerThickness[0] > max_thickness:
            nlayers = math.ceil(materials.layerThickness[0]/float(max_thickness))
            for i in xrange(int(nlayers)):
                newmat.append(Material(k[0],Vhc[0]))
                newthickness.append(materials.layerThickness[0]/float(nlayers))
        # Material should be at least 1cm thick, so if we're here,
        # should give warning and stop. Only warning given for now.
        elif materials.layerThickness[0] < min_thickness*2:
            newthickness = [min_thickness/2., min_thickness/2.]
            newmat = [Material(k[0],Vhc[0]), Material(k[0],Vhc[0])]
            print 'WARNING: a thin (<2cm) single layer element found. May cause error'
        else:
            newthickness = [materials.layerThickness[0]/2., materials.layerThickness[0]/2.]
            newmat = [Material(k[0],Vhc[0]), Material(k[0],Vhc[0])]
    return newmat, newthickness


if __name__ == "__main__":

    DIR_UP_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    epwDir = os.path.join(DIR_UP_PATH,"resources","epw")
    epwFileName = "SGP_Singapore.486980_IWEC.epw"
    uwgParamDir = os.path.join(DIR_UP_PATH,"resources")
    uwgParamFileName = "initialize.uwg"
    uwg = UWG(epwDir, epwFileName, uwgParamDir, uwgParamFileName)
    uwg.read_epw()
    uwg.read_input()
    uwg.hvac_autosize()
    uwg.uwg_main()
