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

import utilities
from material import Material
from simparam import SimParam



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

    # =========================================================================
    # Section 1 - Definitions for constants / other parameters
    # =========================================================================
    #TODO: capitalize for constant covnention
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

    def __init__(self, epwDir, epwFileName, uwgParamDir, uwgParamFileName, destinationDir=None, destinationFile=None):
        self.epwDir = epwDir
        self.epwFileName = epwFileName
        self.uwgParamDir = uwgParamDir
        self.uwgParamFileName = uwgParamFileName
        self.destinationDir = destinationDir
        self.destinationFile = destinationFile

    def __repr__(self):
        return "UWG: {} ".format(self.epwFileName)

    def is_near_zero(self,num,eps=1e-10):
        return abs(float(num)) < eps

    def read_epw(self):
        """Section 2 - Read EPW file
        properties:
            self.newPathName
            self.header     # header data
            self.epwinput   # timestep data for weather
            self.lat        # latitude
            self.lon        # longitude
            self.GMT        # GMT
            self.nSoil      # Number of soil depths
            self.Tsoil      # nSoil x 12 matrix for soil temperture (K)
            self.depth      # nSoil x 1 matrix for soil depth (m)
        """

        # Revise epw file name if not end with epw
        if not self.epwFileName.lower().endswith('.epw'):
            self.epwFileName = self.epwFileName + '.epw'

        # Make dir path to epw file
        climateDataPath = os.path.join(self.epwDir, self.epwFileName)

        # Open epw file and feed csv data to climateDataFile
        try:
            climateDataFile = utilities.read_csv(climateDataPath)
        except OSError:
            # TODO: Figure out why this isn't called when OSError is raised..
            print "Can not find " + climateDataPath

        # Read header lines (1 to 8) from EPW and ensure TMY2 format.
        self.header = climateDataFile[0:8]

        # Read weather data from EPW for each time step in weather file. (lines 8 - end)
        self.epwinput = climateDataFile[8:]

        # Read Lat, Long (line 1 of EPW)
        self.lat = float(self.header[0][6])
        self.lon = float(self.header[0][7])
        self.GMT = float(self.header[0][8])

        # Read in soil temperature data (assumes this is always there)
        # ref: http://bigladdersoftware.com/epx/docs/8-2/auxiliary-programs/epw-csv-format-inout.html
        soilData = self.header[3]
        self.nSoil = int(soilData[1])           # Number of ground temperature depths
        self.Tsoil = utilities.zeros(self.nSoil,12)  # nSoil x 12 matrix for soil temperture (K)
        self.depth = utilities.zeros(self.nSoil,1)   # nSoil x 1 matrix for soil depth (m)

        # Read monthly data for each layer of soil from EPW file
        for i in xrange(self.nSoil):
            self.depth[i][0] = float(soilData[2 + (i*16)])
            # Monthly data
            for j in xrange(12):
                self.Tsoil[i][j] = float(soilData[6 + (i*16) + j]) + 273.15

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

        """

        #TODO: Need to load data from initialize.uwg
        #TODO: Possible take uwg ext, change to py and unpickle/pickle?

        # Run script to generate UCM, UBL, etc.
        nightStart = 18.        # arbitrary values for begin/end hour for night setpoint
        nightEnd = 8.

        """
        simTime = SimParam(dtSim,dtWeather,Month,Day,nDay);
        weather = Weather(climate_data,simTime.timeInitial,simTime.timeFinal);
        forcIP = Forcing(weather.staTemp,weather);
        forc = Forcing;

        % Road (Assume 0.5m of asphalt)
        emis = 0.93;
        asphalt = Material (1.0,1.6e6);
        thickness = 0.05 * ones (ceil(d_road/0.05),1);
        road = Element(alb_road,emis,thickness,[asphalt;asphalt;asphalt;asphalt;asphalt;asphalt;asphalt;asphalt;asphalt;asphalt;],0,293,1);
        road.vegCoverage = min(vegCover/(1-bldDensity),1);

        rural = road;
        rural.vegCoverage = rurVegCover;
        T_init = weather.staTemp(1);
        H_init = weather.staHum(1);

        geoParam = Param(h_ubl1,h_ubl2,h_ref,h_temp,h_wind,c_circ,maxDay,maxNight,...
            latTree,latGrss,albVeg,vegStart,vegEnd,nightStart,nightEnd,windMin,wgmax,c_exch,maxdx,...
            g, cp, vk, r, rv, lv, pi(), sigma, waterDens, lvtt, tt, estt, cl, cpv, b, cm, colburn);
        UBL = UBLDef('C',charLength,weather.staTemp(1),maxdx,geoParam.dayBLHeight,geoParam.nightBLHeight);

        % Define BEM for each DOE type (read the fraction)
        load ('RefDOE.mat');

        % Define building energy models
        k = 0;
        r_glaze = 0;
        SHGC = 0;
        alb_wall = 0;
        area = bld*charLength^2*bldDensity*bldHeight/h_floor;  % building floor area
        for i = 1:16
            for j = 1:3
                if bld(i,j) > 0
                    k = k + 1;
                    BEM(k) = refBEM(i,j,zone);
                    BEM(k).frac = bld(i,j);
                    BEM(k).fl_area = area(i,j);
                    r_glaze = r_glaze + BEM(k).frac * BEM(k).building.glazingRatio;
                    SHGC = SHGC + BEM(k).frac * BEM(k).building.shgc;
                    alb_wall = alb_wall + BEM(k).frac * BEM(k).wall.albedo;
                    BEM(k).Qocc = BEM(k).Qocc;
                    Sch(k) = Schedule(i,j,zone);
                end
            end
        end

        UCM = UCMDef(bldHeight,bldDensity,verToHor,treeCoverage,...
            sensAnth,latAnth,T_init,H_init,weather.staUmod(1),geoParam,r_glaze,SHGC,alb_wall,road);
        UCM.h_mix = h_mix;

        % Reference site class (also include VDM)
        RSM = RSMDef(lat,lon,GMT,h_obs,weather.staTemp(1),weather.staPres(1),geoParam);
        USM = RSMDef(lat,lon,GMT,bldHeight/10,weather.staTemp(1),weather.staPres(1),geoParam);

        # Note:
        # Looks to me like road layer is fixed at 0.5 (init) or user-specified depth and
        # we subtract ground height until that reacehs 0.5 lenght or user-specified length.
        # I.e not sure if soil depth at whcih deep T is taken is dependant on road depth.

        % For .m file, assume the soil depth is close to one of the ground
        % soil depth specified in EPW (0.5, 1.0, 2.0)
        for i = 1:n_soil
            if sum(road.layerThickness) <= depth(i)
                soilindex1 = i;
                break;
            end
        end

        % Same for rural road
        for i = 1:n_soil
            if sum(rural.layerThickness) <= depth(i)
                soilindex2 = i;
                break;
            end
        end
        """

    """
    % =========================================================================
    % Section 6 - HVAC Autosizing (unlimited cooling & heating)
    % =========================================================================
    for j = 1:numel(BEM)
        if autosize
            BEM(j).building.coolCap = 9999;
            BEM(j).building.heatCap = 9999;
        end
    end

    % =========================================================================
    % Section 7 - UWG main section
    % =========================================================================
    N = simTime.days * 24;
    n = 0;
    ph = simTime.dt/3600;       % per hour

    % Data dump variables
    time = transpose(1:1:simTime.days*24);
    WeatherData (N,1) = Forcing;
    UCMData (N,1) = UCMDef;
    UBLData (N,1) = UBLDef;
    RSMData (N,1) = RSMDef;
    USMData (N,1) = RSMDef;

    bTemp = zeros (N,numel(BEM));
    bRHum = zeros (N,numel(BEM));
    bPelec = zeros (N,numel(BEM));
    bQgas = zeros (N,numel(BEM));
    bPequip = zeros (N,numel(BEM));
    bPlight = zeros (N,numel(BEM));
    bQocc = zeros (N,numel(BEM));
    bFluxMass = zeros (N,numel(BEM));
    bFluxRoof = zeros(N,numel(BEM));
    bFluxWall = zeros (N,numel(BEM));
    bFluxSolar = zeros (N,numel(BEM));
    bFluxWindow = zeros (N,numel(BEM));
    bFluxInfil = zeros (N,numel(BEM));
    bFluxVent = zeros (N,numel(BEM));
    bCoolConsump = zeros (N,numel(BEM));
    bHeatConsump = zeros (N,numel(BEM));
    bCoolDemand = zeros (N,numel(BEM));
    bHeatDemand = zeros (N,numel(BEM));
    bTwallext = zeros (N,numel(BEM));
    bTroofext = zeros (N,numel(BEM));
    bTwallin = zeros (N,numel(BEM));
    bTroofin = zeros (N,numel(BEM));
    bTmassin = zeros (N,numel(BEM));
    bCOP = zeros(N,numel(BEM));
    bVent = zeros (N,numel(BEM));

    for it=1:(simTime.nt-1)

        % Update water temperature (estimated)
        if n_soil == 0
            forc.deepTemp = mean([forcIP.temp]);            % for BUBBLE/CAPITOUL/Singapore only
            forc.waterTemp = mean([forcIP.temp]) - 10;      % for BUBBLE/CAPITOUL/Singapore only
        else
            forc.deepTemp = Tsoil(soilindex1,simTime.month);
            forc.waterTemp = Tsoil(3,simTime.month);
        end

        % There's probably a better way to update the weather...
        simTime = UpdateDate(simTime);
        forc.infra = forcIP.infra(ceil(it*ph));
        forc.wind = max(forcIP.wind(ceil(it*ph)),geoParam.windMin);
        forc.uDir = forcIP.uDir(ceil(it*ph));
        forc.hum = forcIP.hum(ceil(it*ph));
        forc.pres = forcIP.pres(ceil(it*ph));
        forc.temp = forcIP.temp(ceil(it*ph));
        forc.rHum = forcIP.rHum(ceil(it*ph));
        forc.prec = forcIP.prec(ceil(it*ph));
        forc.dir = forcIP.dir(ceil(it*ph));
        forc.dif = forcIP.dif(ceil(it*ph));
        UCM.canHum = forc.hum;      % Canyon humidity (absolute) same as rural

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

        elseif strcmp(ext,'.xml')

            for i = 1:numel(BEM)

                % Schedules not used for .xml interface set to zero
                BEM(i).Elec = 0;
                BEM(i).Light = 0;
                BEM(i).Nocc = 0;
                BEM(i).Qocc = 0;
                BEM(i).SWH = 0;         % not used for .xml interface
                BEM(i).Gas = 0;         % not used for .xml interface

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
"""
if __name__ == "__main__":
    # Run the function.
    epwDir = None#'C:\ladybug'
    epwFileName = None#'USA_MA_Boston-Logan.Intl.AP.725090_TMY3.epw'
    uwgParamDir = None
    uwgParamFileName = None
    UWG(epwDir, epwFileName, uwgParamDir, uwgParamFileName)
"""
