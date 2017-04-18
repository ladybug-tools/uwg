function [new_climate_file] = UWG(CL_EPW_PATH,CL_EPW,CL_XML_PATH,CL_XML,CL_RE_PATH,CL_RE)
    % =========================================================================
    %  THE URBAN WEATHER GENERATOR
    % =========================================================================
    % Orig. Author: B. Bueno & edited by A. Nakano & Lingfu Zhang
    % Last modified by Joseph Yang (joeyang@mit.edu) - May, 2016
    % 
    % Notes
    % a. When compiling, add 'z_meso.txt', 'RefDOE.mat','SchDef.m' to the list of files
    % b. Program description can be found in the following papers & websites
    %   - Joseph Yang's Master Thesis (2016)
    %   - https://github.com/hansukyang/UWG_Matlab
    % =========================================================================

    close all;
    ver = 4.1;
    % 4.1 (beta) updates
    %   - changed EPW output to rural wind speed again
    %   - compatibility with xml re-established
    %   - read in 'initialize.m' file for Matlab 
    
    % =========================================================================
    % Section 1 - Definitions for constants / other parameters
    % =========================================================================

    min_thickness = 0.01;   % Minimum layer thickness (to prevent crashing) (m)
    max_thickness = 0.05;   % Maximum layer thickness (m)
    soilTcond = 1;          % http://web.mit.edu/parmstr/Public/NRCan/nrcc29118.pdf (Figly & Snodgrass)
    soilvolHeat = 2e6;      % http://www.europment.org/library/2013/venice/bypaper/MFHEEF/MFHEEF-21.pdf (average taken from Table 1)
    soil = Material(soilTcond,soilvolHeat); % Soil material used for soil-depth padding

    % Physical parameters (moved here from Param.m)
    g = 9.81;               % gravity
    cp = 1004.;             % heat capacity for air (J/kg.K)
    vk = 0.40;              % von karman constant
    r = 287.;               % gas constant
    rv = 461.5;             %
    lv = 2.26e6;            % latent heat of evaporation
    sigma = 5.67e-08 ;      % Stefan Boltzmann constant
    waterDens = 1000;       % water density (kg/m^3)
    lvtt = 2.5008e6;        %
    tt = 273.16;            %
    estt = 611.14;          %
    cl = 4.218e3;           %
    cpv = 1846.1;           %
    b = 9.4;                % Coefficients derived by Louis (1979)
    cm = 7.4;               %
    colburn = (0.713/0.621)^(2./3.); % (Pr/Sc)^(2/3) for Colburn analogy in water evaporation
    
    % Site-specific parameters    
    wgmax = 0.005;          % maximum film water depth on horizontal surfaces (m)

    % =========================================================================
    % Section 2 - Read EPW file
    % =========================================================================
    try
        climate_data = strcat(CL_EPW_PATH,'\',CL_EPW);
        fullyScripted = 1;
    catch
        [epwFileName,epwPathName] = uigetfile('.epw','Select Rural EnergyPlus Weather File');
        climate_data = strcat(epwPathName,epwFileName);
    end

    disp(['Rural weather file selected: ',climate_data])
    epwid = fopen(climate_data);
    C = importdata(climate_data, ',', 8);
    
    % Read header lines (1 to 8) from EPW and ensure TMY2 format
    % Note that TMY3 format is not compatible with the current version of UWG
    header = C.textdata(1:8,1);
    TMY3 = strfind(header(1), 'TMY3');
    if ~isempty(TMY3{1})
        disp('UWG unable to run: UWG requires TMY2 format for weather data');
        new_climate_file = 0;
        return;
    end

    % Read Lat, Long (line 1 of EPW)
    line1 = strsplit(header{1},',');
    lat = str2double(line1{7});
    lon = str2double(line1{8});
    GMT = str2double(line1{9}); 

    % Read in soil temperature data (assumes this is always there)
    soildata = strsplit(header{4},',');
    n_soil = str2double(soildata{2});
    Tsoil = zeros(n_soil,12);
    depth = zeros(n_soil,1);
    
    % Read monthly data for each layer of soil from EPW file
    for i = 1:n_soil
        depth(i) = str2double(soildata{3 + (i-1)*13});
        % Monthly data
        for j = 1:12
            Tsoil(i,j) = str2double(soildata{3+(i-1)*13+j})+273.15;
        end
    end

    % Read weather data from EPW for each time step in weather file
    i = 1;
    readin = fgetl(epwid);
    while (readin ~= -1)
        epwinput.values(i,:) = textscan(readin, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %*[^\n]', 'delimiter',',','MultipleDelimsAsOne',1);
        i=i+1;
        readin = fgetl(epwid);
    end
    epwinput.values(1:8,:) = []; % Ignore the header lines (first 8)
    fclose all;

    % Save location for the new EPW file
    try
        new_climate_file = strcat(CL_RE_PATH,'\',CL_RE);
        newPathName = CL_RE_PATH;
    catch
        [newFileName_withExt,newPathName] = uiputfile('.epw','Select save location');
        new_climate_file = strcat(newPathName,newFileName_withExt);
        newPathName = newPathName(1:end-1);
    end

    [~,newFileName,~] = fileparts(new_climate_file);
    disp(['Save location selected: ',new_climate_file]);
    
    % =========================================================================
    % Section 3 - Read Input File (xlsm, XML, .m, file)
    % =========================================================================

    try
        xml_location = strcat(CL_XML_PATH,'\',CL_XML);
    catch
        [FileName,PathName] = uigetfile('*.xml;*.m;*.xlsm','Select Urban Parameter Input file');
        xml_location = strcat(PathName,FileName);
    end
    
    % Input files for UWG - note that soil layer buffering and 
    % layer thickness control are only performed for XML. (should update) 
    [~,~,ext] = fileparts(xml_location);
    if strcmp(ext,'.xlsm')      % Excel input
        
        % Create simulation & weather class
        [num, ~, ~] = xlsread(xml_location,1,'N24:N34');
        Month = num(1);         % starting month (1-12)
        Day = num(2);           % starting day (1-31)
        nDay = num(3);          % number of days
        dtSim = num(4);         % simulation time step (s)
        
        dtWeather = num(5);     % seconds (s) 
        autosize = num(6);      % autosize HVAC (1 or 0)
        sensOcc = num(7);       % Sensible heat from occupant
        LatFOcc = num(8);       % Latent heat fraction from occupant (normally 0.3)
        RadFOcc = num(9);       % Radiant heat fraction from occupant (normally 0.2)
        RadFEquip = num(10);    % Radiant heat fraction from equipment (normally 0.5)
        RadFLight = num(11);    % Radiant heat fraction from light (normally 0.7)
        
        simTime = SimParam(dtSim,dtWeather,Month,Day,nDay);
        weather = Weather(climate_data,simTime.timeInitial,simTime.timeFinal);
        forcIP = Forcing(weather.staTemp,weather); 
        forc = Forcing;
        
        [~,txt,~] = xlsread(xml_location,1,'R32:R34');
        writeMAT = txt(1);
        writeEPW = txt(2);
        writeXLS = txt(3);

        % Urban microclimate parameters
        [num, ~, ~] = xlsread(xml_location,1,'D4:D32');
        h_ubl1 = num(1);        % ubl height - day (m)
        h_ubl2 = num(2);        % ubl height - night (m)
        h_ref = num(3);         % inversion height
        h_temp = num(4);        % temperature height
        h_wind = num(5);        % wind height
        c_circ = num(6);        % circulation coefficient
        c_exch = num(7);        % exchange coefficient
        maxDay = num(8);        % max day threshhold
        maxNight = num(9);      % max night threshhold
        windMin = num(10);      % min wind speed (m/s)
        h_obs = num(11);        % rural average obstacle height

        % Urban characteristics
        bldHeight = num(12);    % average building height (m)
        h_mix = num(13);        % mixing height (m)
        bldDensity = num(14);   % building density (0-1)
        verToHor = num(15);     % building aspect ratio
        charLength = num(16);   % characteristic length (m)
        maxdx = num(17);        % Max Dx (m)
        alb_road = num(18);     % road albedo
        d_road = num(19);       % road pavement thickness
        sensAnth = num(20);     % non-building sens heat (W/m^2)
        latAnth = num(21);      % non-building lat heat (W/m^2)
        
        % Vegetatin parameters
        vegCover = num(22);     % urban area veg coverage ratio
        treeCoverage = num(23); % urban area tree coverage ratio
        vegStart = num(24);     % vegetation start month
        vegEnd = num(25);       % vegetation end month
        albVeg = num(26);       % Vegetation albedo
        latGrss = num(27);      % latent fraction of grass
        latTree = num(28);      % latent fraction of tree
        rurVegCover = num(28);  % rural vegetation cover

        nightStart = 18;        % arbitrary values (not used for XLSM)
        nightEnd = 8;
        geoParam = Param(h_ubl1,h_ubl2,h_ref,h_temp,h_wind,c_circ,maxDay,maxNight,...
            latTree,latGrss,albVeg,vegStart,vegEnd,nightStart,nightEnd,windMin,wgmax,c_exch,maxdx,...
            g, cp, vk, r, rv, lv, pi(), sigma, waterDens, lvtt, tt, estt, cl, cpv, b, cm, colburn);
        UBL = UBLDef('C',charLength,weather.staTemp(1),maxdx,geoParam.dayBLHeight,geoParam.nightBLHeight); 

        % Traffic schedule
        [SchTraffic, ~, ~] = xlsread(xml_location,1,'H4:K27');
        SchTraffic = transpose(SchTraffic);

        % Define BEM for each DOE type (read the fraction)
        load ('RefDOE.mat');
 
        [zone, ~, ~] = xlsread(xml_location,1,'AA3');
        [num, ~, ~] = xlsread(xml_location,1,'S4:U19');
        [area, ~, ~ ] = xlsread(xml_location,1,'P4:R19');
        
        % Road (Assume 0.5m of asphalt)
        emis = 0.93;
        d_road = 0.5;
        asphalt = Material (1.0,1.6e6);
        thickness = 0.05 * ones (ceil(d_road/0.05),1);
        road = Element(alb_road,emis,thickness,[asphalt;asphalt;asphalt;asphalt;asphalt;asphalt;asphalt;asphalt;asphalt;asphalt;],0,293,1);
        road.vegCoverage = min(vegCover/(1-bldDensity),1);

        % Define building energy models
        k = 0;
        r_glaze = 0;
        SHGC = 0;
        alb_wall = 0;
        for i = 1:16
            for j = 1:3
                if num(i,j) > 0
                    k = k + 1;
                    BEM(k) = refBEM(i,j,zone);
                    BEM(k).frac = num(i,j);
                    BEM(k).fl_area = area(i,j);
                    r_glaze = r_glaze + BEM(k).frac * BEM(k).building.glazingRatio;
                    SHGC = SHGC + BEM(k).frac * BEM(k).building.shgc;
                    alb_wall = alb_wall + BEM(k).frac * BEM(k).wall.albedo;
                    BEM(k).Qocc = BEM(k).Qocc;
                    Sch(k) = Schedule(i,j,zone);      
                end
            end 
        end
        
        % Reference site class (also include VDM)
        RSM = RSMDef(lat,lon,GMT,h_obs,weather.staTemp(1),weather.staPres(1),geoParam);
        USM = RSMDef(lat,lon,GMT,bldHeight/10,weather.staTemp(1),weather.staPres(1),geoParam);

        % Create UCM class (use road characteristics from BEM)
        rural = road;
        rural.vegCoverage = rurVegCover;
        T_init = weather.staTemp(1);
        H_init = weather.staHum(1);
        UCM = UCMDef(bldHeight,bldDensity,verToHor,treeCoverage,...
            sensAnth,latAnth,T_init,H_init,weather.staUmod(1),geoParam,r_glaze,SHGC,alb_wall,road); 
        UCM.h_mix = h_mix;
        
        % Misc. stuff
        soilindex1 = 1;
        soilindex2 = 1;

    elseif strcmp(ext,'.m')
        % Run matlab script to generate UCM, UBL, etc.
        run(xml_location);
        nightStart = 18;        % arbitrary values (not used for XLSM)
        nightEnd = 8;

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

    elseif strcmp(ext,'.xml')
        % Some numbers not specified in XML
        maxdx = 500;            % maximum discretization length for the UBL model (m)
        circCoeff = 1.2;        %
        dayThreshold = 200;     %
        nightThreshold = 50;    %
        windMin = 1;            % minimum wind speed (m s-1)
        c_exch = 1;
        autosize = 1;

        % Process XML file to generate class elements
        xml_input = xml_read(xml_location);
        sim_dt = 300;           % Simulation time step (s)
        weather_dt = 3600;      % Weather data time step (EPW) (s)

        disp(['Urban Parameter file selected: ',xml_location]);

        % Re-naming for file readability
        xmlTyp(1) = xml_input.typology1;
        xmlTyp(2) = xml_input.typology2;
        xmlTyp(3) = xml_input.typology3;
        xmlTyp(4) = xml_input.typology4;
        building(1) = xmlTyp(1).building;
        building(2) = xmlTyp(2).building;
        building(3) = xmlTyp(3).building;
        building(4) = xmlTyp(4).building;
        xmlParam = xml_input.parameter;
        xmlUCM = xml_input.urbanArea;
        xmlRSite = xml_input.referenceSite;

        % Simulation paramters
        simTime = SimParam(sim_dt,weather_dt,xmlParam.simuStartMonth,...
            xmlParam.simuStartDay, xmlParam.simuDuration);
        weather = Weather(climate_data,simTime.timeInitial,simTime.timeFinal);
        
        forcIP = Forcing(weather.staTemp,weather); 
        forc = Forcing;
        
        nightStart = mean([building.nightSetStart]);
        nightEnd = mean([building.nightSetEnd]);
        geoParam = Param(xmlUCM.daytimeBLHeight,xmlUCM.nighttimeBLHeight,...
            xmlUCM.refHeight,xmlParam.tempHeight,xmlParam.windHeight,...
            circCoeff,dayThreshold,nightThreshold,...
            xmlUCM.treeLatent,xmlUCM.grassLatent,xmlUCM.vegAlbedo,xmlUCM.vegStart,xmlUCM.vegEnd,...
            nightStart,nightEnd,windMin,wgmax,c_exch,maxdx,...
            g, cp, vk, r, rv, lv, pi(), sigma, waterDens, lvtt, tt, estt, cl, cpv, b, cm, colburn);
        
        % Define Road Element & buffer to match ground temperature depth
        urbanRoad = xmlUCM.urbanRoad;
        [roadMat, newthickness] = procMat(urbanRoad.materials,max_thickness,min_thickness);
        for i = 1:n_soil
            if sum(newthickness) <= depth(i)
                while(sum(newthickness)<depth(i))
                    newthickness = [newthickness; max_thickness];
                    roadMat = [roadMat soil];
                end
                soilindex1 = i;
                break;
            end
        end
        road = Element(urbanRoad.albedo,urbanRoad.emissivity,newthickness,roadMat,...
            urbanRoad.vegetationCoverage,urbanRoad.initialTemperature + 273.15,urbanRoad.inclination);

        % Define Rural Element
        ruralRoad = xmlRSite.ruralRoad;
        [ruralMat, newthickness] = procMat(ruralRoad.materials,max_thickness,min_thickness);
        for i = 1:n_soil
            if sum(newthickness) <= depth(i)
                while(sum(newthickness)<depth(i))
                    newthickness = [newthickness; max_thickness];
                    ruralMat = [ruralMat soil];
                end
                soilindex2 = i;
                break;
            end
        end
        rural = Element(ruralRoad.albedo,ruralRoad.emissivity,newthickness,...
            ruralMat,ruralRoad.vegetationCoverage,ruralRoad.initialTemperature + 273.15,ruralRoad.inclination);

        frac_total = 0;
        for i = 1:4
            % Define Wall element
            [wallMat, newthickness] = procMat(xmlTyp(i).construction.wall.materials,max_thickness,min_thickness);
            xwall = xmlTyp(i).construction.wall;
            wall = Element(xwall.albedo,xwall.emissivity,newthickness,wallMat,...
                xwall.vegetationCoverage,xwall.initialTemperature + 273.15,xwall.inclination);

            % Define Roof element
            [roofMat, newthickness] = procMat(xmlTyp(i).construction.roof.materials,max_thickness,min_thickness);
            xroof = xmlTyp(i).construction.roof;
            roof = Element(xroof.albedo,xroof.emissivity,newthickness,roofMat,...
                xroof.vegetationCoverage,xroof.initialTemperature + 273.15,xroof.inclination);

            % Define Mass element
            [massMat, newthickness] = procMat(xmlTyp(i).construction.mass.materials,max_thickness,min_thickness);
            xmass = xmlTyp(i).construction.mass;
            mass = Element(xmass.albedo,xmass.emissivity,newthickness,massMat,...
                xmass.vegetationCoverage,xmass.initialTemperature + 273.15,xmass.inclination);

            % Define building typology
            typology = Building(building(i).floorHeight,...
                building(i).nightInternalGains,...
                building(i).dayInternalGains,...
                building(i).radiantFraction,...
                building(i).latentFraction,...
                building(i).infiltration,...
                building(i).ventilation,...
                xmlTyp(i).construction.glazing.glazingRatio,...
                xmlTyp(i).construction.glazing.windowUvalue,...
                xmlTyp(i).construction.glazing.windowSHGC,...
                building(i).coolingSystemType,...
                building(i).coolingCOP,...
                building(i).daytimeCoolingSetPoint + 273.15,...
                building(i).nighttimeCoolingSetPoint + 273.15,...
                building(i).daytimeHeatingSetPoint + 273.15,...
                building(i).nighttimeHeatingSetPoint + 273.15,...
                building(i).coolingCapacity,...
                building(i).heatingEfficiency,...
                building(i).initialT + 273.15);

            % Define Urban Configuration [building,mass,wall,roof]
            BEM(i) = BEMDef(typology,mass,wall,roof,xmlTyp(i).dist/100);
            
            % If only one typology is defined, break
            frac_total = frac_total + xmlTyp(i).dist/100;
            if frac_total >= 1
                break;
            end
        end 

        % Reference site class (also include VDM)
        RSM = RSMDef(lat,lon,GMT,xmlRSite.averageObstacleHeight,...
            weather.staTemp(1),weather.staPres(1),geoParam);


        % Define average typology for UCM (weighted by floor space %)
        SHGC = 0;
        alb_wall = 0;
        r_glaze = 0;
        for i = 1:numel(BEM)
            SHGC = SHGC + BEM(i).frac*BEM(i).building.shgc;
            alb_wall = alb_wall + BEM(i).frac*BEM(i).wall.albedo;
            r_glaze = r_glaze + BEM(i).frac*BEM(i).building.glazingRatio;
        end
        
        T_init = weather.staTemp(1);
        H_init = weather.staHum(1);

        % Define UCM class
        h_bld = xmlUCM.averageBuildingHeight;
        dens = xmlUCM.siteCoverageRatio;
        Qtraffic = xmlUCM.nonBldgLatentAnthropogenicHeat;
        Ltraffic = xmlUCM.nonBldgSensibleHeat;
        VtoH = xmlUCM.facadeToSiteRatio;
        UCM = UCMDef(h_bld,dens,VtoH,xmlUCM.treeCoverage,Ltraffic,Qtraffic,...
            T_init,H_init,weather.staUmod(1),geoParam,r_glaze,SHGC,alb_wall,road); 
        UCM.h_mix = 1;
        UBL = UBLDef('C',xmlUCM.charLength,weather.staTemp(1),geoParam.maxdx,geoParam.dayBLHeight,geoParam.nightBLHeight); 
        USM = RSMDef(lat,lon,GMT,h_bld/10,weather.staTemp(1),weather.staPres(1),geoParam);

    else
        return;
    end  

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

    % =========================================================================
    % Section 8 - Writing new EPW file
    % =========================================================================
    if strcmp('Yes',writeEPW)
        disp('Calculating new Temperature and humidity values')
        for iJ = 1:numel(UCMData)
            epwinput.values{iJ+simTime.timeInitial-8,7}{1,1} = num2str(UCMData(iJ).canTemp- 273.15,'%0.1f'); % dry bulb temperature  [°C]
            epwinput.values{iJ+simTime.timeInitial-8,8}{1,1} = num2str(UCMData(iJ).Tdp,'%0.1f'); % dew point temperature [°C]
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

    % =========================================================================
    % Section 9 - Clean up & write data to Excel/Mat file
    % =========================================================================

    if strcmp('Yes',writeMAT)
        save ('UWGdata.mat','RSMData','UCMData','UBLData','WeatherData','USMData','time');
    end
    
    if strcmp('Yes',writeXLS)
        
        % Open Excel Automation server
        file = 'UWGoutput.xlsx'; % This must be full path name
        Excel = actxserver('Excel.Application');
        Workbooks = Excel.Workbooks;
    
        % Open Excel file
        Workbook=Workbooks.Open(file);
        Excel.Visible=1;

        output = 'UWGoutput.xlsx';
        disp(['Writing output file: ',output]);

        T_rur = transpose([WeatherData.temp])-273.15;
        T_can = transpose([UCMData.canTemp])-273.15;
        T_dp = transpose([UCMData.Tdp]);
        T_ubl = transpose([UBLData.ublTemp])-273.15;
        T_road = transpose([UCMData.roadTemp])-273.15;
        T_wall = transpose([UCMData.wallTemp])-273.15;
        T_roof = transpose([UCMData.roofTemp])-273.15;
        Rh_can = transpose([UCMData.canRHum]);
        Rh_rur = transpose([WeatherData.rHum]);
        wind_can = transpose([UCMData.canWind]);
        wind_rur = transpose([WeatherData.wind]);
        SolRRoad = transpose([UCMData.SolRecRoad]);
        SolRWall = transpose([UCMData.SolRecWall]);
        SolRRoof = transpose([UCMData.SolRecRoof]);
        UHI = T_can - T_rur;
        Q_wall = transpose([UCMData.Q_wall]);
        Q_road = transpose([UCMData.Q_road]);
        Q_roof = transpose([UCMData.Q_roof]);
        Q_window = transpose([UCMData.Q_window]);
        Q_vent = transpose([UCMData.Q_vent]);
        Q_ubl = transpose([UCMData.Q_ubl]);
        Q_hvac = transpose([UCMData.Q_hvac]);
        Q_traffic = transpose([UCMData.sensAnthrop]);
        P_elec = transpose([UCMData.ElecTotal]);
        P_gas = transpose([UCMData.GasTotal]);
        Q_tree = transpose([UCMData.treeSensHeat]*UCM.treeCoverage);
        Q_urb = Q_wall + Q_road + Q_vent + Q_window + Q_hvac + Q_traffic + Q_tree + Q_roof;

        % Clear data in sheets & write heading
        blank = repmat(char(0),[9000 27]);
        Sheets = Excel.ActiveWorkBook.Sheets;

        for i = 2:(4+numel(BEM))
            sheet = get(Sheets, 'Item', i);
            invoke(sheet, 'Activate');
            Activesheet = Excel.Activesheet;
            ActivesheetRange = get(Activesheet,'Range','A1:AA9000');
            set(ActivesheetRange, 'Value', blank);
        end        
        heading = {'TIME','T_RUR','T_UCL','RH_RUR','RH_CAN','Wind_RUR','Wind_CAN','UHI','T_DP','T_UBL ','T_ROAD', 'T_WALL','T_ROOF', 'SolRRoad', 'SolRWall', 'SolRRoof','Q_WALL','Q_ROAD','Q_ROOF','Q_WNDW','Q_VENT','Q_UBL','Q_HVAC','Q_TRAF','P_ELC ','P_GAS','Q_URB'};
        
        % Sheet 2: Average data
        sheet = get(Sheets, 'Item', 2);
        invoke(sheet, 'Activate');
        Activesheet = Excel.Activesheet;
        ActivesheetRange = get(Activesheet,'Range','A1:AA1');
        set(ActivesheetRange, 'Value', heading);
        ActivesheetRange = get(Activesheet,'Range','A27:AA27');
        set(ActivesheetRange, 'Value', heading);
        
        hours = transpose(1:1:24);
        hT_rur = zeros (24,1);
        hT_can = zeros (24,1);
        hT_dp = zeros (24,1);
        hT_ubl = zeros (24,1);
        hT_road = zeros (24,1);
        hT_wall = zeros (24,1);
        hT_roof = zeros (24,1);
        hRh_can = zeros (24,1);
        hRh_rur = zeros (24,1);
        hwind_can = zeros (24,1);
        hwind_rur = zeros (24,1);
        hSolRRoad = zeros (24,1);
        hSolRWall = zeros (24,1);
        hSolRRoof = zeros (24,1);
        hQ_wall = zeros (24,1);
        hQ_roof = zeros (24,1);
        hQ_window = zeros (24,1);
        hQ_road = zeros (24,1);
        hQ_vent = zeros (24,1);
        hQ_ubl = zeros (24,1);
        hQ_hvac = zeros (24,1);
        hQ_traffic = zeros(24,1);
        hP_elec = zeros (24,1);
        hP_gas = zeros (24,1);
        hQ_urb = zeros (24,1);

        days = simTime.days;        % Number of days in the simulation
        for i = 1:N
            hour = mod(i,24);
            days = simTime.days;
            if hour == 0
                hour = 24;
            end
            hT_rur (hour) = hT_rur (hour) + T_rur(i)/days;
            hT_can (hour) = hT_can (hour) + T_can(i)/days;
            hT_dp (hour) = hT_dp (hour) + T_dp(i)/days;
            hT_ubl (hour) = hT_ubl (hour) + T_ubl(i)/days;
            hT_road (hour) = hT_road (hour) + T_road(i)/days;        
            hT_wall (hour) = hT_wall (hour) + T_wall(i)/days;        
            hT_roof (hour) = hT_roof (hour) + T_roof(i)/days;    
            hRh_can (hour) = hRh_can (hour) + Rh_can(i)/days;
            hRh_rur (hour) = hRh_rur (hour) + Rh_rur(i)/days;
            hwind_can (hour) = hwind_can (hour) + wind_can(i)/days;
            hwind_rur(hour) = hwind_rur (hour) + wind_rur(i)/days;
            hSolRRoad (hour) = hSolRRoad(hour) + SolRRoad(i)/days;
            hSolRWall (hour) = hSolRWall(hour) + SolRWall(i)/days;
            hSolRRoof (hour) = hSolRRoof(hour) + SolRRoof(i)/days;
            hQ_wall (hour) = hQ_wall (hour) + Q_wall(i)/days;
            hQ_roof (hour) = hQ_roof (hour) + Q_roof(i)/days;
            hQ_window (hour) = hQ_window (hour) + Q_window(i)/days;
            hQ_road (hour) = hQ_road (hour) + Q_road(i)/days;
            hQ_vent (hour) = hQ_vent (hour) + Q_vent(i)/days;
            hQ_ubl (hour) = hQ_ubl (hour) + Q_ubl(i)/days;
            hQ_hvac (hour) = hQ_hvac (hour) + Q_hvac(i)/days;
            hQ_traffic (hour) = hQ_traffic (hour) + Q_traffic(i)/days;
            hP_elec (hour) = hP_elec(hour) + P_elec(i)/days;
            hP_gas (hour) = hP_gas(hour) + P_gas(i)/days;
            hQ_urb (hour) = hQ_urb(hour) + Q_urb(i)/days;
        end
        hUHI = hT_can - hT_rur;

        % Sheet 2: UCM data (hourly average & hourly data)
        data = [hours,hT_rur,hT_can,hRh_rur,hRh_can,hwind_rur,hwind_can,hUHI,hT_dp,hT_ubl,hT_road,hT_wall,hT_roof,hSolRRoad,hSolRWall,hSolRRoof,hQ_wall,hQ_road,hQ_roof,hQ_window,hQ_vent,hQ_ubl,hQ_hvac,hQ_traffic,hP_elec,hP_gas,hQ_urb];
        ActivesheetRange = get(Activesheet,'Range','A2:AA25');
        set(ActivesheetRange, 'Value', data);

        data = [time,T_rur,T_can,Rh_rur,Rh_can,wind_rur,wind_can,UHI,T_dp,T_ubl,T_road,T_wall,T_roof,SolRRoad,SolRWall,SolRRoof,Q_wall,Q_road,Q_roof,Q_window,Q_vent,Q_ubl,Q_hvac,Q_traffic,P_elec,P_gas,Q_urb];
        datarange = strcat('A28:AA',int2str(numel(time)+27));
        ActivesheetRange = get(Activesheet,'Range',datarange);
        set(ActivesheetRange, 'Value', data);
        
        % Sheet 3 and onwards: Building data (hourly average & hourly data)
        hbTemp = zeros (24,numel(BEM));
        hbRHum = zeros (24,numel(BEM));
        hbPelec = zeros (24,numel(BEM));
        hbQgas = zeros (24,numel(BEM));
        hbPequip = zeros (24,numel(BEM));
        hbPlight = zeros (24,numel(BEM));
        hbQocc = zeros (24,numel(BEM));
        hbFluxMass = zeros (24,numel(BEM));
        hbFluxRoof = zeros(24,numel(BEM));
        hbFluxWall = zeros (24,numel(BEM));
        hbFluxSolar = zeros (24,numel(BEM));
        hbFluxWindow = zeros (24,numel(BEM));
        hbFluxInfil = zeros (24,numel(BEM));
        hbFluxVent = zeros (24,numel(BEM));
        hbCoolConsump = zeros (24,numel(BEM));
        hbHeatConsump = zeros (24,numel(BEM));
        hbCOP = zeros(24,numel(BEM));
        hbCoolDemand = zeros(24,numel(BEM));
        hbTwallin = zeros(24,numel(BEM));
        hbTroofin = zeros(24,numel(BEM));
        hbTmassin = zeros(24,numel(BEM));

        for i = 1:N
            hour = mod(i,24);
            if hour == 0
                hour = 24;
            end
            hbTemp (hour,:) = hbTemp (hour,:) + bTemp(i,:)/days;
            hbRHum (hour,:) = hbRHum (hour,:) + bRHum(i,:)/days;
            hbPelec (hour,:) = hbPelec (hour,:) + bPelec(i,:)/days;
            hbQgas (hour,:) = hbQgas (hour,:) + bQgas(i,:)/days;
            hbPequip (hour,:) = hbPequip (hour,:) + bPequip(i,:)/days;
            hbPlight (hour,:) = hbPlight (hour,:) + bPlight(i,:)/days;
            hbQocc (hour,:) = hbQocc (hour,:) + bQocc(i,:)/days;
            hbFluxMass (hour,:) = hbFluxMass (hour,:) + bFluxMass(i,:)/days;
            hbFluxRoof (hour,:) = hbFluxRoof (hour,:) + bFluxRoof(i,:)/days;
            hbFluxWall (hour,:) = hbFluxWall (hour,:) + bFluxWall(i,:)/days;
            hbFluxSolar (hour,:) = hbFluxSolar (hour,:) + bFluxSolar(i,:)/days;
            hbFluxWindow (hour,:) = hbFluxWindow (hour,:) + bFluxWindow(i,:)/days;
            hbFluxInfil (hour,:) = hbFluxInfil (hour,:) + bFluxInfil(i,:)/days;
            hbFluxVent (hour,:) = hbFluxVent (hour,:) + bFluxVent(i,:)/days;
            hbCoolConsump (hour,:) = hbCoolConsump (hour,:) + bCoolConsump(i,:)/days;
            hbHeatConsump (hour,:) = hbHeatConsump (hour,:) + bHeatConsump(i,:)/days;
            hbCOP (hour,:) = hbCOP (hour,:) + bCOP(i,:)/days;
            hbCoolDemand (hour,:) = hbCoolDemand (hour,:) + bCoolDemand(i,:)/days;
            hbTwallin (hour,:) = hbTwallin (hour,:) + bTwallin(i,:)/days;
            hbTroofin (hour,:) = hbTroofin (hour,:) + bTroofin(i,:)/days;
            hbTmassin (hour,:) = hbTmassin (hour,:) + bTmassin(i,:)/days;
        end
        hbTemp = hbTemp - 273.15;
        bTemp = bTemp -273.15;
        hbTwallin = hbTwallin - 273.15;
        bTwallin = bTwallin - 273.15;
        hbTroofin = hbTroofin - 273.15;
        bTroofin = bTroofin - 273.15;
        bTmassin = bTmassin - 273.15;
        hbTmassin = hbTmassin - 273.15;

        heading = {'TIME','T_indoor','Rhum_ind','T_wall','T_ceil','T_mass','Pelec','Qgas','Pequip','Plight', 'Qocc', 'Qmass','Qroof','Qwall','Qsolar','Qwindow','Qinfil','Qvent','Pcool','QheatDmd','QCoolDmd','COP'};
        for i = 1:numel(BEM)
            
            % Write headings
            sheet = get(Sheets, 'Item', i+2);
            invoke(sheet, 'Activate');
            Activesheet = Excel.Activesheet;
            ActivesheetRange = get(Activesheet,'Range','A1:V1');
            set(ActivesheetRange, 'Value', heading);
            ActivesheetRange = get(Activesheet,'Range','A27:V27');
            set(ActivesheetRange, 'Value', heading);

            % Write 24-hour average data
            data = [hours,hbTemp(:,i),hbRHum(:,i),hbTwallin(:,i),hbTroofin(:,i),hbTmassin(:,i),hbPelec(:,i),hbQgas(:,i),hbPequip(:,i),hbPlight(:,i),hbQocc(:,i),hbFluxMass(:,i),hbFluxRoof(:,i),hbFluxWall(:,i),hbFluxSolar(:,i),...
                hbFluxWindow(:,i),hbFluxInfil(:,i),hbFluxVent(:,i),hbCoolConsump(:,i),hbHeatConsump(:,i),hbCoolDemand(:,i),hbCOP(:,i)];
            ActivesheetRange = get(Activesheet,'Range','A2:V25');
            set(ActivesheetRange, 'Value', data);
            
            % Write hourly datas
            data = [time,bTemp(:,i),bRHum(:,i),bTwallin(:,i),bTroofin(:,i),bTmassin(:,i),bPelec(:,i),bQgas(:,i),bPequip(:,i),bPlight(:,i),bQocc(:,i),bFluxMass(:,i),bFluxRoof(:,i),bFluxWall(:,i),bFluxSolar(:,i),...
                bFluxWindow(:,i),bFluxInfil(:,i),bFluxVent(:,i),bCoolConsump(:,i),bHeatConsump(:,i),bCoolDemand(:,i),bCOP(:,i)];
            datarange = strcat('A28:V',int2str(numel(time)+27));
            ActivesheetRange = get(Activesheet,'Range',datarange);
            set(ActivesheetRange, 'Value', data);
        end
        
        % Write UWG version & disclaimer
        release = strcat('Urban Weather Generator v',num2str(ver));
        disclaimer = 'DISCLAIMER: The data is provided as is without warranty of any kind, either expressed or implied.';

        sheet = get(Sheets, 'Item', 1);
        invoke(sheet, 'Activate');
        Activesheet = Excel.Activesheet;
        ActivesheetRange = get(Activesheet,'Range','A1:A1');
        set(ActivesheetRange, 'Value', release);
        ActivesheetRange = get(Activesheet,'Range','A2:A2');
        set(ActivesheetRange, 'Value', disclaimer);
        
        Save(Workbook);
        delete(Excel);
        clear Excel;

    end

    fclose all;

    if fullyScripted
        disp('Inputs scripted, supressing pop-up notification...');
    else
        h = msgbox('Urban Weather Generation Complete',strcat('UWG',num2str(ver)),'help');
    end

end

function [newmat, newthickness] = procMat(materials,max_thickness,min_thickness)
    % Pocesses material layer so that a material with single
    % layer thickness is divided into two and material layer that is too
    % thick is subdivided
    
    newmat = [];
    newthickness = [];

    k = materials.thermalConductivity;
    Vhc = materials.volumetricHeatCapacity;
    if numel(materials.thickness)>1
        for j = 1:numel(materials.thickness)
            % Break up each layer that's more than 5cm thick
            if materials.thickness(j) > max_thickness
                nlayers = ceil(materials.thickness(j)/max_thickness);
                for l = 1:nlayers
                    newmat = [newmat Material(k{j},Vhc{j})];
                    newthickness = [newthickness; materials.thickness(j)/nlayers];
                end
                
            % Material that's less then min_thickness is not added.
            elseif materials.thickness(j) < min_thickness
%                 newmat = [newmat Material(k{j},Vhc{j})];
%                 newthickness = [newthickness; min_thickness];
                disp('WARNING: Material layer found too thin (<1cm), ignored');
            else
                newmat = [newmat Material(k{j},Vhc{j})];
                newthickness = [newthickness; materials.thickness(j)];
            end
        end
    else
        % Divide single layer into two (UWG assumes at least 2 layers)
        if materials.thickness > max_thickness
            nlayers = ceil(materials.thickness/max_thickness);
            for l = 1:nlayers
                newmat = [newmat Material(k,Vhc)];
                newthickness = [newthickness; materials.thickness/nlayers];
            end
            
        % Material should be at least 1cm thick, so if we're here, 
        % should give warning and stop. Only warning given for now.
        elseif materials.thickness < min_thickness*2
            newthickness = [min_thickness/2; min_thickness/2];
            newmat = [Material(k,Vhc) Material(k,Vhc)];
            disp('WARNING: a thin (<2cm) single layer element found');
            disp('May cause error');

        else
            newthickness = [materials.thickness/2; materials.thickness/2];
            newmat = [Material(k,Vhc) Material(k,Vhc)];
        end
    end
end

function [Tdb, w, phi, h, Tdp, v] = Psychrometrics (Tdb_in, w_in, P)
    % Modified version of Psychometrics by Tea Zakula
    % MIT Building Technology Lab
    % Tdb (dry bulb temperature) and Tdp(dew point temperature) in C
    % w (humidity ratio) in kg/kg of dry air
    % phi (relative humidity) in %
    % h (enthalpy) in J/kg of dry air
    % v (specific volume) in m3/kg of dry air
    % P (Atmospheric Station Pressure) in Pa

    c_air = 1006;   %J/kg, value from ASHRAE Fundamentals
    hlg = 2501000;  %J/kg, value from ASHRAE Fundamentals
    cw  = 1860;     %J/kg, value from ASHRAE Fundamentals
    P = P/1000;     % convert from Pa to kPa

    Tdb = Tdb_in - 273.15;
    w = w_in;

    % phi calculation from Tdb and w
    Pw = w*P/(0.621945+w);              % Partial pressure of water wapor
    Pws = Saturation_pressure(Tdb);
    phi = Pw/Pws*100;

    h = c_air*Tdb+w*(hlg+cw*Tdb);       % Enthalpy
    v = 0.287042*(Tdb+273.15)*(1+1.607858*w)/P; % Specific volume 

    % Dew point 
    pw = (P*w)/(0.621945+w); % water vapor partial pressure in kPa
    alpha = log(pw);
    Tdp = 6.54 + 14.526*alpha+0.7389*(alpha^2)+0.09486*(alpha^3)+0.4569*(pw^0.1984); % valid for Tdp between 0 C and 93 C
end

function [Pws] = Saturation_pressure(Tdb)
    T = Tdb+273.15;
    Pws = exp(-(5.8002206e3)/T+1.3914993+-(4.8640239e-2)*T+(4.1764768e-5)*(T^2)-(1.4452093e-8)*(T^3)+6.5459673*log(T)); %in Pa
    Pws = Pws/1000; % in kPa
end