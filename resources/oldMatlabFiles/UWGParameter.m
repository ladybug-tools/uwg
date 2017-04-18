% -------------------------------------------------------------------------
% Element/Material Definitions 
% -------------------------------------------------------------------------

% [conductivity (W m-1 K-1), Vol heat capacity (J m-3 K-1)]
bldMat = Material(0.67,1.2e6);      % material (concrete? reference?)
roadMat = Material(1.0,1.6e6);      % material (asphalt? reference?)

% Define & build base elements
% [albedo, emissivity, thicknesses (m)(outer layer first),
%  materials, vegetation coverage, initial temperature (K),
%  inclination (horizontal - 1, vertical - 0) ]
wall = Element(0.2,0.9,[0.01;0.05;0.1;0.05;0.01],...
    [bldMat;bldMat;bldMat;bldMat;bldMat],0.,300.,0);
roof = Element(0.2,0.9,[0.01;0.05;0.1;0.05;0.01],...
    [bldMat;bldMat;bldMat;bldMat;bldMat],0.,300.,1);
road = Element(0.5,0.95,[0.05;0.1;0.1;0.5;0.5],...
    [roadMat;roadMat;roadMat;roadMat;roadMat],0.2,300.,1);
rural = Element(0.1,0.95,[0.05;0.1;0.1;0.5;0.5],...
    [roadMat;roadMat;roadMat;roadMat;roadMat],0.73,300.,1);
mass = Element(0.7,0.9,[0.05;0.05],[bldMat;bldMat],0.,300.,0);

% -------------------------------------------------------------------------
% Simulation Parameters
% -------------------------------------------------------------------------
cityName = 'Singapore';   % For plot/printing
LAT = 1.37;
LON = 103.98;
ELEV = 0.1;
dtSim = 300;              % Sim time step
dtWeather = 3600;         % Weather data time-step
monthName = 'July';       % For plot/printing
MONTH = 7;                % Begin month
DAY = 30;                 % Begin day of the month
NUM_DAYS = 7;             % Number of days of simulation
autosize = 0;             % Autosize HVAC 
CityBlock (8,3) = Block(NUM_DAYS * 24,wall,roof,mass,road);

% Create simulation class (SimParam.m)
simTime = SimParam(dtSim,dtWeather,MONTH,DAY,NUM_DAYS);

% Read Rural weather data (EPW file - http://apps1.eere.energy.gov/)
climate_data = char('data/rural_weather_data_changi.epw');
weather = Weather(climate_data,simTime.timeInitial,simTime.timeFinal);

% Building definitions
% Residential building with AC
res_wAC = Building(3.0,... % floorHeight
    4.0,...               % nighttime internal heat gains (W m-2 floor)
    4.0,...               % daytime internal heat gains (W m-2 floor)
    0.2,...               % radiant fraction of internal gains
    0.2,...               % latent fraction of internal gains
    0.5,...               % Infiltration (ACH)
    0.0,...               % Ventilation (ACH)
    0.3,...               % glazing ratio
    2.715,...             % window U-value (W m-2 K)
    0.75,...              % window solar heat gain coefficient
    'AIR',...             % cooling condensation system type {'AIR','WATER'}
    2.5,...               % COP of the cooling system
    1.0,...               % fraction of waste heat released into the canyon
    297.,...              % daytime indoor cooling set-point (K)
    297.,...              % nighttime indoor cooling set-point (K)
    293.,...              % daytime indoor heating set-point (K)
    293.,...              % nighttime indoor heating set-point (K)
    225.,...              % rated cooling system capacity (W m-2 bld)
    0.9,...               % heating system efficiency (-)
    300.);                % intial indoor temp (K)

% -------------------------------------------------------------------------
% Urban Area Definitions (%% need to re-do this!)
% -------------------------------------------------------------------------

% Define Reference (RSMDef(lat,lon,height,initialTemp,initialPres,Param))
RSM = RSMDef(LAT,LON,ELEV,weather.staTemp(1),weather.staPres(1),Param);

T_init = weather.staTemp(1);
Hum_init = weather.staHum(1);
Wind_init = weather.staUmod(1);

UCM = UCMDef(bldHeight,bldDensity,verToHor,treeCoverage,sensAnthrop,latAnthrop,... 
    T_init,Hum_init,Wind_init,r_glaze,SHGC,alb_wall,road,rural);
UBL = UBLDef('C',1000.,weather.staTemp(1),Param.maxdx);
