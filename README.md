# Urban Weather Generator

[![Build Status](https://travis-ci.org/ladybug-tools/uwg.svg?branch=master)](https://travis-ci.org/ladybug-tools/uwg)

The Urban Weather Generator (uwg) is a Python application for modeling the [urban heat island effect](https://en.wikipedia.org/wiki/Urban_heat_island). Specifically, it morphs rural [EnergyPlus weather (.epw) files](http://www.ladybug.tools/epwmap/) to reflect average conditions within the urban canyon using a range of properties including:

* Building geometry (including building height, ground coverage, window:wall area, and facade:site area)
* Building use (including program type, HVAC systems, and occupancy/equipment scheduling)
* Cooling system heat rejection to the outdoors (for Summer)
* Indoor heat leakage to the outdoors (for Winter)
* Urban materials (including the thermal mass, albedo and emissivity of roads, walls, and roofs)
* Anthropogenic heat from traffic (including traffic schedules)
* Vegetation coverage (both trees and shrubs)
* Atmospheric heat transfer from urban boundary and canopy layers

The [original Urban Weather Generator](http://urbanmicroclimate.scripts.mit.edu/uwg.php) was developed by Bruno Bueno for [his PhD thesis at MIT](https://dspace.mit.edu/handle/1721.1/59107).  Since this time, it has been validated 3 times and has been [enhanced by Aiko Nakano](https://dspace.mit.edu/handle/1721.1/108779).  In 2016, Joseph Yang also [improved the engine and added a range of building templates](https://dspace.mit.edu/handle/1721.1/107347).

This repository is a Python translation of the original [MATLAB Urban Weather Generator](https://github.com/hansukyang/UWG_Matlab).

# Example
Here is a Python example that shows how to create and run an Urban Weather Generator object. The example script is available [at resources/uwg_example.py](https://github.com/ladybug-tools/uwg/blob/master/resources/uwg_example.py). Run it through your command prompt in the main uwg directory with the following: ```python -m resources.uwg_example```

```python
from uwg import uwg

# Define the .epw, .uwg filenames to create an uwg object.
# uwg will look for the .epw file in the uwg/resources/epw folder,
# and the .uwg file in the uwg/resources/parameters folder.
epw_filename = "SGP_Singapore.486980_IWEC.epw"      # .epw file name
param_filename = "initialize_singapore.uwg"         # .uwg file name

# Initialize the UWG object and run the simulation
uwg_ = uwg(epw_filename, param_filename)
uwg_.run()
```

Here is the sample .uwg file used in the simulation above. The .uwg file is a a required input where the local building, urban, and geographic features are defined. These features are then used in the simulation to morph the .epw file. This file is available [at resources/initialize_singapore.uwg](https://github.com/ladybug-tools/uwg/blob/master/resources/initialize_singapore.uwg).

```
# =================================================
# REQUIRED PARAMETERS
# =================================================

# Urban characteristics
bldHeight,10,     # average building height (m)
bldDensity,0.5,   # urban area building plan density (0-1)
verToHor,0.8,     # urban area vertical to horizontal ratio
h_mix,1,           # fraction of building HVAC waste heat set to the street canyon [as opposed to the roof]
charLength,1000,  # dimension of a square that encompasses the whole neighborhood [aka. characteristic length] (m)
albRoad,0.1,      # road albedo (0 - 1)
dRoad,0.5,        # road pavement thickness (m)
kRoad,1,          # road pavement conductivity (W/m K)
cRoad,1600000,    # road volumetric heat capacity (J/m^3 K)
sensAnth,20,      # non-building sensible heat at street level [aka. heat from cars, pedestrians, street cooking, etc. ] (W/m^2)
latAnth,2,        # non-building latent heat (W/m^2) (currently not used)

# Climate Zone (Eg. City)   Zone number
# 1A(Miami)                     1
# 2A(Houston)                   2
# 2B(Phoenix)                   3
# 3A(Atlanta)                   4
# 3B-CA(Los Angeles)            5
# 3B(Las Vegas)                 6
# 3C(San Francisco)             7
# 4A(Baltimore)                 8
# 4B(Albuquerque)               9
# 4C(Seattle)                   10
# 5A(Chicago)                   11
# 5B(Boulder)                   12
# 6A(Minneapolis)               13
# 6B(Helena)                    14
# 7(Duluth)                     15
# 8(Fairbanks)                  16

zone,1,

# Vegetation parameters
vegCover,0.2,     # Fraction of the urban ground covered in grass/shrubs only (0-1)
treeCoverage,0.1, # Fraction of the urban ground covered in trees (0-1)
vegStart,4,       # The month in which vegetation starts to evapotranspire (leaves are out)
vegEnd,10,        # The month in which vegetation stops evapotranspiring (leaves fall)
albVeg,0.25,      # Vegetation albedo
latGrss,0.4,      # Fraction of the heat absorbed by grass that is latent (goes to evaporating water)
latTree,0.6,      # Fraction of the heat absorbed by trees that is latent (goes to evaporating water)
rurVegCover,0.9,  # Fraction of the rural ground covered by vegetation

# Traffic schedule [1 to 24 hour],
SchTraffic,
0.2,0.2,0.2,0.2,0.2,0.4,0.7,0.9,0.9,0.6,0.6,0.6,0.6,0.6,0.7,0.8,0.9,0.9,0.8,0.8,0.7,0.3,0.2,0.2, # Weekday
0.2,0.2,0.2,0.2,0.2,0.3,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.6,0.7,0.7,0.7,0.7,0.5,0.4,0.3,0.2,0.2, # Saturday
0.2,0.2,0.2,0.2,0.2,0.3,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.3,0.3,0.2,0.2, # Sunday

# Fraction of building stock for each DOE Building type (pre-80's build, 80's-present build, new)
# Note that sum(bld) must be equal to 1
bld,
0,0,0,    # FullServiceRestaurant
0,0,0,    # Hospital
0,0,0,    # LargeHotel
0,0.4,0,  # LargeOffice
0,0,0,    # MediumOffice
0,0.6,0,  # MidRiseApartment
0,0,0,    # OutPatient
0,0,0,    # PrimarySchool
0,0,0,    # QuickServiceRestaurant
0,0,0,    # SecondarySchool
0,0,0,    # SmallHotel
0,0,0,    # SmallOffice
0,0,0,    # Stand-aloneRetail
0,0,0,    # StripMall
0,0,0,    # SuperMarket
0,0,0,    # Warehouse

# =================================================
# OPTIONAL URBAN PARAMETERS
# =================================================

albRoof,0.5,  # roof albedo (0 - 1)
vegRoof,0.1,  # Fraction of the roofs covered in grass/shrubs (0-1)
glzR,0.5,     # Glazing Ratio. If not provided, all buildings are assumed to have 40% glazing ratio
hvac,0,       # HVAC TYPE; 0 = Fully Conditioned (21C-24C); 1 = Mixed Mode Natural Ventilation (19C-29C + windows open >22C); 2 = Unconditioned (windows open >22C)

# =================================================
# OPTIONAL PARAMETERS FOR SIMULATION CONTROL,
# =================================================

# Simulation parameters,
Month,1,        # starting month (1-12)
Day,1,          # starting day (1-31)
nDay,31,        # number of days to run simultion
dtSim,300,      # simulation time step (s)
dtWeather,3600, # weather time step (s)

autosize,0,     # autosize HVAC (1 for yes; 0 for no)
sensOcc,100,    # Sensible heat per occupant (W)
LatFOcc,0.3,    # Latent heat fraction from occupant (normally 0.3)
RadFOcc,0.2,    # Radiant heat fraction from occupant (normally 0.2)
RadFEquip,0.5,  # Radiant heat fraction from equipment (normally 0.5)
RadFLight,0.7,  # Radiant heat fraction from light (normally 0.7)

#Urban climate parameters
h_ubl1,1000,    # ubl height - day (m)
h_ubl2,80,      # ubl height - night (m)
h_ref,150,      # inversion height (m)
h_temp,2,       # temperature height (m)
h_wind,10,      # wind height (m)
c_circ,1.2,     # circulation coefficient (default = 1.2 per Bruno (2012))
c_exch,1,       # exchange coefficient (default = 1; ref Bruno (2014))
maxDay,150,     # max day threshold (W/m^2)
maxNight,20,    # max night threshold (W/m^2)
windMin,1,      # min wind speed (m/s)
h_obs,0.1,      # rural average obstacle height (m)
```
