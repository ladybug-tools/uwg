# Urban Weather Generator

[![Build Status](https://travis-ci.org/saeranv/urbanWeatherGen.svg?branch=master)](https://travis-ci.org/saeranv/urbanWeatherGen)

The Urban Weather Generator (urbanWeatherGen) is a Python application for modeling the [urban heat island effect](https://en.wikipedia.org/wiki/Urban_heat_island).  Specifically, it morphs rural [EnergyPlus weather (.epw) files](http://www.ladybug.tools/epwmap/) to reflect average conditions within the urban canyon using a range of properties including:

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

# Quickstart Example
Here is a Python example that shows how to create and run an Urban Weather Generator object. The quickstart example file is available [at resources/quickstart.py](https://github.com/ladybug-tools/urbanWeatherGen/blob/master/resources/quickstart.py). Run it through your command prompt in the main urbanWeatherGen directory with the following: ```python -m resources.quickstart```

```python
from UWG import UWG

# Define the .epw, .uwg filenames to create an UWG object.
# UWG will look for the .epw file int the UWG/resources/epw folder,
# and the .uwg file in the UWG/resources/parameters folder.
epw_filename = "SGP_Singapore.486980_IWEC.epw"      # .epw file name
param_filename = "initialize_singapore.uwg"         # .uwg file name

# Initialize the UWG object and run the simulation
uwg = UWG(epw_filename, param_filename)
uwg.run()
```
