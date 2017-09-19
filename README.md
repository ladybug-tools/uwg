# Urban Weather Generator

[![Build Status](https://travis-ci.org/saeranv/UWG_Python.svg?branch=master)](https://travis-ci.org/saeranv/UWG_Python)

The Urban Weather Generator (urbanWeatherGen) is a Python application for modeling the [urban heat island effect](https://en.wikipedia.org/wiki/Urban_heat_island).  Specifically, it morphs rural [EnergyPlus weather files (.epw)](http://www.ladybug.tools/epwmap/) to reflect average conditions within the urban canyon using a range of properties including:

* Building geometry (including building height, ground coverage, and facade:site area)
* Urban materials (including the mass and reflectivity of roads, walls, and roofs)
* Anthropogenic heat from traffic (including traffic schedules)
* Cooling system heat rejection to the outdoors (for Summer)
* Indoor heat leakage to the outdoors (for Winter)
* Vegetation coverage (both trees and shrubs)


The [original Urban Weather Generator](http://urbanmicroclimate.scripts.mit.edu/uwg.php) was developed by Bruno Bueno for [his PhD thesis at MIT](https://dspace.mit.edu/handle/1721.1/59107).  Since this time, it ha been validated 3 times and has been [improved by Aiko Nakano](https://dspace.mit.edu/handle/1721.1/108779).  In 2016, Joseph Yang also [improved the engine and added a range of building templates](https://dspace.mit.edu/handle/1721.1/107347).

This repository is a Python translation of the original [MATLAB Urban Weather Generator](https://github.com/hansukyang/UWG_Matlab).
