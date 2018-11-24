import cPickle
import os

def get_bem_pickle(uwg_data, refdir, fname):

    # Pickle objects, protocol 1 b/c binary file
    pkl_file_path = os.path.join(refdir,fname+".pkl")
    with open(pkl_file_path, "wb") as outf:
        cPickle.dump(uwg_data.BEM, outf)

    return pkl_file_path

def get_uwg_file(uwg_object, refdir, fname):
    uwg_file_path = os.path.join(refdir,fname+".uwg")
    f = open(uwg_file_path, "w")

    f.write("# =================================================\n")
    f.write("# REQUIRED PARAMETERS\n")
    f.write("# =================================================\n")
    f.write("\n")
    f.write("# Urban characteristics\n")
    f.write("bldHeight,{},\n".format(uwg_object.bldHeight))
    f.write("bldDensity,{},\n".format(uwg_object.bldDensity))
    f.write("verToHor,{},\n".format(uwg_object.verToHor))
    f.write("h_mix,{},\n".format(uwg_object.h_mix))
    f.write("charLength,{},\n".format(uwg_object.charLength))  # dimension of a square that encompasses the whole neighborhood [aka. characteristic length] (m)
    f.write("albRoad,{},\n".format(uwg_object.alb_road))      # road albedo (0 - 1)
    f.write("dRoad,{},\n".format(uwg_object.d_road))        # road pavement thickness (m)
    f.write("kRoad,{},\n".format(uwg_object.kRoad))         # road pavement conductivity (W/m K)
    f.write("cRoad,{},\n".format(uwg_object.cRoad))    # road volumetric heat capacity (J/m^3 K)
    f.write("sensAnth,{},\n".format(uwg_object.sensAnth))      # non-building sensible heat at street level [aka. heat from cars, pedestrians, street cooking, etc. ] (W/m^2)
    f.write("latAnth,{},\n".format(uwg_object.latAnth))        # non-building latent heat (W/m^2) (currently not used)
    f.write("\n")
    f.write("zone,{},\n".format(uwg_object.zone + 1))
    f.write("\n")
    f.write("# Vegetation parameters\n")
    f.write("vegCover,{},\n".format(uwg_object.vegCover))     # Fraction of the urban ground covered in grass/shrubs only (0-1)
    f.write("treeCoverage,{},\n".format(uwg_object.treeCoverage)) # Fraction of the urban ground covered in trees (0-1)
    f.write("vegStart,{},\n".format(uwg_object.vegStart))       # The month in which vegetation starts to evapotranspire (leaves are out)
    f.write("vegEnd,{},\n".format(uwg_object.vegEnd))        # The month in which vegetation stops evapotranspiring (leaves fall)
    f.write("albVeg,{},\n".format(uwg_object.albVeg))      # Vegetation albedo
    f.write("rurVegCover,{},\n".format(uwg_object.rurVegCover))  # Fraction of the rural ground covered by vegetation
    f.write("latGrss,{},\n".format(uwg_object.latGrss))      # Fraction of the heat absorbed by grass that is latent. Used in UWG only to calculate sensible heat fraction.
    f.write("latTree,{},\n".format(uwg_object.latTree))      # Fraction of the heat absorbed by trees that is latent. Used in UWG only to calculate sensible heat fraction.
    f.write("\n")
    f.write("# Traffic schedule [1 to 24 hour],\n")
    f.write("SchTraffic,\n")
    for i in range(3):
        for j in range(24):
            f.write("{},".format(uwg_object.SchTraffic[i][j]))
        f.write("\n")
    f.write("\n")
    f.write("# Fraction of building stock for each DOE Building type (pre-80's build, 80's-present build, new)\n")
    f.write("# Note that sum(bld) must be equal to 1\n")
    f.write("bld,\n")
    for i in range(16):
        for j in range(3):
            f.write("{},".format(uwg_object.bld[i][j]))
        f.write("\n")
    f.write("\n")
    f.write("# =================================================\n")
    f.write("# OPTIONAL URBAN PARAMETERS\n")
    f.write("# =================================================\n")
    f.write("# If not provided, optional parameters are taken from corresponding DOE Reference building\n")
    f.write("albRoof,{},\n".format(uwg_object.albRoof if uwg_object.albRoof else ""))  # roof albedo (0 - 1)
    f.write("vegRoof,{},\n".format(uwg_object.vegRoof if uwg_object.vegRoof else ""))  # Fraction of the roofs covered in grass/shrubs (0 - 1)
    f.write("glzR,{},\n".format(uwg_object.glzR if uwg_object.glzR else ""))     # Glazing Ratio (0 - 1)
    f.write("SHGC,{},\n".format(uwg_object.SHGC if uwg_object.SHGC else ""))     # Solar Heat Gain Coefficient (0 - 1)
    f.write("albWall,{},\n".format(uwg_object.albWall if uwg_object.albWall else ""))  # wall albedo (0 - 1)
    f.write("\n")
    f.write("# =================================================\n")
    f.write("# OPTIONAL PARAMETERS FOR SIMULATION CONTROL,\n")
    f.write("# =================================================\n")
    f.write("\n")
    f.write("# Simulation parameters,\n")
    f.write("Month,{},\n".format(uwg_object.Month))        # starting month (1-12)
    f.write("Day,{},\n".format(uwg_object.Day))          # starting day (1-31)
    f.write("nDay,{},\n".format(uwg_object.nDay))        # number of days to run simultion
    f.write("dtSim,{},\n".format(uwg_object.dtSim))      # simulation time step (s)
    f.write("dtWeather,{},\n".format(uwg_object.dtWeather)), # weather time step (s)
    f.write("\n")
    f.write("# HVAC system and internal loads\n")
    f.write("autosize,{},\n".format(uwg_object.autosize))     # autosize HVAC (1 for yes; 0 for no)
    f.write("sensOcc,{},\n".format(uwg_object.sensOcc))    # Sensible heat per occupant (W)
    f.write("LatFOcc,{},\n".format(uwg_object.LatFOcc))    # Latent heat fraction from occupant (normally 0.3)
    f.write("RadFOcc,{},\n".format(uwg_object.RadFOcc))    # Radiant heat fraction from occupant (normally 0.2)
    f.write("RadFEquip,{},\n".format(uwg_object.RadFEquip))  # Radiant heat fraction from equipment (normally 0.5)
    f.write("RadFLight,{},\n".format(uwg_object.RadFLight))  # Radiant heat fraction from light (normally 0.7)
    f.write("\n")
    f.write("#Urban climate parameters\n")
    f.write("h_ubl1,{},\n".format(uwg_object.h_ubl1))    # ubl height - day (m)
    f.write("h_ubl2,{},\n".format(uwg_object.h_ubl2))      # ubl height - night (m)
    f.write("h_ref,{},\n".format(uwg_object.h_ref))      # inversion height (m)
    f.write("h_temp,{},\n".format(uwg_object.h_temp))       # temperature height (m)
    f.write("h_wind,{},\n".format(uwg_object.h_wind))     # wind height (m)
    f.write("c_circ,{},\n".format(uwg_object.c_circ))   # circulation coefficient (default = 1.2 per Bruno (2012))
    f.write("c_exch,{},\n".format(uwg_object.c_exch))      # exchange coefficient (default = 1; ref Bruno (2014))
    f.write("maxDay,{},\n".format(uwg_object.maxDay))     # max day threshold (W/m^2)
    f.write("maxNight,{},\n".format(uwg_object.maxNight))    # max night threshold (W/m^2)
    f.write("windMin,{},\n".format(uwg_object.windMin))      # min wind speed (m/s)
    f.write("h_obs,{},\n".format(uwg_object.h_obs))      # rural average obstacle height (m)

    f.close()

    return uwg_file_path

if x and y and z and run:

    _uwg_object = x
    _dir = y
    _fname = z

    print "Remember not to run_ the 'RunUWG' component. Just write the file!"
    uwg_file = get_uwg_file(_uwg_object, _dir, _fname)
    bem_file = get_bem_pickle(_uwg_object, _dir, _fname)

    fout = [uwg_file, bem_file]

    print fout
