import os
import re

fname = "ucm_properties.txt" #"bemdef_properties.txt"#"building_properties_python.txt"#"building_properties.txt" # insert your file path here
curr_dir = os.path.abspath(os.path.dirname(__file__))
fpath = os.path.join(curr_dir,fname)

def building_properties_matlab():
    f = open(fpath,'r')
    L = f.readlines()
    f.close()


    blst = L[:25]
    flst = L[26:51]
    fprintf = [L[78]] * len(blst)

    for i in xrange(len(flst)):
        #print i,blst[i]
        name = blst[i][3:]
        name = "".join(name.split()) # removes all white space from string

        fprintf[i] = fprintf[i][:14]+str(i)+fprintf[i][14:]

        fprintf[i] = fprintf[i].replace('uValue', name)
        fprintf[i] = "".join(fprintf[i].split())
        #print fprintf[i]
        #print "fclose(fileID{});".format(i)
        print r"fprintf(fileID{},'\n');".format(i)
        #flst[i] = flst[i][:6]+str(i)+flst[i][6:]
        #flst[i] = flst[i].replace('uValue', name)
        #flst[i] = "".join(flst[i].split())
        #print flst[i]

def building_properties_python():
    """
    matlab_uval = float(file_uval.next())
    assert refDOE[bldType][bldEra][climateZone].uValue == pytest.approx(matlab_uval, abs=1e-15),\
        "type={},era={},czone={}".format(bldType+1, bldEra+1, climateZone+1)
    """
    f = open(fpath,'r')
    L = f.readlines()
    f.close()

    blst = L[:25]

    a = "elif bpropid == 'XXX':\n"
    b = "\tassert refDOE[bldType][bldEra][climateZone].floorHeight == pytest.approx(matlab_ref_value, abs=1e-15),\n"
    c = "\t\t'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)\n"

    for i in xrange(len(blst)):
        name = "".join(blst[i].split()) # removes all white space from string
        #print "'{}',".format(name)
        a_ = a.replace('XXX', name)
        b_ = b.replace("floorHeight", name)

        s = a_ + b_ + c
        #print s

        fp = "matlab_uval_path = os.path.join(self.DIR_MATLAB_PATH,'matlab_ref_uValue_.txt')"
        fp = fp.replace("uval",name)
        fp = fp.replace("uValue",name)
        #print fp

def bemdef_properties():
    f = open(fpath,'r')
    L = f.readlines()
    f.close()

    blst = L[:8]
    for i in xrange(len(blst)):
        name = blst[i][5:].split(" ")[0]
        #print "'{}'".format(name)

    bemlst=["mass","wall","roof"]
    elemlst =[
    'albedo',
    'emissivity',
    'layerThickness',
    'layerThermalCond',
    'layerVolHeat',
    'vegCoverage',
    'layerTemp',
    'horizontal'
    ]

    #export to file.txt
    fp_ = "fprintf(fileIDZZ,'%.16f\n',refBEM(i,j,k).YY.XX);"
    for bi in xrange(len(bemlst)):
        for ei in xrange(len(elemlst)):
            boundary_prop = elemlst[ei]
            boundary = bemlst[bi]

            fp = fp_.replace("ZZ",boundary + "_" + boundary_prop)
            fp = fp.replace("YY",boundary)
            fp = fp.replace("XX",boundary_prop)
            fp = "".join(fp.split())
            #print fp

    #generate file write
    fid_ = "fileIDXX=fopen('C:/saeran/master/git/UWG_Python/tests/matlab_ref/matlab_readDOE/matlab_ref_bemdef_floorHeight_.txt','w');"
    for bi in xrange(len(bemlst)):
        for ei in xrange(len(elemlst)):
            boundary_prop = elemlst[ei]
            boundary = bemlst[bi]
            fid = fid_.replace("floorHeight_", boundary + "_" + boundary_prop)
            fid = fid.replace("XX", boundary + "_" + boundary_prop)
            #print fid

    #close file
    fid_ = "fclose(fileIDXX);"
    for bi in xrange(len(bemlst)):
        for ei in xrange(len(elemlst)):
            boundary_prop = elemlst[ei]
            boundary = bemlst[bi]
            fid = fid_.replace("XX", boundary + "_" + boundary_prop)
            #print fid

    #make assert statement in Python
    for bi in xrange(len(bemlst)):
        for ei in xrange(len(elemlst)):
            bemid = bemlst[bi] + "_" + elemlst[ei]
            a = "elif bemid == '{}':\n".format(bemid)
            if "layer" in bemid:
                b = "\tassert refBEM[bldType][bldEra][climateZone].{}.{}[0] == pytest.approx(matlab_ref_value, abs=1e-15),\n".format(bemlst[bi],elemlst[ei])
            else:
                b = "\tassert refBEM[bldType][bldEra][climateZone].{}.{} == pytest.approx(matlab_ref_value, abs=1e-15),\n".format(bemlst[bi],elemlst[ei])

            c = "\t\t'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)\n"
            print a+b+c

def building_properties_schedule():
    f = open(fpath,'r')
    L = f.readlines()
    f.close()

    blst = L[:13]
    for i in xrange(len(blst)):
        name = blst[i].split(" ")[0]
        #print "'{}',".format(name)

    schlst =[
    'Elec',
    'Light',
    'Gas',
    'Occ',
    'Cool',
    'Heat',
    'SWH',
    'Qelec',
    'Qlight',
    'Nocc',
    'Qgas',
    'Vent',
    'Vswh',
    ]

    #export to file.txt
    fp_ = "fprintf(fileIDZZ,'%.16f',refBEM(i,j,k).XX);"
    for bi in xrange(len(schlst)):
        sch = schlst[bi]

        fp = fp_.replace("ZZ",sch)
        fp = fp.replace("XX",sch)
        fp = "".join(fp.split())
        #print fp

    #generate file write
    fid_ = "fileIDXX=fopen('C:/saeran/master/git/UWG_Python/tests/matlab_ref/matlab_readDOE/matlab_ref_sch_ZZ_.txt','w');"
    for bi in xrange(len(schlst)):
        sch = schlst[bi]
        fid = fid_.replace("ZZ", sch)
        fid = fid.replace("XX", sch)
        #print fid

    #close file
    fid_ = "fclose(fileIDXX);"
    for bi in xrange(len(schlst)):
        sch = schlst[bi]
        fid = fid_.replace("XX", sch)
        #print fid

    #make assert statement in Python
    for si in xrange(len(schlst)):
        a = "elif schid == '{}':\n".format(schlst[si])

        b = "\tassert Schedule[bldType][bldEra][climateZone].{} == pytest.approx(matlab_ref_value, abs=1e-15),\n".format(schlst[si])

        c = "\t\t'btype={},era={},czone={}'.format(bldType+1, bldEra+1, climateZone+1)\n"

        print a+b+c

def simparam_properties():
    f = open(fpath,'r')
    L = f.readlines()
    f.close()

    simlst = L[:16]
    for i in xrange(len(simlst)):
        name = simlst[i].split(" ")[0]
        name = "".join(name.split()) # removes all white space from string
        #print "{}".format(name[17:-1])

    simlst = L[17:33]
    for i in xrange(len(simlst)):
        name = "".join(simlst[i].split())
        # %fprintf(fileID, '%.16f\n',horSol);
        print "fprintf(fileID, {}, obj.{});".format(r".16f\n", name)

def ucm_properties():
    f = open(fpath,'r')
    L = f.readlines()
    f.close()

    proplst = L[:31]
    for i in xrange(len(proplst)):
        name = proplst[i].split(" ")[0]
        name = "".join(name.split()) # removes all white space from string
        a = "{}".format(name[:])
        #print "fprintf(fileID, '{}', obj.{});".format(r".16f\n", a)
        print "self.uwg.UCM.{},".format(a)
if __name__ == "__main__":

    if fname=="":
        pass
    elif fname=="building_properties_matlab.txt":
        building_properties_matlab()
    elif fname=="building_properties_python.txt":
        building_properties_python()
    elif fname=="bemdef_properties.txt":
        bemdef_properties()
    elif fname=="building_properties_schedule.txt":
        building_properties_schedule()
    elif fname=="simparam_properties.txt":
        simparam_properties()
    elif fname=="ucm_properties.txt":
        ucm_properties()
