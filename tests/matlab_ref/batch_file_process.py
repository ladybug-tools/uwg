import os
import re

fname = "batch_input.txt"
curr_dir = os.path.abspath(os.path.dirname(__file__))
fpath = os.path.join(curr_dir,fname)

def batch_file_process():
    f = open(fpath,'r')
    L = f.readlines()
    f.close()

    # Clean raw data
    prefix2delete = ""
    prefix2add = ""
    objectref = ""
    methodref = ""
    i = 0
    while i < len(L):
        line = L[i].split(" ")[0]
        line = "".join(line.split()) # removes all white space from string

        if "prefix2delete" in line:
            prefix2delete = line.split("=")[-1]
            L.pop(i)
        elif "prefix2add" in line:
            prefix2add = line.split("=")[-1]
            L.pop(i)
        elif "objectref" in line:
            objectref = line.split("=")[-1]
            L.pop(i)
        elif "methodref" in line:
            methodref = line.split("=")[-1]
            L.pop(i)
        elif line == "":
            L.pop(i)
        else:
            line = line.split(",")[0] if "," in line else line

            L[i] = line
            i += 1

    # For uwg_python_val
    #for i in xrange(len(L)):
    #    L[i] = L[i].split(prefix2delete)[-1] if prefix2delete!="" else L[i]
    #    print '{prefix}.{valuename},'.format(valuename=L[i], prefix=prefix2add)

    # create three outputs
    print "%%{"
    print "fileID = fopen('..\\UWG_Python\\tests\\matlab_ref\\matlab_{}\\matlab_ref_{}.txt','w');".format(objectref,objectref+"_"+methodref)
    print "format long;"
    for i in xrange(len(L)):
        L[i] = L[i].split(prefix2delete)[-1] if prefix2delete!="" else L[i]
        print 'fprintf(fileID, "{value}", {prefix}.{valuename});'.format(value=r"%.16f\n", valuename=L[i], prefix=prefix2add)
    print "fclose(fileID);"
    print "%}%"


if __name__ == "__main__":
    batch_file_process()

"""
ARCHIVE_RSM_VDM:
self.uwg.rural.aeroCond        # Convection coef (refL uwg, eq.12)
self.uwg.rural.waterStorage    # thickness of water film (m) (only for horizontal surfaces)
self.uwg.rural.solAbs          # solar radiation absorbed (W m-2)
self.uwg.rural.lat             # surface latent heat flux (W m-2)
self.uwg.rural.sens            # surface sensible heat flux (W m-2)
self.uwg.rural.flux             # external surface heat flux (W m-2)"

prefix2delete=self.uwg.rural.
prefix2add=RSM
objectref=rsmdef
methodref=surfflux
"""
