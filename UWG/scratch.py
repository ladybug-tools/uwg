import pprint
import sys

pp = lambda t: pprint.pprint(t)
pl = lambda t: map(lambda t_: pp(t_), t)


def text_parse(s="self.densityProfC[0],self.z, self.dz"):
    s = str(s) if type(s) != type("") else s
    #remove white space and split
    sp = lambda t: "".join(t.split()).split(",")
    # process for matlab
    um = lambda t: map(lambda t_: "fprintf('{0}-%.16f\n',{1});".format(t_[5:7],\
        t_.replace("self","obj").replace("[","(").replace("]",")").replace("0","1")), sp(t))

    up = lambda t: map(lambda t_: "&print '{}-',{}&".format(t_[5:7],t_), sp(t))


    pl(um(s))
    print ''
    pl(up(s))

def scratch():
    for iz in xrange(0,17,1):
        print iz+1
        for izz in xrange(iz,0,-1):
            #for izz in reversed(range(iz)[1:]):
            print iz+1, '-', izz+1
        print '--------'

if __name__ == "__main__":
    scratch()
