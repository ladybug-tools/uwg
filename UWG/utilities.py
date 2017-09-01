"""Collection of useful methods."""

from csv import reader as csv_reader

def zeros(h, w):
    """create a matrix of zeros.

    Args:
        h: Height of the matrix.
        w: Width of the matrix.
    """
    return [[0 for x in range(w)] for y in range(h)]

def read_csv(file_name_):
    file_ = open(file_name_,"r")
    gen_ = csv_reader(file_, delimiter=",")
    list_ = map(lambda r: r,gen_)
    file_.close()
    return list_

def str2fl(x):
    #Recurses through lists and converts lists of string to float
    def helper_to_fl(s_):
        if s_ == "":
            return "null"
        if "," in s_:
            s_ = s_.replace(",","")
        try:
            return float(s_)
        except:
            return (s_)
    fl_lst = []
    if isinstance(x[0], basestring):
        return map(lambda s: helper_to_fl(s), x)
    elif type(x[0]) == type([]):
        for xi in xrange(len(x)):
            fl_lst.append(str2fl(x[xi]))
        return fl_lst
    else:
        print 'Fail to convert to list of floats; type error {a} is {b}'.format(a=x[0], b=type(x[0]))
        return False

def vector_times_scalar(vector,scalar):
    new_vector = []
    for i in xrange(len(vector)):
        coordinate = vector[i]
        new_vector.append(coordinate * scalar)
    return new_vector

def vector_add_scalar(vector,scalar):
    # based on matrix + scalar op in matlab
    new_vector = []
    for i in xrange(len(vector)):
        coordinate = vector[i]
        new_vector.append(coordinate + scalar)
    return new_vector
