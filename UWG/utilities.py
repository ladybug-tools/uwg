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
