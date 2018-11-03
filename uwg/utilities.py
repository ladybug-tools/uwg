"""Collection of useful methods."""
import os
from csv import reader as csv_reader
import sys

try:
    range = xrange
except NameError:
    pass


def zeros(h, w):
    """create a (h x w) matrix of zeros.

    Args:
        h: Height of the matrix.
        w: Width of the matrix.
    """
    return [[0 for x in range(w)] for y in range(h)]


def read_csv(file_name_):
    # open csv file and read
    if os.path.exists(file_name_):
        if sys.version_info[0] >= 3:
            file_ = open(file_name_, "r", errors='ignore')
        else:
            file_ = open(file_name_, "r")

        gen_ = csv_reader(file_, delimiter=",")
        L = [r for r in gen_]
        file_.close()
        return L
    else:
        raise Exception("File name: '{}' does not exist.".format(file_name_))


def is_near_zero(num, eps=1e-10):
    return abs(num) < eps


def str2fl(x):
    """Recurses through lists and converts lists of string to float

    Args:
        x: string or list of strings
    """
    def helper_to_fl(s_):
        """ deals with odd string imports converts to float"""
        if s_ == "":
            return "null"
        elif "," in s_:
            s_ = s_.replace(",", "")

        try:
            return float(s_)
        except:
            return (s_)

    fl_lst = []
    if isinstance(x[0], str):                # Check if list of strings, then sent to conversion
        for xi in range(len(x)):
            fl_lst.append(helper_to_fl(x[xi]))
    elif isinstance(x[0], list):                    # Check if list of lists, then recurse
        for xi in range(len(x)):
            fl_lst.append(str2fl(x[xi]))
    else:
        return False

    return fl_lst
