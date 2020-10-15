"""Collection of useful methods."""
from csv import reader as csv_reader
import sys

try:
    range = xrange
except NameError:
    pass


try:
    import math
    INFPOS = math.inf
    INFNEG = -1 * math.inf
except AttributeError:
    # python 2
    INFPOS = float('inf')
    INFNEG = float('-inf')


# DOE References
REF_BLDTYPE = ('fullservicerestaurant', 'hospital', 'largehotel', 'largeoffice',
               'medoffice', 'midriseapartment', 'outpatient', 'primaryschool',
               'quickservicerestaurant', 'secondaryschool', 'smallhotel',
               'smalloffice', 'standaloneretail', 'stripmall', 'supermarket',
               'warehouse')
REF_BUILTERA = ('pre80', 'pst80', 'new')
REF_ZONETYPE = ('1A', '2A', '2B', '3A', '3B-CA', '3B', '3C', '4A', '4B', '4C', '5A',
                '5B', '6A', '6B', '7', '8')
REF_BLDTYPE_SET = {'fullservicerestaurant', 'hospital', 'largehotel', 'largeoffice',
                   'medoffice', 'midriseapartment', 'outpatient', 'primaryschool',
                   'quickservicerestaurant', 'secondaryschool', 'smallhotel',
                   'smalloffice', 'standaloneretail', 'stripmall', 'supermarket',
                   'warehouse'}
REF_ZONETYPE_SET = {'1A', '2A', '2B', '3A', '3B-CA', '3B', '3C', '4A', '4B', '4C', '5A',
                    '5B', '6A', '6B', '7', '8'}
REF_BUILTERA_SET = {'pre80', 'pst80', 'new'}


def is_near_zero(num, eps=1e-10):
    return abs(float(num)) < eps


def read_csv(file_path):
    """Open csv file and read.

    Args:
        file_path: Text string for file path.

    Returns:
        List of file lines as str type.
    """
    if sys.version_info[0] >= 3:
        file_ = open(file_path, "r", errors='ignore')
    else:
        file_ = open(file_path, "r")

    gen_ = csv_reader(file_, delimiter=",")
    L = [r for r in gen_]
    file_.close()
    return L


def str2fl(x):
    """Recurses through lists and converts lists of string to float

    Args:
        x: string or list of strings
    """
    def helper_to_fl(s_):
        """Deals with odd string imports converts to float"""
        if s_ == "":
            return "null"
        elif "," in s_:
            s_ = s_.replace(",", "")

        try:
            return float(s_)
        except (ValueError, TypeError):
            return (s_)

    fl_lst = []
    if isinstance(x[0], str):
        # Check if list of strings then conversion
        for xi in range(len(x)):
            fl_lst.append(helper_to_fl(x[xi]))
    elif isinstance(x[0], list):
        # Check if list of lists, then recurse
        for xi in range(len(x)):
            fl_lst.append(str2fl(x[xi]))
    else:
        return False

    return fl_lst


def float_in_range(value, mi=INFNEG, ma=INFPOS, input_name=''):
    """Check a float value to be between minimum and maximum."""
    assert mi <= value <= ma, 'Input number {} must be between {} and {}. ' \
        'Got {}'.format(input_name, mi, ma, value)
    return value


def float_in_range_excl(value, mi=INFNEG, ma=INFPOS, input_name=''):
    """Check a float value to be greater than minimum and less than maximum."""
    assert mi < value < ma, 'Input number {} must be greater than {} ' \
        'and less than {}. Got {}'.format(input_name, mi, ma, value)
    return value


def float_in_range_excl_incl(value, mi=INFNEG, ma=INFPOS, input_name=''):
    """Check a float value to be greater than minimum and less than/equal to maximum."""
    assert mi < value <= ma, 'Input number {} must be greater than {} and less than ' \
        'or equal to {}. Got {}'.format(input_name, mi, ma, value)
    return value


def float_in_range_incl_excl(value, mi=INFNEG, ma=INFPOS, input_name=''):
    """Check a float value to be greater than/equal to minimum and less than maximum."""
    assert mi <= value < ma, 'Input number {} must be greater than or equal to {} ' \
        'and less than {}. Got {}'.format(input_name, mi, ma, value)
    return value


def int_in_range(value, mi=INFNEG, ma=INFPOS, input_name=''):
    """Check an integer value to be between minimum and maximum."""
    number = int(value)
    assert mi <= number <= ma, 'Input integer {} must be between {} and {}. ' \
        'Got {}.'.format(input_name, mi, ma, value)
    return number


def float_positive(value, input_name=''):
    """Check a float value to be positive."""
    return float_in_range(value, 0, INFPOS, input_name)


def int_positive(value, input_name=''):
    """Check if an integer value is positive."""
    return int_in_range(value, 0, INFPOS, input_name)
