"""Urban Weather Generator Library.
This module is only if you want to run UWG from cli
"""

import sys
import os
import UWG


#TODO: revise with better args to UWG in default
#TODO: probably better way to get sys.arg
if len(sys.argv) > 5:
    epw_dir = sys.argv[1]
    epw_file_name = sys.argv[2]
    uwg_param_dir = sys.arg[3]
    uwg_param_file_name = sys.arg[4]
else:
    du = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    epw_dir = os.path.join(du,"UWG\\data\\epw")
    epw_file_name = "SGP_Singapore.486980_IWEC.epw"
    uwg_param_dir = None
    uwg_param_file_name = None

UWG.UWG(epw_dir, epw_file_name, uwg_param_dir, uwg_param_file_name)
