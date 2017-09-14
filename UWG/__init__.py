"""Urban Weather Generator Library."""

import sys
import os

__all__ = [
    "UWG",
    "utilities",
    "material",
    "element",
    "building",
    "BEMDef",
    "forcing",
    "param",
    "psychrometrics",
    "schdef",
    "simparam",
    "UCMDef",
    "urbflux",
    "weather",
    ]

if __name__ == '__main__':
    # path to the real library
    path_to_up_dir = os.path.join(os.path.dirname(__file__),'..')
    sys.path.insert(0, path_to_up_dir)
