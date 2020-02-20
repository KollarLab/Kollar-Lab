"""
This file should be imported in any example file trying to use the GraphCodes or DrawCodes functions.
It adds the root directory to the python path for the sessiob
Functions inside of folders can be imported with import {foldername}.{functionname}

Date: 2-20-2020
@author: Martin Ritter
"""

import os
import sys

rootDir  = os.path.dirname(os.path.realpath(__file__))

if not rootDir in sys.path:
    sys.path.append(rootDir)