"""
This file shouold be imported in any example file trying to use the GraphCodes or DrawCodes functions.
It adds the appropriate paths to the python interpreter.
Whenever a new folder is made with a new category of functions in the main repository, its path should be added here
Functions inside of folders can be imported with import {foldername}.{functionname} or just import {functionname}

Date: 2-20-2020
@author: Martin Ritter
"""

import os
import sys

rootDir  = os.path.dirname(os.path.realpath(__file__))
drawDir  = os.path.join(rootDir, 'DrawCodes')
graphDir = os.path.join(rootDir, 'GraphCodes')

if not drawDir in sys.path:
    sys.path.append(drawDir)

if not graphDir in sys.path:
    sys.path.append(graphDir)
