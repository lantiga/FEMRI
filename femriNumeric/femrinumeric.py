""" This module loads all the classes from the femriNumeric library into its
namespace.  This is a required module."""

import os
import vtk

__all__ = []

for item in __all__:
        exec('import '+item)

if os.name == 'posix':
        from libfemriNumericPython import *
else:
        from femriNumericPython import *
    
