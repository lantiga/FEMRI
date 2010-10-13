""" This module loads all the classes from the femriSymbolic library into its
namespace.  This is a required module."""

import os
import vtk

__all__ = []

for item in __all__:
        exec('import '+item)

if os.name == 'posix':
        from libfemriSymbolicPython import *
else:
        from femriSymbolicPython import *
    
