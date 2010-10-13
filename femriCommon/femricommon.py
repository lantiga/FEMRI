""" This module loads all the classes from the vtkBvg library into its
namespace.  This is a required module."""

import os
import vtk

__all__ = []

for item in __all__:
        exec('import '+item)

if os.name == 'posix':
        from libfemriCommonPython import *
else:
        from femriCommonPython import *
    
