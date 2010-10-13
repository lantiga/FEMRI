#!/usr/bin/env python

## Program:   femriemtpykspace.py
## Module:    $RCSfile: femriemptykspace.py,v $
## Language:  Python
## Date:      $Date: 2008/11/03 17:01:14 $
## Version:   $Revision: 1.1 $

##   Copyright (c) Luca Antiga, David Steinman. All rights reserved.
##   See LICENCE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

import sys
import math
import vtk
import femricommon
import femrikspace

from vmtk import pypes

femriemptykspace = 'femriEmptyKSpace'

class femriEmptyKSpace(femrikspace.femriKSpace):

    def __init__(self):

        femrikspace.femriKSpace.__init__(self)

        self.DCValue = 0.0

        self.SetScriptName('femriemptykspace')
        self.SetInputMembers([
            ['DCValue','dcvalue','float',1]
            ])
        self.SetOutputMembers([
            ])

    def Execute(self):

        kSpacing = [1.0/self.FOV[0], 1.0/self.FOV[1], 1.0/self.FOV[2]]
      
        origin = [self.FOVCenter[0] - self.FOV[0]/2.0, self.FOVCenter[1] - self.FOV[1]/2.0, self.FOVCenter[2] - self.FOV[2]/2.0]
        
        if self.KSpaceDimensionality == 2:
            origin[2] = 0
        
        kSpace = vtk.vtkImageData()
        kSpace.SetScalarTypeToDouble()
        kSpace.SetNumberOfScalarComponents(2)
        kSpace.SetDimensions(self.MatrixSize)
        kSpace.SetSpacing(kSpacing)
        kSpace.SetOrigin(origin)
        kSpace.AllocateScalars()
        kSpace.GetPointData().GetScalars().FillComponent(0,0.0)
        kSpace.GetPointData().GetScalars().FillComponent(1,0.0)
        kSpace.GetPointData().GetScalars().SetComponent(0,0,self.DCValue)
  
        self.KSpace = self.ComputeKSpaceOperation(kSpace)


if __name__=='__main__':
    main = pypes.pypeMain()
    main.Arguments = sys.argv
    main.Execute()
