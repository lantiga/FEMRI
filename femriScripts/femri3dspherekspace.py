#!/usr/bin/env python

## Program:   femri3dspherekspace.py
## Module:    $RCSfile: femri3dspherekspace.py,v $
## Language:  Python
## Date:      $Date: 2007/07/07 15:38:01 $
## Version:   $Revision: 1.2 $

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

femri3dspherekspace = 'femri3DSphereKSpace'

class femri3DSphereKSpace(femrikspace.femriKSpace):

    def __init__(self):

        femrikspace.femriKSpace.__init__(self)

        self.Center = [0.0,0.0,0.0]
        self.Radius = 1.0

        self.KSpaceDimensionality = 3

        self.SetScriptName('femri3dspherekspace')
        self.SetInputMembers([
            ['Center','center','float',3],
            ['Radius','radius','float',1,'(0.0,)']
            ])
        self.SetOutputMembers([
            ])

    def Execute(self):

        origin = [self.FOVCenter[0] - self.FOV[0]/2.0, self.FOVCenter[1] - self.FOV[1]/2.0, self.FOVCenter[2] - self.FOV[2]/2.0]
      
        kSpaceGenerator = femricommon.vtkfemri3DSphereKSpaceGenerator()
        kSpaceGenerator.SetCenter(self.Center)
        kSpaceGenerator.SetRadius(self.Radius)
        kSpaceGenerator.SetFOV(self.FOV)
        kSpaceGenerator.SetOrigin(origin)
        kSpaceGenerator.SetMatrix(self.MatrixSize)
        kSpaceGenerator.Update()

        self.KSpace = self.ComputeKSpaceOperation(kSpaceGenerator.GetOutput()) 

if __name__=='__main__':
    main = pypes.pypeMain()
    main.Arguments = sys.argv
    main.Execute()
