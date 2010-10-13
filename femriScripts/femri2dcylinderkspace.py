#!/usr/bin/env python

## Program:   femri2dcylinderkspace.py
## Module:    $RCSfile: femri2dcylinderkspace.py,v $
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

femri2dcylinderkspace = 'femri2DCylinderKSpace'

class femri2DCylinderKSpace(femrikspace.femriKSpace):

    def __init__(self):

        femrikspace.femriKSpace.__init__(self)

        self.SliceThickness = 1.0
        self.Radius = 1.0
        self.Center = [0.0,0.0,0.0]
        self.RotationAngle = 0.0
        self.TiltingAngle = 0.0

        self.KSpaceDimensionality = 2

        self.SetScriptName('femri2dcylinderkspace')
        self.SetInputMembers([
            ['SliceThickness','slicethickness','float',1,'(0.0,)'],
            ['Radius','radius','float',1,'(0.0,)'],
            ['Center','center','float',3],
            ['RotationAngle','rotation','float',1],
            ['TiltingAngle','tilting','float',1]
            ])
        self.SetOutputMembers([
            ])

    def Execute(self):

        origin = [self.FOVCenter[0] - self.FOV[0]/2.0, self.FOVCenter[1] - self.FOV[1]/2.0, 0]
      
        kSpaceGenerator = femricommon.vtkfemri2DCylinderKSpaceGenerator()
        kSpaceGenerator.SetRadius(self.Radius)
        kSpaceGenerator.SetTheta(math.radians(self.RotationAngle))
        kSpaceGenerator.SetCenter(self.Center)
        kSpaceGenerator.SetTiltingAngle(math.radians(self.TiltingAngle))
        kSpaceGenerator.SetFOV(self.FOV)
        kSpaceGenerator.SetSliceThickness(self.SliceThickness)
        kSpaceGenerator.SetOrigin(origin)
        kSpaceGenerator.SetMatrix(self.MatrixSize)
        kSpaceGenerator.Update()

        self.KSpace = self.ComputeKSpaceOperation(kSpaceGenerator.GetOutput()) 

if __name__=='__main__':
    main = pypes.pypeMain()
    main.Arguments = sys.argv
    main.Execute()
