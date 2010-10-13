#!/usr/bin/env python

## Program:   femrikspacezeropadding.py
## Module:    $RCSfile: femrikspaceerror.py,v $
## Language:  Python
## Date:      $Date: 2008/11/03 17:01:14 $
## Version:   $Revision: 1.1 $

##   Copyright (c) Luca Antiga, David Steinman. All rights reserved.
##   See LICENCE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

import sys
import vtk
from vmtk import pypes
import femricommon

femrikspaceerror = 'femriKSpaceError'

class femriKSpaceError(pypes.pypeScript):

    def __init__(self):

        pypes.pypeScript.__init__(self)

        self.KSpace = None
        self.KSpace2 = None
        self.ErrorImage = None
        
        self.SetScriptName('femrikspaceerror')
        self.SetInputMembers([
            ['KSpace','kspace','vtkImageData',1,'','','vmtkimagereader'],
            ['KSpace2','kspace2','vtkImageData',1,'','','vmtkimagereader']
            ])
        self.SetOutputMembers([
            ['ErrorImage','o','vtkImageData',1,'','','vmtkimagewriter']
            ])

    def Execute(self):

        if not self.KSpace:
            self.PrintError('Error: No input KSpace.')

        imageMathematics = vtk.vtkImageMathematics()
        imageMathematics.SetInput(self.KSpace)
        imageMathematics.SetInput2(self.KSpace2)
        imageMathematics.SetOperationToSubtract()
        imageMathematics.Update()

        imageMagnitude = vtk.vtkImageMagnitude()
        imageMagnitude.SetInput(imageMathematics.GetOutput())
        imageMagnitude.Update()

        self.ErrorImage = imageMagnitude.GetOutput()


if __name__=='__main__':
    main = pypes.pypeMain()
    main.Arguments = sys.argv
    main.Execute()
