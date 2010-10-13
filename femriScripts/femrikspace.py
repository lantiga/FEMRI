#!/usr/bin/env python

## Program:   femrikspace.py
## Module:    $RCSfile: femrikspace.py,v $
## Language:  Python
## Date:      $Date: 2008/11/03 17:00:31 $
## Version:   $Revision: 1.3 $

##   Copyright (c) Luca Antiga, David Steinman. All rights reserved.
##   See LICENCE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

import sys
import vtk
from vmtk import pypes

femrikspace = 'femriKSpace'

class femriKSpace(pypes.pypeScript):

    def __init__(self):

        pypes.pypeScript.__init__(self)

        self.KSpace = None
        self.KSpaceOperation = 'none'
        
        self.FOV = [1.0, 1.0, 1.0]
        self.FOVCenter = [0.0, 0.0, 0.0]
        self.MatrixSize = [16, 16, 16]
        self.KSpaceDimensionality = 3
        
        self.AcquireSymmetricKSpace = 1

        self.MValue = 1.0

        self.SetScriptName('femrikspace')
        self.SetInputMembers([
            ['KSpace','kspace','vtkImageData',1,'','','vmtkimagereader'],
            ['KSpaceOperation','operation','str',1,'["add","subtract","none"]'],
            ['FOV','fov','float',3],
            ['FOVCenter','fovcenter','float',3],
            ['MatrixSize','matrix','int',3,'(0,)'],
            ['KSpaceDimensionality','dimensionality','int',1,'(2,3)'],
            ['MValue','mvalue','float',1],
            ['AcquireSymmetricKSpace','symmetric','bool',1]
            ])
        self.SetOutputMembers([
            ['KSpace','okspace','vtkImageData',1,'','','vmtkimagewriter'],
            ['KSpaceDimensionality','dimensionality','int',1],
            ['FOV','ofov','float',3],
            ['FOVCenter','ofovcenter','float',3],
            ['MatrixSize','omatrix','int',3],
            ['KSpaceDimensionality','odimensionality','int',1]
            ])

    def ComputeKSpaceOperation(self,kSpace):

        kSpaceMValue = vtk.vtkImageMathematics()
        kSpaceMValue.SetInput(kSpace)
        kSpaceMValue.SetOperationToMultiplyByK()
        kSpaceMValue.SetConstantK(self.MValue)
        kSpaceMValue.Update()

        kSpace = kSpaceMValue.GetOutput()

        if not self.KSpace:
            return kSpace
 
        if self.KSpaceOperation == 'none':
            return kSpace
       
        kSpaceMaths = vtk.vtkImageMathematics()
        if self.KSpaceOperation == 'add':
            kSpaceMaths.SetOperationToAdd()
        elif self.KSpaceOperation == 'subtract':
            kSpaceMaths.SetOperationToSubtract()
        else:
            self.PrintError('KSpaceOperation not supported.')
            return kSpace
            
        kSpaceMaths.SetInput1(self.KSpace)
        kSpaceMaths.SetInput2(kSpace)
        kSpaceMaths.Update()

        return kSpaceMaths.GetOutput()
        
    def Execute(self):
 
        pass 
  

if __name__=='__main__':
    main = pypes.pypeMain()
    main.Arguments = sys.argv
    main.Execute()
