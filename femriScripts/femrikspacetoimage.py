#!/usr/bin/env python

## Program:   femrikspacetoimage.py
## Module:    $RCSfile: femrikspacetoimage.py,v $
## Language:  Python
## Date:      $Date: 2007/07/07 15:38:01 $
## Version:   $Revision: 1.2 $

##   Copyright (c) Luca Antiga, David Steinman. All rights reserved.
##   See LICENCE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

import sys
import vtk
from vmtk import pypes

femrikspacetoimage = 'femriKSpaceToImage'

class femriKSpaceToImage(pypes.pypeScript):

    def __init__(self):

        pypes.pypeScript.__init__(self)

        self.Image = None
        self.KSpace = None

        self.KSpaceDimensionality = 2

        self.SetScriptName('femrikspacetoimage')
        self.SetInputMembers([
            ['KSpace','kspace','vtkImageData',1,'','','vmtkimagereader'],
            ['KSpaceDimensionality','dimensionality','int',1,'(2,3)']
            ])
        self.SetOutputMembers([
            ['Image','o','vtkImageData',1,'','','vmtkimagewriter']
            ])

    def Execute(self):

        if not self.KSpace:
            self.PrintError('Error: no KSpace.')

        ifft = vtk.vtkImageRFFT()
        ifft.SetInput(self.KSpace)
        ifft.SetDimensionality(self.KSpaceDimensionality)
        ifft.Update()

        ifftMagnitude = vtk.vtkImageMagnitude()
        ifftMagnitude.SetInput(ifft.GetOutput())
        ifftMagnitude.Update()

        origin = self.KSpace.GetOrigin()
        kspacing = self.KSpace.GetSpacing()
        dimensions = self.KSpace.GetDimensions()
        spacing = [1.0/(dimensions[0]*kspacing[0]),1.0/(dimensions[1]*kspacing[1]),1.0/(dimensions[2]*kspacing[2])]
  
        imageInformation = vtk.vtkImageChangeInformation()
        imageInformation.SetInput(ifftMagnitude.GetOutput())
        imageInformation.SetOutputSpacing(spacing)
        imageInformation.SetOutputOrigin(origin)
        imageInformation.Update()
        
        self.Image = imageInformation.GetOutput()
  
       
if __name__=='__main__':
    main = pypes.pypeMain()
    main.Arguments = sys.argv
    main.Execute()
