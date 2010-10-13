#!/usr/bin/env python

## Program:   femriimagetokspace.py
## Module:    $RCSfile: femriimagetokspace.py,v $
## Language:  Python
## Date:      $Date: 2008/10/28 16:46:59 $
## Version:   $Revision: 1.1 $

##   Copyright (c) Luca Antiga, David Steinman. All rights reserved.
##   See LICENCE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

import sys
import vtk
from vmtk import pypes

femriimagetokspace = 'femriImageToKSpace'

class femriImageToKSpace(pypes.pypeScript):

    def __init__(self):

        pypes.pypeScript.__init__(self)

        self.Image = None
        self.KSpace = None

        self.KSpaceDimensionality = 2

        self.SetScriptName('femriimagetokspace')
        self.SetInputMembers([
            ['Image','o','vtkImageData',1,'','','vmtkimagereader'],
            ['KSpaceDimensionality','dimensionality','int',1,'(2,3)']
            ])
        self.SetOutputMembers([
            ['KSpace','kspace','vtkImageData',1,'','','vmtkimagewriter']
            ])

    def Execute(self):

        if not self.Image:
            self.PrintError('Error: no Image.')

        fft = vtk.vtkImageFFT()
        fft.SetInput(self.Image)
        fft.SetDimensionality(self.KSpaceDimensionality)
        fft.Update()

        origin = self.Image.GetOrigin()
        spacing = self.Image.GetSpacing()
        dimensions = self.Image.GetDimensions()
        kspacing = [1.0/(dimensions[0]*spacing[0]),1.0/(dimensions[1]*spacing[1]),1.0/(dimensions[2]*spacing[2])]
  
        imageInformation = vtk.vtkImageChangeInformation()
        imageInformation.SetInput(fft.GetOutput())
        imageInformation.SetOutputSpacing(kspacing)
        imageInformation.SetOutputOrigin(origin)
        imageInformation.Update()
        
        self.KSpace = imageInformation.GetOutput()
  
       
if __name__=='__main__':
    main = pypes.pypeMain()
    main.Arguments = sys.argv
    main.Execute()
