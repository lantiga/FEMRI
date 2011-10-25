#!/usr/bin/env python

## Program:   femrikspacecenter.py
## Module:    $RCSfile: femrikspacecenter.py,v $
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
import femricommon

femrikspacecenter = 'femriKSpaceCenter'

class femriKSpaceCenter(pypes.pypeScript):

    def __init__(self):

        pypes.pypeScript.__init__(self)

        self.KSpace = None
        
        self.SetScriptName('femrikspacecenter')
        self.SetInputMembers([
            ['KSpace','kspace','vtkImageData',1,'','','vmtkimagereader']
            ])
        self.SetOutputMembers([
            ['KSpace','okspace','vtkImageData',1,'','','vmtkimagewriter']
            ])

    def Execute(self):

        if not self.KSpace:
            self.PrintError('Error: No input KSpace.')

        kSpaceCenter = vtk.vtkImageFourierCenter()
        kSpaceCenter.SetInput(self.KSpace)
        kSpaceCenter.Update()
  
        self.KSpace = kSpaceCenter.GetOutput()
  

if __name__=='__main__':
    main = pypes.pypeMain()
    main.Arguments = sys.argv
    main.Execute()
