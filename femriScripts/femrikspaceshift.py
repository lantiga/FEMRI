#!/usr/bin/env python

## Program:   femrikspaceshift.py
## Module:    $RCSfile: femrikspaceshift.py,v $
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

femrikspaceshift = 'femriKSpaceShift'

class femriKSpaceShift(pypes.pypeScript):

    def __init__(self):

        pypes.pypeScript.__init__(self)

        self.KSpace = None
        
        self.Translation = [0.0, 0.0, 0.0]
        
        self.SetScriptName('femrikspaceshift')
        self.SetInputMembers([
            ['KSpace','kspace','vtkImageData',1,'','','vmtkimagereader'],
            ['Translation','translation','float',3]
            ])
        self.SetOutputMembers([
            ['KSpace','okspace','vtkImageData',1,'','','vmtkimagewriter']
            ])

    def Execute(self):

        if not self.KSpace:
            self.PrintError('Error: No input KSpace.')

        kSpaceShift = femricommon.vtkfemriKSpaceShift()
        kSpaceShift.SetInput(self.KSpace)
        kSpaceShift.SetTranslation(self.Translation)
        kSpaceShift.Update()
  
        self.KSpace = kSpaceShift.GetOutput()
  

if __name__=='__main__':
    main = pypes.pypeMain()
    main.Arguments = sys.argv
    main.Execute()
