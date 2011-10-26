#!/usr/bin/env python

## Program:   femrisymbolicsimulator.py
## Module:    $RCSfile: femrinumericsurfacekspace.py,v $
## Language:  Python
## Date:      $Date: 2008/11/04 11:23:41 $
## Version:   $Revision: 1.3 $

##   Copyright (c) Luca Antiga, David Steinman. All rights reserved.
##   See LICENCE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

import sys
import vtk
import femrinumeric
import femrikspace

from vmtk import pypes

femrinumericsurfacekspace = 'femriNumericSurfaceKSpace'

class femriNumericSurfaceKSpace(femrikspace.femriKSpace):

    def __init__(self):

        femrikspace.femriKSpace.__init__(self)

        self.Surface = None
        
        self.MagnetizationValue = 1.0

        self.UseExactAlgorithm = True

        self.QuadratureOrder = 1
        self.NumberOfSubdivisions = 0
        self.UseOptimalAlgorithm = 0
        self.UseMonolithic = 0
        self.ErrorThreshold = 1E-4

        self.CellNormalsArrayName = ''

        self.SetScriptName('femrinumericsurfacekspace')
        self.SetInputMembers([
            ['Surface','i','vtkPolyData',1,'','','vmtksurfacereader'],
            ['MagnetizationValue','magnetizationvalue','float',1],
            ['UseExactAlgorithm','useexact','bool',1],
            ['QuadratureOrder','qorder','int',1,'(0,)'],
            ['UseOptimalAlgorithm','useoptimal','bool',1],
            ['UseMonolithic','usemonolithic','bool',1],
            ['CellNormalsArrayName','cellnormals','str',1],
            ['ErrorThreshold','error','float',1,'(0.0,)']
            ])
        self.SetOutputMembers([
            ])

    def AcquireKSpaceExact(self,surface,origin,spacing):
        if self.UseMonolithic:
            kSpaceAcquisition = femrinumeric.vtkfemriPolyDataExactKSpaceGeneratorMonolithic()
        else:
            kSpaceAcquisition = femrinumeric.vtkfemriPolyDataExactKSpaceGenerator()
        kSpaceAcquisition.SetInput(self.Surface)
        kSpaceAcquisition.SetKSpaceDimensionality(self.KSpaceDimensionality)
        kSpaceAcquisition.SetMatrix(self.MatrixSize)
        kSpaceAcquisition.SetFOV(self.FOV)
        kSpaceAcquisition.SetOrigin(origin)
        kSpaceAcquisition.SetAcquireSymmetricKSpace(self.AcquireSymmetricKSpace)
        kSpaceAcquisition.SetMagnetizationValue(self.MagnetizationValue)
        if self.CellNormalsArrayName:
            kSpaceAcquisition.SetCellNormalsArrayName(self.CellNormalsArrayName)
        kSpaceAcquisition.Update()

        return kSpaceAcquisition.GetOutput()

    def AcquireKSpaceNumeric(self,surface,origin,spacing):
        
        kSpaceAcquisition = femrinumeric.vtkfemriPolyDataNumericKSpaceGenerator()
        kSpaceAcquisition.SetInput(self.Surface)
        kSpaceAcquisition.SetKSpaceDimensionality(self.KSpaceDimensionality)
        kSpaceAcquisition.SetMatrix(self.MatrixSize)
        kSpaceAcquisition.SetFOV(self.FOV)
        kSpaceAcquisition.SetOrigin(origin)
        kSpaceAcquisition.SetQuadratureOrder(self.QuadratureOrder)
        kSpaceAcquisition.SetAcquireSymmetricKSpace(self.AcquireSymmetricKSpace)
        kSpaceAcquisition.SetMagnetizationValue(self.MagnetizationValue)
        kSpaceAcquisition.SetUseOptimalAlgorithm(self.UseOptimalAlgorithm)
        kSpaceAcquisition.SetErrorThreshold(self.ErrorThreshold)
        kSpaceAcquisition.Update()

        self.PrintLog("Max quadrature order used: %d" % kSpaceAcquisition.GetMaximumQuadratureOrderUsed())
        self.PrintLog("No. of Gauss point evaluations: %d" % kSpaceAcquisition.GetNumberOfGaussPointEvaluations())
        self.PrintLog("Avg no. of Gauss point evaluations per kSpace location: %d" % (kSpaceAcquisition.GetNumberOfGaussPointEvaluations()/kSpaceAcquisition.GetOutput().GetNumberOfPoints()))

        return kSpaceAcquisition.GetOutput()
        
    def Execute(self):

        if self.Surface == None:
            self.PrintError('Error: no Surface.')

        acquiredSurface = self.Surface
        acquiredKSpace = None
        
        origin = [self.FOVCenter[0] - self.FOV[0]/2.0,
                  self.FOVCenter[1] - self.FOV[1]/2.0,
                  self.FOVCenter[2] - self.FOV[2]/2.0]

        spacing = [self.FOV[0] / self.MatrixSize[0],
                   self.FOV[1] / self.MatrixSize[1],
                   self.FOV[2] / self.MatrixSize[2]]

        if self.UseExactAlgorithm:
            acquiredKSpace = self.AcquireKSpaceExact(self.Surface,origin,spacing)
        else:
            acquiredKSpace = self.AcquireKSpaceNumeric(self.Surface,origin,spacing)

        self.KSpace = self.ComputeKSpaceOperation(acquiredKSpace)

if __name__=='__main__':
    main = pypes.pypeMain()
    main.Arguments = sys.argv
    main.Execute()
