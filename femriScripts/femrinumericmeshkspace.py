#!/usr/bin/env python

## Program:   femrisymbolicsimulator.py
## Module:    $RCSfile: femrinumericmeshkspace.py,v $
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
import femrinumeric
import femrikspace

from vmtk import pypes

femrinumericmeshkspace = 'femriNumericMeshKSpace'

class femriNumericMeshKSpace(femrikspace.femriKSpace):

    def __init__(self):

        femrikspace.femriKSpace.__init__(self)

        self.Mesh = None
        
        self.UniformMagnetization = 0
        self.MagnetizationValue = 1.0
        self.MagnetizationArrayName = ''

        self.SliceSelection = 0
        self.SliceProfile = 'ideal'
        self.SliceThickness = 1.0
        self.SliceSpacing = 1.0

        self.QuadratureOrder = 1
        self.NumberOfSubdivisions = 0
        self.UseOptimalAlgorithm = 0
        self.ErrorThreshold = 1E-4
        
        self.SetScriptName('femrinumericmeshkspace')
        self.SetInputMembers([
            ['Mesh','i','vtkUnstructuredGrid',1,'','','vmtkmeshreader'],
            ['MagnetizationArrayName','magnetizationarray','str',1],
            ['UniformMagnetization','uniformmagnetization','bool',1],
            ['MagnetizationValue','magnetizationvalue','float',1],
            ['SliceSelection','sliceselection','bool',1],
            ['SliceProfile','sliceprofile','str',1,'["ideal","trapezoidal","quadratic"]'],
            ['SliceThickness','slicethickness','float',1,'(0.0,)'],
            ['SliceSpacing','slicespacing','float',1,'(0.0,)'],
            ['QuadratureOrder','qorder','int',1,'(0,)'],
            ['NumberOfSubdivisions','subdivisions','int',1,'(0,)'],
            ['UseOptimalAlgorithm','useoptimal','bool',1],
            ['ErrorThreshold','error','float',1,'(0.0,)']
            ])
        self.SetOutputMembers([
            ])

    def AcquireKSpace(self,mesh,origin,spacing):
        
        kSpaceAcquisition = femrinumeric.vtkfemriUnstructuredGridNumericKSpaceGenerator()
        kSpaceAcquisition.SetInput(self.Mesh)
        kSpaceAcquisition.SetKSpaceDimensionality(self.KSpaceDimensionality)
        kSpaceAcquisition.SetMatrix(self.MatrixSize)
        kSpaceAcquisition.SetFOV(self.FOV)
        kSpaceAcquisition.SetOrigin(origin)
        kSpaceAcquisition.SetQuadratureOrder(self.QuadratureOrder)
        kSpaceAcquisition.SetNumberOfSubdivisions(self.NumberOfSubdivisions)
        kSpaceAcquisition.SetUseOptimalAlgorithm(self.UseOptimalAlgorithm)
        kSpaceAcquisition.SetErrorThreshold(self.ErrorThreshold)
        kSpaceAcquisition.SetAcquireSymmetricKSpace(self.AcquireSymmetricKSpace)
        kSpaceAcquisition.SetUniformMagnetization(self.UniformMagnetization)
        kSpaceAcquisition.SetMagnetizationValue(self.MagnetizationValue)
        kSpaceAcquisition.SetMagnetizationArrayName(self.MagnetizationArrayName)
        kSpaceAcquisition.SetSliceSelection(self.SliceSelection)
        kSpaceAcquisition.SetSliceThickness(self.SliceThickness)
        kSpaceAcquisition.SetSliceOrigin(origin[2])
        if self.SliceProfile == 'ideal':
            kSpaceAcquisition.SetSliceProfileToIdeal()
        elif self.SliceProfile == 'trapezoidal':
            kSpaceAcquisition.SetSliceProfileToTrapezoidal()
        elif self.SliceProfile == 'quadratic':
            kSpaceAcquisition.SetSliceProfileToQuadratic()
        else:
            self.PrintError('Invalid sliceprofile')
        kSpaceAcquisition.Update()

        self.PrintLog("Max quadrature order used: %d" % kSpaceAcquisition.GetMaximumQuadratureOrderUsed())
        self.PrintLog("No. of Gauss point evaluations: %d" % kSpaceAcquisition.GetNumberOfGaussPointEvaluations())
        self.PrintLog("Avg no. of Gauss point evaluations per kSpace location: %d" % (kSpaceAcquisition.GetNumberOfGaussPointEvaluations()/kSpaceAcquisition.GetOutput().GetNumberOfPoints()))

        return kSpaceAcquisition.GetOutput()
        
    def Execute(self):

        if self.Mesh == None:
            self.PrintError('Error: no Mesh.')

        acquiredMesh = self.Mesh
        acquiredKSpace = None
        
        if (self.KSpaceDimensionality == 3) or not self.SliceSelection:

            origin = [self.FOVCenter[0] - self.FOV[0]/2.0,
                      self.FOVCenter[1] - self.FOV[1]/2.0,
                      self.FOVCenter[2] - self.FOV[2]/2.0]

            spacing = [self.FOV[0] / self.MatrixSize[0],
                       self.FOV[1] / self.MatrixSize[1],
                       self.FOV[2] / self.MatrixSize[2]]

            acquiredKSpace = self.AcquireKSpace(self.Mesh,origin,spacing)

        elif self.KSpaceDimensionality == 2:
 
            spacing = [self.FOV[0] / self.MatrixSize[0],
                       self.FOV[1] / self.MatrixSize[1],
                       self.FOV[2] / self.MatrixSize[2]]
           
            kSpaceAppend = vtk.vtkImageAppend()
            kSpaceAppend.SetAppendAxis(2)

            sliceLocations = []
            sliceLocation = self.FOVCenter[2] - self.FOV[2]/2.0
            while (sliceLocation <= self.FOVCenter[2] + self.FOV[2]/2.0):
                sliceLocations.append(sliceLocation)
                sliceLocation += self.SliceSpacing
            
            bounds = self.Mesh.GetBounds()
           
            for sliceLocation in sliceLocations:
       
                self.PrintLog("Processing slice at " + str(sliceLocation))

                origin = [self.FOVCenter[0] - self.FOV[0]/2.0,
                          self.FOVCenter[1] - self.FOV[1]/2.0,
                          sliceLocation]

                sliceKSpace = self.AcquireKSpace(self.Mesh,origin,spacing)
                kSpaceAppend.AddInput(sliceKSpace)
               
            kSpaceAppend.Update()
            
            acquiredKSpace = kSpaceAppend.GetOutput()

        self.KSpace = self.ComputeKSpaceOperation(acquiredKSpace)

if __name__=='__main__':
    main = pypes.pypeMain()
    main.Arguments = sys.argv
    main.Execute()
