#!/usr/bin/env python

## Program:   femrisymbolicmeshkspace.py
## Module:    $RCSfile: femrisymbolicmeshkspace.py,v $
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
import femrisymbolic
import femrikspace

from vmtk import pypes

femrisymbolicmeshkspace = 'femriSymbolicMeshKSpace'

class femriSymbolicMeshKSpace(femrikspace.femriKSpace):

    def __init__(self):

        femrikspace.femriKSpace.__init__(self)

        self.Mesh = None

        self.SliceThickness = 1.0
        self.MagnetizationFunction = '1'

        self.SliceModel = 0

        self.SetScriptName('femrisymbolicmeshkspace')
        self.SetInputMembers([
            ['Mesh','i','vtkUnstructuredGrid',1,'','','vmtkmeshreader'],
            ['SliceThickness','slicethickness','float',1,'(0.0,)'],
            ['MagnetizationFunction','magnetizationfunction','str',1],
            ['SliceModel','slicemodel','bool',1]
            ])
        self.SetOutputMembers([
            ])

##    def AngleBetweenNormals(self, normal0, normal1):
##
##        sum = [0.0, 0.0, 0.0]
##        sum[0] = normal0[0] + normal1[0];
##        sum[1] = normal0[1] + normal1[1];
##        sum[2] = normal0[2] + normal1[2];
##        sumNorm = vtk.vtkMath.Norm(sum);
##                                                                                                                                     
##        difference = [0.0, 0.0, 0.0]
##        difference[0] = normal0[0] - normal1[0];
##        difference[1] = normal0[1] - normal1[1];
##        difference[2] = normal0[2] - normal1[2];
##        differenceNorm = vtk.vtkMath.Norm(difference);
##        
##        return 2.0 * math.atan2(differenceNorm,sumNorm);

    def AcquireKSpace(self,mesh,origin,spacing):
        
        kSpaceAcquisition = femrisymbolic.vtkfemriUnstructuredGridKSpaceGenerator()
        kSpaceAcquisition.SetInput(mesh)
        kSpaceAcquisition.SetFunctionString(self.MagnetizationFunction)
        kSpaceAcquisition.SetKSpaceDimensionality(self.KSpaceDimensionality)
        kSpaceAcquisition.SetMatrix(self.MatrixSize)
        kSpaceAcquisition.SetFOV(self.FOV)
        kSpaceAcquisition.SetOrigin(origin)
        kSpaceAcquisition.UseElementIntegralCacheOn()
        kSpaceAcquisition.AcquireSymmetricKSpaceOn()
##         kSpaceAcquisition.UseElementIntegralCacheOff()
##         kSpaceAcquisition.AcquireSymmetricKSpaceOff()
        kSpaceAcquisition.Update()

        return kSpaceAcquisition.GetOutput()

    def Execute(self):

        if self.Mesh == None:
            self.PrintError('Error: no Mesh.')

##        if (self.FOVNormal != [0.0, 0.0, 1.0]):
##                    
##            translation = [-self.FOVCenter[0], -self.FOVCenter[1], -self.FOVCenter[2]]
##    
##            referenceNormal = [0.0, 0.0, 1.0]
##            rotationAxis = [0.0, 0.0, 0.0]
##            vtk.vtkMath.Normalize(self.FOVNormal)
##            vtk.vtkMath.Cross(self.FOVNormal,referenceNormal,rotationAxis)
##            angle = self.AngleBetweenNormals(referenceNormal,self.FOVNormal) / vtk.vtkMath.Pi() * 180.0
##        
##            transform = vtk.vtkTransform()
##            transform.PostMultiply()
##            transform.Translate(translation)
##            transform.RotateWXYZ(angle,rotationAxis)
##            transform.Translate(self.FOVCenter)
##    
##            transformFilter = vtk.vtkTransformFilter()
##            transformFilter.SetInput(self.Mesh)
##            transformFilter.SetTransform(transform)
##            transformFilter.Update()
##            
##            acquiredMesh = transformFilter.GetOutput()

        if (self.KSpaceDimensionality == 3) or not self.SliceModel:

            origin = [self.FOVCenter[0] - self.FOV[0]/2.0,
                      self.FOVCenter[1] - self.FOV[1]/2.0,
                      self.FOVCenter[2] - self.FOV[2]/2.0]

            spacing = [self.FOV[0] / self.MatrixSize[0],
                       self.FOV[1] / self.MatrixSize[1],
                       self.FOV[2] / self.MatrixSize[2]]

            self.KSpace = self.AcquireKSpace(self.Mesh,origin,spacing)

        elif self.KSpaceDimensionality == 2:
            
            kSpaceAppend = vtk.vtkImageAppend()
            kSpaceAppend.SetAppendAxis(2)

            sliceLocations = []
            sliceLocation = self.FOVCenter[2] - self.FOV[2]/2.0
            while (sliceLocation < self.FOVCenter[2] + self.FOV[2]/2.0):
                sliceLocations.append(sliceLocation)
                sliceLocation += self.SliceSpacing

            spacing = [self.FOV[0] / self.MatrixSize[0],
                       self.FOV[1] / self.MatrixSize[1],
                       self.FOV[2] / self.MatrixSize[2]]
            
            bounds = self.Mesh.GetBounds()
           
            for sliceLocation in sliceLocations:

                self.PrintLog("Processing slice at" + float(sliceLocation))

                origin = [self.FOVCenter[0] - self.FOV[0]/2.0,
                          self.FOVCenter[1] - self.FOV[1]/2.0,
                          sliceLocation]
                
                clipper1 = vtk.vtkClipDataSet()
                clipper1.SetInput(self.Mesh)
                clipper1.InsideOutOff()
                
                plane1 = vtk.vtkPlane()
                plane1.SetNormal(0.0,0.0,1.0)
                plane1.SetOrigin(0.0,0.0,sliceLocation - self.SliceThickness / 2.0)

                clipper1.SetClipFunction(plane1)
                clipper1.Update()

                clipper2 = vtk.vtkClipDataSet()
                clipper2.SetInput(clipper1.GetOutput())
                clipper2.InsideOutOn()
                
                plane2 = vtk.vtkPlane()
                plane2.SetNormal(0.0,0.0,1.0)
                plane2.SetOrigin(0.0,0.0,sliceLocation + self.SliceThickness / 2.0)

                clipper2.SetClipFunction(plane2)
                clipper2.Update()

                clipper2Bounds = clipper2.GetOutput().GetBounds()

                cleaner = vtk.vtkExtractUnstructuredGrid()
                cleaner.SetInput(clipper2.GetOutput())
                cleaner.ExtentClippingOn()
                cleaner.SetExtent(clipper2Bounds[0],clipper2Bounds[1],
                                  clipper2Bounds[2],clipper2Bounds[3],
                                  sliceLocation-self.SliceThickness/2.0,sliceLocation+self.SliceThickness/2.0)
                cleaner.Update()

                tetraFilter = vtk.vtkDataSetTriangleFilter()
                tetraFilter.SetInput(cleaner.GetOutput())
                tetraFilter.Update()

                sliceMesh = tetraFilter.GetOutput()

                self.PrintLog("Number of integration elements:" + int(sliceMesh.GetNumberOfCells()))

                sliceKSpace = self.AcquireKSpace(sliceMesh,origin,spacing)

                kSpaceAppend.AddInput(sliceKSpace)
                
            kSpaceAppend.Update()
            
            self.KSpace = self.ComputeKSpaceOperation(kSpaceAppend.GetOutput())


if __name__=='__main__':
    main = pypes.pypeMain()
    main.Arguments = sys.argv
    main.Execute()
