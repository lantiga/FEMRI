#!/usr/bin/env python

## Program:   femrithicknessmeasurement.py
## Module:    $RCSfile: femrithicknessmeasurement.py,v $
## Language:  Python
## Date:      $Date: 2008/10/28 16:46:59 $
## Version:   $Revision: 1.3 $

##   Copyright (c) Luca Antiga, David Steinman. All rights reserved.
##   See LICENCE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

import sys
import math
import vtk
from vmtk import vtkvmtk

from vmtk import pypes

femrithicknessmeasurement = 'femriThicknessMeasurement'


class femri2DThickness(object):

    def __init__(self):
    
        self.Contour = None
        self.Center = [0.0, 0.0, 0.0]
        self.TiltingAngle = 0.0
        self.RotationAngle = 0.0
        self.Locations = [0.0, 0.0]
        self.Thickness = 0.0
        self.Thickness3D = 0.0

    def LineFunction(self,point,theta,origin=[0.0,0.0]):

        normal = [math.cos(theta),math.sin(theta)]
        position = [point[0]-origin[0],point[1]-origin[1]]
        return normal[0]*position[0]+normal[1]*position[1]
 
    def FindLocations(self,points):
        pass
        
    def Execute(self):

        plane = vtk.vtkPlane()
        plane.SetOrigin(self.Center)
        plane.SetNormal(-math.sin(self.RotationAngle),math.cos(self.RotationAngle),0.0)
        cutter = vtk.vtkCutter()
        cutter.SetInput(self.Contour)
        cutter.SetCutFunction(plane)
        cutter.SetValue(0,0.0)
        cutter.Update()
        cutPoints = cutter.GetOutput().GetPoints()
        if not cutPoints or cutPoints.GetNumberOfPoints()==0:
            self.Locations = [0.0, 0.0]
            self.Thickness = 0.0
            self.Thickness3D = 0.0
            return
        self.Locations = self.FindLocations(cutPoints)
        self.Thickness = self.Locations[1] - self.Locations[0]
        self.Thickness3D = self.Thickness * math.cos(self.TiltingAngle)

  
class femri2DPlaneThickness(femri2DThickness):
    
    def __init__(self):
          femri2DThickness.__init__(self)

    def FindLocations(self,points):
        lineFunctions = []
        for i in range(points.GetNumberOfPoints()):
            point = points.GetPoint(i)
            lineFunctions.append(self.LineFunction(point,self.RotationAngle,self.Center))
        if len(lineFunctions) < 2:
            return [0.0,0.0]
        lineFunctions.sort()
        return [lineFunctions[0], lineFunctions[-1]]

        
class femri2DCylinderThickness(femri2DThickness):
    
    def __init__(self):
          femri2DThickness.__init__(self)

    def FindLocations(self,points):
        lineFunctions = []
        for i in range(points.GetNumberOfPoints()):
            point = points.GetPoint(i)
            lineFunctions.append(self.LineFunction(point,self.RotationAngle,self.Center))
        if len(lineFunctions) < 2:
            return [0.0,0.0]
        lineFunctions.sort()
        return [lineFunctions[0], lineFunctions[-1]]


class femri2DHollowCylinderThickness(femri2DThickness):
    
    def __init__(self):
          femri2DThickness.__init__(self)
          self.InnerArea = 0.0
          self.OuterArea = 0.0

    def FindLocations(self,points):
        lineFunctions = []
        for i in range(points.GetNumberOfPoints()):
            point = points.GetPoint(i)
            lineFunctions.append(self.LineFunction(point,self.RotationAngle,self.Center))
        lineFunctions = [value for value in lineFunctions if value > 0.0]
        if len(lineFunctions) < 2:
            return [0.0,0.0]
        lineFunctions.sort()
        return [lineFunctions[0], lineFunctions[-1]]

    def ComputeAreas(self):
        stripper = vtk.vtkStripper()
        stripper.SetInput(self.Contour)
        stripper.SetMaximumLength(100000)
        stripper.Update()
        contour = stripper.GetOutput()
        if contour.GetNumberOfCells() != 2:
            self.InnerArea = 0.0
            self.OuterArea = 0.0
            return
        areas = [0.0,0.0]
        for n in range(2):
            points = contour.GetCell(n).GetPoints()
            numberOfPoints = points.GetNumberOfPoints()
            for i in range(numberOfPoints):
                point0 = points.GetPoint(i)
                point1 = points.GetPoint((i+1)%numberOfPoints)
                areas[n] += vtk.vtkTriangle.TriangleArea(self.Center,point0,point1)
        self.InnerArea = min(areas)
        self.OuterArea = max(areas)

       
class femri3DCylinderThickness(object):

    def __init__(self):
    
        self.Contour = None
        self.Center = [0.0, 0.0, 0.0]
        self.TiltingAngle = 0.0
        self.RotationAngle = 0.0
        self.Locations = [0.0, 0.0]
        self.Thickness = 0.0
        self.Thickness3D = 0.0

    def LineFunction(self,point,theta,origin=[0.0,0.0]):

        normal = [math.cos(theta),math.sin(theta)]
        position = [point[0]-origin[0],point[1]-origin[1]]
        return normal[0]*position[0]+normal[1]*position[1]
 
    def FindLocations(self,points):
        return [0.0,0.0]
        
    def ComputeMeanRadius(self,line):
        numberOfPoints = line.GetNumberOfPoints()
        if numberOfPoints == 0:
            return 0.0
        meanRadius = 0.0
        weightSum = 0.0
        for pointId in range(numberOfPoints):
            point = line.GetPoint(pointId)
            weight =  0.5 * math.sqrt(vtk.vtkMath.Distance2BetweenPoints(point,line.GetPoint((pointId+numberOfPoints-1)%numberOfPoints))) + \
                      0.5 * math.sqrt(vtk.vtkMath.Distance2BetweenPoints(point,line.GetPoint((pointId+1)%numberOfPoints)))
            meanRadius += weight * math.sqrt(vtk.vtkMath.Distance2BetweenPoints(self.Center,point))
            weightSum += weight
        meanRadius /= weightSum
        return meanRadius
  
    def Execute(self):

        normal = [math.cos(self.RotationAngle) * math.sin(self.TiltingAngle), 
                  math.sin(self.RotationAngle) * math.sin(self.TiltingAngle), 
                  math.cos(self.TiltingAngle)]
        
        plane = vtk.vtkPlane()
        plane.SetOrigin(self.Center)
        plane.SetNormal(normal)
        cutter = vtk.vtkCutter()
        cutter.SetInput(self.Contour)
        cutter.SetCutFunction(plane)
        cutter.SetValue(0,0.0)
        cutter.Update()
        cleaner = vtk.vtkCleanPolyData()
        cleaner.SetInput(cutter.GetOutput())
        cleaner.Update()
        stripper = vtk.vtkStripper()
        stripper.SetInput(cleaner.GetOutput())
        stripper.Update()
        cutLine = stripper.GetOutput()
        self.Locations = [0.0, 0.0]
        self.Thickness = 2.0 * self.ComputeMeanRadius(cutLine)
        self.Thickness3D = self.Thickness
#        self.Contour = cutLine


class femriThicknessMeasurement(pypes.pypeScript):

    def __init__(self):

        pypes.pypeScript.__init__(self)

        self.Image = None
        self.LevelSets = None
        self.FeatureImage = None
        self.Contour = None
  
        self.Shape = 'thickplane2d'

        self.GaussFiltering = 0
        self.StandardDeviations = [1.0,1.0,1.0]
        self.RadiusFactors = [10.0,10.0,10.0]
  
        self.Center = [0.0,0.0,0.0]
        self.RotationAngle = 0.0
        self.TiltingAngle = 0.0

        self.Method = 'levelsets'

        self.FeatureImageType = 'gradient'
        self.UpwindFactor = 1.0

        self.CurvatureScaling = 0.0
        self.LevelSetsIterations = 1000

        self.FWHMBackground = None
        self.FWHMRatio = 0.5
        self.FWHMRegion = 'image'
        self.FWHMLevel = None

        self.Smoothing = 0
        self.SmoothingIterations = 30
        self.SmoothingPassBand = 0.1
        
        self.Thickness = 0.0
        self.Thickness3D = 0.0
        self.Locations = [0.0,0.0]

        self.SetScriptName('femrithicknessmeasurement')
        self.SetInputMembers([
            ['Image','i','vtkImageData',1,'','','vmtkimagereader'],
            ['Shape','shape','str',1],
            ['Center','center','float',3],
            ['RotationAngle','rotation','float',1],
            ['TiltingAngle','tilting','float',1],
            ['Method','method','str',1,'["levelsets","fwhm"]'],
            ['FeatureImageType','featureimagetype','str',1],
            ['UpwindFactor','upwindfactor','float',1,'(0.0,1.0)'],
            ['GaussFiltering','gauss','bool',1],
            ['StandardDeviations','sigmas','float',3,'(0.0,)'],
            ['RadiusFactors','radiusfactors','float',3,'(0.0,)'],
            ['LevelSetsIterations','iterations','int',1,'(0,)'],
            ['CurvatureScaling','curvature','float',1],
            ['FWHMBackground','fwhmbackground','float',1],
            ['FWHMRatio','fwhmratio','float',1],
            ['FWHMRegion','fwhmregion','str',1,'["image","midline"]'],
            ['Smoothing','smoothing','bool',1],
            ['SmoothingIterations','smoothingiterations','int',1,'(0,)'],
            ['SmoothingPassBand','smoothingpassband','float',1,'(0.0,)']
            ])
        self.SetOutputMembers([
            ['LevelSets','levelsets','vtkImageData',1,'','','vmtkimagewriter'],
            ['FeatureImage','featureimage','vtkImageData',1,'','','vmtkimagewriter'],
            ['Contour','contour','vtkPolyData',1,'','','vmtksurfacewriter'],
            ['Thickness','thickness','float',1],
            ['Thickness3D','thickness3d','float',1],
            ['Locations','locations','float',2],
            ['FWHMLevel','fwhmlevel','float',1]
            ])

    def Execute(self):

        if self.GaussFiltering:
            gauss = vtk.vtkImageGaussianSmooth()
            gauss.SetInput(self.Image)
            gauss.SetStandardDeviations(self.StandardDeviations)
            gauss.SetRadiusFactors(self.RadiusFactors)
            if self.Shape.find('2d') != -1:
                gauss.SetDimensionality(2)
            elif self.Shape.find('3d') != -1:
                gauss.SetDimensionality(3)
            else:
                gauss.SetDimensionality(3)
            gauss.Update()
            self.Image = gauss.GetOutput()

        scalarRange = [0.0,0.0]
        if self.FWHMRegion == 'image':
            scalarRange = self.Image.GetScalarRange()
        elif self.FWHMRegion == 'midline':
            extent = self.Image.GetWholeExtent()
            newYExtent = extent[2]+(extent[3]-extent[2])/2
            clip = vtk.vtkImageClip()
            clip.SetInput(self.Image)
            clip.SetOutputWholeExtent(extent[0],extent[1],newYExtent,newYExtent,extent[4],extent[5])
            clip.ClipDataOn()
            clip.Update()
            scalarRange = clip.GetOutput().GetScalarRange()

        self.FWHMLevel = (scalarRange[1] - scalarRange[0]) * self.FWHMRatio + scalarRange[0]
        if self.FWHMBackground != None:
            self.FWHMLevel = (scalarRange[1] - self.FWHMBackground) * self.FWHMRatio + self.FWHMBackground

        if self.Method == 'levelsets':
            if self.Shape.find('2d') != -1:
                gradientMagnitude = vtk.vtkImageGradientMagnitude()
                gradientMagnitude.SetDimensionality(2)
            elif self.Shape.find('3d') != -1:
                if self.FeatureImageType == 'gradient':
                    gradientMagnitude = vtkvmtk.vtkvmtkGradientMagnitudeImageFilter()
                elif self.FeatureImageType == 'upwind':
                    gradientMagnitude = vtkvmtk.vtkvmtkUpwindGradientMagnitudeImageFilter()
                    gradientMagnitude.SetUpwindFactor(self.UpwindFactor)
                else:
                    self.PrintError('Unsupported feature image type: choices are "gradient", "upwind".')
                    return
            else:
                gradientMagnitude = vtk.vtkImageGradientMagnitude()
            gradientMagnitude.SetInput(self.Image)
            gradientMagnitude.Update()
            boundedReciprocal = vtkvmtk.vtkvmtkBoundedReciprocalImageFilter()
            boundedReciprocal.SetInput(gradientMagnitude.GetOutput())
            boundedReciprocal.Update()
            self.FeatureImage = vtk.vtkImageData()
            self.FeatureImage.DeepCopy(boundedReciprocal.GetOutput())
            levelSetsFilter = None
            if self.Shape.find('2d') != -1:
                levelSetsFilter = vtkvmtk.vtkvmtkGeodesicActiveContourLevelSet2DImageFilter()
            elif self.Shape.find('3d') != -1:
                levelSetsFilter = vtkvmtk.vtkvmtkGeodesicActiveContourLevelSetImageFilter()
            else:
                levelSetsFilter = vtkvmtk.vtkvmtkGeodesicActiveContourLevelSetImageFilter()
            levelSetsFilter.SetInput(self.Image)
            levelSetsFilter.SetFeatureImage(boundedReciprocal.GetOutput())
            levelSetsFilter.SetDerivativeSigma(0.0)
            levelSetsFilter.SetAutoGenerateSpeedAdvection(1)
            levelSetsFilter.SetNumberOfIterations(self.LevelSetsIterations)
            levelSetsFilter.SetPropagationScaling(0.0)
            levelSetsFilter.SetCurvatureScaling(self.CurvatureScaling)
            levelSetsFilter.SetAdvectionScaling(1.0)
            levelSetsFilter.SetIsoSurfaceValue(self.FWHMLevel)
            levelSetsFilter.SetInterpolateSurfaceLocation(1)
            levelSetsFilter.SetMaximumRMSError(1E-10)
            levelSetsFilter.Update()

            self.LevelSets = vtk.vtkImageData()
            self.LevelSets.DeepCopy(levelSetsFilter.GetOutput())

            contourFilter = vtk.vtkMarchingContourFilter()
            contourFilter.SetInput(self.LevelSets)
            contourFilter.SetValue(0,0.0)
            contourFilter.Update()

            self.Contour = contourFilter.GetOutput()

        elif self.Method == 'fwhm':
            contourFilter = vtk.vtkMarchingContourFilter()
            contourFilter.SetInput(self.Image)
            contourFilter.SetValue(0,self.FWHMLevel)
            contourFilter.Update()

            self.Contour = contourFilter.GetOutput()
        else:
            self.PrintError('Unsupported method: choices are "levelsets", "fwhm".')
            return


        if self.Smoothing:
            smoothingFilter = vtk.vtkWindowedSincPolyDataFilter()
            smoothingFilter.SetInput(self.Contour)
            smoothingFilter.SetNumberOfIterations(self.SmoothingIterations)
            smoothingFilter.SetPassBand(self.SmoothingPassBand)
            smoothingFilter.Update()
            self.Contour = smoothingFilter.GetOutput()

        measurementFilter = None

        if self.Shape is 'thickplane2d':
            measurementFilter = femri2DPlaneThickness()
        elif self.Shape is 'cylinder2d':
            measurementFilter = femri2DCylinderThickness()
        elif self.Shape is 'hollowcylinder2d':
            measurementFilter = femri2DHollowCylinderThickness()
        elif self.Shape is 'cylinder3d':
            measurementFilter = femri3DCylinderThickness()
        else:
            self.PrintError('Unsupported shape: choices are "thickplane2d", "cylinder2d", "hollowcylinder2d", "cylinder3d".')
            return

        measurementFilter.Contour = self.Contour
        measurementFilter.Center = self.Center
        measurementFilter.TiltingAngle = math.radians(self.TiltingAngle)
        measurementFilter.RotationAngle = math.radians(self.RotationAngle)
        measurementFilter.Execute()

        if self.Shape is 'hollowcylinder2d':
            measurementFilter.ComputeAreas()
            self.InnerArea = measurementFilter.InnerArea
            self.OuterArea = measurementFilter.OuterArea

        self.Thickness = measurementFilter.Thickness
        self.Thickness3D = measurementFilter.Thickness3D
        self.Locations = measurementFilter.Locations

        self.Contour = measurementFilter.Contour
  
        self.OutputText('\n')
        self.OutputText('Thickness: ' + str(self.Thickness) + '\n')
        self.OutputText('Thickness3D: ' + str(self.Thickness3D) + '\n')
        self.OutputText('Locations: ' + str(self.Locations) + '\n')
        if self.Shape is 'hollowcylinder2d':
            self.OutputText('InnerArea: ' + str(self.InnerArea) + '\n')
            self.OutputText('OuterArea: ' + str(self.OuterArea) + '\n')
        self.OutputText('\n')


if __name__=='__main__':
    main = pypes.pypeMain()
    main.Arguments = sys.argv
    main.Execute()

