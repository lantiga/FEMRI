/*=========================================================================

  Program:   FEMRI
  Module:    $RCSfile: vtkfemriPolyDataExactKSpaceGeneratorMonolithic.cxx,v $
  Language:  C++
  Date:      $Date: 2008/11/03 17:00:30 $
  Version:   $Revision: 1.2 $

  Copyright (c) Luca Antiga, David Steinman. All rights reserved.
  See LICENCE file for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm 
  for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "vtkfemriPolyDataExactKSpaceGeneratorMonolithic.h"
#include "vtkfemriKSpaceShift.h"
#include "vtkImageData.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkObjectFactory.h"
#include "vtkPolyData.h"
#include "vtkCell.h"
#include "vtkCellData.h"
#include "vtkTriangle.h"

vtkStandardNewMacro(vtkfemriPolyDataExactKSpaceGeneratorMonolithic);
vtkCxxRevisionMacro(vtkfemriPolyDataExactKSpaceGeneratorMonolithic, "$Revision: 1.2 $");

vtkfemriPolyDataExactKSpaceGeneratorMonolithic::vtkfemriPolyDataExactKSpaceGeneratorMonolithic()
{
  this->Matrix[0] = this->Matrix[1] = this->Matrix[2] = 0;
  this->FOV[0] = this->FOV[1] = this->FOV[2] = 0.0;
  this->Origin[0] = this->Origin[1] = this->Origin[2] = 0.0;
 
  this->KSpaceDimensionality = 3;
  this->AcquireSymmetricKSpace = 0;
  this->MagnetizationValue = 1.0;
//  this->NormalizeKSpace = 0;
  this->SetNumberOfInputPorts(1);
}

vtkfemriPolyDataExactKSpaceGeneratorMonolithic::~vtkfemriPolyDataExactKSpaceGeneratorMonolithic()
{
}

int vtkfemriPolyDataExactKSpaceGeneratorMonolithic::FillInputPortInformation(int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
  return 1;
}

int vtkfemriPolyDataExactKSpaceGeneratorMonolithic::RequestInformation (
  vtkInformation * vtkNotUsed(request),
  vtkInformationVector** vtkNotUsed( inputVector ),
  vtkInformationVector *outputVector)
{
  vtkInformation* outInfo = outputVector->GetInformationObject(0);

  int wholeExtent[6];
  wholeExtent[0] = wholeExtent[2] = wholeExtent[4] = 0;
  wholeExtent[1] = this->Matrix[0] - 1;
  wholeExtent[3] = this->Matrix[1] - 1;
  wholeExtent[5] = this->Matrix[2] - 1;
  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),wholeExtent,6);
  
  double spacing[3];
  spacing[0] = 1.0 / this->FOV[0];
  spacing[1] = 1.0 / this->FOV[1];
  spacing[2] = 1.0 / this->FOV[2];
  outInfo->Set(vtkDataObject::SPACING(),spacing,3);
 
  outInfo->Set(vtkDataObject::ORIGIN(),this->Origin,3);

  vtkDataObject::SetPointDataActiveScalarInfo(outInfo, VTK_DOUBLE, 2);
  return 1;
}

double vtkfemriPolyDataExactKSpaceGeneratorMonolithic::ComputeVoxelVolume()
{
  double voxelVolume = 1.0;
  switch (this->KSpaceDimensionality)
    {
    case 2:
      voxelVolume = this->FOV[0]*this->FOV[1] / double(this->Matrix[0]*this->Matrix[1]);
      break;
    case 3:
      voxelVolume = this->FOV[0]*this->FOV[1]*this->FOV[2] / double(this->Matrix[0]*this->Matrix[1]*this->Matrix[2]);
      break;
    default:
      vtkErrorMacro(<< "Error: unsupported KSpace dimensionality");
      return 0.0;
    }
  return voxelVolume;
}

int vtkfemriPolyDataExactKSpaceGeneratorMonolithic::RequestData(
    vtkInformation* vtkNotUsed( request ),
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkImageData *output = vtkImageData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
 
  int updateExtent[6];
  outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(),updateExtent);
  output->SetExtent(updateExtent);
  output->AllocateScalars();

  vtkDoubleArray* newScalars = vtkDoubleArray::SafeDownCast(output->GetPointData()->GetScalars());
  newScalars->FillComponent(0,0.0);
  newScalars->FillComponent(1,0.0);
  int numberOfPoints = newScalars->GetNumberOfTuples();
  int percentage = 0;
  cout<<"Progress "<<10*percentage<<"%"<<endl;

  this->Initialize();

  int kSpaceDimensions[3];
  double kSpaceSpacing[3];
  int extent[6];
  
  output->GetWholeExtent(extent);
  output->GetSpacing(kSpaceSpacing);
  output->GetDimensions(kSpaceDimensions);

  int ijk[3];
  int kdx[3];
  double frequency[3];

  double volume = this->ComputeVolume();
  double voxelVolume = this->ComputeVoxelVolume();

  vtkPolyData* input = vtkPolyData::SafeDownCast(this->GetInput());

  vtkPoints* points = input->GetPoints();
  vtkDataArray* normals = input->GetCellData()->GetNormals();

  if (!normals)
    {
    vtkErrorMacro("No surface normals found. Please provide a surface with outwards normals defined on cells.");
    }

  double value[2];
  value[0] = 0.0;
  value[1] = 0.0;

  double pi = vtkMath::Pi();
  double twoPi = 2.0 * pi;
  double pi3 = pi * pi * pi;

  vtkIdType npts, *pts;

  int numberOfCells = input->GetNumberOfCells();

  input->BuildCells();

  int i, j;
  for (i=0; i<numberOfCells; i++)
    {
    input->GetCellPoints(i,npts,pts); 

    if (npts != 3)
      {
      continue;
      }

    double x1[3], x2[3], x3[3];
    points->GetPoint(pts[0],x1);
    points->GetPoint(pts[1],x2);
    points->GetPoint(pts[2],x3);

    double A = vtkTriangle::TriangleArea(x1,x2,x3);

    if (frequency[0] == 0.0 && frequency[1] == 0.0 && frequency[2] == 0.0)
      {
      value[0] += this->MagnetizationValue * A;
      value[1] += 0.0;
      continue;
      }
 
    //double normal[3];
    //vtkTriangle::ComputeNormal(x1,x2,x3,normal);

    double triangleNormal[3];
    normals->GetTuple(i,triangleNormal);

    //bool pointsReversed = false;

    //if (vtkMath::Dot(normal,triangleNormal) < 0.0) {
    //  pointsReversed = true;
    //  normal[0] *= -1.0;
    //  normal[1] *= -1.0;
    //  normal[2] *= -1.0;
    //}

    //double kdotnt = vtkMath::Dot(frequency,normal);

    //double cross[3];
    //vtkTriangle::ComputeNormal(x1,x2,x3,cross);

    //double normalDot = vtkMath::Dot(triangleNormal,cross);
    //bool pointsReversed = normalDot < 0 ? true : false;

    double kdotnt = vtkMath::Dot(frequency,triangleNormal);

    double cross[3];
    vtkTriangle::ComputeNormal(x1,x2,x3,cross);

    double normalDot = vtkMath::Dot(triangleNormal,cross);
    bool pointsReversed = normalDot < 0 ? true : false;
//

    if (pointsReversed)
      {
      points->GetPoint(pts[0],x1);
      points->GetPoint(pts[2],x2);
      points->GetPoint(pts[1],x3);
      }

    double x1x3[3], x2x3[3], x1x2[3];
    this->Vector(x1,x3,x1x3);
    this->Vector(x2,x3,x2x3);
    this->Vector(x1,x2,x1x2);

    //TODO: Acquire symmetric KSpace
    for (ijk[2]=extent[4]; ijk[2]<=extent[5]; ijk[2]++)
      {
      kdx[2] = ijk[2] <= kSpaceDimensions[2]/2 ? ijk[2] : ijk[2] - kSpaceDimensions[2];
      frequency[2] = double(kdx[2]) * kSpaceSpacing[2];
      for (ijk[1]=extent[2]; ijk[1]<=extent[3]; ijk[1]++)
        {
        kdx[1] = ijk[1] <= kSpaceDimensions[1]/2 ? ijk[1] : ijk[1] - kSpaceDimensions[1];
        frequency[1] = double(kdx[1]) * kSpaceSpacing[1];
        for (ijk[0]=extent[0]; ijk[0]<=extent[1]; ijk[0]++)
          {
          kdx[0] = ijk[0] <= kSpaceDimensions[0]/2 ? ijk[0] : ijk[0] - kSpaceDimensions[0];
          frequency[0] = double(kdx[0]) * kSpaceSpacing[0];
          double fourierValue[2];

          double k[3];
          k[0] = frequency[0];
          k[1] = frequency[1];
          k[2] = frequency[2];

          double k2 = this->Dot(k,k);
          double phi_e = this->Dot(k,x3);
          double k_r = this->Dot(k,x1x3);
          double k_s = this->Dot(k,x2x3);
          double k_rs = this->Dot(k,x1x2);
          //double k_dot_n = this->Dot(k,normal);
          double k_dot_n = this->Dot(k,triangleNormal);

          const double tol = 1E-8;

          if (fabs(k_s) < tol) 
            {
            if (fabs(k_r) < tol)
              {
              value[0] += k_dot_n / (twoPi * k2) * A * sin(twoPi * phi_e);
              value[1] += k_dot_n / (twoPi * k2) * A * cos(twoPi * phi_e);
              continue;
              }
            value[0] += k_dot_n / (twoPi * k2) * A * ( cos(twoPi*phi_e)/(pi*k_r) + sin(twoPi*phi_e)/(twoPi*twoPi*k_r*k_r) - cos(twoPi*k_r)*sin(twoPi*phi_e)/(twoPi*twoPi*k_r*k_r) - cos(twoPi*phi_e)*sin(twoPi*k_r)/(twoPi*twoPi*k_r*k_r));
            value[1] += k_dot_n / (twoPi * k2) * A * (-sin(twoPi*phi_e)/(pi*k_r) + cos(twoPi*phi_e)/(twoPi*twoPi*k_r*k_r) + sin(twoPi*k_r)*sin(twoPi*phi_e)/(twoPi*twoPi*k_r*k_r) - cos(twoPi*phi_e)*cos(twoPi*k_r)/(twoPi*twoPi*k_r*k_r));
            continue;
            }

          if (fabs(k_r) < tol) 
            {
            value[0] += k_dot_n / (twoPi * k2) * A * ( cos(twoPi*phi_e)/(pi*k_s) + sin(twoPi*phi_e)/(twoPi*twoPi*k_s*k_s) - cos(twoPi*k_s)*sin(twoPi*phi_e)/(twoPi*twoPi*k_s*k_s) - cos(twoPi*phi_e)*sin(twoPi*k_s)/(twoPi*twoPi*k_s*k_s));
            value[1] += k_dot_n / (twoPi * k2) * A * (-sin(twoPi*phi_e)/(pi*k_s) + cos(twoPi*phi_e)/(twoPi*twoPi*k_s*k_s) + sin(twoPi*k_s)*sin(twoPi*phi_e)/(twoPi*twoPi*k_s*k_s) - cos(twoPi*phi_e)*cos(twoPi*k_s)/(twoPi*twoPi*k_s*k_s));
            continue;
            }

          if (fabs(k_rs) < tol) 
            {
            value[0] += k_dot_n / (twoPi * k2) * A * (-sin(twoPi*phi_e)/(twoPi*twoPi*k_s*k_s) + sin(twoPi*k_s)*sin(twoPi*phi_e)/(pi*k_s) + cos(twoPi*k_s)*sin(twoPi*phi_e)/(twoPi*twoPi*k_s*k_s) + cos(twoPi*phi_e)*sin(twoPi*k_s)/(twoPi*twoPi*k_s*k_s) - cos(twoPi*k_s)*cos(twoPi*phi_e)/(pi*k_s));
            value[1] += k_dot_n / (twoPi * k2) * A * (-cos(twoPi*phi_e)/(twoPi*twoPi*k_s*k_s) + cos(twoPi*k_s)*sin(twoPi*phi_e)/(pi*k_s) - sin(twoPi*k_s)*sin(twoPi*phi_e)/(twoPi*twoPi*k_s*k_s) + cos(twoPi*phi_e)*cos(twoPi*k_s)/(twoPi*twoPi*k_s*k_s) + sin(twoPi*k_s)*cos(twoPi*phi_e)/(pi*k_s));
            continue;
            }

          value[0] = - 1.0 / (4.0 * pi3 * k2) * A * k_dot_n / k_s * ((sin(twoPi * (phi_e + k_r)) - sin(twoPi * (phi_e + k_s))) / k_rs - (sin(twoPi * (phi_e + k_r)) - sin(twoPi * (phi_e))) / k_r );
          value[1] = - 1.0 / (4.0 * pi3 * k2) * A * k_dot_n / k_s * ((cos(twoPi * (phi_e + k_r)) - cos(twoPi * (phi_e + k_s))) / k_rs - (cos(twoPi * (phi_e + k_r)) - cos(twoPi * (phi_e))) / k_r );

          fourierValue[0] = value[0];
          fourierValue[1] = value[1];

          fourierValue[0] *= volume / voxelVolume;
          fourierValue[1] *= volume / voxelVolume;
          double translation[3];
          translation[0] = -this->Origin[0];
          translation[1] = -this->Origin[1];
          translation[2] = 0.0; 
          if (this->KSpaceDimensionality == 3)
            {
            translation[2] = -this->Origin[2];
            }
          double phaseShiftedFourierValue[2];
          vtkfemriKSpaceShift::ShiftPhase(fourierValue,frequency,translation,phaseShiftedFourierValue);
          vtkIdType pointId = output->ComputePointId(ijk);

          double accumulatedFourierValue[2];
          newScalars->GetTuple(pointId,accumulatedFourierValue);
          accumulatedFourierValue[0] += phaseShiftedFourierValue[0];
          accumulatedFourierValue[1] += phaseShiftedFourierValue[1];

          newScalars->SetTuple(pointId,accumulatedFourierValue);
          }
        }
      }

    if (10*i/numberOfCells != percentage)
      {
      percentage = 10*i/numberOfCells;
      cout<<"Progress "<<10*percentage<<"%"<<endl;
      }
    }
#if 0
  if (this->NormalizeKSpace)
    {
    double value[2];
    double magnitudeSum = 0;
    int numberOfSamples = newScalars->GetNumberOfTuples();
    for (int i=0; i<numberOfSamples; i++)
      {
      newScalars->GetTuple(i,value);
      magnitudeSum += sqrt(value[0]*value[0] + value[1]*value[1]);
      }
    for (int i=0; i<numberOfSamples; i++)
      {
      newScalars->GetTuple(i,value);
      value[0] /= magnitudeSum;
      value[1] /= magnitudeSum;
      newScalars->SetTuple(i,value);
      }
    }
#endif
 
  return 1;             
}

void vtkfemriPolyDataExactKSpaceGeneratorMonolithic::PrintSelf(ostream& os, vtkIndent indent)
{
  Superclass::PrintSelf(os,indent);
}
