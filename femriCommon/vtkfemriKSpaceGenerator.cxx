/*=========================================================================

  Program:   FEMRI
  Module:    $RCSfile: vtkfemriKSpaceGenerator.cxx,v $
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

#include "vtkfemriKSpaceGenerator.h"
#include "vtkfemriKSpaceShift.h"
#include "vtkImageData.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkObjectFactory.h"

vtkCxxRevisionMacro(vtkfemriKSpaceGenerator, "$Revision: 1.2 $");

vtkfemriKSpaceGenerator::vtkfemriKSpaceGenerator()
{
  this->Matrix[0] = this->Matrix[1] = this->Matrix[2] = 0;
  this->FOV[0] = this->FOV[1] = this->FOV[2] = 0.0;
  this->Origin[0] = this->Origin[1] = this->Origin[2] = 0.0;
 
  this->KSpaceDimensionality = 3;
  this->AcquireSymmetricKSpace = 0;
//  this->NormalizeKSpace = 0;
  this->SetNumberOfInputPorts(0);
}

vtkfemriKSpaceGenerator::~vtkfemriKSpaceGenerator()
{
}

int vtkfemriKSpaceGenerator::RequestInformation (
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

double vtkfemriKSpaceGenerator::ComputeVoxelVolume()
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

int vtkfemriKSpaceGenerator::RequestData(
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
        this->EvaluateFourierFunction(frequency,fourierValue);
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
        newScalars->SetComponent(pointId,0,phaseShiftedFourierValue[0]);
        newScalars->SetComponent(pointId,1,phaseShiftedFourierValue[1]);
        if (10*pointId/numberOfPoints != percentage)
          {
          percentage = 10*pointId/numberOfPoints;
          cout<<"Progress "<<10*percentage<<"%"<<endl;
          }
        }
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

void vtkfemriKSpaceGenerator::PrintSelf(ostream& os, vtkIndent indent)
{
  Superclass::PrintSelf(os,indent);
}
