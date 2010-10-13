/*=========================================================================

  Program:   FEMRI
  Module:    $RCSfile: vtkfemriKSpaceShift.cxx,v $
  Language:  C++
  Date:      $Date: 2007/03/19 13:31:25 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Luca Antiga, David Steinman. All rights reserved.
  See LICENCE file for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm 
  for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "vtkfemriKSpaceShift.h"
#include "vtkImageData.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkObjectFactory.h"

vtkStandardNewMacro(vtkfemriKSpaceShift);
vtkCxxRevisionMacro(vtkfemriKSpaceShift, "$Revision: 1.1.1.1 $");

vtkfemriKSpaceShift::vtkfemriKSpaceShift()
{
  this->Translation[0] = this->Translation[1] = this->Translation[2] = 0.0;
  this->SetNumberOfInputPorts(1);
  //TODO: force working with number of components = 2
}

vtkfemriKSpaceShift::~vtkfemriKSpaceShift()
{
}


int vtkfemriKSpaceShift::RequestData(
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
  
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkImageData *input = vtkImageData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkDataArray* inScalars = input->GetPointData()->GetScalars();
  
  int kSpaceDimensions[3];
  double kSpaceSpacing[3];
  int extent[6];
  
  output->GetWholeExtent(extent);
  output->GetSpacing(kSpaceSpacing);
  output->GetDimensions(kSpaceDimensions);

  int ijk[3];
  int kdx[3];
  double frequency[3];
 
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
        vtkIdType pointId = input->ComputePointId(ijk);
        double fourierValue[2];
        fourierValue[0] = inScalars->GetComponent(pointId,0);
        fourierValue[1] = inScalars->GetComponent(pointId,1);
        double phaseShiftedFourierValue[2];
        this->ShiftPhase(fourierValue,frequency,this->Translation,phaseShiftedFourierValue);
        newScalars->SetComponent(pointId,0,phaseShiftedFourierValue[0]);
        newScalars->SetComponent(pointId,1,phaseShiftedFourierValue[1]);
        }
      }
    }
  
  return 1;             
}

void vtkfemriKSpaceShift::ShiftPhase(double value[2], double frequency[3], double translation[3], double shiftedValue[2])
{
  double phaseShift = -2.0 * vtkMath::Pi() * (frequency[0] * translation[0] + frequency[1] * translation[1] + frequency[2] * translation[2]);
  double shiftComplex[2];
  shiftComplex[0] = cos(phaseShift);
  shiftComplex[1] = sin(phaseShift);
  shiftedValue[0] = value[0] * shiftComplex[0] - value[1] * shiftComplex[1];
  shiftedValue[1] = value[0] * shiftComplex[1] + value[1] * shiftComplex[0];
}

void vtkfemriKSpaceShift::PrintSelf(ostream& os, vtkIndent indent)
{
  Superclass::PrintSelf(os,indent);
}
