/*=========================================================================

  Program:   femri
  Module:    $RCSfile: vtkfemriUnstructuredGridKSpaceGenerator.cxx,v $
  Language:  C++
  Date:      $Date: 2007/03/19 13:31:27 $
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

#include "vtkfemriUnstructuredGridKSpaceGenerator.h"
#include "vtkfemriUnstructuredGridFourierIntegrator.h"
#include "vtkUnstructuredGrid.h"
#include "vtkImageData.h"
#include "vtkPointData.h"
#include "vtkDoubleArray.h"
#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkObjectFactory.h"

#include "ginac.h"

using namespace GiNaC;
using namespace std;

vtkStandardNewMacro(vtkfemriUnstructuredGridKSpaceGenerator);
vtkCxxRevisionMacro(vtkfemriUnstructuredGridKSpaceGenerator, "$Revision: 1.1.1.1 $");

vtkfemriUnstructuredGridKSpaceGenerator::vtkfemriUnstructuredGridKSpaceGenerator()
{
  this->FunctionString = NULL;

  this->FOV[0] = this->FOV[1] = this->FOV[2] = 0.0;
  this->Origin[0] = this->Origin[1] = this->Origin[2] = 0.0;

  this->KSpaceDimensionality = 3;

  this->Matrix[0] = 16;
  this->Matrix[1] = 16;
  this->Matrix[2] = 16;

//   this->UsePreComputedSymbolicIntegral = 0;
  this->UseElementIntegralCache = 0;

  this->AcquireSymmetricKSpace = 0;

  Digits=40;
}

vtkfemriUnstructuredGridKSpaceGenerator::~vtkfemriUnstructuredGridKSpaceGenerator()
{
  if (this->FunctionString)
    {
    delete[] this->FunctionString;
    this->FunctionString = NULL;
    }
}

int vtkfemriUnstructuredGridKSpaceGenerator::FillInputPortInformation(
  int vtkNotUsed( port ), vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
  return 1;
}

int vtkfemriUnstructuredGridKSpaceGenerator::RequestInformation (
  vtkInformation * vtkNotUsed(request),
  vtkInformationVector ** vtkNotUsed( inputVector ),
  vtkInformationVector *outputVector)
{
  vtkInformation* outInfo = outputVector->GetInformationObject(0);

  vtkIdType wholeExtent[6];
  wholeExtent[0] = wholeExtent[2] = wholeExtent[4] = 0;
  wholeExtent[1] = this->Matrix[0]-1;
  wholeExtent[3] = this->Matrix[1]-1;
  wholeExtent[5] = this->Matrix[2]-1;
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

int vtkfemriUnstructuredGridKSpaceGenerator::RequestData(
  vtkInformation* vtkNotUsed( request ),
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
  vtkUnstructuredGrid *input = vtkUnstructuredGrid::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkImageData *output = vtkImageData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  ex K = symbolic_matrix(3,1,"k");
  ex X = symbolic_matrix(3,1,"x");

  try
    {
    this->Function = ex(this->FunctionString,lst(X[0],X[1],X[2],K[0],K[1],K[2]));
    //TODO: allow to select encoding axes
    switch (this->KSpaceDimensionality)
      {
      case 1:
        this->Function *= ex("exp(-2*Pi*I*(k0*x0))",lst(X[0],X[1],X[2],K[0],K[1],K[2]));
        break;
      case 2:
        this->Function *= ex("exp(-2*Pi*I*(k0*x0+k1*x1))",lst(X[0],X[1],X[2],K[0],K[1],K[2]));
        break;
      case 3:
        this->Function *= ex("exp(-2*Pi*I*(k0*x0+k1*x1+k2*x2))",lst(X[0],X[1],X[2],K[0],K[1],K[2]));
        break;
      default:
        cerr << "Unsupported dimensionality." << endl;
        break;
      }
    cout<<"Integrand: "<<this->Function<<endl;
    }
  catch (exception& p)
    {
    cerr << p.what() << endl;
    }

  int updateExtent[6];
  outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(),updateExtent);
  
  output->SetExtent(updateExtent);
  output->AllocateScalars();

  vtkDoubleArray* newScalars;
  newScalars = vtkDoubleArray::SafeDownCast(output->GetPointData()->GetScalars());

  if (input->GetNumberOfCells()==0)
    {
    newScalars->FillComponent(0,0.0);
    newScalars->FillComponent(1,0.0);
    return 1;
    }

  vtkfemriUnstructuredGridFourierIntegrator* meshIntegrator = vtkfemriUnstructuredGridFourierIntegrator::New();
  meshIntegrator->SetMesh(input);
  meshIntegrator->SetFunction(this->Function);
  meshIntegrator->SetX(X);
  meshIntegrator->SetK(K);

  vtkIdType extent[6];
  vtkIdType kSpaceDimensions[3];
  double kSpaceSpacing[3];
  double frequency[3];
  vtkIdType ijk[3];
  vtkIdType pointId;
  vtkIdType kdx[3];

  output->GetWholeExtent(extent);
  output->GetSpacing(kSpaceSpacing);
  output->GetDimensions(kSpaceDimensions);

  double phaseShift[3];
  phaseShift[0] = - this->Origin[0];
  phaseShift[1] = - this->Origin[1];
  phaseShift[2] = - this->Origin[2];

  double fourierIntegralRe, fourierIntegralIm;

//   if (this->UsePreComputedSymbolicIntegral)
//     {
//     meshIntegrator->PreComputeSymbolicIntegral();
//     meshIntegrator->SetUsePreComputedSymbolicIntegral(this->UsePreComputedSymbolicIntegral);
//     }

  if (this->UseElementIntegralCache)
    {
    meshIntegrator->BuildElementIntegralCache();
    meshIntegrator->SetUseElementIntegralCache(this->UseElementIntegralCache);
    }

  time_t start, end;
  time(&start);
  for (ijk[2]=extent[4]; ijk[2]<=extent[5]; ijk[2]++)
    {
    kdx[2] = ijk[2] <= kSpaceDimensions[2]/2 ? ijk[2] : ijk[2] - kSpaceDimensions[2];
    frequency[2] = double(kdx[2]) * kSpaceSpacing[2];

    // TODO: partially evaluate integral with frequency[2] (if not symmetry)

    for (ijk[1]=extent[2]; ijk[1]<=extent[3]; ijk[1]++)
      {
      kdx[1] = ijk[1] <= kSpaceDimensions[1]/2 ? ijk[1] : ijk[1] - kSpaceDimensions[1];
      frequency[1] = double(kdx[1]) * kSpaceSpacing[1];

      // TODO: partially evaluate integral with frequency[1] (if not symmetry)

      for (ijk[0]=extent[0]; ijk[0]<=extent[1]; ijk[0]++)
        {
        kdx[0] = ijk[0] <= kSpaceDimensions[0]/2 ? ijk[0] : ijk[0] - kSpaceDimensions[0];
        frequency[0] = double(kdx[0]) * kSpaceSpacing[0];

        if ((this->AcquireSymmetricKSpace) && (this->KSpaceDimensionality == 2) && (ijk[1]> kSpaceDimensions[1]/2))
          {
          vtkIdType sijk[3];
          sijk[0] = ijk[0] == 0 ? 0 : kSpaceDimensions[0] - ijk[0];
          sijk[1] = ijk[1] == 0 ? 0 : kSpaceDimensions[1] - ijk[1];
          sijk[2] = ijk[2];
          vtkIdType symmetricPointId = output->ComputePointId(sijk);
          double symmFourierIntegralRe = newScalars->GetComponent(symmetricPointId,0);
          double symmFourierIntegralIm = newScalars->GetComponent(symmetricPointId,1);
          fourierIntegralRe = symmFourierIntegralRe;
          fourierIntegralIm = -symmFourierIntegralIm;
          }
        else if ((this->AcquireSymmetricKSpace) && (this->KSpaceDimensionality == 3) && (ijk[2] > kSpaceDimensions[2]/2))
          {
          vtkIdType sijk[3];
          sijk[0] = ijk[0] == 0 ? 0 : kSpaceDimensions[0] - ijk[0];
          sijk[1] = ijk[1] == 0 ? 0 : kSpaceDimensions[1] - ijk[1];
          sijk[2] = ijk[2] == 0 ? 0 : kSpaceDimensions[2] - ijk[2];
          vtkIdType symmetricPointId = output->ComputePointId(sijk);
          double symmFourierIntegralRe = newScalars->GetComponent(symmetricPointId,0);
          double symmFourierIntegralIm = newScalars->GetComponent(symmetricPointId,1);
          fourierIntegralRe = symmFourierIntegralRe;
          fourierIntegralIm = -symmFourierIntegralIm;
          }
        else
          {
          // TODO: only substitute frequency[0] (if not symmetry)

          numeric fourierIntegral = meshIntegrator->ComputeIntegral(lst(K[0]==frequency[0], K[1]==frequency[1], K[2]==frequency[2]));
          numeric phaseShiftPhase = - 2 * ex_to<numeric>(PiEvalf()) * (frequency[0] * phaseShift[0] + frequency[1] * phaseShift[1] + frequency[2] * phaseShift[2]);
          numeric phaseShiftComplex = cos(phaseShiftPhase) + I * sin(phaseShiftPhase);
          numeric phaseShiftedFourierIntegral = fourierIntegral * phaseShiftComplex;
          fourierIntegralRe = phaseShiftedFourierIntegral.real().to_double();
          fourierIntegralIm = phaseShiftedFourierIntegral.imag().to_double();
          }

        pointId = output->ComputePointId(ijk);

        newScalars->SetComponent(pointId,0,fourierIntegralRe);
        newScalars->SetComponent(pointId,1,fourierIntegralIm);
        if (frequency[0] == 0.0 && frequency[1] == 0.0 && frequency[2] == 0.0)
          {
          cout << "DC: " << fourierIntegralRe << " + " << fourierIntegralIm << "*I" << endl;
          }
	else
	  {
          cout<<"Row: "<<ijk[1]<<" Column: "<<ijk[0]<<" Frequency: ("<<frequency[1]<<", "<<frequency[0]<<") "<<endl;
	  }
        }
      cout<<"Row: "<<ijk[1]<<" Frequency: "<<frequency[1]<<endl;
      }
    }
  time(&end);
  double elapsedTime;
  elapsedTime = difftime(end,start);
  cout<<"Computation time: "<<elapsedTime<<endl;
  cout<<"NumberOfIntegralComputations: "<< meshIntegrator->GetNumberOfComputations() << endl;
  cout<<"NumberOfIntegralEvaluations: "<< meshIntegrator->GetNumberOfEvaluations() << endl;

  meshIntegrator->Delete();

  return 1;
}
