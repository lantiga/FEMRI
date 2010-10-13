/*=========================================================================

  Program:   femri
  Module:    $RCSfile: vtkfemriPolyDataNumericKSpaceGenerator.cxx,v $
  Language:  C++
  Date:      $Date: 2008/11/04 14:38:31 $
  Version:   $Revision: 1.7 $

  Copyright (c) Luca Antiga, David Steinman. All rights reserved.
  See LICENCE file for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm 
  for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "vtkfemriPolyDataNumericKSpaceGenerator.h"
#include "vtkfemriOptimalQuadratureOrderCalculator.h"
#include "vtkPolyData.h"
#include "vtkCell.h"
#include "vtkPointData.h"
#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkObjectFactory.h"
#include "vtkfemriGaussQuadrature.h"
#include "vtkTriangle.h"
#include "vtkQuadraticTriangle.h"
#include "vtkQuad.h"
#include "vtkQuadraticQuad.h"
#include "vtkPolyDataNormals.h"

#include "vtkImageData.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkObjectFactory.h"


vtkStandardNewMacro(vtkfemriPolyDataNumericKSpaceGenerator);
vtkCxxRevisionMacro(vtkfemriPolyDataNumericKSpaceGenerator, "$Revision: 1.7 $");

vtkfemriPolyDataNumericKSpaceGenerator::vtkfemriPolyDataNumericKSpaceGenerator()
{
  this->OptimalQuadratureOrderCalculator = NULL;

  this->UseOptimalAlgorithm = 0;
  this->ErrorThreshold = 1E-4;

  this->MagnetizationValue = 1.0;
  this->QuadratureOrder = 1;

  this->NumberOfGaussPointEvaluations = 0;
  this->MaximumQuadratureOrderUsed = 0;

  this->SetNumberOfInputPorts(1);
}

vtkfemriPolyDataNumericKSpaceGenerator::~vtkfemriPolyDataNumericKSpaceGenerator()
{
  if (this->OptimalQuadratureOrderCalculator)
    {
    this->OptimalQuadratureOrderCalculator->Delete();
    this->OptimalQuadratureOrderCalculator = NULL;
    }
}

int vtkfemriPolyDataNumericKSpaceGenerator::FillInputPortInformation(int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  return 1;
}

void vtkfemriPolyDataNumericKSpaceGenerator::Initialize()
{
  if (this->OptimalQuadratureOrderCalculator)
    {
    this->OptimalQuadratureOrderCalculator->Delete();
    this->OptimalQuadratureOrderCalculator = NULL;
    }
  if (this->UseOptimalAlgorithm)
    {
    vtkImageData* output = vtkImageData::SafeDownCast(this->GetOutput());
    int kSpaceDimensions[3];
    output->GetDimensions(kSpaceDimensions);
    double kSpaceSpacing[3];
    output->GetSpacing(kSpaceSpacing);
    double maxFrequency[3];
    maxFrequency[0] = kSpaceDimensions[0]/2 * kSpaceSpacing[0];
    maxFrequency[1] = kSpaceDimensions[1]/2 * kSpaceSpacing[1];
    maxFrequency[2] = kSpaceDimensions[2]/2 * kSpaceSpacing[2];
    vtkDataSet* input = vtkDataSet::SafeDownCast(this->GetInput());
    double scaledErrorThreshold = this->ErrorThreshold * this->ComputeVoxelVolume();
    this->OptimalQuadratureOrderCalculator = vtkfemriOptimalQuadratureOrderCalculator::New();
    this->OptimalQuadratureOrderCalculator->SetDataSet(input);
//    this->OptimalQuadratureOrderCalculator->SetErrorThreshold(this->ErrorThreshold);
    this->OptimalQuadratureOrderCalculator->SetErrorThreshold(scaledErrorThreshold);
    this->OptimalQuadratureOrderCalculator->SetMaxFrequency(maxFrequency);
    this->OptimalQuadratureOrderCalculator->SetFrequencySpacing(kSpaceSpacing);
    this->OptimalQuadratureOrderCalculator->SetCyclesPerElementResolution(0.01);
    this->OptimalQuadratureOrderCalculator->Initialize();
    }

  this->NumberOfGaussPointEvaluations = 0;
  this->MaximumQuadratureOrderUsed = 0;
}

void vtkfemriPolyDataNumericKSpaceGenerator::EvaluateFourierFunction(double frequency[3], double value[2])
{
  vtkDataSet* input = vtkDataSet::SafeDownCast(this->GetInput());

  vtkDataArray* normals = input->GetPointData()->GetNormals();

  if (!normals)
    {
    vtkErrorMacro("No surface normals found. Please provide a surface with outwards normals defined on it.");
    }

  vtkfemriGaussQuadrature* gaussQuadrature = vtkfemriGaussQuadrature::New();

  value[0] = 0.0;
  value[1] = 0.0;

  int numberOfCells = input->GetNumberOfCells();
  int i;
  for (i=0; i<numberOfCells; i++)
    {
    vtkCell* cell = input->GetCell(i);
    
    if (cell->GetCellDimension() != 2)
      {
      continue;
      }
 
    if (this->UseOptimalAlgorithm)
      {
      int optimalOrder = this->OptimalQuadratureOrderCalculator->ComputeOptimalQuadratureOrder(i,frequency);
//      if (optimalOrder < this->OptimalQuadratureOrderCalculator->GetMaxQuadratureOrder())
//        {
//        switch (cell->GetCellType())
//          {
//          case VTK_QUADRATIC_QUAD:
//          case VTK_QUADRATIC_TRIANGLE:
//              optimalOrder += 1;
//              break;
//          }
//        }
//      else
//        {
//        vtkWarningMacro("Cannot increase quadrature order to account for quadratic element Jacobian because order would exceed maximum allowable order.");
//        }
      if (optimalOrder > this->MaximumQuadratureOrderUsed)
        {
        this->MaximumQuadratureOrderUsed = optimalOrder;
        }
      gaussQuadrature->SetOrder(optimalOrder);
      }
    else
      {
      gaussQuadrature->SetOrder(this->QuadratureOrder);
      }

    gaussQuadrature->Initialize(cell->GetCellType());

    int numberOfCellPoints = cell->GetNumberOfPoints();
 
    double twoPi = 2.0 * vtkMath::Pi();
    int subId = 0;
    double quadraturePCoords[3], quadraturePoint[3];
    double quadratureWeight;
    double* weights = new double[numberOfCellPoints];
    int numberOfQuadraturePoints = 0;
    bool preComputedQuadratureRule = false;
//    if (this->QuadratureOrder == 1 && cell->GetCellType() == VTK_TRIANGLE)
//      {
//      numberOfQuadraturePoints = 1;
//      quadraturePCoords[0] = 0.33333333333333333333333333333333;
//      quadraturePCoords[1] = 0.33333333333333333333333333333333;
//      quadraturePCoords[2] = 0.0;
//      quadratureWeight = 0.5;
//      preComputedQuadratureRule = true;
//      }
//    else
      {
      numberOfQuadraturePoints = gaussQuadrature->GetNumberOfQuadraturePoints();
      }
    int q, j, k;
    for (q=0; q<numberOfQuadraturePoints; q++)
      {
      if (!preComputedQuadratureRule)
        {
        gaussQuadrature->GetQuadraturePoint(q,quadraturePCoords);
        quadratureWeight = gaussQuadrature->GetQuadratureWeight(q);
        }
      this->EvaluateCellLocation(cell,quadraturePCoords,quadraturePoint,weights);
      double normal[3], pointNormal[3];
      normal[0] = normal[1] = normal[2] = 0.0;
      for (k=0; k<numberOfCellPoints; k++)
        {
        normals->GetTuple(cell->GetPointId(k),pointNormal);
        normal[0] += weights[k] * pointNormal[0];
        normal[1] += weights[k] * pointNormal[1];
        normal[2] += weights[k] * pointNormal[2];
        }
      double jacobian = this->ComputeJacobian(cell,quadraturePCoords);

      double kdotx = vtkMath::Dot(quadraturePoint,frequency);
      double kdotn = vtkMath::Dot(normal,frequency);
      double twoPik2 = twoPi * (frequency[0] * frequency[0] + frequency[1] * frequency[1] + frequency[2] * frequency[2]);

      if (frequency[0] == 0.0 && frequency[1] == 0.0 && frequency[2] == 0.0)
        {
        value[0] += this->MagnetizationValue * jacobian * quadratureWeight * normal[2] * quadraturePoint[2];
        value[1] += 0.0;
        }
      else
        {
        value[0] += this->MagnetizationValue * jacobian * quadratureWeight * kdotn / twoPik2 * sin(twoPi * kdotx);
        value[1] += this->MagnetizationValue * jacobian * quadratureWeight * kdotn / twoPik2 * cos(twoPi * kdotx);
        }
      }

    this->NumberOfGaussPointEvaluations += numberOfQuadraturePoints;
    delete[] weights;
    }

  gaussQuadrature->Delete();
}

#if 0
int vtkfemriPolyDataNumericKSpaceGenerator::RequestData(
    vtkInformation* vtkNotUsed( request ),
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector)
{
  vtkPolyData* input = vtkPolyData::SafeDownCast(this->GetInput());

  vtkDataArray* normals = input->GetPointData()->GetNormals();

  if (!normals)
    {
    vtkErrorMacro("No surface normals found. Please provide a surface with outwards normals defined on it.");
    }

  vtkfemriGaussQuadrature* gaussQuadrature = vtkfemriGaussQuadrature::New();

  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkImageData *output = vtkImageData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
 
  vtkIdType updateExtent[6];
  outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(),updateExtent);
  output->SetExtent(updateExtent);
  output->AllocateScalars();

  vtkDoubleArray* newScalars = vtkDoubleArray::SafeDownCast(output->GetPointData()->GetScalars());
  newScalars->FillComponent(0,0.0);
  newScalars->FillComponent(1,0.0);
 
  vtkIdType kSpaceDimensions[3];
  double kSpaceSpacing[3];
  vtkIdType extent[6];
  
  output->GetWholeExtent(extent);
  output->GetSpacing(kSpaceSpacing);
  output->GetDimensions(kSpaceDimensions);

  vtkIdType ijk[3];
  vtkIdType kdx[3];
  double frequency[3];

  int progress = 0;
  cout<<"Progress: "<<progress<<"%"<<endl;
  int numberOfCells = input->GetNumberOfCells();
  int i;
  for (i=0; i<numberOfCells; i++)
    {
    vtkCell* cell = input->GetCell(i);

    int currentProgress = (int)((float)i/numberOfCells*100);
    if (progress != currentProgress)
      {
      progress = currentProgress;
      cout<<"Progress: "<<progress<<"%"<<endl;
      }
 
    if (cell->GetCellDimension() != 2)
      {
      continue;
      }

    if (this->QuadratureOrder > 1)
      {
      gaussQuadrature->SetOrder(this->QuadratureOrder);
      gaussQuadrature->Initialize(cell->GetCellType());
      }
 
    int numberOfCellPoints = cell->GetNumberOfPoints();

    double twoPi = 2.0 * vtkMath::Pi();
    int subId = 0;
    double quadraturePCoords[3], quadraturePoint[3];
    double quadratureWeight;
    double* weights = new double[numberOfCellPoints];
    int numberOfQuadraturePoints = 0;
    if (this->QuadratureOrder > 1)
      {
      numberOfQuadraturePoints = gaussQuadrature->GetNumberOfQuadraturePoints();
      }
    else
      {
      numberOfQuadraturePoints = 1;
      quadraturePCoords[0] = 0.33333333333333333333333333333333;
      quadraturePCoords[1] = 0.33333333333333333333333333333333;
      quadraturePCoords[2] = 0.0;
      quadratureWeight = 0.5;
      }
    int q, j, k;
    for (q=0; q<numberOfQuadraturePoints; q++)
      {
      if (this->QuadratureOrder > 1)
        {
        gaussQuadrature->GetQuadraturePoint(q,quadraturePCoords);
        quadratureWeight = gaussQuadrature->GetQuadratureWeight(q);
        }
      EvaluateCellLocation(cell,quadraturePCoords,quadraturePoint,weights);
      double normal[3], pointNormal[3];
      normal[0] = normal[1] = normal[2] = 0.0;
      for (k=0; k<numberOfCellPoints; k++)
        {
        normals->GetTuple(cell->GetPointId(k),pointNormal);
        normal[0] += weights[k] * pointNormal[0];
        normal[1] += weights[k] * pointNormal[1];
        normal[2] += weights[k] * pointNormal[2];
        }
      double jacobian = this->ComputeJacobian(cell,quadraturePCoords);

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

            double kdotx = vtkMath::Dot(quadraturePoint,frequency);
            double kdotn = vtkMath::Dot(normal,frequency);
            double twoPik2 = twoPi * (frequency[0] * frequency[0] + frequency[1] * frequency[1] + frequency[2] * frequency[2]);

            vtkIdType pointId = output->ComputePointId(ijk);

            double fourierValue[2];
            fourierValue[0] = newScalars->GetComponent(pointId,0);
            fourierValue[1] = newScalars->GetComponent(pointId,1);

            if (frequency[0] == 0.0 && frequency[1] == 0.0 && frequency[2] == 0.0)
              {
              fourierValue[0] += this->MagnetizationValue;
//              fourierValue[0] += 1.0;
              fourierValue[1] += 0.0;
              }
            else
              {
              fourierValue[0] += this->MagnetizationValue * jacobian * quadratureWeight * kdotn / twoPik2 * sin(twoPi * kdotx);
              fourierValue[1] += this->MagnetizationValue * jacobian * quadratureWeight * kdotn / twoPik2 * cos(twoPi * kdotx);
              }

            newScalars->SetComponent(pointId,0,fourierValue[0]);
            newScalars->SetComponent(pointId,1,fourierValue[1]);
            }
          }
        }
      }
    delete[] weights;
    this->NumberOfGaussPointEvaluations += numberOfQuadraturePoints;
    } 

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

        double translation[3];
        translation[0] = -this->Origin[0];
        translation[1] = -this->Origin[1];
        translation[2] = -this->Origin[2];
        if (this->KSpaceDimensionality == 2)
          {
          translation[2] = 0.0; 
          }

        double fourierValue[2];
        vtkIdType pointId = output->ComputePointId(ijk);
        fourierValue[0] = newScalars->GetComponent(pointId,0);
        fourierValue[1] = newScalars->GetComponent(pointId,1);

        double phaseShiftedFourierValue[2];
        vtkfemriKSpaceShift::ShiftPhase(fourierValue,frequency,translation,phaseShiftedFourierValue);
        newScalars->SetComponent(pointId,0,phaseShiftedFourierValue[0]);
        newScalars->SetComponent(pointId,1,phaseShiftedFourierValue[1]);
        }
      }
    }    

  gaussQuadrature->Delete();
  return 1;             
}
#endif

void vtkfemriPolyDataNumericKSpaceGenerator::EvaluateCellLocation(vtkCell* cell, double pcoords[3], double x[3], double* weights)
{
  double pt0[3], pt1[3], pt2[3];
  double u3;
  int subId = 0;
  int i;
  switch (cell->GetCellType())
  {
    case VTK_TRIANGLE:
      cell->Points->GetPoint(0, pt0);
      cell->Points->GetPoint(1, pt1);
      cell->Points->GetPoint(2, pt2);
      u3 = 1.0 - pcoords[0] - pcoords[1];
      for (i=0; i<3; i++)
        {
        x[i] = pt0[i]*u3 + pt1[i]*pcoords[0] + pt2[i]*pcoords[1];
        }
      weights[0] = u3;
      weights[1] = pcoords[0];
      weights[2] = pcoords[1]; 
      break;
    default:
      cell->EvaluateLocation(subId,pcoords,x,weights);
  }
}

void vtkfemriPolyDataNumericKSpaceGenerator::CellInterpolationDerivs(double* pcoords, double* derivs)
{
  //TODO: not needed, no virtual functions are called (so far. They will be starting with the next VTK release).
}

double vtkfemriPolyDataNumericKSpaceGenerator::ComputeJacobian(vtkCell* cell, double pcoords[3])
{
  double jacobian = 0.0;
  
  int cellDimension = cell->GetCellDimension();
  
  if (cellDimension != 2)
  {
    return 0.0;
  }

  int numberOfCellPoints = cell->GetNumberOfPoints();
  double* derivs = new double[2*numberOfCellPoints];
  
  switch (cell->GetCellType())
  {
    case VTK_TRIANGLE:
//      vtkTriangle::SafeDownCast(cell)->InterpolationDerivs(pcoords,derivs);
      vtkTriangle::InterpolationDerivs(pcoords,derivs);
      break;
    case VTK_QUADRATIC_TRIANGLE:
//      vtkQuadraticTriangle::SafeDownCast(cell)->InterpolationDerivs(pcoords,derivs);
      vtkQuadraticTriangle::InterpolationDerivs(pcoords,derivs);
      break;
    case VTK_QUAD:
      vtkQuad::InterpolationDerivs(pcoords,derivs);
      break;
    case VTK_QUADRATIC_QUAD:
      vtkQuadraticQuad::InterpolationDerivs(pcoords,derivs);
      break;
    default:
      vtkErrorMacro("Error: unsupported cell type.");
      return 0.0;
  }

  int i, j;

  double jacobianMatrixTr[2][3];
  for (i=0; i<3; i++)
  {
    jacobianMatrixTr[0][i] = jacobianMatrixTr[1][i] = 0.0;
  }

  double x[3];
  for (j=0; j<numberOfCellPoints; j++)
  {
    cell->GetPoints()->GetPoint(j,x);
    for (i=0; i<3; i++)
    {
      jacobianMatrixTr[0][i] += x[i] * derivs[j];
      jacobianMatrixTr[1][i] += x[i] * derivs[numberOfCellPoints+j];
    }
  }
  delete[] derivs;

  double jacobianMatrixSquared[2][2];
  jacobianMatrixSquared[0][0] = vtkMath::Dot(jacobianMatrixTr[0],jacobianMatrixTr[0]);
  jacobianMatrixSquared[0][1] = vtkMath::Dot(jacobianMatrixTr[0],jacobianMatrixTr[1]);
  jacobianMatrixSquared[1][0] = vtkMath::Dot(jacobianMatrixTr[1],jacobianMatrixTr[0]);
  jacobianMatrixSquared[1][1] = vtkMath::Dot(jacobianMatrixTr[1],jacobianMatrixTr[1]);

  double jacobianSquared = vtkMath::Determinant2x2(jacobianMatrixSquared[0],jacobianMatrixSquared[1]);

  if (jacobianSquared < 0.0)
  {
    jacobianSquared = fabs(jacobianSquared);
  }

  jacobian = sqrt(jacobianSquared);
 
  return jacobian;
}

