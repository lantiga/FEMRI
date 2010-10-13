/*=========================================================================

  Program:   femri
  Module:    $RCSfile: vtkfemriOptimalQuadratureOrderCalculator.cxx,v $
  Language:  C++
  Date:      $Date: 2008/11/03 17:00:30 $
  Version:   $Revision: 1.3 $

  Copyright (c) Luca Antiga, David Steinman. All rights reserved.
  See LICENCE file for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm 
  for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "vtkfemriOptimalQuadratureOrderCalculator.h"
#include "vtkfemriGaussQuadrature.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkCell.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"

vtkStandardNewMacro(vtkfemriOptimalQuadratureOrderCalculator);
vtkCxxRevisionMacro(vtkfemriOptimalQuadratureOrderCalculator, "$Revision: 1.3 $");

vtkfemriOptimalQuadratureOrderCalculator::vtkfemriOptimalQuadratureOrderCalculator()
{
  this->ErrorThreshold = 1E-4;
  this->MaxQuadratureOrder = 43;
  this->MaxFrequency[0] = 0.0;
  this->MaxFrequency[1] = 0.0;
  this->MaxFrequency[2] = 0.0;
  this->FrequencySpacing[0] = 0.0;
  this->FrequencySpacing[1] = 0.0;
  this->FrequencySpacing[2] = 0.0;
  this->CyclesPerElementResolution = 1.0;
  this->DataSet = NULL;
  this->DesignTable = NULL;
}

vtkfemriOptimalQuadratureOrderCalculator::~vtkfemriOptimalQuadratureOrderCalculator()
{
  if (this->DesignTable)
    {
    this->DesignTable->Delete();
    this->DesignTable = NULL;
    }
  if (this->DataSet)
    {
    this->DataSet->Delete();
    this->DataSet = NULL;
    }
}

void vtkfemriOptimalQuadratureOrderCalculator::Initialize()
{
  if (this->DataSet == NULL)
    {
    vtkErrorMacro(<<"Error: DataSet not set!");
    return;
    }

  double cyclesPerElement;
  double maxCyclesPerElement = 0.0;
  vtkIdType numberOfCells = this->DataSet->GetNumberOfCells();
  int KIJK[3];
  KIJK[0] = vtkMath::Floor(this->MaxFrequency[0] / this->FrequencySpacing[0]) + 1;
  KIJK[1] = vtkMath::Floor(this->MaxFrequency[1] / this->FrequencySpacing[1]) + 1;
  KIJK[2] = vtkMath::Floor(this->MaxFrequency[2] / this->FrequencySpacing[2]) + 1;
  double frequency[3];
  int kijk[3];
  for (kijk[2]=-KIJK[2]; kijk[2]<KIJK[2]; kijk[2]++)
    {
    for (kijk[1]=-KIJK[1]; kijk[1]<KIJK[1]; kijk[1]++)
      {
      for (kijk[0]=-KIJK[0]; kijk[0]<KIJK[0]; kijk[0]++)
        {
        frequency[0] = kijk[0] * this->FrequencySpacing[0];
        frequency[1] = kijk[1] * this->FrequencySpacing[1];
        frequency[2] = kijk[2] * this->FrequencySpacing[2];
        for (int i=0; i<numberOfCells; i++)
          {
          cyclesPerElement = this->ComputeCyclesPerElement(i,frequency);
          if (cyclesPerElement > maxCyclesPerElement)
            {
            maxCyclesPerElement = cyclesPerElement;
            }
          }
        }
      }
    }

  this->BuildDesignTable(maxCyclesPerElement);
}

void vtkfemriOptimalQuadratureOrderCalculator::BuildDesignTable(double maxCyclesPerElement)
{
  if (this->DesignTable)
    {
    this->DesignTable->Delete();
    this->DesignTable = NULL;
    }
//  int cyclesPerElementSteps = vtkMath::Floor(maxCyclesPerElement/this->CyclesPerElementResolution)+1;
  int cyclesPerElementSteps = (int)ceil(maxCyclesPerElement/this->CyclesPerElementResolution);
  this->DesignTable = vtkImageData::New();
  this->DesignTable->SetScalarTypeToDouble();
  this->DesignTable->SetDimensions(cyclesPerElementSteps,this->MaxQuadratureOrder+1,1);  
  this->DesignTable->SetSpacing(this->CyclesPerElementResolution,1.0,1.0);
  this->DesignTable->AllocateScalars();
  vtkDataArray* errorScalars = this->DesignTable->GetPointData()->GetScalars();
  errorScalars->FillComponent(0,VTK_DOUBLE_MIN);

  vtkfemriGaussQuadrature* quadratureRule = vtkfemriGaussQuadrature::New();

  double twoPi = 2.0 * vtkMath::Pi();

  int ijk[3];
  ijk[0] = ijk[1] = ijk[2] = 0;

  bool sufficientMaxOrder = true;

  for (int order=0; order<this->MaxQuadratureOrder+1; order++)
    {
    quadratureRule->SetOrder(order);
    quadratureRule->Initialize(VTK_LINE);

    ijk[1] = order;
    double maxRowError = 0.0;
    for (int cycleStep=0; cycleStep<cyclesPerElementSteps; cycleStep++)
      {
      ijk[0] = cycleStep;
      double cycles = cycleStep * this->CyclesPerElementResolution;
      vtkIdType id = this->DesignTable->ComputePointId(ijk);
      double exactIntegral = 0.0;
      if (cycles == 0)
        {
        exactIntegral = 1.0;
        }
      else
        {
        exactIntegral = 1.0 / (twoPi*cycles) * sqrt(2.0 - 2.0 * cos(twoPi*cycles));
        }
      double quadratureIntegralReal, quadratureIntegralImag;
      quadratureIntegralReal = quadratureIntegralImag = 0.0;
      int numberOfQuadraturePoints = quadratureRule->GetNumberOfQuadraturePoints();
      for (int q=0; q<numberOfQuadraturePoints; q++)
        {
        double quadraturePoint = quadratureRule->GetQuadraturePoint(q)[0];
        quadratureIntegralReal += quadratureRule->GetQuadratureWeight(q) * cos(twoPi*cycles*quadraturePoint);
        quadratureIntegralImag += -1.0 * quadratureRule->GetQuadratureWeight(q) * sin(twoPi*cycles*quadraturePoint);
        }
      double quadratureIntegral = sqrt(quadratureIntegralReal*quadratureIntegralReal + quadratureIntegralImag*quadratureIntegralImag);
      double error = fabs(exactIntegral - quadratureIntegral);
      errorScalars->SetComponent(id,0,error);
      if (error > maxRowError)
        {
        maxRowError = error;
        }
      }
    if (order == this->MaxQuadratureOrder && maxRowError > this->ErrorThreshold)
      {
      sufficientMaxOrder = false;
      }
    }

  quadratureRule->Delete();

  if (!sufficientMaxOrder)
    {
    vtkWarningMacro(<< "MaxQuadratureOrder is not sufficient to provide the desired error bound over the whole mesh. Please choose a larger MaxQuadratureOrder.");
    }
}

double vtkfemriOptimalQuadratureOrderCalculator::ComputeCyclesPerElement(vtkIdType cellId, double frequency[3])
{
  vtkCell* cell = this->DataSet->GetCell(cellId);
  int numberOfEdges = cell->GetNumberOfEdges();
  double maxCycles = 0.0;
  int i, j;
  for (i=0; i<numberOfEdges; i++)
  {
    vtkCell* edge = cell->GetEdge(i);
    int edgePoints = edge->GetNumberOfPoints();
    double edgeLength = 0.0;
    for (j=0; j<edgePoints-1; j++)
    {
      double point0[3], point1[3];
      double edgeVector[3];
      edge->GetPoints()->GetPoint(j,point0);
      edge->GetPoints()->GetPoint(j+1,point1);
      edgeVector[0] = point1[0] - point0[0];
      edgeVector[1] = point1[1] - point0[1];
      edgeVector[2] = point1[2] - point0[2];
      edgeLength += vtkMath::Norm(edgeVector);
    }
    for (j=0; j<edgePoints-1; j++)
    {
      double point0[3], point1[3];
      double edgeVector[3];
      edge->GetPoints()->GetPoint(j,point0);
      edge->GetPoints()->GetPoint(j+1,point1);
      edgeVector[0] = point1[0] - point0[0];
      edgeVector[1] = point1[1] - point0[1];
      edgeVector[2] = point1[2] - point0[2];
      double subEdgeLength = vtkMath::Norm(edgeVector);
      double cycles = fabs(vtkMath::Dot(frequency,edgeVector)) / subEdgeLength * edgeLength;
//      cout<<cycles<<" "<<frequency[0]<<" "<<frequency[1]<<" "<<frequency[2]<<" "<< <<endl;
      if (cycles > maxCycles)
      {
        maxCycles = cycles;
      }
    }
  }
  return maxCycles;
}

int vtkfemriOptimalQuadratureOrderCalculator::ComputeOptimalQuadratureOrder(vtkIdType cellId, double frequency[3])
{
  double cyclesPerElement = this->ComputeCyclesPerElement(cellId,frequency);

//  int cyclesPerElementStep = vtkMath::Floor(cyclesPerElement/this->CyclesPerElementResolution)+1;
  int cyclesPerElementStep = (int)ceil(cyclesPerElement/this->CyclesPerElementResolution);

  int ijk[3];
  ijk[0] = ijk[1] = ijk[2] = 0;
  ijk[0] = cyclesPerElementStep;

  if (ijk[0] > this->DesignTable->GetDimensions()[0]-1)
    {
    ijk[0]--;
    }

  int optimalQuadratureOrder = this->MaxQuadratureOrder;

  vtkIdType id;
  double error;
  vtkDataArray* errorScalars = this->DesignTable->GetPointData()->GetScalars();
  for (int i=this->MaxQuadratureOrder; i>=0; i--)
    {
    ijk[1] = i;
    id = this->DesignTable->ComputePointId(ijk);
    error = errorScalars->GetComponent(id,0);
    if (error > this->ErrorThreshold)
      {
      if (i < this->MaxQuadratureOrder)
        {
        optimalQuadratureOrder = i+1;
        }
      break;
      }
    else
      {
      optimalQuadratureOrder = i;
      }
    }

//  cout<<optimalQuadratureOrder<<" "<<cyclesPerElement<<" "<<error<<" "<<id<<" "<<endl;
  return optimalQuadratureOrder;
}
#if 0
int vtkfemriOptimalQuadratureOrderCalculator::ComputeOptimalQuadratureOrder(vtkIdType cellId, double frequency[3])
{
  double cyclesPerElement = this->ComputeCyclesPerElement(cellId,frequency);

  int cyclesPerElementStep = vtkMath::Floor(cyclesPerElement/this->CyclesPerElementResolution);

  int ijk0[3], ijk1[3];
  ijk0[0] = ijk0[1] = ijk0[2] = 0;
  ijk1[0] = ijk1[1] = ijk1[2] = 0;
  ijk0[0] = cyclesPerElementStep;
  ijk1[0] = cyclesPerElementStep+1;

  if (ijk1[0] > this->DesignTable->GetDimensions()[0]-1)
    {
    ijk1[0]--;
    }

  int optimalQuadratureOrder = this->MaxQuadratureOrder;

  vtkIdType id0, id1;
  double error0, error1, error;
  vtkDataArray* errorScalars = this->DesignTable->GetPointData()->GetScalars();
  for (int i=this->MaxQuadratureOrder; i>=0; i--)
    {
    ijk0[1] = i;
    ijk1[1] = i;
    id0 = this->DesignTable->ComputePointId(ijk0);
    id1 = this->DesignTable->ComputePointId(ijk1);
    error0 = errorScalars->GetComponent(id0,0);
    error1 = errorScalars->GetComponent(id1,0);
    error = (error1 - error0) * (cyclesPerElement - cyclesPerElementStep*this->CyclesPerElementResolution) + error0;
    if (error > this->ErrorThreshold)
      {
      if (i < this->MaxQuadratureOrder)
        {
        optimalQuadratureOrder = i+1;
        }
      break;
      }
    else
      {
      optimalQuadratureOrder = i;
      }
    }

//  cout<<cyclesPerElement<<" "<<optimalQuadratureOrder<<" "<<error<<" "<<error0<<" "<<error1<<" "<<id0<<" "<<id1<<" "<<this->DesignTable->GetDimensions()[0]<<endl;
  return optimalQuadratureOrder;
}
#endif
#if 0
void vtkfemriOptimalQuadratureOrderCalculator::InitDesignCurves()
{
  const unsigned int MAX_MAP_QUAD_ORDERS = 20; 
  
  // array of supported discrete quadrature order values (must be same for all curves)  
  unsigned int qarr[MAX_MAP_QUAD_ORDERS] = {1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39};
  
  // arrays of corresponding cycles/element values for error E=1e-4, 1e-3, etc.
  // these are the reference values
  
  // 1D curves
  double ce2[MAX_MAP_QUAD_ORDERS] = {0.26, 0.77, 1.34, 1.93, 2.53, 3.13, 3.91, 4.71, 4.97, 5.59, 6.21, 6.93, 7.73, 8.41, 8.70, 9.33, 9.98, 10.00, 10.00, 10.00};
  double ce3[MAX_MAP_QUAD_ORDERS] = {0.08, 0.42, 0.87, 1.36, 1.89, 2.43, 2.98, 3.54, 4.11, 4.69, 5.27, 5.85, 6.44, 7.03, 7.62, 8.22, 8.81, 9.41, 10.00, 10.00};
  double ce4[MAX_MAP_QUAD_ORDERS] = {0.03, 0.24, 0.58, 1.00, 1.45, 1.94, 2.44, 2.96, 3.49, 4.03, 4.58, 5.13, 5.69, 6.25, 6.82, 7.39, 7.96, 8.54, 9.12, 9.70};
  double ce5[MAX_MAP_QUAD_ORDERS] = {0.01, 0.13, 0.39, 0.74, 1.14, 1.57, 2.03, 2.51, 3.00, 3.51, 4.03, 4.55, 5.08, 5.62, 6.16, 6.71, 7.26, 7.82, 8.38, 8.94};
  double ce6[MAX_MAP_QUAD_ORDERS] = {0.01, 0.08, 0.27, 0.55, 0.89, 1.28, 1.70, 2.14, 2.60, 3.08, 3.57, 4.06, 4.57, 5.08, 5.60, 6.13, 6.67, 7.20, 7.73, 8.28};
  
  // vectors to hold the cycles/element and qorder values
  vector<double> cycElemE2;
  vector<double> cycElemE3;
  vector<double> cycElemE4;
  vector<double> cycElemE5;
  vector<double> cycElemE6;
  vector<unsigned int> QOrders;
  
  // initialize vectors with values
  for(int i=0; i < MAX_MAP_QUAD_ORDERS; i++)
  {
    cycElemE2.push_back(ce2[i]);
    cycElemE3.push_back(ce3[i]);
    cycElemE4.push_back(ce4[i]);
    cycElemE5.push_back(ce5[i]);
    cycElemE6.push_back(ce6[i]);
    QOrders.push_back(qarr[i]);
  }
  
  // create design curve objects from vectors and corresponding error values
  this->DesignCurves->AddDesignCurve(1e-2,cycElemE2,QOrders);
  this->DesignCurves->AddDesignCurve(1e-3,cycElemE3,QOrders);
  this->DesignCurves->AddDesignCurve(1e-4,cycElemE4,QOrders);
  this->DesignCurves->AddDesignCurve(1e-5,cycElemE5,QOrders);
  this->DesignCurves->AddDesignCurve(1e-6,cycElemE6,QOrders);
}

double vtkfemriOptimalQuadratureOrderCalculator::ComputeCyclesPerElement(vtkCell* cell, double frequency[3])
{
  int numberOfEdges = cell->GetNumberOfEdges();
 
  double maxCycles = 0.0;
  int i, j;
  for (i=0; i<numberOfEdges; i++)
  {
    vtkCell* edge = cell->GetEdge(i);
    int edgePoints = edge->GetNumberOfPoints();
    for (j=0; j<edgePoints-1; j++)
    {
      double point0[3], point1[3];
      double edgeVector[3];
      edge->GetPoints()->GetPoint(j,point0);
      edge->GetPoints()->GetPoint(j+1,point1);
      edgeVector[0] = point1[0] - point0[0];
      edgeVector[1] = point1[1] - point0[1];
      edgeVector[2] = point1[2] - point0[2];
      double cycles = fabs(vtkMath::Dot(frequency,edgeVector));
      if (cycles > maxCycles)
      {
        maxCycles = cycles;
      }
    }
  }
  return maxCycles;
}

int vtkfemriOptimalQuadratureOrderCalculator::ComputeOptimalQuadratureOrder(vtkCell* cell, double frequency[3])
{
  double cycles = this->ComputeCyclesPerElement(cell,frequency);
  int qorder = -1;  // default -1 as cycles/element was not found in interpolated design curve
  
  // initialize an empty design curve to interpolate
  femriIntegrationDesignCurve currCurve = femriIntegrationDesignCurve(this->ErrorThreshold);
  this->DesignCurves->Interpolate(&currCurve);
  
  for (unsigned int i = 0; i < currCurve.GetNumEntries(); i++)
  {
    if ((cycles >= 0) && (cycles <= currCurve.GetCycles(i)))
      qorder = currCurve.GetQOrder(i);
  }
  // default to 42 if cycles/element is out of range of curve
  if ((cycles >= 0) && (cycles > currCurve.GetCycles(currCurve.GetNumEntries() - 1)))
    qorder = 42;
  
  if (qorder == -1)  // should never reach this, but just in case
  {
    vtkErrorMacro("ERROR: optimal order not found.");
    return -1;
  }
  
  return qorder;
}
#endif

