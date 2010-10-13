/*=========================================================================

  Program:   femri
  Module:    $RCSfile: vtkfemriUnstructuredGridNumericKSpaceGenerator.cxx,v $
  Language:  C++
  Date:      $Date: 2008/11/04 15:46:07 $
  Version:   $Revision: 1.13 $

  Copyright (c) Luca Antiga, David Steinman. All rights reserved.
  See LICENCE file for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm 
  for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "vtkfemriUnstructuredGridNumericKSpaceGenerator.h"
#include "vtkfemriOptimalQuadratureOrderCalculator.h"
#include "vtkUnstructuredGrid.h"
#include "vtkCell.h"
#include "vtkPointData.h"
#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkObjectFactory.h"
#include "vtkfemriGaussQuadrature.h"
#include "vtkHexahedron.h"
#include "vtkQuadraticHexahedron.h"
#include "vtkWedge.h"
#include "vtkQuadraticWedge.h"
#if VTK_MAJOR_VERSION > 5 || (VTK_MAJOR_VERSION == 5 && VTK_MINOR_VERSION > 0)
#include "vtkBiQuadraticQuadraticWedge.h"
#endif
#include "vtkTetra.h"
#include "vtkQuadraticTetra.h"

vtkStandardNewMacro(vtkfemriUnstructuredGridNumericKSpaceGenerator);
vtkCxxRevisionMacro(vtkfemriUnstructuredGridNumericKSpaceGenerator, "$Revision: 1.13 $");

vtkfemriUnstructuredGridNumericKSpaceGenerator::vtkfemriUnstructuredGridNumericKSpaceGenerator()
{
  this->OptimalQuadratureOrderCalculator = NULL;

  this->UseOptimalAlgorithm = 1;
  this->ErrorThreshold = 1E-4;

  this->NumberOfSubdivisions = 0;
  this->QuadratureOrder = 4;

  this->SliceSelection = 0;
  this->SliceThickness = 1.0;
  this->SliceOrigin = 1.0;
  this->SliceProfile = FEMRI_IDEAL_SLICE_PROFILE;

  this->UniformMagnetization = 1;
  this->MagnetizationValue = 1.0;
  this->MagnetizationArrayName = NULL;

  this->NumberOfGaussPointEvaluations = 0;
  this->MaximumQuadratureOrderUsed = 0;

  // Initialize default design curves and interpolator object
//  this->InitDesignCurves();  

  this->SetNumberOfInputPorts(1);
}
#if 0
int vtkfemriUnstructuredGridNumericKSpaceGenerator::InitDesignCurves()
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
  
  // 3D curves
  //double ce4[MAX_MAP_QUAD_ORDERS] = {0.03, 0.21, 0.6, 1.08, 1.59, 2.16, 2.73, 3.33, 3.93, 4.53, 5.13, 5.73, 6.33, 6.93, 7.53, 8.13, 8.73, 9.33, 9.93, 10.53};
  //double ce3[MAX_MAP_QUAD_ORDERS] = {0.06, 0.36, 0.87, 1.47, 2.07, 2.7, 3.33, 3.96, 4.62, 5.28, 5.94, 6.6, 7.26, 7.92, 8.58, 9.24, 9.90, 10.56, 11.22, 11.88}; 

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
  femriIntegrationDesignCurve *CE2 = new femriIntegrationDesignCurve(1e-2, cycElemE2, QOrders);
  femriIntegrationDesignCurve *CE3 = new femriIntegrationDesignCurve(1e-3, cycElemE3, QOrders);
  femriIntegrationDesignCurve *CE4 = new femriIntegrationDesignCurve(1e-4, cycElemE4, QOrders);
  femriIntegrationDesignCurve *CE5 = new femriIntegrationDesignCurve(1e-5, cycElemE5, QOrders);
  femriIntegrationDesignCurve *CE6 = new femriIntegrationDesignCurve(1e-6, cycElemE6, QOrders);
    
  // initialize interpolator object holding all the design curves
  this->designCurves = new femriIntegrationDesignCurveInterpolator();
    
  // add references
  this->designCurves->AddDesignCurve(CE2);
  this->designCurves->AddDesignCurve(CE3);
  this->designCurves->AddDesignCurve(CE4);
  this->designCurves->AddDesignCurve(CE5);
  this->designCurves->AddDesignCurve(CE6);
  
  return 1;  
}
#endif

vtkfemriUnstructuredGridNumericKSpaceGenerator::~vtkfemriUnstructuredGridNumericKSpaceGenerator()
{
  if (this->MagnetizationArrayName)
    {
    delete[] this->MagnetizationArrayName;
    this->MagnetizationArrayName = NULL;
    }
  if (this->OptimalQuadratureOrderCalculator)
    {
    this->OptimalQuadratureOrderCalculator->Delete();
    this->OptimalQuadratureOrderCalculator = NULL;
    }
}

int vtkfemriUnstructuredGridNumericKSpaceGenerator::FillInputPortInformation(int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
  return 1;
}

void vtkfemriUnstructuredGridNumericKSpaceGenerator::Initialize()
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

void vtkfemriUnstructuredGridNumericKSpaceGenerator::EvaluateFourierFunction(double frequency[3], double value[2])
{
  vtkUnstructuredGrid* input = vtkUnstructuredGrid::SafeDownCast(this->GetInput());
  vtkImageData* output = vtkImageData::SafeDownCast(this->GetOutput());

  vtkDataArray* magnetizationArray = NULL;
  if (!this->UniformMagnetization)
  {
    if (!this->MagnetizationArrayName) 
    {
      vtkErrorMacro("Error: MagnetizationArrayName not specified.");
      return;
    }

    magnetizationArray = input->GetPointData()->GetArray(this->MagnetizationArrayName);
    if (!magnetizationArray)
    {
      vtkErrorMacro("Error: MagnetizationArray with name specified does not exist.");
      return;
    }
  }
  
  vtkfemriGaussQuadrature* gaussQuadrature = vtkfemriGaussQuadrature::New();

  value[0] = 0.0;
  value[1] = 0.0;

  int numberOfCells = input->GetNumberOfCells();
  int i;
  for (i=0; i<numberOfCells; i++)
  {
    vtkCell* cell = input->GetCell(i);
    
    if (cell->GetCellDimension() != 3)
    {
      continue;
    }
  
    if (this->UseOptimalAlgorithm)
    {
      int optimalOrder = this->OptimalQuadratureOrderCalculator->ComputeOptimalQuadratureOrder(i,frequency);
      switch (cell->GetCellType())
        {
        case VTK_QUADRATIC_HEXAHEDRON:
        case VTK_QUADRATIC_WEDGE:
#if VTK_MAJOR_VERSION > 5 || (VTK_MAJOR_VERSION == 5 && VTK_MINOR_VERSION > 0)
        case VTK_BIQUADRATIC_QUADRATIC_WEDGE:
#endif
        case VTK_QUADRATIC_TETRA:
          if (optimalOrder < this->OptimalQuadratureOrderCalculator->GetMaxQuadratureOrder())
            {
            optimalOrder += 1;
            }
          else
            {
            vtkWarningMacro("Cannot increase quadrature order to account for quadratic element Jacobian because order would exceed maximum allowable order.");
            }
          break;
        }

      if (!this->UniformMagnetization)
        {
        if (optimalOrder < this->OptimalQuadratureOrderCalculator->GetMaxQuadratureOrder())
          {
          optimalOrder += 1;
          }
        else
          {
          vtkWarningMacro("Cannot increase quadrature order to account for non-constant magnetization because order would exceed maximum allowable order.");
          }
        switch (cell->GetCellType())
          {
          case VTK_QUADRATIC_HEXAHEDRON:
          case VTK_QUADRATIC_WEDGE:
#if VTK_MAJOR_VERSION > 5 || (VTK_MAJOR_VERSION == 5 && VTK_MINOR_VERSION > 0)
          case VTK_BIQUADRATIC_QUADRATIC_WEDGE:
#endif
          case VTK_QUADRATIC_TETRA:
            if (optimalOrder < this->OptimalQuadratureOrderCalculator->GetMaxQuadratureOrder())
              {
              optimalOrder += 1;
              }
            else
              {
              vtkWarningMacro("Cannot increase quadrature order to account for non-linear magnetization because order would exceed maximum allowable order.");
              }
            break;
          }
        }
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
    int q, j;
    for (q=0; q<gaussQuadrature->GetNumberOfQuadraturePoints(); q++)
    {
      gaussQuadrature->GetQuadraturePoint(q,quadraturePCoords);
      quadratureWeight = gaussQuadrature->GetQuadratureWeight(q);
      cell->EvaluateLocation(subId,quadraturePCoords,quadraturePoint,weights);
      double kdotx = vtkMath::Dot(quadraturePoint,frequency);
      double complexExpRe = cos(-twoPi * kdotx);
      double complexExpIm = sin(-twoPi * kdotx);
      double magnetization = 0.0;
      if (this->UniformMagnetization)
      {
        magnetization = this->MagnetizationValue;
      }
      else
      {
        for (j=0; j<numberOfCellPoints; j++)
        {
          magnetization += magnetizationArray->GetComponent(cell->GetPointId(j),0) * weights[j];
        }
      }
      if (this->SliceSelection)
      {
        double sliceMagnetization = 0.0;
        switch (this->SliceProfile)
        {
          case FEMRI_IDEAL_SLICE_PROFILE:
            sliceMagnetization = this->IdealSliceProfile(quadraturePoint[2]);
            break;
          case FEMRI_TRAPEZOIDAL_SLICE_PROFILE:
            sliceMagnetization = this->TrapezoidalSliceProfile(quadraturePoint[2]);
            break;
          case FEMRI_QUADRATIC_SLICE_PROFILE:
            sliceMagnetization = this->QuadraticSliceProfile(quadraturePoint[2]);
            break;
          default:
            vtkErrorMacro("Invalid slice profile."); 
        }
        magnetization *= sliceMagnetization;
      }
      double jacobian = this->ComputeJacobian(cell,quadraturePCoords);
      value[0] += magnetization * jacobian * quadratureWeight * complexExpRe;
      value[1] += magnetization * jacobian * quadratureWeight * complexExpIm;
    }
    this->NumberOfGaussPointEvaluations += gaussQuadrature->GetNumberOfQuadraturePoints();
    delete[] weights;
  }
  
  gaussQuadrature->Delete();
}

double vtkfemriUnstructuredGridNumericKSpaceGenerator::ComputeJacobian(vtkCell* cell, double pcoords[3])
{
  int numberOfCellPoints = cell->GetNumberOfPoints();
  double* derivs = new double[3*numberOfCellPoints];
  switch (cell->GetCellType())
  {
    case VTK_HEXAHEDRON:
      vtkHexahedron::SafeDownCast(cell)->InterpolationDerivs(pcoords,derivs);
      break;
    case VTK_QUADRATIC_HEXAHEDRON:
      vtkQuadraticHexahedron::SafeDownCast(cell)->InterpolationDerivs(pcoords,derivs);
      break;
    case VTK_WEDGE:
      vtkWedge::SafeDownCast(cell)->InterpolationDerivs(pcoords,derivs);
      break;
    case VTK_QUADRATIC_WEDGE:
      vtkQuadraticWedge::SafeDownCast(cell)->InterpolationDerivs(pcoords,derivs);
      break;
#if VTK_MAJOR_VERSION > 5 || (VTK_MAJOR_VERSION == 5 && VTK_MINOR_VERSION > 0)
    case VTK_BIQUADRATIC_QUADRATIC_WEDGE:
      vtkBiQuadraticQuadraticWedge::SafeDownCast(cell)->InterpolationDerivs(pcoords,derivs);
      break;
#endif
    case VTK_TETRA:
#if VTK_MAJOR_VERSION > 5 || (VTK_MAJOR_VERSION == 5 && VTK_MINOR_VERSION > 0)
      vtkTetra::SafeDownCast(cell)->InterpolationDerivs(pcoords,derivs);
#else
      vtkTetra::SafeDownCast(cell)->InterpolationDerivs(derivs);
#endif
      break;
    case VTK_QUADRATIC_TETRA:
      vtkQuadraticTetra::SafeDownCast(cell)->InterpolationDerivs(pcoords,derivs);
      break;
    default:
      vtkErrorMacro("Error: unsupported cell type.");
      delete[] derivs;
      return 0.0;
  }

  int i, j;
  
  double jacobianMatrix[3][3];
  for (i=0; i < 3; i++)
  {
    jacobianMatrix[0][i] = jacobianMatrix[1][i] = jacobianMatrix[2][i] = 0.0;
  }
  double x[3];

  for (j=0; j<numberOfCellPoints; j++)
  {
    cell->GetPoints()->GetPoint(j,x);
    for (i=0; i<3; i++)
    {
      jacobianMatrix[0][i] += x[i] * derivs[j];
      jacobianMatrix[1][i] += x[i] * derivs[numberOfCellPoints+j];
      jacobianMatrix[2][i] += x[i] * derivs[2*numberOfCellPoints+j];
    }
  }

  delete[] derivs;
 
  double jacobian = vtkMath::Determinant3x3(jacobianMatrix);

  if (jacobian < 0.0)
  {
    jacobian = fabs(jacobian);
  }

  return jacobian;
}

#if 0
int vtkfemriUnstructuredGridNumericKSpaceGenerator::ComputeOptimalOrder(vtkCell* cell, double frequency[3])
{
  double cycles = this->ComputeCyclesPerElement(cell,frequency);
  int qorder = -1;  // default -1 as cycles/element was not found in interpolated design curve
  
  // initialize an empty design curve to interpolate
  femriIntegrationDesignCurve currCurve = femriIntegrationDesignCurve(this->ErrorThreshold);
  this->designCurves->Interpolate(&currCurve);
  
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

double vtkfemriUnstructuredGridNumericKSpaceGenerator::ComputeCyclesPerElement(vtkCell* cell, double frequency[3])
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
#endif

double vtkfemriUnstructuredGridNumericKSpaceGenerator::TrapezoidalSliceProfile(double z)
{
  double mag = 0.0;
  double slope;
  
  double z_slice_bound_upper = this->SliceOrigin + 0.5 * this->SliceThickness;
  double z_slice_bound_lower = this->SliceOrigin - 0.5 * this->SliceThickness;
  
  double z_tail_bound_upper = z_slice_bound_upper + 0.3 * 0.5 * this->SliceThickness;
  double z_tail_bound_lower = z_slice_bound_lower - 0.3 * 0.5 * this->SliceThickness;

  if ((z <= z_slice_bound_upper) && (z >= z_slice_bound_lower))
  {
    mag = 1.0;
  }
  else if ((z > z_slice_bound_upper) && (z <= z_tail_bound_upper))
  {
    slope = -1/(z_tail_bound_upper - z_slice_bound_upper);
    mag = slope*(z - z_slice_bound_upper) + 1.0;  
  }
  else if ((z < z_slice_bound_lower) && (z >= z_tail_bound_lower))
  {
    slope = 1/(z_slice_bound_lower - z_tail_bound_lower);
    mag = slope*(z - z_tail_bound_lower);  
  }
  else 
  { 
    mag = 0.0; 
  }
 
  return mag;
}

double vtkfemriUnstructuredGridNumericKSpaceGenerator::IdealSliceProfile(double z)
{
  double mag = 0.0;

  double z_slice_bound_upper = this->SliceOrigin + 0.5 * this->SliceThickness;
  double z_slice_bound_lower = this->SliceOrigin - 0.5 * this->SliceThickness;
  
  if ((z <= z_slice_bound_upper) && (z >= z_slice_bound_lower))
  {
    mag = 1.0;
  }
  else 
  { 
    mag = 0.0;
  }
  
  return mag;
}

double vtkfemriUnstructuredGridNumericKSpaceGenerator::QuadraticSliceProfile(double z)
{
  double mag = 0;

  double z_slice_bound_upper = this->SliceOrigin + 0.5 * this->SliceThickness;
  double z_slice_bound_lower = this->SliceOrigin - 0.5 * this->SliceThickness;
  
  if ((z <= z_slice_bound_upper) && (z >= z_slice_bound_lower))
  {
    mag = -0.7 * (mag - this->SliceOrigin) * (mag - this->SliceOrigin) + 1.0;
  }
  else 
  { 
    mag = 0.0; 
  }
  
  return mag;
}
