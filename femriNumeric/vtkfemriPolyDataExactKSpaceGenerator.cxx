/*=========================================================================

  Program:   femri
  Module:    $RCSfile: vtkfemriPolyDataExactKSpaceGenerator.cxx,v $
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

#include "vtkfemriPolyDataExactKSpaceGenerator.h"
#include "vtkPolyData.h"
#include "vtkCell.h"
#include "vtkCellData.h"
#include "vtkTriangle.h"
#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkObjectFactory.h"

vtkStandardNewMacro(vtkfemriPolyDataExactKSpaceGenerator);
vtkCxxRevisionMacro(vtkfemriPolyDataExactKSpaceGenerator, "$Revision: 1.7 $");

vtkfemriPolyDataExactKSpaceGenerator::vtkfemriPolyDataExactKSpaceGenerator()
{
  this->MagnetizationValue = 1.0;

  this->SetNumberOfInputPorts(1);
}

vtkfemriPolyDataExactKSpaceGenerator::~vtkfemriPolyDataExactKSpaceGenerator()
{
}

int vtkfemriPolyDataExactKSpaceGenerator::FillInputPortInformation(int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
  return 1;
}

void vtkfemriPolyDataExactKSpaceGenerator::Initialize()
{
}

void vtkfemriPolyDataExactKSpaceGenerator::EvaluateFourierFunction(double frequency[3], double value[2])
{
  vtkDataSet* input = vtkDataSet::SafeDownCast(this->GetInput());

  vtkDataArray* normals = input->GetCellData()->GetNormals();

  if (!normals)
    {
    vtkErrorMacro("No surface normals found. Please provide a surface with outwards normals defined on cells.");
    }

  value[0] = 0.0;
  value[1] = 0.0;

  double pi = vtkMath::Pi();
  double twoPi = 2.0 * pi;
  double pi3 = pi * pi * pi;

  double k[3];
  k[0] = frequency[0];
  k[1] = frequency[1];
  k[2] = frequency[2];

  double k2 = this->Dot(k,k);

  int numberOfCells = input->GetNumberOfCells();
  int i, j;
  for (i=0; i<numberOfCells; i++)
    {
    vtkTriangle* triangle = vtkTriangle::SafeDownCast(input->GetCell(i));
    
    if (triangle == NULL)
      {
      //vtkErrorMacro("Error: cell not a vtkTriangle, skipping.");
      continue;
      }
 
    if (frequency[0] == 0.0 && frequency[1] == 0.0 && frequency[2] == 0.0)
      {
      value[0] += this->MagnetizationValue * triangle->ComputeArea();
      value[1] += 0.0;
      continue;
      }

    double A = triangle->ComputeArea();
 
    double triangleNormal[3];
    normals->GetTuple(i,triangleNormal);

    double kdotnt = vtkMath::Dot(frequency,triangleNormal);

    double x1[3], x2[3], x3[3];
    triangle->GetPoints()->GetPoint(0,x1);
    triangle->GetPoints()->GetPoint(1,x2);
    triangle->GetPoints()->GetPoint(2,x3);

    double cross[3];
    vtkTriangle::ComputeNormal(x1,x2,x3,cross);

    double normalDot = vtkMath::Dot(triangleNormal,cross);
    bool pointsReversed = normalDot < 0 ? true : false;

    if (pointsReversed)
      {
      triangle->GetPoints()->GetPoint(0,x1);
      triangle->GetPoints()->GetPoint(2,x2);
      triangle->GetPoints()->GetPoint(1,x3);
      }

    double phi_e = this->Dot(k,x3);
    double x1x3[3], x2x3[3], x1x2[3];
    this->Vector(x1,x3,x1x3);
    this->Vector(x2,x3,x2x3);
    this->Vector(x1,x2,x1x2);
    double k_r = this->Dot(k,x1x3);
    double k_s = this->Dot(k,x2x3);
    double k_rs = this->Dot(k,x1x2);
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

    value[0] += - 1.0 / (4.0 * pi3 * k2) * A * k_dot_n / k_s * ((sin(twoPi * (phi_e + k_r)) - sin(twoPi * (phi_e + k_s))) / k_rs - (sin(twoPi * (phi_e + k_r)) - sin(twoPi * (phi_e))) / k_r );
    value[1] += - 1.0 / (4.0 * pi3 * k2) * A * k_dot_n / k_s * ((cos(twoPi * (phi_e + k_r)) - cos(twoPi * (phi_e + k_s))) / k_rs - (cos(twoPi * (phi_e + k_r)) - cos(twoPi * (phi_e))) / k_r );
    }
}

