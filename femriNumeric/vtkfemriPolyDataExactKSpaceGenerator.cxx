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
  this->CellNormalsArrayName = NULL;

  this->SetNumberOfInputPorts(1);
}

vtkfemriPolyDataExactKSpaceGenerator::~vtkfemriPolyDataExactKSpaceGenerator()
{
  if (this->CellNormalsArrayName)
    {
    delete[] this->CellNormalsArrayName;
    this->CellNormalsArrayName = NULL;
    }
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

  if (this->CellNormalsArrayName)
    {
    normals = input->GetCellData()->GetArray(this->CellNormalsArrayName);
    }

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
      vtkErrorMacro("Error: cell not a vtkTriangle, skipping.");
      continue;
      }
 
    double A = triangle->ComputeArea();
 
    double triangleNormal[3];
    normals->GetTuple(i,triangleNormal);

    double x1[3], x2[3], x3[3];
    triangle->GetPoints()->GetPoint(0,x1);
    triangle->GetPoints()->GetPoint(1,x2);
    triangle->GetPoints()->GetPoint(2,x3);

    if (frequency[0] == 0.0 && frequency[1] == 0.0 && frequency[2] == 0.0)
      {
      value[0] += this->MagnetizationValue * A * triangleNormal[2] * (x1[2]+x2[2]+x3[2]) / 3.0;
      value[1] += 0.0;
      continue;
      }

    double kdotnt = vtkMath::Dot(frequency,triangleNormal);

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

    double kr = k_r;
    double ks = k_s;
    double kdotn = k_dot_n;
    double phi = phi_e;

    //if (fabs(kdotn) < tol) 
    //  {
    //  continue;
    //  }

    if (fabs(k_s) < tol) 
      {
      if (fabs(k_r) < tol)
        {
        value[0] += A*kdotn*sin(2*pi*phi)/(2*pi*k2);
        value[1] += A*kdotn*cos(2*pi*phi)/(2*pi*k2);
        continue;
        }
      value[0] += A*kdotn*cos(2*pi*phi)/(2*pi*pi*k2*kr) + A*kdotn*sin(2*pi*phi)/(4*pi*pi*pi*k2*kr*kr) - A*kdotn*cos(2*pi*kr)*sin(2*pi*phi)/(4*pi*pi*pi*k2*kr*kr) - A*kdotn*cos(2*pi*phi)*sin(2*pi*kr)/(4*pi*pi*pi*k2*kr*kr);
      value[1] += -A*kdotn*sin(2*pi*phi)/(2*pi*pi*k2*kr) + A*kdotn*cos(2*pi*phi)/(4*pi*pi*pi*k2*kr*kr) - A*kdotn*cos(2*pi*kr)*cos(2*pi*phi)/(4*pi*pi*pi*k2*kr*kr) + A*kdotn*sin(2*pi*kr)*sin(2*pi*phi)/(4*pi*pi*pi*k2*kr*kr);
      continue;
      }

    if (fabs(k_r) < tol) 
      {
      value[0] += A*kdotn*cos(2*pi*phi)/(2*pi*pi*k2*ks) + A*kdotn*sin(2*pi*phi)/(4*pi*pi*pi*k2*ks*ks) - A*kdotn*cos(2*pi*ks)*sin(2*pi*phi)/(4*pi*pi*pi*k2*ks*ks) - A*kdotn*cos(2*pi*phi)*sin(2*pi*ks)/(4*pi*pi*pi*k2*ks*ks);
      value[1] += -A*kdotn*sin(2*pi*phi)/(2*pi*pi*k2*ks) + A*kdotn*cos(2*pi*phi)/(4*pi*pi*pi*k2*ks*ks) - A*kdotn*cos(2*pi*ks)*cos(2*pi*phi)/(4*pi*pi*pi*k2*ks*ks) + A*kdotn*sin(2*pi*ks)*sin(2*pi*phi)/(4*pi*pi*pi*k2*ks*ks);
      continue;
      }

    if (fabs(k_rs) < tol) 
      {
      value[0] += -A*kdotn*sin(2*pi*phi)/(4*pi*pi*pi*k2*ks*ks) + A*kdotn*sin(2*pi*ks)*sin(2*pi*phi)/(2*pi*pi*k2*ks) - A*kdotn*cos(2*pi*ks)*cos(2*pi*phi)/(2*pi*pi*k2*ks) + A*kdotn*cos(2*pi*ks)*sin(2*pi*phi)/(4*pi*pi*pi*k2*ks*ks) + A*kdotn*cos(2*pi*phi)*sin(2*pi*ks)/(4*pi*pi*pi*k2*ks*ks);
      value[1] += -A*kdotn*cos(2*pi*phi)/(4*pi*pi*pi*k2*ks*ks) + A*kdotn*cos(2*pi*ks)*sin(2*pi*phi)/(2*pi*pi*k2*ks) + A*kdotn*cos(2*pi*phi)*sin(2*pi*ks)/(2*pi*pi*k2*ks) - A*kdotn*sin(2*pi*ks)*sin(2*pi*phi)/(4*pi*pi*pi*k2*ks*ks) + A*kdotn*cos(2*pi*ks)*cos(2*pi*phi)/(4*pi*pi*pi*k2*ks*ks);
      continue;
      }

    value[0] += -A*kdotn*sin(2*pi*phi)/(4*pi*pi*pi*k2*kr*ks) - A*kdotn*cos(2*pi*kr)*sin(2*pi*phi)/(4*pi*pi*pi*k2*ks*(kr - ks)) - A*kdotn*cos(2*pi*phi)*sin(2*pi*kr)/(4*pi*pi*pi*k2*ks*(kr - ks)) + A*kdotn*cos(2*pi*kr)*sin(2*pi*phi)/(4*pi*pi*pi*k2*kr*ks) + A*kdotn*cos(2*pi*phi)*sin(2*pi*kr)/(4*pi*pi*pi*k2*kr*ks) + A*kdotn*cos(2*pi*ks)*sin(2*pi*phi)/(4*pi*pi*pi*k2*ks*(kr - ks)) + A*kdotn*cos(2*pi*phi)*sin(2*pi*ks)/(4*pi*pi*pi*k2*ks*(kr - ks));
    value[1] += -A*kdotn*cos(2*pi*phi)/(4*pi*pi*pi*k2*kr*ks) - A*kdotn*sin(2*pi*kr)*sin(2*pi*phi)/(4*pi*pi*pi*k2*kr*ks) - A*kdotn*cos(2*pi*kr)*cos(2*pi*phi)/(4*pi*pi*pi*k2*ks*(kr - ks)) - A*kdotn*sin(2*pi*ks)*sin(2*pi*phi)/(4*pi*pi*pi*k2*ks*(kr - ks)) + A*kdotn*cos(2*pi*kr)*cos(2*pi*phi)/(4*pi*pi*pi*k2*kr*ks) + A*kdotn*cos(2*pi*ks)*cos(2*pi*phi)/(4*pi*pi*pi*k2*ks*(kr - ks)) + A*kdotn*sin(2*pi*kr)*sin(2*pi*phi)/(4*pi*pi*pi*k2*ks*(kr - ks));
    }
}

