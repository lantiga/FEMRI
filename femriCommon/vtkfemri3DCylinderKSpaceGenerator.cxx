/*=========================================================================

  Program:   FEMRI
  Module:    $RCSfile: vtkfemri3DCylinderKSpaceGenerator.cxx,v $
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

#include "vtkfemri3DCylinderKSpaceGenerator.h"
#include "vtkfemri2DCylinderKSpaceGenerator.h"
#include "vtkfemri2DBoxKSpaceGenerator.h"
#include "vtkfemriKSpaceShift.h"
#include "vtkImageData.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkObjectFactory.h"

vtkStandardNewMacro(vtkfemri3DCylinderKSpaceGenerator);
vtkCxxRevisionMacro(vtkfemri3DCylinderKSpaceGenerator, "$Revision: 1.3 $");

vtkfemri3DCylinderKSpaceGenerator::vtkfemri3DCylinderKSpaceGenerator()
{
  this->Center[0] = this->Center[1] = this->Center[2] = 0.0;
  this->Radius = 0.0;
  this->Radii[0] = this->Radii[1] = 0.0;
  this->TiltingAngle = 0.0;
  this->Theta = 0.0;
  this->Length = 0.0;
  this->KSpaceDimensionality = 3;
}

vtkfemri3DCylinderKSpaceGenerator::~vtkfemri3DCylinderKSpaceGenerator()
{
}

void vtkfemri3DCylinderKSpaceGenerator::EvaluateFourierFunction(double frequency[3], double value[2])
{
  value[0] = value[1] = 0.0;

  double radius = this->Radius;
  double ratio = 1.0;
  if (this->Radii[0] > 0.0 && this->Radii[1] > 0.0)
  {
    radius = this->Radii[0];
    ratio = this->Radii[1] / this->Radii[0];
  }

  double ws =  cos(this->Theta) * frequency[0] + sin(this->Theta) * frequency[1];
  double wt = -sin(this->Theta) * frequency[0] + cos(this->Theta) * frequency[1];
  double wss = cos(this->TiltingAngle) * ws - sin(this->TiltingAngle) * frequency[2];
  double wu = sin(this->TiltingAngle) * ws + cos(this->TiltingAngle) * frequency[2];
//  value[0] = 4.0 * radius * radius * vtkfemri2DCylinderKSpaceGenerator::JincPi(2.0*radius * sqrt(wss * wss + ratio * ratio * wt * wt)) * this->Length * vtkfemri2DBoxKSpaceGenerator::SincPi(this->Length * wu);
  value[0] = vtkfemri2DCylinderKSpaceGenerator::JincPi(2.0*radius * sqrt(wss * wss + ratio * ratio * wt * wt)) * vtkfemri2DBoxKSpaceGenerator::SincPi(this->Length * wu);
  if ((this->Center[0] != 0.0) || (this->Center[1] != 0.0) || (this->Center[2] != 0.0))
    {
    double translation[3];
    translation[0] = this->Center[0];
    translation[1] = this->Center[1];
    translation[2] = this->Center[2];
    double shiftedValue[3];
    vtkfemriKSpaceShift::ShiftPhase(value,frequency,translation,shiftedValue);
    value[0] = shiftedValue[0];
    value[1] = shiftedValue[1];
    }
}

void vtkfemri3DCylinderKSpaceGenerator::PrintSelf(ostream& os, vtkIndent indent)
{
  Superclass::PrintSelf(os,indent);
}
