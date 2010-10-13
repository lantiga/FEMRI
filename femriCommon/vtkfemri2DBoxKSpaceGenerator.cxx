/*=========================================================================

  Program:   FEMRI
  Module:    $RCSfile: vtkfemri2DBoxKSpaceGenerator.cxx,v $
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

#include "vtkfemri2DBoxKSpaceGenerator.h"
#include "vtkfemriKSpaceShift.h"
#include "vtkImageData.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkObjectFactory.h"

vtkStandardNewMacro(vtkfemri2DBoxKSpaceGenerator);
vtkCxxRevisionMacro(vtkfemri2DBoxKSpaceGenerator, "$Revision: 1.2 $");

vtkfemri2DBoxKSpaceGenerator::vtkfemri2DBoxKSpaceGenerator()
{
  this->Center[0] = this->Center[1] = this->Center[2] = 0.0;
  this->Thickness = 0.0;
  this->Length = 0.0;
  this->TiltingAngle = 0.0;
  this->Theta = 0.0;
  this->KSpaceDimensionality = 2;
}

vtkfemri2DBoxKSpaceGenerator::~vtkfemri2DBoxKSpaceGenerator()
{
}

void vtkfemri2DBoxKSpaceGenerator::EvaluateFourierFunction(double frequency[3], double value[2])
{
  value[0] = value[1] = 0.0;
  double ws =  cos(this->Theta) * frequency[0] + sin(this->Theta) * frequency[1];
  double wt = -sin(this->Theta) * frequency[0] + cos(this->Theta) * frequency[1];
//  value[0] = this->Thickness / cos(this->TiltingAngle) * vtkfemri2DBoxKSpaceGenerator::SincPi(this->Thickness / cos(this->TiltingAngle) * ws) * this->Length * vtkfemri2DBoxKSpaceGenerator::SincPi(this->Length * wt);
  value[0] = vtkfemri2DBoxKSpaceGenerator::SincPi(this->Thickness / cos(this->TiltingAngle) * ws) * vtkfemri2DBoxKSpaceGenerator::SincPi(this->Length * wt);
  if (this->SliceThickness > 0.0 && this->TiltingAngle > 0.0)
    {
    double sigma = this->SliceThickness * tan(this->TiltingAngle);
//    value[0] *= sigma * vtkfemri2DBoxKSpaceGenerator::SincPi(sigma * ws);
    value[0] *= vtkfemri2DBoxKSpaceGenerator::SincPi(sigma * ws);
    }
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

void vtkfemri2DBoxKSpaceGenerator::PrintSelf(ostream& os, vtkIndent indent)
{
  Superclass::PrintSelf(os,indent);
}
