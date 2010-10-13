/*=========================================================================

  Program:   FEMRI
  Module:    $RCSfile: vtkfemri3DBoxKSpaceGenerator.cxx,v $
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

#include "vtkfemri3DBoxKSpaceGenerator.h"
#include "vtkfemri2DBoxKSpaceGenerator.h"
#include "vtkfemriKSpaceShift.h"
#include "vtkImageData.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkObjectFactory.h"

vtkStandardNewMacro(vtkfemri3DBoxKSpaceGenerator);
vtkCxxRevisionMacro(vtkfemri3DBoxKSpaceGenerator, "$Revision: 1.2 $");

vtkfemri3DBoxKSpaceGenerator::vtkfemri3DBoxKSpaceGenerator()
{
  this->Center[0] = this->Center[1] = this->Center[2] = 0.0;
  this->Lengths[0] = this->Lengths[1] = this->Lengths[2] = 0.0;
  this->TiltingAngle = 0.0;
  this->Theta = 0.0;
  this->KSpaceDimensionality = 3;
}

vtkfemri3DBoxKSpaceGenerator::~vtkfemri3DBoxKSpaceGenerator()
{
}

void vtkfemri3DBoxKSpaceGenerator::EvaluateFourierFunction(double frequency[3], double value[2])
{
  value[0] = value[1] = 0.0;
  double ws =  cos(this->Theta) * frequency[0] + sin(this->Theta) * frequency[1];
  double wt = -sin(this->Theta) * frequency[0] + cos(this->Theta) * frequency[1];
  double wss = cos(this->TiltingAngle) * ws - sin(this->TiltingAngle) * frequency[2];
  double wu = sin(this->TiltingAngle) * ws + cos(this->TiltingAngle) * frequency[2];
//  value[0]  = this->Lengths[0] * vtkfemri2DBoxKSpaceGenerator::SincPi(this->Lengths[0] * wss);
//  value[0] *= this->Lengths[1] * vtkfemri2DBoxKSpaceGenerator::SincPi(this->Lengths[1] * wt);
//  value[0] *= this->Lengths[2] * vtkfemri2DBoxKSpaceGenerator::SincPi(this->Lengths[2] * wu);
  value[0]  = vtkfemri2DBoxKSpaceGenerator::SincPi(this->Lengths[0] * wss) * vtkfemri2DBoxKSpaceGenerator::SincPi(this->Lengths[1] * wt) * vtkfemri2DBoxKSpaceGenerator::SincPi(this->Lengths[2] * wu);
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

void vtkfemri3DBoxKSpaceGenerator::PrintSelf(ostream& os, vtkIndent indent)
{
  Superclass::PrintSelf(os,indent);
}
