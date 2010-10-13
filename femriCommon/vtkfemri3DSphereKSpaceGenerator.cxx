/*=========================================================================

  Program:   FEMRI
  Module:    $RCSfile: vtkfemri3DSphereKSpaceGenerator.cxx,v $
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

#include "vtkfemri3DSphereKSpaceGenerator.h"
#include "vtkfemriKSpaceShift.h"
#include "vtkImageData.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkObjectFactory.h"

vtkStandardNewMacro(vtkfemri3DSphereKSpaceGenerator);
vtkCxxRevisionMacro(vtkfemri3DSphereKSpaceGenerator, "$Revision: 1.2 $");

vtkfemri3DSphereKSpaceGenerator::vtkfemri3DSphereKSpaceGenerator()
{
  this->Center[0] = this->Center[1] = this->Center[2] = 0.0;
  this->Radius = 0.0;
  this->KSpaceDimensionality = 3;
}

vtkfemri3DSphereKSpaceGenerator::~vtkfemri3DSphereKSpaceGenerator()
{
}

void vtkfemri3DSphereKSpaceGenerator::EvaluateFourierFunction(double frequency[3], double value[2])
{
  value[0] = value[1] = 0.0;
  double rho = 2.0 * vtkMath::Pi() * sqrt(frequency[0]*frequency[0] + frequency[1]*frequency[1] + frequency[2]*frequency[2]) * this->Radius;
  if (rho == 0.0)
    {
    value[0] = 1.0;
    }
  else
    {
    value[0] = 3.0 * (sin(rho) - rho * cos(rho)) / (rho*rho*rho);
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

void vtkfemri3DSphereKSpaceGenerator::PrintSelf(ostream& os, vtkIndent indent)
{
  Superclass::PrintSelf(os,indent);
}
