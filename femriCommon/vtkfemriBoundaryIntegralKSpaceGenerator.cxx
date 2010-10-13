/*=========================================================================

  Program:   FEMRI
  Module:    $RCSfile: vtkfemriBoundaryIntegralKSpaceGenerator.cxx,v $
  Language:  C++
  Date:      $Date: 2007/03/19 13:31:25 $
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

#include "vtkfemriBoundaryIntegralKSpaceGenerator.h"
#include "vtkImageData.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkObjectFactory.h"

vtkStandardNewMacro(vtkfemriBoundaryIntegralKSpaceGenerator);
vtkCxxRevisionMacro(vtkfemriBoundaryIntegralKSpaceGenerator, "$Revision: 1.1.1.1 $");

vtkfemriBoundaryIntegralKSpaceGenerator::vtkfemriBoundaryIntegralKSpaceGenerator()
{
  this->BoundarySurface = NULL;
}

vtkfemriBoundaryIntegralKSpaceGenerator::~vtkfemriBoundaryIntegralKSpaceGenerator()
{
  if (this->BoundarySurface)
  {
    this->BoundarySurface->Delete();
    this->BoundarySurface = NULL;
  }
}

void vtkfemriBoundaryIntegralKSpaceGenerator::EvaluateFourierFunction(double frequency[3], double value[2])
{
  value[0] = value[1] = 0.0;
}

void vtkfemriBoundaryIntegralKSpaceGenerator::PrintSelf(ostream& os, vtkIndent indent)
{
  Superclass::PrintSelf(os,indent);
}
