/*=========================================================================

  Program:   FEMRI
  Module:    $RCSfile: vtkfemri3DCylinderKSpaceGenerator.h,v $
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

// .NAME vtkfemri3DCylinderKSpaceGenerator - ..
// .SECTION Description
// ..

#ifndef __vtkfemri3DCylinderKSpaceGenerator_h
#define __vtkfemri3DCylinderKSpaceGenerator_h

#include "vtkfemriKSpaceGenerator.h"
#include "vtkfemriCommonWin32Header.h"

#include "vtkMath.h"

class VTK_FEMRI_COMMON_EXPORT vtkfemri3DCylinderKSpaceGenerator : public vtkfemriKSpaceGenerator
{
public:
  vtkTypeRevisionMacro(vtkfemri3DCylinderKSpaceGenerator,vtkfemriKSpaceGenerator);
  void PrintSelf(ostream& os, vtkIndent indent);
  static vtkfemri3DCylinderKSpaceGenerator* New();

  vtkSetVectorMacro(Center,double,3);
  vtkGetVectorMacro(Center,double,3);

  vtkSetMacro(Radius,double);
  vtkGetMacro(Radius,double);
 
  vtkSetVectorMacro(Radii,double,2);
  vtkGetVectorMacro(Radii,double,2);
  
  vtkSetMacro(TiltingAngle,double);
  vtkGetMacro(TiltingAngle,double);

  vtkSetMacro(Theta,double);
  vtkGetMacro(Theta,double);

  vtkSetMacro(Length,double);
  vtkGetMacro(Length,double);

  virtual void EvaluateFourierFunction(double frequency[3], double value[2]);
  
protected:
  vtkfemri3DCylinderKSpaceGenerator();
  ~vtkfemri3DCylinderKSpaceGenerator();

  virtual double ComputeVolume()
  {
    return  vtkMath::Pi() * this->Radii[0] * this->Radii[1] * this->Length;
  }

  double Center[3];
  double Radius;
  double Radii[2];
  double TiltingAngle;
  double Theta;
  double Length;

private:
  vtkfemri3DCylinderKSpaceGenerator(const vtkfemri3DCylinderKSpaceGenerator&);  // Not implemented.
  void operator=(const vtkfemri3DCylinderKSpaceGenerator&);  // Not implemented.
};

#endif
