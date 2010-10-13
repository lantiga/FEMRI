/*=========================================================================

  Program:   FEMRI
  Module:    $RCSfile: vtkfemri2DCylinderKSpaceGenerator.h,v $
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

// .NAME vtkfemri2DCylinderKSpaceGenerator - ..
// .SECTION Description
// ..

#ifndef __vtkfemri2DCylinderKSpaceGenerator_h
#define __vtkfemri2DCylinderKSpaceGenerator_h

#include "vtkfemriKSpaceGenerator.h"
#include "vtkfemriCommonWin32Header.h"

#include "vtkMath.h"

class VTK_FEMRI_COMMON_EXPORT vtkfemri2DCylinderKSpaceGenerator : public vtkfemriKSpaceGenerator
{
public:
  vtkTypeRevisionMacro(vtkfemri2DCylinderKSpaceGenerator,vtkfemriKSpaceGenerator);
  void PrintSelf(ostream& os, vtkIndent indent);
  static vtkfemri2DCylinderKSpaceGenerator* New();

  vtkSetVectorMacro(Center,double,3);
  vtkGetVectorMacro(Center,double,3);

  vtkSetMacro(Radius,double);
  vtkGetMacro(Radius,double);
  
  vtkSetMacro(TiltingAngle,double);
  vtkGetMacro(TiltingAngle,double);

  vtkSetMacro(Theta,double);
  vtkGetMacro(Theta,double);

  vtkSetMacro(SliceThickness,double);
  vtkGetMacro(SliceThickness,double);

  static double Jinc(double x)
    {
    if (x==0.0)
      {
      return 1.0;
      }
    return 2.0 * j1(x) / x;
    }
 
  static double JincPi(double x)
    {
    if (x==0.0)
      {
      return 1.0;
      }
    return 2.0 * j1(vtkMath::Pi()*x) / (vtkMath::Pi()*x);
    }
 
  virtual void EvaluateFourierFunction(double frequency[3], double value[2]);
  
protected:
  vtkfemri2DCylinderKSpaceGenerator();
  ~vtkfemri2DCylinderKSpaceGenerator();

  virtual double ComputeVolume()
  {
    return vtkMath::Pi() * this->Radius * this->Radius / cos(this->TiltingAngle);
  }

  double Center[3];
  double Radius;
  double TiltingAngle;
  double Theta;
  double SliceThickness;

private:
  vtkfemri2DCylinderKSpaceGenerator(const vtkfemri2DCylinderKSpaceGenerator&);  // Not implemented.
  void operator=(const vtkfemri2DCylinderKSpaceGenerator&);  // Not implemented.
};

#endif
