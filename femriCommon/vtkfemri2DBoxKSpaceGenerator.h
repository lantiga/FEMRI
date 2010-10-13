/*=========================================================================

  Program:   FEMRI
  Module:    $RCSfile: vtkfemri2DBoxKSpaceGenerator.h,v $
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

// .NAME vtkfemri2DBoxKSpaceGenerator - ..
// .SECTION Description
// ..

#ifndef __vtkfemri2DBoxKSpaceGenerator_h
#define __vtkfemri2DBoxKSpaceGenerator_h

#include "vtkfemriKSpaceGenerator.h"
#include "vtkfemriCommonWin32Header.h"

#include "vtkMath.h"

class VTK_FEMRI_COMMON_EXPORT vtkfemri2DBoxKSpaceGenerator : public vtkfemriKSpaceGenerator
{
public:
  vtkTypeRevisionMacro(vtkfemri2DBoxKSpaceGenerator,vtkfemriKSpaceGenerator);
  void PrintSelf(ostream& os, vtkIndent indent);
  static vtkfemri2DBoxKSpaceGenerator* New();

  vtkSetVectorMacro(Center,double,3);
  vtkGetVectorMacro(Center,double,3);

  vtkSetMacro(Thickness,double);
  vtkGetMacro(Thickness,double);
 
  vtkSetMacro(Length,double);
  vtkGetMacro(Length,double);
  
  vtkSetMacro(TiltingAngle,double);
  vtkGetMacro(TiltingAngle,double);

  vtkSetMacro(Theta,double);
  vtkGetMacro(Theta,double);

  vtkSetMacro(SliceThickness,double);
  vtkGetMacro(SliceThickness,double);

  static double SincPi(double x)
    {
    if (x==0.0)
      {
      return 1.0;
      }
    return sin(vtkMath::Pi()*x) / (vtkMath::Pi()*x);
    }

  virtual void EvaluateFourierFunction(double frequency[3], double value[2]);
  
protected:
  vtkfemri2DBoxKSpaceGenerator();
  ~vtkfemri2DBoxKSpaceGenerator();

  virtual double ComputeVolume()
  {
    return this->Thickness / cos(this->TiltingAngle) * this->Length;
  }

  double Center[3];
  double Thickness;
  double Length;
  double TiltingAngle;
  double Theta;
  double SliceThickness;

private:
  vtkfemri2DBoxKSpaceGenerator(const vtkfemri2DBoxKSpaceGenerator&);  // Not implemented.
  void operator=(const vtkfemri2DBoxKSpaceGenerator&);  // Not implemented.
};

#endif
