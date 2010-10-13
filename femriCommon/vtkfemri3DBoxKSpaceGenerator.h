/*=========================================================================

  Program:   FEMRI
  Module:    $RCSfile: vtkfemri3DBoxKSpaceGenerator.h,v $
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

// .NAME vtkfemri3DBoxKSpaceGenerator - ..
// .SECTION Description
// ..

#ifndef __vtkfemri3DBoxKSpaceGenerator_h
#define __vtkfemri3DBoxKSpaceGenerator_h

#include "vtkfemriKSpaceGenerator.h"
#include "vtkfemriCommonWin32Header.h"

#include "vtkMath.h"

class VTK_FEMRI_COMMON_EXPORT vtkfemri3DBoxKSpaceGenerator : public vtkfemriKSpaceGenerator
{
public:
  vtkTypeRevisionMacro(vtkfemri3DBoxKSpaceGenerator,vtkfemriKSpaceGenerator);
  void PrintSelf(ostream& os, vtkIndent indent);
  static vtkfemri3DBoxKSpaceGenerator* New();

  vtkSetVectorMacro(Center,double,3);
  vtkGetVectorMacro(Center,double,3);

  vtkSetVectorMacro(Lengths,double,3);
  vtkGetVectorMacro(Lengths,double,3);
  
  vtkSetMacro(TiltingAngle,double);
  vtkGetMacro(TiltingAngle,double);

  vtkSetMacro(Theta,double);
  vtkGetMacro(Theta,double);

  virtual void EvaluateFourierFunction(double frequency[3], double value[2]);
  
protected:
  vtkfemri3DBoxKSpaceGenerator();
  ~vtkfemri3DBoxKSpaceGenerator();

  virtual double ComputeVolume()
  {
    return this->Lengths[0] * this->Lengths[1] * this->Lengths[2];
  }

  double Center[3];
  double TiltingAngle;
  double Theta;
  double Lengths[3];

private:
  vtkfemri3DBoxKSpaceGenerator(const vtkfemri3DBoxKSpaceGenerator&);  // Not implemented.
  void operator=(const vtkfemri3DBoxKSpaceGenerator&);  // Not implemented.
};

#endif
