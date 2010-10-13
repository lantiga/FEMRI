/*=========================================================================

  Program:   FEMRI
  Module:    $RCSfile: vtkfemri3DSphereKSpaceGenerator.h,v $
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

// .NAME vtkfemri3DSphereKSpaceGenerator - ..
// .SECTION Description
// ..

#ifndef __vtkfemri3DSphereKSpaceGenerator_h
#define __vtkfemri3DSphereKSpaceGenerator_h

#include "vtkfemriKSpaceGenerator.h"
#include "vtkfemriCommonWin32Header.h"

#include "vtkMath.h"

class VTK_FEMRI_COMMON_EXPORT vtkfemri3DSphereKSpaceGenerator : public vtkfemriKSpaceGenerator
{
public:
  vtkTypeRevisionMacro(vtkfemri3DSphereKSpaceGenerator,vtkfemriKSpaceGenerator);
  void PrintSelf(ostream& os, vtkIndent indent);
  static vtkfemri3DSphereKSpaceGenerator* New();

  vtkSetVectorMacro(Center,double,3);
  vtkGetVectorMacro(Center,double,3);

  vtkSetMacro(Radius,double);
  vtkGetMacro(Radius,double);
  
  virtual void EvaluateFourierFunction(double frequency[3], double value[2]);
  
protected:
  vtkfemri3DSphereKSpaceGenerator();
  ~vtkfemri3DSphereKSpaceGenerator();

  virtual double ComputeVolume()
  {
    return 4.0 / 3.0 * vtkMath::Pi() * this->Radius * this->Radius * this->Radius;
  }

  double Center[3];
  double Radius;

private:
  vtkfemri3DSphereKSpaceGenerator(const vtkfemri3DSphereKSpaceGenerator&);  // Not implemented.
  void operator=(const vtkfemri3DSphereKSpaceGenerator&);  // Not implemented.
};

#endif
