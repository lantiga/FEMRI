/*=========================================================================

  Program:   FEMRI
  Module:    $RCSfile: vtkfemriBoundaryIntegralKSpaceGenerator.h,v $
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

// .NAME vtkfemriBoundaryIntegralKSpaceGenerator - ..
// .SECTION Description
// ..

#ifndef __vtkfemriBoundaryIntegralKSpaceGenerator_h
#define __vtkfemriBoundaryIntegralKSpaceGenerator_h

#include "vtkfemriKSpaceGenerator.h"
#include "vtkfemriCommonWin32Header.h"
#include "vtkPolyData.h"

#include "vtkMath.h"

class VTK_FEMRI_COMMON_EXPORT vtkfemriBoundaryIntegralKSpaceGenerator : public vtkfemriKSpaceGenerator
{
public:
  vtkTypeRevisionMacro(vtkfemriBoundaryIntegralKSpaceGenerator,vtkfemriKSpaceGenerator);
  void PrintSelf(ostream& os, vtkIndent indent);
  static vtkfemriBoundaryIntegralKSpaceGenerator* New();

  virtual void EvaluateFourierFunction(double frequency[3], double value[2]);
  
  vtkSetObjectMacro(BoundarySurface,vtkPolyData);
  vtkGetObjectMacro(BoundarySurface,vtkPolyData);
  
protected:
  vtkfemriBoundaryIntegralKSpaceGenerator();
  ~vtkfemriBoundaryIntegralKSpaceGenerator();

  virtual double ComputeVolume()
  {
    return 1.0;
  }

  vtkPolyData* BoundarySurface;

private:
  vtkfemriBoundaryIntegralKSpaceGenerator(const vtkfemriBoundaryIntegralKSpaceGenerator&);  // Not implemented.
  void operator=(const vtkfemriBoundaryIntegralKSpaceGenerator&);  // Not implemented.
};

#endif
