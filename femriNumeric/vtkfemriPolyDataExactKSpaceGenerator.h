/*=========================================================================

  Program:   femri
  Module:    vtkfemriPolyDataExactKSpaceGenerator.h
  Language:  C++
  Date:      $Date: 2006/04/06 16:46:43 $
  Version:   $Revision: 1.4 $

  Copyright (c) Luca Antiga, David Steinman. All rights reserved.
  See LICENCE file for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm 
  for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

// .NAME vtkfemriPolyDataExactKSpaceGenerator - ..
// .SECTION Description
// ..

#ifndef __vtkfemriPolyDataExactKSpaceGenerator_h
#define __vtkfemriPolyDataExactKSpaceGenerator_h

#include "vtkfemriKSpaceGenerator.h"
#include "vtkfemriNumericWin32Header.h"

class vtkPolyData;
class vtkCell;

class VTK_FEMRI_NUMERIC_EXPORT vtkfemriPolyDataExactKSpaceGenerator : public vtkfemriKSpaceGenerator
{
public:
  vtkTypeRevisionMacro(vtkfemriPolyDataExactKSpaceGenerator,vtkImageAlgorithm);
  static vtkfemriPolyDataExactKSpaceGenerator* New();

  vtkSetMacro(MagnetizationValue,double);
  vtkGetMacro(MagnetizationValue,double);  

  virtual void Initialize();

  virtual void EvaluateFourierFunction(double frequency[3], double value[2]);

  static double Dot(double a[3], double b[3]) 
  { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

  static double Vector(double a[3], double b[3], double c[3]) 
  { c[0] = a[0] - b[0]; c[1] = a[1] - b[1]; c[2] = a[2] - b[2]; }

protected:
  vtkfemriPolyDataExactKSpaceGenerator();
  ~vtkfemriPolyDataExactKSpaceGenerator();
  
  virtual int FillInputPortInformation(int, vtkInformation*);

  virtual double ComputeVolume()
  {
    return 1.0;
  }

  double MagnetizationValue;

private:
  vtkfemriPolyDataExactKSpaceGenerator(const vtkfemriPolyDataExactKSpaceGenerator&);  // Not implemented.
  void operator=(const vtkfemriPolyDataExactKSpaceGenerator&);  // Not implemented.
};

#endif
