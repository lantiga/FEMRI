/*=========================================================================

  Program:   FEMRI
  Module:    $RCSfile: vtkfemriPolyDataExactKSpaceGeneratorMonolithic.h,v $
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

// .NAME vtkfemriPolyDataExactKSpaceGeneratorMonolithic - ..
// .SECTION Description
// ..

#ifndef __vtkfemriPolyDataExactKSpaceGeneratorMonolithic_h
#define __vtkfemriPolyDataExactKSpaceGeneratorMonolithic_h

#include "vtkImageAlgorithm.h"
#include "vtkfemriNumericWin32Header.h"

class vtkPolyData;
class vtkCell;

class VTK_FEMRI_NUMERIC_EXPORT vtkfemriPolyDataExactKSpaceGeneratorMonolithic : public vtkImageAlgorithm
{
public:
  vtkTypeRevisionMacro(vtkfemriPolyDataExactKSpaceGeneratorMonolithic,vtkImageAlgorithm);
  static vtkfemriPolyDataExactKSpaceGeneratorMonolithic* New();

  void PrintSelf(ostream& os, vtkIndent indent);

  vtkSetMacro(KSpaceDimensionality,int);
  vtkGetMacro(KSpaceDimensionality,int);

  vtkSetVectorMacro(Matrix,int,3);
  vtkGetVectorMacro(Matrix,int,3);

  vtkSetVectorMacro(FOV,double,3);
  vtkGetVectorMacro(FOV,double,3);

  vtkSetVectorMacro(Origin,double,3);
  vtkGetVectorMacro(Origin,double,3);

  vtkSetMacro(AcquireSymmetricKSpace,int);
  vtkGetMacro(AcquireSymmetricKSpace,int);
  vtkBooleanMacro(AcquireSymmetricKSpace,int);

  vtkSetMacro(NormalizeKSpace,int);
  vtkGetMacro(NormalizeKSpace,int);
  vtkBooleanMacro(NormalizeKSpace,int);

  vtkSetMacro(MagnetizationValue,double);
  vtkGetMacro(MagnetizationValue,double);  

  virtual void Initialize() {}

  static double Dot(double a[3], double b[3]) 
  { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

  static double Vector(double a[3], double b[3], double c[3]) 
  { c[0] = a[0] - b[0]; c[1] = a[1] - b[1]; c[2] = a[2] - b[2]; }
 
protected:
  vtkfemriPolyDataExactKSpaceGeneratorMonolithic();
  ~vtkfemriPolyDataExactKSpaceGeneratorMonolithic();

  virtual double ComputeVoxelVolume();

  virtual int RequestInformation (vtkInformation *, vtkInformationVector**, vtkInformationVector *);
  virtual int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

  virtual int FillInputPortInformation(int, vtkInformation*);

  virtual double ComputeVolume()
  {
    return 1.0;
  }

  int Matrix[3];
  double FOV[3];
  double Origin[3];
  int KSpaceDimensionality;

  int AcquireSymmetricKSpace;
  int NormalizeKSpace;
  double MagnetizationValue;

private:
  vtkfemriPolyDataExactKSpaceGeneratorMonolithic(const vtkfemriPolyDataExactKSpaceGeneratorMonolithic&);  // Not implemented.
  void operator=(const vtkfemriPolyDataExactKSpaceGeneratorMonolithic&);  // Not implemented.
};

#endif
