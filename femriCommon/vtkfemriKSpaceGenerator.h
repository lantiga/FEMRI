/*=========================================================================

  Program:   FEMRI
  Module:    $RCSfile: vtkfemriKSpaceGenerator.h,v $
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

// .NAME vtkfemriKSpaceGenerator - ..
// .SECTION Description
// ..

#ifndef __vtkfemriKSpaceGenerator_h
#define __vtkfemriKSpaceGenerator_h

#include "vtkImageAlgorithm.h"
#include "vtkfemriCommonWin32Header.h"

class VTK_FEMRI_COMMON_EXPORT vtkfemriKSpaceGenerator : public vtkImageAlgorithm
{
public:
  vtkTypeRevisionMacro(vtkfemriKSpaceGenerator,vtkImageAlgorithm);
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

  virtual void Initialize() {}

  virtual void EvaluateFourierFunction(double frequency[3], double value[2]) = 0;
  
protected:
  vtkfemriKSpaceGenerator();
  ~vtkfemriKSpaceGenerator();

  virtual double ComputeVolume() = 0;
  virtual double ComputeVoxelVolume();

  virtual int RequestInformation (vtkInformation *, vtkInformationVector**, vtkInformationVector *);
  virtual int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

  int Matrix[3];
  double FOV[3];
  double Origin[3];
  int KSpaceDimensionality;

  int AcquireSymmetricKSpace;
  int NormalizeKSpace;

private:
  vtkfemriKSpaceGenerator(const vtkfemriKSpaceGenerator&);  // Not implemented.
  void operator=(const vtkfemriKSpaceGenerator&);  // Not implemented.
};

#endif
