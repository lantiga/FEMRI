/*=========================================================================

  Program:   femri
  Module:    $RCSfile: vtkfemriUnstructuredGridKSpaceGenerator.h,v $
  Language:  C++
  Date:      $Date: 2007/03/19 13:31:27 $
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

// .NAME vtkfemriUnstructuredGridToImageFilter - ..
// .SECTION Description
// ..
// .SECTION Thanks
// This class was developed by Luca Antiga, PhD \n
// Imaging Research Labs - Robarts Research Institute \n
// Bioengineering Dept - Mario Negri Institute for Pharmacological Research \n 
// email: lantiga@imaging.robarts.ca - homepage: http://www.imaging.robarts.ca/~lantiga

#ifndef __vtkfemriUnstructuredGridKSpaceGenerator_h
#define __vtkfemriUnstructuredGridKSpaceGenerator_h

#include "vtkImageAlgorithm.h"
#include "vtkfemriSymbolicWin32Header.h"
#include "ginac.h"

class vtkUnstructuredGrid;

class VTK_FEMRI_SYMBOLIC_EXPORT vtkfemriUnstructuredGridKSpaceGenerator : public vtkImageAlgorithm
{
public:
  vtkTypeRevisionMacro(vtkfemriUnstructuredGridKSpaceGenerator,vtkImageAlgorithm);
  static vtkfemriUnstructuredGridKSpaceGenerator* New();

  vtkSetStringMacro(FunctionString);
  vtkGetStringMacro(FunctionString);

  vtkSetMacro(KSpaceDimensionality,int);
  vtkGetMacro(KSpaceDimensionality,int);

  vtkSetVector3Macro(Matrix,int);
  vtkGetVector3Macro(Matrix,int);

  vtkSetVector3Macro(FOV,double);
  vtkGetVector3Macro(FOV,double);

  vtkSetVector3Macro(Origin,double);
  vtkGetVector3Macro(Origin,double);

//   vtkSetMacro(UsePreComputedSymbolicIntegral,int);
//   vtkGetMacro(UsePreComputedSymbolicIntegral,int);
//   vtkBooleanMacro(UsePreComputedSymbolicIntegral,int);

  vtkSetMacro(UseElementIntegralCache,int);
  vtkGetMacro(UseElementIntegralCache,int);
  vtkBooleanMacro(UseElementIntegralCache,int);

  vtkSetMacro(AcquireSymmetricKSpace,int);
  vtkGetMacro(AcquireSymmetricKSpace,int);
  vtkBooleanMacro(AcquireSymmetricKSpace,int);

protected:
  vtkfemriUnstructuredGridKSpaceGenerator();
  ~vtkfemriUnstructuredGridKSpaceGenerator();

  virtual int FillInputPortInformation(int, vtkInformation*);
  virtual int RequestInformation (vtkInformation *, 
                                  vtkInformationVector **, 
                                  vtkInformationVector *);
  virtual int RequestData (vtkInformation *, 
                           vtkInformationVector **, 
                           vtkInformationVector *);

  char* FunctionString;

  int AcquireSymmetricKSpace;

  int KSpaceDimensionality;
  int Matrix[3];

  double FOV[3];
  double Origin[3];

//BTX
  GiNaC::ex Function;
//ETX

//   int UsePreComputedSymbolicIntegral;
  int UseElementIntegralCache;

private:
  vtkfemriUnstructuredGridKSpaceGenerator(const vtkfemriUnstructuredGridKSpaceGenerator&);  // Not implemented.
  void operator=(const vtkfemriUnstructuredGridKSpaceGenerator&);  // Not implemented.
};

#endif
