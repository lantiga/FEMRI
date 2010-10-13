/*=========================================================================

  Program:   FEMRI
  Module:    $RCSfile: vtkfemriKSpaceShift.h,v $
  Language:  C++
  Date:      $Date: 2007/03/19 13:31:25 $
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

// .NAME vtkfemriKSpaceShift - ..
// .SECTION Description
// ..

#ifndef __vtkfemriKSpaceShift_h
#define __vtkfemriKSpaceShift_h

#include "vtkImageAlgorithm.h"
#include "vtkfemriCommonWin32Header.h"

class VTK_FEMRI_COMMON_EXPORT vtkfemriKSpaceShift : public vtkImageAlgorithm
{
public:
  vtkTypeRevisionMacro(vtkfemriKSpaceShift,vtkImageAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  static vtkfemriKSpaceShift* New();

  vtkSetVectorMacro(Translation,double,3);
  vtkGetVectorMacro(Translation,double,3);

  static void ShiftPhase(double value[2], double frequency[3], double translation[3], double shiftedValue[2]);
  
protected:
  vtkfemriKSpaceShift();
  ~vtkfemriKSpaceShift();

  virtual int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

  double Translation[3];

private:
  vtkfemriKSpaceShift(const vtkfemriKSpaceShift&);  // Not implemented.
  void operator=(const vtkfemriKSpaceShift&);  // Not implemented.
};

#endif
