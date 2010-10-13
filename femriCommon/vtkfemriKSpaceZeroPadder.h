/*=========================================================================

  Program:   FEMRI
  Module:    $RCSfile: vtkfemriKSpaceZeroPadder.h,v $
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

// .NAME vtkfemriKSpaceZeroPadder - ..
// .SECTION Description
// ..

#ifndef __vtkfemriKSpaceZeroPadder_h
#define __vtkfemriKSpaceZeroPadder_h

#include "vtkSimpleImageToImageFilter.h"
#include "vtkfemriCommonWin32Header.h"

class vtkImageData;

class VTK_FEMRI_COMMON_EXPORT vtkfemriKSpaceZeroPadder: public vtkSimpleImageToImageFilter
{
public:
  vtkTypeRevisionMacro(vtkfemriKSpaceZeroPadder,vtkSimpleImageToImageFilter);
  void PrintSelf(ostream& os, vtkIndent indent);
               
  static vtkfemriKSpaceZeroPadder *New();
        
  vtkSetVectorMacro(PadSize, int, 3);
  vtkGetVectorMacro(PadSize, int, 3);

  protected:
  vtkfemriKSpaceZeroPadder();
  ~vtkfemriKSpaceZeroPadder() {};     
  
  virtual int RequestInformation (vtkInformation *, vtkInformationVector**, vtkInformationVector *);
  virtual void SimpleExecute(vtkImageData* input, vtkImageData* output);

  int PadSize[3];

private:
  vtkfemriKSpaceZeroPadder(const vtkfemriKSpaceZeroPadder&);  // Not implemented.
  void operator=(const vtkfemriKSpaceZeroPadder&);  // Not implemented.
};

#endif
