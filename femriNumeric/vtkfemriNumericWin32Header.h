/*=========================================================================

  Program:   FEMRI
  Module:    $RCSfile: vtkfemriNumericWin32Header.h,v $
  Language:  C++
  Date:      $Date: 2007/03/19 13:31:26 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) 1993-2002 Ken Martin, Will Schroeder, Bill Lorensen 
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkfemriNumericWin32Header - manage Windows system differences
// .SECTION Description
// The vtkfemriNumericWin32Header captures some system differences between Unix
// and Windows operating systems. 

#ifndef __vtkfemriNumericWin32Header_h
#define __vtkfemriNumericWin32Header_h

#include <vtkfemriNumericConfigure.h>

#if defined(WIN32) && !defined(FEMRI_STATIC)
#if defined(vtkfemriNumeric_EXPORTS)
#define VTK_FEMRI_NUMERIC_EXPORT __declspec( dllexport ) 
#else
#define VTK_FEMRI_NUMERIC_EXPORT __declspec( dllimport ) 
#endif
#else
#define VTK_FEMRI_NUMERIC_EXPORT
#endif

#endif
