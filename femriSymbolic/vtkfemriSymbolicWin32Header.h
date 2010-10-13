/*=========================================================================

  Program:   FEMRI
  Module:    $RCSfile: vtkfemriSymbolicWin32Header.h,v $
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
// .NAME vtkfemriSymbolicWin32Header - manage Windows system differences
// .SECTION Description
// The vtkfemriSymbolicWin32Header captures some system differences between Unix
// and Windows operating systems. 

#ifndef __vtkfemriSymbolicWin32Header_h
#define __vtkfemriSymbolicWin32Header_h

#include <vtkfemriSymbolicConfigure.h>

#if defined(WIN32) && !defined(FEMRI_STATIC)
#if defined(vtkfemriSymbolic_EXPORTS)
#define VTK_FEMRI_SYMBOLIC_EXPORT __declspec( dllexport ) 
#else
#define VTK_FEMRI_SYMBOLIC_EXPORT __declspec( dllimport ) 
#endif
#else
#define VTK_FEMRI_SYMBOLIC_EXPORT
#endif

#endif
