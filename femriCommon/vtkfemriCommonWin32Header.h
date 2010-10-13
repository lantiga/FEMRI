/*=========================================================================

  Program:   FEMRI
  Module:    $RCSfile: vtkfemriCommonWin32Header.h,v $
  Language:  C++
  Date:      $Date: 2007/03/19 13:31:25 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) 1993-2002 Ken Martin, Will Schroeder, Bill Lorensen 
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkfemriCommonWin32Header - manage Windows system differences
// .SECTION Description
// The vtkfemriCommonWin32Header captures some system differences between Unix
// and Windows operating systems. 

#ifndef __vtkfemriCommonWin32Header_h
#define __vtkfemriCommonWin32Header_h

#include <vtkfemriCommonConfigure.h>

#if defined(WIN32) && !defined(FEMRI_STATIC)
#if defined(vtkfemriCommon_EXPORTS)
#define VTK_FEMRI_COMMON_EXPORT __declspec( dllexport ) 
#else
#define VTK_FEMRI_COMMON_EXPORT __declspec( dllimport ) 
#endif
#else
#define VTK_FEMRI_COMMON_EXPORT
#endif

#endif
