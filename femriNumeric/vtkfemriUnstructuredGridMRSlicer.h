/*=========================================================================

  Program:   femri
  Module:    $RCSfile: vtkfemriUnstructuredGridMRSlicer.h,v $
  Language:  C++
  Date:      $Date: 2007/03/19 13:31:26 $
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

// .NAME vtkfemriUnstructuredGridMRSlicer - ..
// .SECTION Description
// ..
// .SECTION Thanks
// This class was developed by Paul Simedrea\n
// Imaging Research Labs - Robarts Research Institute \n
// email: simedrea@imaging.robarts.ca - homepage: http://www.imaging.robarts.ca/~simedrea

#ifndef __vtkfemriUnstructuredGridMRSlicer_h
#define __vtkfemriUnstructuredGridMRSlicer_h

#include "vtkUnstructuredGridToUnstructuredGridFilter.h"
#include "vtkfemriNumericWin32Header.h"

class vtkUnstructuredGrid;

class VTK_FEMRI_NUMERIC_EXPORT vtkfemriUnstructuredGridMRSlicer : public vtkUnstructuredGridToUnstructuredGridFilter
{
public:
  vtkTypeRevisionMacro(vtkfemriUnstructuredGridMRSlicer,vtkUnstructuredGridToUnstructuredGridFilter);
  static vtkfemriUnstructuredGridMRSlicer* New();
  
  vtkSetMacro(SliceProfile,int);
  vtkGetMacro(SliceProfile,int);
  
  vtkSetMacro(SliceThickness,double);
  vtkGetMacro(SliceThickness,double);
    
  // perhaps for eventual multi-slice acquisition? (Or can run mesh again through code)  

  //vtkSetMacro(SliceSpacing,int);
  //vtkGetMacro(SliceSpacing,int);

  vtkSetVector3Macro(SliceOrigin,double);
  vtkGetVector3Macro(SliceOrigin,double);

  vtkSetVector3Macro(SliceOrientation,double);
  vtkGetVector3Macro(SliceOrientation,double);

  vtkSetMacro(Residual, double);
  vtkGetMacro(Residual, double);
  
protected:
  vtkfemriUnstructuredGridMRSlicer();
  ~vtkfemriUnstructuredGridMRSlicer();

  void ExecuteInformation();
  void ExecuteData(vtkDataObject *outp);
  
  // function for subdividing tets near the slice boundary
  void SubdivideTetsNearSliceBoundary();
  
  // function applying magnetization profile after subdivision
  void AssignSliceProfile();

  // trapezoid slice profile
  double TrapezoidSlice(double z);

  double TrapezoidResidual(double z_min, double z_max);


  int NumSliceBoundarySubdivisions;
  int SliceProfile;
  
  double Residual;
  double SliceThickness;  
  double SliceOrigin[3];
  double SliceOrientation[3];
    
private:
  vtkfemriUnstructuredGridMRSlicer(const vtkfemriUnstructuredGridMRSlicer&);  // Not implemented.
  void operator=(const vtkfemriUnstructuredGridMRSlicer&);  // Not implemented.
};

#endif
