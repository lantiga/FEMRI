/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkSubdivideSelectTetra.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkSubdivideSelectTetra - subdivide one tetrahedron into twelve for every tetra
// .SECTION Description
// This filter subdivides tetrahedra in an unstructured grid into twelve tetrahedra.


#ifndef __vtkSubdivideSelectTetra_h
#define __vtkSubdivideSelectTetra_h

#include "vtkUnstructuredGridToUnstructuredGridFilter.h"
#include "vtkIdList.h"

class VTK_GRAPHICS_EXPORT vtkSubdivideSelectTetra : public vtkUnstructuredGridToUnstructuredGridFilter
{
public:
  static vtkSubdivideSelectTetra *New();
  vtkTypeRevisionMacro(vtkSubdivideSelectTetra,vtkUnstructuredGridToUnstructuredGridFilter);
  void PrintSelf(ostream& os, vtkIndent indent);
 
  vtkIdList * GetCellsToSubdivideIdList() { return this->IdList; }
  void SetCellsToSubdivideIdList(vtkIdList *ids) { this->IdList = ids; }   

protected:
  vtkSubdivideSelectTetra();
  ~vtkSubdivideSelectTetra() {};

  void Execute();

private:
  vtkSubdivideSelectTetra(const vtkSubdivideSelectTetra&);  // Not implemented.
  void operator=(const vtkSubdivideSelectTetra&);  // Not implemented.
  
  vtkIdList *IdList;
  
};

#endif


