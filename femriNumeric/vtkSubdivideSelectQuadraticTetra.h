/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkSubdivideSelectQuadraticTetra.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkSubdivideSelectQuadraticTetra - subdivide one tetrahedron into twelve for every tetra
// .SECTION Description
// This filter subdivides tetrahedra in an unstructured grid into twelve tetrahedra.


#ifndef __vtkSubdivideSelectQuadraticTetra_h
#define __vtkSubdivideSelectQuadraticTetra_h

#include "vtkUnstructuredGridToUnstructuredGridFilter.h"
#include "vtkIdList.h"

class VTK_GRAPHICS_EXPORT vtkSubdivideSelectQuadraticTetra : public vtkUnstructuredGridToUnstructuredGridFilter
{
public:
  static vtkSubdivideSelectQuadraticTetra *New();
  vtkTypeRevisionMacro(vtkSubdivideSelectQuadraticTetra,vtkUnstructuredGridToUnstructuredGridFilter);
  void PrintSelf(ostream& os, vtkIndent indent);
 
  vtkIdList * GetCellsToSubdivideIdList() { return this->IdList; }
  void SetCellsToSubdivideIdList(vtkIdList *ids) { this->IdList = ids; }   
  
  vtkSetMacro(SubdivideAllCells, int);
  vtkGetMacro(SubdivideAllCells, int);
  vtkBooleanMacro(SubdivideAllCells, int);

protected:
  vtkSubdivideSelectQuadraticTetra();
  ~vtkSubdivideSelectQuadraticTetra() {};

  int SubdivideAllCells;  

  void Execute();

private:
  vtkSubdivideSelectQuadraticTetra(const vtkSubdivideSelectQuadraticTetra&);  // Not implemented.
  void operator=(const vtkSubdivideSelectQuadraticTetra&);  // Not implemented.
  
  vtkIdList *IdList;
  
};

#endif


