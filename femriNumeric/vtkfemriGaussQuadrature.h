/*=========================================================================

  Program:   femri
  Module:    $RCSfile: vtkfemriGaussQuadrature.h,v $
  Language:  C++
  Date:      $Date: 2007/07/07 15:38:00 $
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

#ifndef __vtkfemriGaussQuadrature_h
#define __vtkfemriGaussQuadrature_h

#include "vtkObject.h"
#include "vtkfemriNumericWin32Header.h"

#include "vtkDoubleArray.h"

class VTK_FEMRI_NUMERIC_EXPORT vtkfemriGaussQuadrature : public vtkObject
{
public:
  vtkTypeRevisionMacro(vtkfemriGaussQuadrature,vtkObject);
  static vtkfemriGaussQuadrature* New();

  vtkGetObjectMacro(QuadraturePoints,vtkDoubleArray);
  vtkGetObjectMacro(QuadratureWeights,vtkDoubleArray);

  vtkSetMacro(Order,int);
  vtkGetMacro(Order,int);
  
  int GetNumberOfQuadraturePoints()
  {
    return this->QuadraturePoints->GetNumberOfTuples();
  }
  
  double* GetQuadraturePoint(vtkIdType id)
  {
    return this->QuadraturePoints->GetTuple(id);
  }
 
  void GetQuadraturePoint(vtkIdType id, double* quadraturePoint)
  {
    this->QuadraturePoints->GetTuple(id,quadraturePoint);
  }
  
  double GetQuadratureWeight(vtkIdType id)
  {
    return this->QuadratureWeights->GetValue(id);
  }
 
  void Initialize(vtkIdType cellType);
 
  void Initialize1DGauss();
  void Initialize1DJacobi(int alpha, int beta);
  void ScaleTo01();
 
protected:
  vtkfemriGaussQuadrature();
  ~vtkfemriGaussQuadrature();

  void TensorProductQuad(vtkfemriGaussQuadrature* q1D);
  void TensorProductTriangle(vtkfemriGaussQuadrature* gauss1D, vtkfemriGaussQuadrature* jacA1D);
  
  void TensorProductHexahedron(vtkfemriGaussQuadrature* q1D);
  void TensorProductWedge(vtkfemriGaussQuadrature* q1D, vtkfemriGaussQuadrature* q2D);
  void TensorProductTetra(vtkfemriGaussQuadrature* gauss1D, vtkfemriGaussQuadrature* jacA1D, vtkfemriGaussQuadrature* jacB1D);
  
  vtkDoubleArray* QuadraturePoints;
  vtkDoubleArray* QuadratureWeights;
 
  int Order;
  int QuadratureType;
  vtkIdType CellType;
  int PreviousOrder;

private:  
  vtkfemriGaussQuadrature(const vtkfemriGaussQuadrature&);  // Not implemented.
  void operator=(const vtkfemriGaussQuadrature&);  // Not implemented.

};

#endif
