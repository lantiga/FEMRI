/*=========================================================================

  Program:   femri
  Module:    $RCSfile: vtkfemriOptimalQuadratureOrderCalculator.h,v $
  Language:  C++
  Date:      $Date: 2008/11/03 17:00:30 $
  Version:   $Revision: 1.3 $

  Copyright (c) Luca Antiga, David Steinman. All rights reserved.
  See LICENCE file for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm 
  for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __vtkfemriOptimalQuadratureOrderCalculator_h
#define __vtkfemriOptimalQuadratureOrderCalculator_h

#include "vtkImageData.h"
#include "vtkCell.h"
#include "vtkObject.h"
#include "vtkfemriNumericWin32Header.h"

class VTK_FEMRI_NUMERIC_EXPORT vtkfemriOptimalQuadratureOrderCalculator : public vtkObject
{
public:
  vtkTypeRevisionMacro(vtkfemriOptimalQuadratureOrderCalculator,vtkObject);
  static vtkfemriOptimalQuadratureOrderCalculator* New();

  void Initialize();

  vtkSetVector3Macro(MaxFrequency,double);
  vtkGetVector3Macro(MaxFrequency,double);
 
  vtkSetVector3Macro(FrequencySpacing,double);
  vtkGetVector3Macro(FrequencySpacing,double);
 
  vtkSetMacro(ErrorThreshold,double);
  vtkGetMacro(ErrorThreshold,double);
 
  vtkSetMacro(MaxQuadratureOrder,int);
  vtkGetMacro(MaxQuadratureOrder,int);

  vtkSetMacro(CyclesPerElementResolution,double);
  vtkGetMacro(CyclesPerElementResolution,double);
 
  vtkGetObjectMacro(DesignTable,vtkImageData);
  vtkSetObjectMacro(DesignTable,vtkImageData);

  vtkGetObjectMacro(DataSet,vtkDataSet);
  vtkSetObjectMacro(DataSet,vtkDataSet);

  int ComputeOptimalQuadratureOrder(vtkIdType cellId, double frequency[3]);

  void BuildDesignTable(double maxCyclesPerElement);

protected:
  vtkfemriOptimalQuadratureOrderCalculator();
  ~vtkfemriOptimalQuadratureOrderCalculator();

  double ComputeCyclesPerElement(vtkIdType cellId, double frequency[3]);

  double MaxFrequency[3];
  double FrequencySpacing[3];
  double ErrorThreshold;
  int MaxQuadratureOrder;
  double CyclesPerElementResolution;

  vtkDataSet* DataSet;
  vtkImageData* DesignTable;

private:  
  vtkfemriOptimalQuadratureOrderCalculator(const vtkfemriOptimalQuadratureOrderCalculator&);  // Not implemented.
  void operator=(const vtkfemriOptimalQuadratureOrderCalculator&);  // Not implemented.

};

#endif
