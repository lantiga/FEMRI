/*=========================================================================

  Program:   femri
  Module:    vtkfemriPolyDataNumericKSpaceGenerator.h
  Language:  C++
  Date:      $Date: 2006/04/06 16:46:43 $
  Version:   $Revision: 1.4 $

  Copyright (c) Luca Antiga, David Steinman. All rights reserved.
  See LICENCE file for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm 
  for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

// .NAME vtkfemriPolyDataNumericKSpaceGenerator - ..
// .SECTION Description
// ..

#ifndef __vtkfemriPolyDataNumericKSpaceGenerator_h
#define __vtkfemriPolyDataNumericKSpaceGenerator_h

#include "vtkfemriKSpaceGenerator.h"
#include "vtkfemriNumericWin32Header.h"

class vtkPolyData;
class vtkCell;
class vtkfemriOptimalQuadratureOrderCalculator;

class VTK_FEMRI_NUMERIC_EXPORT vtkfemriPolyDataNumericKSpaceGenerator : public vtkfemriKSpaceGenerator
{
public:
  vtkTypeRevisionMacro(vtkfemriPolyDataNumericKSpaceGenerator,vtkImageAlgorithm);
  static vtkfemriPolyDataNumericKSpaceGenerator* New();

  vtkSetMacro(QuadratureOrder,int);
  vtkGetMacro(QuadratureOrder,int);

  vtkSetMacro(MagnetizationValue,double);
  vtkGetMacro(MagnetizationValue,double);  

  vtkSetMacro(UseOptimalAlgorithm,int);
  vtkGetMacro(UseOptimalAlgorithm,int);  
  vtkBooleanMacro(UseOptimalAlgorithm,int);
  
  vtkSetMacro(ErrorThreshold,double);
  vtkGetMacro(ErrorThreshold,double); 

  vtkGetMacro(NumberOfGaussPointEvaluations,int);
  vtkGetMacro(MaximumQuadratureOrderUsed,int);

  virtual void Initialize();

  virtual void EvaluateFourierFunction(double frequency[3], double value[2]);

protected:
  vtkfemriPolyDataNumericKSpaceGenerator();
  ~vtkfemriPolyDataNumericKSpaceGenerator();
  
  virtual int FillInputPortInformation(int, vtkInformation*);

//  virtual int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

  virtual double ComputeVolume()
  {
    return 1.0;
  }

  double ComputeJacobian(vtkCell* cell, double pcoords[3]);
  void EvaluateCellLocation(vtkCell* cell, double pcoords[3], double x[3], double* weights);
  void CellInterpolationDerivs(double* pcoords, double* derivs);

  double MagnetizationValue;
  int QuadratureOrder;

  int UseOptimalAlgorithm;
  double ErrorThreshold;

  int NumberOfGaussPointEvaluations;
  int MaximumQuadratureOrderUsed;

  vtkfemriOptimalQuadratureOrderCalculator* OptimalQuadratureOrderCalculator;

private:
  vtkfemriPolyDataNumericKSpaceGenerator(const vtkfemriPolyDataNumericKSpaceGenerator&);  // Not implemented.
  void operator=(const vtkfemriPolyDataNumericKSpaceGenerator&);  // Not implemented.
};

#endif
