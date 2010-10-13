/*=========================================================================

  Program:   femri
  Module:    $RCSfile: vtkfemriUnstructuredGridNumericKSpaceGenerator.h,v $
  Language:  C++
  Date:      $Date: 2008/11/04 11:23:41 $
  Version:   $Revision: 1.5 $

  Copyright (c) Luca Antiga, David Steinman. All rights reserved.
  See LICENCE file for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm 
  for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

// .NAME vtkfemriUnstructuredGridToImageFilter - ..
// .SECTION Description
// ..
// .SECTION Thanks
// This class was developed by Luca Antiga, PhD \n
// Imaging Research Labs - Robarts Research Institute \n
// Bioengineering Dept - Mario Negri Institute for Pharmacological Research \n 
// email: lantiga@imaging.robarts.ca - homepage: http://www.imaging.robarts.ca/~lantiga

#ifndef __vtkfemriUnstructuredGridNumericKSpaceGenerator_h
#define __vtkfemriUnstructuredGridNumericKSpaceGenerator_h

#include "vtkfemriKSpaceGenerator.h"
#include "vtkfemriNumericWin32Header.h"

class vtkUnstructuredGrid;
class vtkCell;
class vtkfemriOptimalQuadratureOrderCalculator;

class VTK_FEMRI_NUMERIC_EXPORT vtkfemriUnstructuredGridNumericKSpaceGenerator : public vtkfemriKSpaceGenerator
{
public:
  vtkTypeRevisionMacro(vtkfemriUnstructuredGridNumericKSpaceGenerator,vtkfemriKSpaceGenerator);
  static vtkfemriUnstructuredGridNumericKSpaceGenerator* New();

  vtkSetMacro(QuadratureOrder,int);
  vtkGetMacro(QuadratureOrder,int);

  vtkSetMacro(NumberOfSubdivisions,int);
  vtkGetMacro(NumberOfSubdivisions,int);

  vtkSetMacro(UseOptimalAlgorithm,int);
  vtkGetMacro(UseOptimalAlgorithm,int);  
  vtkBooleanMacro(UseOptimalAlgorithm,int);
  
  vtkSetMacro(ErrorThreshold,double);
  vtkGetMacro(ErrorThreshold,double); 

  vtkSetMacro(UniformMagnetization,int);
  vtkGetMacro(UniformMagnetization,int);  
  vtkBooleanMacro(UniformMagnetization,int);

  vtkSetMacro(MagnetizationValue,double);
  vtkGetMacro(MagnetizationValue,double);  

  vtkSetStringMacro(MagnetizationArrayName);
  vtkGetStringMacro(MagnetizationArrayName);

  vtkGetMacro(NumberOfGaussPointEvaluations,int);
  vtkGetMacro(MaximumQuadratureOrderUsed,int);

  vtkSetMacro(SliceSelection,int);
  vtkGetMacro(SliceSelection,int);  
  vtkBooleanMacro(SliceSelection,int);

  vtkSetMacro(SliceProfile,int);
  vtkGetMacro(SliceProfile,int);
  void SetSliceProfileToIdeal()
  { this->SliceProfile = FEMRI_IDEAL_SLICE_PROFILE; }
  void SetSliceProfileToTrapezoidal()
  { this->SliceProfile = FEMRI_TRAPEZOIDAL_SLICE_PROFILE; }
  void SetSliceProfileToQuadratic()
  { this->SliceProfile = FEMRI_QUADRATIC_SLICE_PROFILE; }

  vtkSetMacro(SliceThickness,double);
  vtkGetMacro(SliceThickness,double);

  vtkSetMacro(SliceOrigin,double);
  vtkGetMacro(SliceOrigin,double);

  virtual void Initialize();

  virtual void EvaluateFourierFunction(double frequency[3], double value[2]);

protected:
  vtkfemriUnstructuredGridNumericKSpaceGenerator();
  ~vtkfemriUnstructuredGridNumericKSpaceGenerator();
  
  virtual int FillInputPortInformation(int, vtkInformation*);

//  int ComputeOptimalOrder(vtkCell* cell, double frequency[3]);
//  double ComputeCyclesPerElement(vtkCell* cell, double frequency[3]);

  virtual double ComputeVolume()
  {
    return 1.0;
  }

  double ComputeJacobian(vtkCell* cell, double pcoords[3]);
  
  double TrapezoidalSliceProfile(double z); 
  double IdealSliceProfile(double z); 
  double QuadraticSliceProfile(double z); 

  int UseOptimalAlgorithm;
  double ErrorThreshold;

  int QuadratureOrder;
  int NumberOfSubdivisions;

  int SliceSelection;
  int SliceProfile;
  double SliceThickness; 
  double SliceOrigin;

  int UniformMagnetization;
  double MagnetizationValue;
  char* MagnetizationArrayName;

  int NumberOfGaussPointEvaluations;
  int MaximumQuadratureOrderUsed;

//BTX
  enum
  {
    FEMRI_IDEAL_SLICE_PROFILE,
    FEMRI_TRAPEZOIDAL_SLICE_PROFILE,
    FEMRI_QUADRATIC_SLICE_PROFILE
  };
//ETX
  
//  int InitDesignCurves(void);
//  femriIntegrationDesignCurveInterpolator *designCurves;

  vtkfemriOptimalQuadratureOrderCalculator* OptimalQuadratureOrderCalculator;
  
private:
  vtkfemriUnstructuredGridNumericKSpaceGenerator(const vtkfemriUnstructuredGridNumericKSpaceGenerator&);  // Not implemented.
  void operator=(const vtkfemriUnstructuredGridNumericKSpaceGenerator&);  // Not implemented.
};

#endif
