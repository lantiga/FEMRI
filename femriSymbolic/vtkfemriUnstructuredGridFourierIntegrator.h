/*=========================================================================

  Program:   femri
  Module:    $RCSfile: vtkfemriUnstructuredGridFourierIntegrator.h,v $
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

// .NAME vtkfemriUnstructuredGridFourierIntegrator - ..
// .SECTION Description
// ..
// .SECTION Thanks
// This class was developed by Luca Antiga, PhD \n
// Imaging Research Labs - Robarts Research Institute \n
// Bioengineering Dept - Mario Negri Institute for Pharmacological Research \n 
// email: lantiga@imaging.robarts.ca - homepage: http://www.imaging.robarts.ca/~lantiga

#ifndef __vtkfemriUnstructuredGridFourierIntegrator_h
#define __vtkfemriUnstructuredGridFourierIntegrator_h

#include "vtkObject.h"
#include "vtkfemriSymbolicWin32Header.h"
#include "ginac.h"

class vtkUnstructuredGrid;

class VTK_FEMRI_SYMBOLIC_EXPORT vtkfemriUnstructuredGridFourierIntegrator : public vtkObject
{
public:
  vtkTypeRevisionMacro(vtkfemriUnstructuredGridFourierIntegrator,vtkObject);
  static vtkfemriUnstructuredGridFourierIntegrator* New();

  virtual void SetMesh(vtkUnstructuredGrid* mesh);
  vtkGetObjectMacro(Mesh,vtkUnstructuredGrid);

  vtkSetMacro(UseElementIntegralCache,int);
  vtkGetMacro(UseElementIntegralCache,int);
  vtkBooleanMacro(UseElementIntegralCache,int);

  //BTX
  void SetFunction(const GiNaC::ex& function) { this->Function = function; }
  const GiNaC::ex& GetFunction() { return this->Function; }

  GiNaC::numeric ComputeIntegral(const GiNaC::lst& frequency_substitutions);
  void ComputeIntegral(const GiNaC::lst& frequency_substitutions, double& re, double& im)
    {
    GiNaC::numeric fourierIntegral = this->ComputeIntegral(frequency_substitutions);
    re = fourierIntegral.real().to_double();
    im = fourierIntegral.imag().to_double();
    }

  void SetK(const GiNaC::ex& k) { this->K = k; }
  const GiNaC::ex& GetK() { return this->K; }

  void SetX(const GiNaC::ex& x) { this->X = x; }
  const GiNaC::ex& GetX() { return this->X; }

  void BuildElementIntegralCache();
  //ETX

  vtkGetMacro(NumberOfComputations,int);
  vtkGetMacro(NumberOfEvaluations,int);

protected:
  vtkfemriUnstructuredGridFourierIntegrator();
  ~vtkfemriUnstructuredGridFourierIntegrator();

  //BTX

  void GetElementShapeFunctions(int cellType, int npts, GiNaC::matrix& shape_functions);
  void ComputeElementIntegral(int cellType, int npts, int* pts, const GiNaC::ex& integrand, GiNaC::ex& integral);

  void EvaluateElementIntegral(const GiNaC::ex& integral, GiNaC::numeric& numeric_integral);
  void EvaluateCachedElementIntegral(const GiNaC::ex& cached_integral, const GiNaC::lst& frequency_substitutions, GiNaC::numeric& numeric_integral);

//   void SubstituteInCachedTetrahedronIntegral(const GiNaC::ex& cached_integral, const GiNaC::lst& frequency_substitutions, GiNaC::ex& subs_integral);

  GiNaC::ex Function;
  GiNaC::ex K;
  GiNaC::ex X;
  GiNaC::matrix V;
  GiNaC::symbol v0, v1, v2;
  GiNaC::matrix P;

  std::vector<GiNaC::ex> ElementIntegralCache;

  //ETX

  int UseElementIntegralCache;

  vtkUnstructuredGrid* Mesh;
  
  bool FirstRun;

  int NumberOfEvaluations;
  int NumberOfComputations;

private:
  vtkfemriUnstructuredGridFourierIntegrator(const vtkfemriUnstructuredGridFourierIntegrator&);  // Not implemented.
  void operator=(const vtkfemriUnstructuredGridFourierIntegrator&);  // Not implemented.
};

#endif
