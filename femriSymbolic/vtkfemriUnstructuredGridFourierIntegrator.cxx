/*=========================================================================

  Program:   femri
  Module:    $RCSfile: vtkfemriUnstructuredGridFourierIntegrator.cxx,v $
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

#include "vtkfemriUnstructuredGridFourierIntegrator.h"
#include "vtkObjectFactory.h"
#include "vtkUnstructuredGrid.h"
#include "vtkCellArray.h"
#include "vtkUnsignedCharArray.h"

#include "ginac.h"
#include "ginacfemriFourierIntegrators.h"

using namespace GiNaC;
using namespace std;

vtkStandardNewMacro(vtkfemriUnstructuredGridFourierIntegrator);
vtkCxxRevisionMacro(vtkfemriUnstructuredGridFourierIntegrator, "$Revision: 1.1.1.1 $");
vtkCxxSetObjectMacro(vtkfemriUnstructuredGridFourierIntegrator,Mesh,vtkUnstructuredGrid);

vtkfemriUnstructuredGridFourierIntegrator::vtkfemriUnstructuredGridFourierIntegrator()
{
  this->Mesh=NULL; 
  this->Function = numeric(1); 
  this->K = symbolic_matrix(3,1,"k"); 
  this->X = symbolic_matrix(3,1,"x"); 
  this->v0 = symbol("v0");
  this->v1 = symbol("v1");
  this->v2 = symbol("v2");
  this->V = matrix(3,1,lst(v0, v1, v2));

  this->UseElementIntegralCache = 0;

  this->NumberOfEvaluations = 0;
  this->NumberOfComputations = 0;
}

vtkfemriUnstructuredGridFourierIntegrator::~vtkfemriUnstructuredGridFourierIntegrator()
{
  if (this->Mesh) 
    {
    Mesh->Delete();
    }
}

numeric vtkfemriUnstructuredGridFourierIntegrator::ComputeIntegral(const lst& frequency_substitutions)
{
  ex integrand;
  integrand = this->Function.subs(frequency_substitutions,subs_options::no_pattern);

  size_t numberOfCells = this->Mesh->GetNumberOfCells();

  vtkCellArray *meshCellArray;
  vtkUnsignedCharArray *meshCellTypesArray;
  unsigned char cellType;
  int npts, *pts;
  
  meshCellArray = this->Mesh->GetCells();
  meshCellTypesArray = this->Mesh->GetCellTypesArray();

  meshCellArray->InitTraversal();

  numeric numeric_integral_sum(0.0);
  
  vector<ex> ksubs(3);
  ksubs[0] = this->K[0].subs(frequency_substitutions,subs_options::no_pattern);
  ksubs[1] = this->K[1].subs(frequency_substitutions,subs_options::no_pattern);
  ksubs[2] = this->K[2].subs(frequency_substitutions,subs_options::no_pattern);

  for (size_t i=0; i<numberOfCells; i++)
    {
    cellType = meshCellTypesArray->GetValue(i);
    meshCellArray->GetNextCell(npts,pts);
    if (!((cellType == VTK_TETRA) || (cellType == VTK_HEXAHEDRON) || (cellType == VTK_WEDGE) || (cellType == VTK_PYRAMID)))
      {
      continue;
      }
    
    ex integral;
    numeric numeric_integral;
    if (this->UseElementIntegralCache)
      {
      try
        {
        integral = this->ElementIntegralCache[i];
        this->EvaluateCachedElementIntegral(integral,frequency_substitutions,numeric_integral);
        ++this->NumberOfEvaluations;
        }
      catch (exception &p)
        {
        this->ComputeElementIntegral(cellType,npts,pts,integrand,integral);
        this->EvaluateElementIntegral(integral,numeric_integral);
        ++this->NumberOfComputations;
        }
      }
    else
      {
      this->ComputeElementIntegral(cellType,npts,pts,integrand,integral);
      this->EvaluateElementIntegral(integral,numeric_integral);
      ++this->NumberOfComputations;
      }
    numeric_integral_sum += numeric_integral;
    }

  return numeric_integral_sum;
}

// void vtkfemriUnstructuredGridFourierIntegrator::SubstituteInCachedTetrahedronIntegral(const GiNaC::ex& cached_integral, const GiNaC::lst& frequency_substitutions, GiNaC::ex& subs_integral)
// {
//   try
//     {
//     cached_integral.subs(frequency_substitutions,subs_options::no_pattern);
//     }
//   catch (exception& e)
//     {
//     cerr << "Error: " << e.what() << endl;
//     }
// }

void vtkfemriUnstructuredGridFourierIntegrator::BuildElementIntegralCache()
{
  size_t numberOfCells = this->Mesh->GetNumberOfCells();

  vtkCellArray *meshCellArray;
  vtkUnsignedCharArray *meshCellTypesArray;
  unsigned char cellType;
  int npts, *pts;
  
  meshCellArray = this->Mesh->GetCells();
  meshCellTypesArray = this->Mesh->GetCellTypesArray();

  meshCellArray->InitTraversal();

  ex integral;
  
  this->ElementIntegralCache.clear();
  this->ElementIntegralCache.resize(numberOfCells);

  for (size_t i=0; i<numberOfCells; i++)
    {
    cellType = meshCellTypesArray->GetValue(i);
    meshCellArray->GetNextCell(npts,pts);

    if (!((cellType == VTK_TETRA) || (cellType == VTK_HEXAHEDRON) || (cellType == VTK_WEDGE) || (cellType == VTK_PYRAMID)))
      {
      continue;
      }
    
    this->ComputeElementIntegral(cellType,npts,pts,this->Function,integral);
    this->ElementIntegralCache[i] = integral;
    }
}

void vtkfemriUnstructuredGridFourierIntegrator::EvaluateElementIntegral(const ex& integral, numeric& numeric_integral)
{
  try
    {
    numeric_integral = ex_to<numeric>(evalf(integral));
    }
  catch (exception &p) 
    {
    cerr << "Error: " << p.what() << endl;
    }
}

void vtkfemriUnstructuredGridFourierIntegrator::EvaluateCachedElementIntegral(const ex& cached_integral, const lst& frequency_substitutions, numeric& numeric_integral)
{
  try
    {
    numeric_integral = ex_to<numeric>(cached_integral.subs(frequency_substitutions,subs_options::no_pattern).evalf());
    }
  catch (exception& e)
    {
//     cerr << "Error: " << e.what() << endl;
    // Here we just rethrow in order to detect singularities.
    throw e;
    }
}

void vtkfemriUnstructuredGridFourierIntegrator::GetElementShapeFunctions(int cellType, int npts, matrix& shape_functions)
{
  switch (cellType)
    {
    case VTK_TETRA:
      shape_functions = matrix(npts,1,lst(
                                 1-v0-v1-v2,
                                 v0,
                                 v1,
                                 v2));
      break;
    case VTK_HEXAHEDRON:
      shape_functions = matrix(npts,1,lst(
                                 (1-v0)*(1-v1)*(1-v2),
                                    v0 *(1-v1)*(1-v2),
                                    v0 *   v1 *(1-v2),
                                 (1-v0)*   v1 *(1-v2),
                                 (1-v0)*(1-v1)*   v2 ,
                                    v0 *(1-v1)*   v2 ,
                                    v0 *   v1 *   v2 ,
                                 (1-v0)*   v1 *   v2 ));
      break;
//     case VTK_WEDGE:
//       break;
//     case VTK_PYRAMID:
//       break;
//     case VTK_QUADRATIC_TETRA:
//       break;
//     case VTK_QUADRATIC_HEXAHEDRON:
//       break;
//     case VTK_QUADRATIC_WEDGE:
//       break;
//     case VTK_QUADRATIC_PYRAMID:
//       break;
    default:
      throw invalid_argument("Element type not supported.");
    }
}

void vtkfemriUnstructuredGridFourierIntegrator::ComputeElementIntegral(int cellType, int npts, int* pts, const ex& integrand, ex& integral)
{

  matrix shape_functions;
  this->GetElementShapeFunctions(cellType,npts,shape_functions);

  matrix x_coords = matrix(npts,1);
  matrix y_coords = matrix(npts,1);
  matrix z_coords = matrix(npts,1);

  double point[3];

  for (int i=0; i<npts; i++)
    {
    this->Mesh->GetPoint(pts[i],point);
    x_coords(i,0) = numeric(point[0]);
    y_coords(i,0) = numeric(point[1]);
    z_coords(i,0) = numeric(point[2]);
    }

  try
    {
    //TODO: fix element orientation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //eg: if element not left-handed, permute columns in x_coords, y_coords, z_coords before building the transformation matrix, so that element points are reordered

    matrix transformation_matrix = matrix(3, 1);
    transformation_matrix(0,0) = (shape_functions.transpose() * x_coords).evalm();
    transformation_matrix(1,0) = (shape_functions.transpose() * y_coords).evalm();
    transformation_matrix(2,0) = (shape_functions.transpose() * z_coords).evalm();

    lst std_substitution_list = lst(this->X[0] == transformation_matrix[0][0], 
                                    this->X[1] == transformation_matrix[1][0], 
                                    this->X[2] == transformation_matrix[2][0]);

    matrix jacobian_matrix = matrix(3, 3);
    jacobian_matrix(0,0) = transformation_matrix[0][0].diff(v0);
    jacobian_matrix(0,1) = transformation_matrix[1][0].diff(v0);
    jacobian_matrix(0,2) = transformation_matrix[2][0].diff(v0);
    jacobian_matrix(1,0) = transformation_matrix[0][0].diff(v1);
    jacobian_matrix(1,1) = transformation_matrix[1][0].diff(v1);
    jacobian_matrix(1,2) = transformation_matrix[2][0].diff(v1);
    jacobian_matrix(2,0) = transformation_matrix[0][0].diff(v2);
    jacobian_matrix(2,1) = transformation_matrix[1][0].diff(v2);
    jacobian_matrix(2,2) = transformation_matrix[2][0].diff(v2);

//     ex jacobian_determinant = jacobian_matrix.determinant();
    //TODO: the absolute value here is a hack to solve the element left-handed problem (which leads to negative Jacobians)
    ex jacobian_determinant = abs(jacobian_matrix.determinant());

    ex std_integrand = jacobian_determinant * integrand.subs(std_substitution_list,subs_options::no_pattern);

    switch (cellType)
      {
      case VTK_TETRA:
      {
        calc_std_tetrahedron_fourier_integral do_tetrahedron_integral;
        integral = do_tetrahedron_integral(std_integrand,this->v0,this->v1,this->v2);
        break;
      }
      case VTK_HEXAHEDRON:
      {
        calc_std_hexahedron_fourier_integral do_hexahedron_integral;
        integral = do_hexahedron_integral(std_integrand,this->v0,this->v1,this->v2);
        break;
      }
//       case VTK_WEDGE:
//         break;
//       case VTK_PYRAMID:
//         break;
//       case VTK_QUADRATIC_TETRA:
//         break;
//       case VTK_QUADRATIC_HEXAHEDRON:
//         break;
//       case VTK_QUADRATIC_WEDGE:
//         break;
//       case VTK_QUADRATIC_PYRAMID:
//         break;
      default:
        throw invalid_argument("Element type not supported.");
      }
   
    }
  catch (exception &p) 
    {
    cerr << "Error: " << p.what() << endl;
    }

}

