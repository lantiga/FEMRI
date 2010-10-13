/*=========================================================================

  Program:   femri
  Module:    $RCSfile: vtkfemriUnstructuredGridMRSlicer.cxx,v $
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

#include "vtkfemriUnstructuredGridMRSlicer.h"
#include "vtkfemriUnstructuredGridNumericFourierIntegrator.h"
#include "vtkfemriUnstructuredGridToImageFilter.h"
#include "vtkUnstructuredGrid.h"
#include "vtkPointData.h"
#include "vtkDoubleArray.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkSubdivideSelectTetra.h"
#include "vtkSubdivideSelectQuadraticTetra.h"


using namespace std;

vtkStandardNewMacro(vtkfemriUnstructuredGridMRSlicer);
vtkCxxRevisionMacro(vtkfemriUnstructuredGridMRSlicer, "$Revision: 1.1.1.1 $");

vtkfemriUnstructuredGridMRSlicer::vtkfemriUnstructuredGridMRSlicer()
{
    this->SliceProfile = 0; //default to ideal profile;
}

vtkfemriUnstructuredGridMRSlicer::~vtkfemriUnstructuredGridMRSlicer()
{
}

 
double vtkfemriUnstructuredGridMRSlicer::TrapezoidResidual(double z_min, double z_max)
{
    const double dz = 0.01; // discretization step along z (chosen to be very small)
    
    double r2 = 0;
    double r2tmp = 0;
    
    double ML1 = this->TrapezoidSlice(z_min);
    double ML2 = this->TrapezoidSlice(z_max);
     
    double m = (ML2 - ML1)/(z_max - z_min);
    //cout << "m=" << m << " (ML2,ML1)=(" << ML2 << "," << ML1 << ")" << endl;
    
    for (double z = z_min; z <= z_max; z += dz)
    {
        double ML = m*(z - z_min) + ML1;
        
        double M = this->TrapezoidSlice(z);
        //cout << z_min << "," << z_max << " M=" << M << " ML=" << ML << endl;
        r2 += (M - ML)*(M - ML);
        //if (r2tmp > r2) { r2 = r2tmp; } 
    }
    
    return r2;
    
}
  
double vtkfemriUnstructuredGridMRSlicer::TrapezoidSlice(double z)
{
    double mag = 0;
    double slope;
    // constant section of trapezoid
    double z_slice_bound_upper = this->GetSliceOrigin()[2] + this->GetSliceThickness()/2;
    double z_slice_bound_lower = this->GetSliceOrigin()[2] - this->GetSliceThickness()/2;
    
    // upper and lower tail ends of the trapezoidal function
    double z_tail_bound_upper = z_slice_bound_upper + this->GetSliceThickness()/2;
    double z_tail_bound_lower = z_slice_bound_lower - this->GetSliceThickness()/2;

    if ((z <= z_slice_bound_upper) && (z >= z_slice_bound_lower))
    {
        mag = 1.0;
    }
    else if ((z > z_slice_bound_upper) && (z <= z_tail_bound_upper))
    {
        slope = -1/(z_tail_bound_upper - z_slice_bound_upper);
        mag = slope*(z - z_slice_bound_upper) + 1.0;  
    }
    else if ((z < z_slice_bound_lower) && (z >= z_tail_bound_lower))
    {
        slope = 1/(z_slice_bound_lower - z_tail_bound_lower);
        mag = slope*(z - z_tail_bound_lower);  
    }
    else { mag = 0.0; }
   
    return mag;
    
}

void vtkfemriUnstructuredGridMRSlicer::ExecuteInformation()
{

}

void vtkfemriUnstructuredGridMRSlicer::ExecuteData(vtkDataObject *outp)
{
    
    vtkUnstructuredGrid* input = this->GetInput();
    vtkUnstructuredGrid* output = this->GetOutput();
    
    vtkUnstructuredGrid *tmpgrid = vtkUnstructuredGrid::New();    
    
    vtkIdType cellId, pointId, numPoints; 
    vtkIdList *subdiv_list;
  
    vtkIdType numCells = input->GetNumberOfCells();
    vtkIdType nsubdiv = numCells; // initialize to non-zero;
    const double residual_threshold = this->Residual;;
        
    tmpgrid->DeepCopy(input);
    /*
    for (double z=-1.0; z <= 1.0; z=z+0.01) 
    {
        cout << this->TrapezoidSlice(z) << " "; 
    }
      */ 
    
    subdiv_list = vtkIdList::New(); 
    
    // iteratively subdivide cells until residual small enough   
    while (nsubdiv > 0) 
    {   
        nsubdiv = 0;
        subdiv_list->Reset();          
        // figure out which cells need subdividing
        for(cellId=0; cellId < numCells; cellId++)
        {
            double bounds[6];        
            tmpgrid->GetCellBounds(cellId, bounds);      
            //cout << "Computing residual..." << endl; 
            double residual = this->TrapezoidResidual(bounds[4], bounds[5]);
            
            //cout << "bounds=" << bounds[4] << "," << bounds[5] << 
            //cout << " residual=" << residual << endl;
            if (residual > residual_threshold)
            {
                
                subdiv_list->InsertNextId(cellId);
                //cout << cellId << endl;
                nsubdiv++;
            }            
        }
          
          
        // seems like this could be done more elegantly, but enh...  
        if (tmpgrid->GetCellType(0) == VTK_TETRA)
        {
            // subdivide each tagged cell once
            vtkSubdivideSelectTetra *subdiv = vtkSubdivideSelectTetra::New();  
            subdiv->SetCellsToSubdivideIdList(subdiv_list);
            tmpgrid->Update();
            subdiv->SetInput(tmpgrid);
            subdiv->Update();
            // copy output of subdivided mesh to slicer output    
            output->DeepCopy(subdiv->GetOutput()); 
            subdiv->Delete();
        }
        else if (tmpgrid->GetCellType(0) == VTK_QUADRATIC_TETRA)
        {
            // subdivide each tagged cell once
            vtkSubdivideSelectQuadraticTetra *subdiv = vtkSubdivideSelectQuadraticTetra::New();  
            subdiv->SetCellsToSubdivideIdList(subdiv_list);
            tmpgrid->Update();
            subdiv->SetInput(tmpgrid);
            subdiv->Update();
            // copy output of subdivided mesh to slicer output    
            output->DeepCopy(subdiv->GetOutput()); 
            subdiv->Delete();
        }
        
        
        tmpgrid->DeepCopy(output); 
        
        numCells = tmpgrid->GetNumberOfCells();
            
    }
    
    
    numPoints = output->GetNumberOfPoints();
    
    vtkDoubleArray *scalars = vtkDoubleArray::New(); 
    
    // Assign magnetization values to slice-refined mesh
    for (pointId=0; pointId < numPoints; pointId++)
    {
       double mag = this->TrapezoidSlice(output->GetPoint(pointId)[2]); 
       scalars->InsertNextValue(mag);    
    }
    
    output->GetPointData()->SetScalars(scalars);
    
    subdiv_list->Delete();
    tmpgrid->Delete();

}
