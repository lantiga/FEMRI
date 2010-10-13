/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkSubdivideSelectTetra.cxx,v $

  Extended by Paul Simedrea (2006)

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkSubdivideSelectTetra.h"

#include "vtkCellType.h"
#include "vtkGenericCell.h"
#include "vtkMergePoints.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkUnstructuredGrid.h"

vtkCxxRevisionMacro(vtkSubdivideSelectTetra, "$Revision: 1.1.1.1 $");
vtkStandardNewMacro(vtkSubdivideSelectTetra);

// Description:
// Construct with all types of clipping turned off.
vtkSubdivideSelectTetra::vtkSubdivideSelectTetra()
{
}


double tet_jacobian_det(double *p1, double *p2, double *p3, double *p4) {
    double J=0;
    
    cout << "p1=(" << p1[0] << "," << p1[1] << "," << p1[2] << ")"<<endl;
    cout << "p2=(" << p2[0] << "," << p2[1] << "," << p2[2] << ")"<<endl;
    cout << "p3=(" << p3[0] << "," << p3[1] << "," << p3[2] << ")"<<endl;
    cout << "p4=(" << p4[0] << "," << p4[1] << "," << p4[2] << ")"<<endl;
     
    
    double c1 = p2[0] - p1[0]; cout << "c1 = " << p2[0] << "-" << p1[0] << "=" << c1 << endl; 
    double c2 = p3[0] - p1[0]; cout << "c2 = " << p3[0] << "-" << p1[0] << "=" << c1 << endl;
    double c3 = p4[0] - p1[0]; cout << "c3 = " << p4[0] << "-" << p1[0] << "=" << c1 << endl;
   
    double d1 = p2[1] - p1[1]; cout << "d1 = " << p2[1] << "-" << p1[1] << "=" << d1 << endl;
    double d2 = p3[1] - p1[1]; cout << "d2 = " << p3[1] << "-" << p1[1] << "=" << d2 << endl;
    double d3 = p4[1] - p1[1]; cout << "d3 = " << p4[1] << "-" << p1[1] << "=" << d3 << endl;
   
    double e1 = p2[2] - p1[2]; cout << "e1 = " << p2[2] << "-" << p1[2] << "=" << e1 << endl;
    double e2 = p3[2] - p1[2]; cout << "e2 = " << p3[2] << "-" << p1[2] << "=" << e2 << endl;
    double e3 = p4[2] - p1[2]; cout << "e3 = " << p4[2] << "-" << p1[2] << "=" << e3 << endl;
   
    J = c1*(d2*e3 - e2*d3) - c2*(d1*e3 - e1*d3) + c3*(d1*e2 - e1*d2);
    cout << "J= " << J << endl;
    /*  should add check to make sure Jacobian isn't too small (0) */ 
    
    return J;

}


void vtkSubdivideSelectTetra::Execute()
{
  vtkUnstructuredGrid *input=(vtkUnstructuredGrid *)this->GetInput();
  vtkIdType numPts = input->GetNumberOfPoints();
  vtkIdType numCells = input->GetNumberOfCells();
  vtkIdType numsubCells = 0;
  vtkPoints *inPts=input->GetPoints();
  vtkIdType cellId, i;
  vtkIdType pts[4];
  vtkGenericCell *cell;
  vtkPointData *pd = input->GetPointData();
  vtkUnstructuredGrid *output = this->GetOutput();
  vtkPointData *outputPD = output->GetPointData();
  vtkPoints *newPts;
  vtkIdType ptId;

  double weights[4], x0[3], x1[3], x2[3], x3[3], x[3];
  int p0, p1, p2, p3;
  vtkIdType center, e01, e02, e03, e12, e13, e23;
  vtkMergePoints *locator;
  vtkIdType tokenId;
  
   
  
  vtkDebugMacro(<<"Executing mesh subdivide");

    
  // Check if all cells are same type and tetrahedrons
  // Change this possibly later on to accomodate inhomogeneous meshes.   
  
  if (input->GetCellType(0) != VTK_TETRA)
    {  
    vtkErrorMacro(<<"cells are not of type VTK_TETRA.");
    return;        
    }
  if (input->IsHomogeneous() == 0)
   
    {
    vtkErrorMacro(<<"mesh is not homogenous.");
    return;
    }

  
  // Copy original points and point data
  newPts = vtkPoints::New();
  newPts->Allocate(5*numPts,numPts);
  outputPD->InterpolateAllocate(pd,5*numPts,numPts);
  
  output->Allocate(numCells);
  output->SetPoints(newPts);

  locator = vtkMergePoints::New();
  locator->InitPointInsertion (newPts, input->GetBounds());

  // copy points from input to output
  for (ptId=0; ptId < numPts; ptId++)
    {
    // add existing points in the locator so that any copies are later detected.
    locator->InsertUniquePoint(inPts->GetPoint(ptId),tokenId);
    outputPD->CopyData(pd,ptId,ptId);
    }

  cell = vtkGenericCell::New();
  // loop over tetrahedra, generating sixteen new ones for each. This is
  // done by introducing mid-edge nodes and a single mid-tetra node.
  for(cellId=0; cellId < numCells; cellId++)
    {
    input->GetCell(cellId, cell);
    // get tetra points
    cell->Points->GetPoint(0,x0);
    cell->Points->GetPoint(1,x1);
    cell->Points->GetPoint(2,x2);
    cell->Points->GetPoint(3,x3);

    p0 = cell->PointIds->GetId(0);
    p1 = cell->PointIds->GetId(1);
    p2 = cell->PointIds->GetId(2);
    p3 = cell->PointIds->GetId(3);

    //std::cerr << this->IdList->IsId(cellId) << std::endl;
    if (this->IdList->IsId(cellId) >= 0)
    {   
    
        // compute center point
        weights[0] = weights[1] = weights[2] = weights[3] = 0.25;
        for (i=0; i<3; i++)
        {
        x[i] = 0.25*(x0[i] + x1[i] + x2[i] + x3[i]);
        }
        locator->InsertUniquePoint(x,center);
        outputPD->InterpolatePoint(pd, center, cell->PointIds, weights);
    
        // compute edge points
        // edge 0-1
        for (i=0; i<3; i++)
        {
            x[i] = 0.5 * (x1[i] + x0[i]);
        }
        locator->InsertUniquePoint(x,e01);
        outputPD->InterpolateEdge(pd, e01, p0, p1, 0.5);
    
        // edge 1-2
        for (i=0; i<3; i++)
        {
            x[i] = 0.5 * (x2[i] + x1[i]);
        }
        locator->InsertUniquePoint(x,e12);
        outputPD->InterpolateEdge(pd, e12, p1, p2, 0.5);
    
        // edge 2-0
        for (i=0; i<3; i++)
        {
            x[i] = 0.5 * (x2[i] + x0[i]);
        }
        locator->InsertUniquePoint(x,e02);
        outputPD->InterpolateEdge(pd, e02, p2, p0, 0.5);
    
        // edge 0-3
        for (i=0; i<3; i++)
        {
            x[i] = 0.5 * (x3[i] + x0[i]);
        }
        locator->InsertUniquePoint(x,e03);
        outputPD->InterpolateEdge(pd, e03, p0, p3, 0.5);
    
        // edge 1-3
        for (i=0; i<3; i++)
        {
            x[i] = 0.5 * (x3[i] + x1[i]);
        }
        locator->InsertUniquePoint(x,e13);
        outputPD->InterpolateEdge(pd, e13, p1, p3, 0.5);
    
        // edge 2-3
        for (i=0; i<3; i++)
        {
            x[i] = 0.5 * (x3[i] + x2[i]);
        }
        locator->InsertUniquePoint(x,e23);
        outputPD->InterpolateEdge(pd, e23, p2, p3, 0.5);

        // Now create tetrahedra
        // First, four tetra from each vertex
        pts[0] = p0;
        pts[1] = e01;
        pts[2] = e02;
        pts[3] = e03;
        output->InsertNextCell(VTK_TETRA, 4, pts);
//        vtkPoints *locpts = locator->GetPoints();
                           
        pts[0] = e01;
        pts[1] = e12;
        pts[2] = e13;
        pts[3] = p1;
        //cout << "ids= " << pts[0] << "," << pts[1] << "," << pts[2] << "," << pts[3] << endl;  
        output->InsertNextCell(VTK_TETRA, 4, pts);
        
        pts[0] = p2;        
        pts[1] = e02;
        pts[2] = e12;
        pts[3] = e23;           
        output->InsertNextCell(VTK_TETRA, 4, pts);
                
        pts[0] = e03;
        pts[1] = e13;
        pts[2] = e23;
        pts[3] = p3;       
        output->InsertNextCell(VTK_TETRA, 4, pts);


        // Now four tetra from cut-off tetra corners
        pts[0] = e01;
        pts[1] = e02;
        pts[2] = e03;
        pts[3] = center;
        output->InsertNextCell(VTK_TETRA, 4, pts);
        
        pts[0] = center;
        pts[1] = e01;
        pts[2] = e12;
        pts[3] = e13;
        output->InsertNextCell(VTK_TETRA, 4, pts);
        
        pts[0] = center;        
        pts[1] = e23;
        pts[2] = e12;
        pts[3] = e02;
        output->InsertNextCell(VTK_TETRA, 4, pts);
        
        pts[1] = e03;
        pts[2] = e13;
        pts[3] = e23;
        output->InsertNextCell(VTK_TETRA, 4, pts);

        // Now four tetra from triangles on tetra faces
        pts[0] = center;
        pts[1] = e01;
        pts[2] = e02;
        pts[3] = e12;      
        output->InsertNextCell(VTK_TETRA, 4, pts);
        
        pts[1] = e01;
        pts[2] = e13;
        pts[3] = e03;   
        output->InsertNextCell(VTK_TETRA, 4, pts);

        pts[1] = e12;
        pts[2] = e23;
        pts[3] = e13;                     
        output->InsertNextCell(VTK_TETRA, 4, pts);
        
        pts[1] = e03;
        pts[2] = e23;
        pts[3] = e02;
        output->InsertNextCell(VTK_TETRA, 4, pts);
    
        
  
  /*   
        double v1[3];
        double v2[3];
        double v3[3];
        double v4[3];
        
        locpts->GetPoint(pts[0],v1);
        locpts->GetPoint(pts[1],v2);
        locpts->GetPoint(pts[2],v3);
        locpts->GetPoint(pts[3],v4);
        
        double J = tet_jacobian_det(v1, v2, v3, v4);
  */   
        
        numsubCells++;
        
    } // if cell in id list 
    else
    {
        pts[0] = p0;
        pts[1] = p1;
        pts[2] = p2;
        pts[3] = p3;
        output->InsertNextCell(VTK_TETRA, 4, pts);
      
    }     
  } //for all cells
  cell->Delete();

  vtkDebugMacro(<<"Subdivided " << numsubCells << " cells");
  
  locator->Delete();
  newPts->Delete();
  output->Squeeze();
}

void vtkSubdivideSelectTetra::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}


