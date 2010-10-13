/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkSubdivideSelectQuadraticTetra.cxx,v $

  Extended by Paul Simedrea (2006)

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkSubdivideSelectQuadraticTetra.h"

#include "vtkCellType.h"
#include "vtkGenericCell.h"
#include "vtkMergePoints.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkQuadraticEdge.h"

vtkCxxRevisionMacro(vtkSubdivideSelectQuadraticTetra, "$Revision: 1.1.1.1 $");
vtkStandardNewMacro(vtkSubdivideSelectQuadraticTetra);

// Description:
// Construct with all types of clipping turned off.
vtkSubdivideSelectQuadraticTetra::vtkSubdivideSelectQuadraticTetra()
{
    this->SubdivideAllCellsOff();
}



void vtkSubdivideSelectQuadraticTetra::Execute()
{
  vtkUnstructuredGrid *input=(vtkUnstructuredGrid *)this->GetInput();
  vtkIdType numPts = input->GetNumberOfPoints();
  vtkIdType numCells = input->GetNumberOfCells();
  vtkIdType numsubCells = 0;
  vtkPoints *inPts=input->GetPoints();
  vtkIdType cellId;
  vtkIdType pts[10];
  vtkGenericCell *cell;
  vtkPointData *pd = input->GetPointData();
  vtkUnstructuredGrid *output = this->GetOutput();
  vtkPointData *outputPD = output->GetPointData();
  vtkPoints *newPts;
  vtkIdType ptId;

  double weights[10], x0[3], x1[3], x2[3], x3[3], x[3];
  double xm0[3], xm1[3], xm2[3], xm3[3], xm4[3], xm5[3];
  
  
  double pcenter[3], xcenter[3];
  
  int p0, p1, p2, p3, m0, m1, m2, m3, m4, m5;
  vtkIdType center, e01, e02, e03, e12, e13, e23;
  vtkIdType tokenId;
  vtkMergePoints *locator;
  
   
  
  vtkDebugMacro(<<"Executing mesh subdivide");

  // Check if all cells are same type and tetrahedrons
  // Change this possibly later on to accomodate inhomogeneous meshes.   
  if (input->IsHomogeneous() == 0 ||
      input->GetCellType(0) != VTK_QUADRATIC_TETRA)
    {
    vtkErrorMacro(<<"all cells must be quadratic tetrahedra.");
    return;
    }

  // Copy original points and point data
  newPts = vtkPoints::New();
  newPts->Allocate(5*numPts,numPts);
  outputPD->InterpolateAllocate(pd,5*numPts,numPts);
  
  output->Allocate(numCells);
  output->SetPoints(newPts);

  locator = vtkMergePoints::New();
  locator->InitPointInsertion (newPts, input->GetBounds(),5*numPts);

  // copy points from input to output
  for (ptId=0; ptId < numPts; ptId++)
    {
    // add existing points in the locator so that any copies are later detected.
    locator->InsertNextPoint(inPts->GetPoint(ptId));
    outputPD->CopyData(pd,ptId,ptId);
    }

  cell = vtkGenericCell::New();
  // loop over tetrahedra, generating sixteen new ones for each. This is
  // done by introducing mid-edge nodes and a single mid-tetra node.
  
  
  for(cellId=0; cellId < numCells; cellId++)
    {
    
    input->GetCell(cellId, cell);
        
    // get tetra point locations
    cell->Points->GetPoint(0,x0);
    cell->Points->GetPoint(1,x1);
    cell->Points->GetPoint(2,x2);
    cell->Points->GetPoint(3,x3);

    // get tetra midpoint node locations
    cell->Points->GetPoint(4,xm0);
    cell->Points->GetPoint(5,xm1);
    cell->Points->GetPoint(6,xm2);
    cell->Points->GetPoint(7,xm3);
    cell->Points->GetPoint(8,xm4);
    cell->Points->GetPoint(9,xm5);
      
    // get tetra node ids
    p0 = cell->PointIds->GetId(0);
    p1 = cell->PointIds->GetId(1);
    p2 = cell->PointIds->GetId(2);
    p3 = cell->PointIds->GetId(3);

    // get tetra midpoint node ids
    m0 = cell->PointIds->GetId(4);
    m1 = cell->PointIds->GetId(5);
    m2 = cell->PointIds->GetId(6);
    m3 = cell->PointIds->GetId(7);
    m4 = cell->PointIds->GetId(8);
    m5 = cell->PointIds->GetId(9);


    //std::cerr << this->IdList->IsId(cellId) << std::endl;
    if ((this->IdList->IsId(cellId) >= 0) || this->GetSubdivideAllCells())
        {   
        
        double xweights[10];
        double xcenter[3];
        vtkIdType tokenId;
        vtkIdType subId;
        double *pcoords;
        
        double xsm0[3], xsm1[3], xsm2[3], xsm3[3], xsm4[3], xsm5[3];
        double xsm6[3], xsm7[3], xsm8[3], xsm9[3], xsm10[3], xsm11[3];
        double xsm12[3], xsm13[3], xsm14[3], xsm15[3], xsm16[3], xsm17[3];
        double xsm18[3], xsm19[3], xsm20[3], xsm21[3], xsm22[3], xsm23[3];
         
  
        // sub-tetra 0 midpoints
        double pm0[3]; pm0[0] = 0.25; pm0[1] = 0.0;  pm0[2] = 0.0;
        double pm1[3]; pm1[0] = 0.25; pm1[1] = 0.25; pm1[2] = 0.0;
        double pm2[3]; pm2[0] = 0.0;  pm2[1] = 0.25; pm2[2] = 0.0;
        double pm3[3]; pm3[0] = 0.0;  pm3[1] = 0.0;  pm3[2] = 0.25;
        double pm4[3]; pm4[0] = 0.25; pm4[1] = 0.0;  pm4[2] = 0.25;
        double pm5[3]; pm5[0] = 0.0;  pm5[1] = 0.25; pm5[2] = 0.25;
                
        // sub-tetra 1 midpoints
        double pm6[3]; pm6[0] = 0.25; pm6[1] = 0.0;  pm6[2] = 0.5;
        double pm7[3]; pm7[0] = 0.25; pm7[1] = 0.25; pm7[2] = 0.5;
        double pm8[3]; pm8[0] = 0.0;  pm8[1] = 0.25; pm8[2] = 0.5;
        double pm9[3]; pm9[0] = 0.0;  pm9[1] = 0.0;  pm9[2] = 0.75;
        double pm10[3]; pm10[0] = 0.25; pm10[1] = 0.0;  pm10[2] = 0.75;
        double pm11[3]; pm11[0] = 0.0;  pm11[1] = 0.25; pm11[2] = 0.75;
        
        // sub-tetra 2 midpoints
        double pm12[3]; pm12[0] = 0.75; pm12[1] = 0.0;  pm12[2] = 0.0;
        double pm13[3]; pm13[0] = 0.75; pm13[1] = 0.25; pm13[2] = 0.0;
        double pm14[3]; pm14[0] = 0.5;  pm14[1] = 0.25; pm14[2] = 0.0;
        double pm15[3]; pm15[0] = 0.5;  pm15[1] = 0.0;  pm15[2] = 0.25;
        double pm16[3]; pm16[0] = 0.75; pm16[1] = 0.0;  pm16[2] = 0.25;
        double pm17[3]; pm17[0] = 0.5;  pm17[1] = 0.25; pm17[2] = 0.25;
     
        // sub-tetra 3 midpoints
        double pm18[3]; pm18[0] = 0.25; pm18[1] = 0.5;  pm18[2] = 0.0;
        double pm19[3]; pm19[0] = 0.25; pm19[1] = 0.75; pm19[2] = 0.0;
        double pm20[3]; pm20[0] = 0.0;  pm20[1] = 0.75; pm20[2] = 0.0;
        double pm21[3]; pm21[0] = 0.0;  pm21[1] = 0.5;  pm21[2] = 0.25;
        double pm22[3]; pm22[0] = 0.25; pm22[1] = 0.5;  pm22[2] = 0.25;
        double pm23[3]; pm23[0] = 0.0;  pm23[1] = 0.75; pm23[2] = 0.25;
        
       
       
        pcoords = cell->GetParametricCoords();
           
        double pnode0[3]; pnode0[0] = pcoords[0]; pnode0[1] = pcoords[1]; pnode0[2] = pcoords[2];
        double pnode1[3]; pnode1[0] = pcoords[3]; pnode1[1] = pcoords[4]; pnode1[2] = pcoords[5];
        double pnode2[3]; pnode2[0] = pcoords[6]; pnode2[1] = pcoords[7]; pnode2[2] = pcoords[8];
        double pnode3[3]; pnode3[0] = pcoords[9]; pnode3[1] = pcoords[10]; pnode3[2] = pcoords[11];
       
        double pnode4[3]; pnode4[0] = pcoords[12]; pnode4[1] = pcoords[13]; pnode4[2] = pcoords[14];
        double pnode5[3]; pnode5[0] = pcoords[15]; pnode5[1] = pcoords[16]; pnode5[2] = pcoords[17];
        double pnode6[3]; pnode6[0] = pcoords[18]; pnode6[1] = pcoords[19]; pnode6[2] = pcoords[20];
        double pnode7[3]; pnode7[0] = pcoords[21]; pnode7[1] = pcoords[22]; pnode7[2] = pcoords[23];
        double pnode8[3]; pnode8[0] = pcoords[24]; pnode8[1] = pcoords[25]; pnode8[2] = pcoords[26];
        double pnode9[3]; pnode9[0] = pcoords[27]; pnode9[1] = pcoords[28]; pnode9[2] = pcoords[29];
       
       
         // compute center point (10)
        cell->GetParametricCenter(pcenter);
        cell->EvaluateLocation(subId,pcenter,xcenter,xweights);                 
         
        locator->InsertUniquePoint(xcenter,center);
        outputPD->InterpolatePoint(pd, center, cell->PointIds, xweights);
    
        // compute edge sub-midpoints 
               
        // sub-tetra 0 midpoints
        
        // edge sub-midpoint sm0 (11)
        vtkIdType xm0id;
        cell->EvaluateLocation(subId,pm0, xsm0, xweights);
        locator->InsertUniquePoint(xsm0,xm0id);
        outputPD->InterpolatePoint(pd, xm0id, cell->PointIds, xweights);
        
        // edge sub-midpoint sm1 (12)
        vtkIdType xm1id;
        cell->EvaluateLocation(subId,pm1, xsm1, xweights);
        locator->InsertUniquePoint(xsm1,xm1id);
        outputPD->InterpolatePoint(pd, xm1id, cell->PointIds, xweights);
       
        // edge sub-midpoint sm2 (13)
        vtkIdType xm2id;
        cell->EvaluateLocation(subId,pm2, xsm2, xweights);
        locator->InsertUniquePoint(xsm2,xm2id);
        outputPD->InterpolatePoint(pd, xm2id, cell->PointIds, xweights);
       
        // edge sub-midpoint sm3 (14)
        vtkIdType xm3id;
        cell->EvaluateLocation(subId,pm3, xsm3, xweights);
        locator->InsertUniquePoint(xsm3,xm3id);
        outputPD->InterpolatePoint(pd, xm3id, cell->PointIds, xweights);
        
        // edge sub-midpoint sm4 (15)
        vtkIdType xm4id;
        cell->EvaluateLocation(subId,pm4, xsm4, xweights);
        locator->InsertUniquePoint(xsm4,xm4id);
        outputPD->InterpolatePoint(pd, xm4id, cell->PointIds, xweights);
                
        // edge sub-midpoint sm5 (16)
        vtkIdType xm5id;
        cell->EvaluateLocation(subId,pm5, xsm5, xweights);
        locator->InsertUniquePoint(xsm5,xm5id);
        outputPD->InterpolatePoint(pd, xm5id, cell->PointIds, xweights);
        
        
        // sub-tetra 1 midpoints
        
        // edge sub-midpoint sm6 (17)
        vtkIdType xm6id;
        cell->EvaluateLocation(subId,pm6, xsm6, xweights);
        locator->InsertUniquePoint(xsm6,xm6id);
        outputPD->InterpolatePoint(pd, xm6id, cell->PointIds, xweights);
        
        // edge sub-midpoint sm7 (18)
        vtkIdType xm7id;
        cell->EvaluateLocation(subId,pm7, xsm7, xweights);
        locator->InsertUniquePoint(xsm7,xm7id);
        outputPD->InterpolatePoint(pd, xm7id, cell->PointIds, xweights);
       
        // edge sub-midpoint sm8 (19)
        vtkIdType xm8id;
        cell->EvaluateLocation(subId,pm8, xsm8, xweights);
        locator->InsertUniquePoint(xsm8,xm8id);
        outputPD->InterpolatePoint(pd, xm8id, cell->PointIds, xweights);
       
        // edge sub-midpoint sm9 (20)
        vtkIdType xm9id;
        cell->EvaluateLocation(subId,pm9, xsm9, xweights);
        locator->InsertUniquePoint(xsm9,xm9id);
        outputPD->InterpolatePoint(pd, xm9id, cell->PointIds, xweights);
        
        // edge sub-midpoint sm10 (21)
        vtkIdType xm10id;
        cell->EvaluateLocation(subId,pm10, xsm10, xweights);
        locator->InsertUniquePoint(xsm10,xm10id);
        outputPD->InterpolatePoint(pd, xm10id, cell->PointIds, xweights);
                
        // edge sub-midpoint sm11 (22)
        vtkIdType xm11id;
        cell->EvaluateLocation(subId,pm11, xsm11, xweights);
        locator->InsertUniquePoint(xsm11,xm11id);
        outputPD->InterpolatePoint(pd, xm11id, cell->PointIds, xweights);
        
        
        
        // sub-tetra 2 midpoints
        
        // edge sub-midpoint sm12 (23)
        vtkIdType xm12id;
        cell->EvaluateLocation(subId,pm12, xsm12, xweights);
        locator->InsertUniquePoint(xsm12,xm12id);
        outputPD->InterpolatePoint(pd, xm12id, cell->PointIds, xweights);
        
        // edge sub-midpoint sm13 (24)
        vtkIdType xm13id;
        cell->EvaluateLocation(subId,pm13, xsm13, xweights);
        locator->InsertUniquePoint(xsm13,xm13id);
        outputPD->InterpolatePoint(pd, xm13id, cell->PointIds, xweights);
       
        // edge sub-midpoint sm14 (25)
        vtkIdType xm14id;
        cell->EvaluateLocation(subId,pm14, xsm14, xweights);
        locator->InsertUniquePoint(xsm14,xm14id);
        outputPD->InterpolatePoint(pd, xm14id, cell->PointIds, xweights);
       
        // edge sub-midpoint sm15 (26)
        vtkIdType xm15id;
        cell->EvaluateLocation(subId,pm15, xsm15, xweights);
        locator->InsertUniquePoint(xsm15,xm15id);
        outputPD->InterpolatePoint(pd, xm15id, cell->PointIds, xweights);
        
        // edge sub-midpoint sm16 (27)
        vtkIdType xm16id;
        cell->EvaluateLocation(subId,pm16, xsm16, xweights);
        locator->InsertUniquePoint(xsm16,xm16id);
        outputPD->InterpolatePoint(pd, xm16id, cell->PointIds, xweights);
                
        // edge sub-midpoint sm17 (28)
        vtkIdType xm17id;
        cell->EvaluateLocation(subId,pm17, xsm17, xweights);
        locator->InsertUniquePoint(xsm17,xm17id);
        outputPD->InterpolatePoint(pd, xm17id, cell->PointIds, xweights);
        
  
        // sub-tetra 3 midpoints
        
        // edge sub-midpoint sm18 (29)
        vtkIdType xm18id;
        cell->EvaluateLocation(subId,pm18, xsm18, xweights);
        locator->InsertUniquePoint(xsm18,xm18id);
        outputPD->InterpolatePoint(pd, xm18id, cell->PointIds, xweights);
        
        // edge sub-midpoint sm19 (30)
        vtkIdType xm19id;
        cell->EvaluateLocation(subId,pm19, xsm19, xweights);
        locator->InsertUniquePoint(xsm19,xm19id);
        outputPD->InterpolatePoint(pd, xm19id, cell->PointIds, xweights);
       
        // edge sub-midpoint sm20 (31)
        vtkIdType xm20id;
        cell->EvaluateLocation(subId,pm20, xsm20, xweights);
        locator->InsertUniquePoint(xsm20,xm20id);
        outputPD->InterpolatePoint(pd, xm20id, cell->PointIds, xweights);
       
        // edge sub-midpoint sm21 (32)
        vtkIdType xm21id;
        cell->EvaluateLocation(subId,pm21, xsm21, xweights);
        locator->InsertUniquePoint(xsm21,xm21id);
        outputPD->InterpolatePoint(pd, xm21id, cell->PointIds, xweights);
        
        // edge sub-midpoint sm22 (33)
        vtkIdType xm22id;
        cell->EvaluateLocation(subId,pm22, xsm22, xweights);
        locator->InsertUniquePoint(xsm22,xm22id);
        outputPD->InterpolatePoint(pd, xm22id, cell->PointIds, xweights);
                
        // edge sub-midpoint sm23 (34)
        vtkIdType xm23id;
        cell->EvaluateLocation(subId,pm23, xsm23, xweights);
        locator->InsertUniquePoint(xsm23,xm23id);
        outputPD->InterpolatePoint(pd, xm23id, cell->PointIds, xweights);
        
  
        // internal sub-midpoint sm24 (35)
        vtkIdType xm24id;
        double xim0[3];
        for (int i=0; i<3; i++)
        {
            x[i] = 0.5 * (pcenter[i] + pnode4[i]);
        }
        cell->EvaluateLocation(subId,x, xim0, xweights);
        locator->InsertUniquePoint(xim0,xm24id);
        outputPD->InterpolatePoint(pd, xm24id, cell->PointIds, xweights);        

        // internal sub-midpoint sm25 (36)
        vtkIdType xm25id;
        double xim1[3];
        for (int i=0; i<3; i++)
        {
            x[i] = 0.5 * (pcenter[i] + pnode5[i]);
        }
        cell->EvaluateLocation(subId,x, xim1, xweights);
        locator->InsertUniquePoint(xim1,xm25id);
        outputPD->InterpolatePoint(pd, xm25id, cell->PointIds, xweights);
        
        // internal sub-midpoint sm26 (37)
        vtkIdType xm26id;
        double xim2[3];
        for (int i=0; i<3; i++)
        {
            x[i] = 0.5 * (pcenter[i] + pnode6[i]);
        }
        cell->EvaluateLocation(subId,x, xim2, xweights);
        locator->InsertUniquePoint(xim2,xm26id);
        outputPD->InterpolatePoint(pd, xm26id, cell->PointIds, xweights);
        
        // internal sub-midpoint sm27 (38)
        vtkIdType xm27id;
        double xim3[3];
        for (int i=0; i<3; i++)
        {
            x[i] = 0.5 * (pcenter[i] + pnode7[i]);
        }
        cell->EvaluateLocation(subId,x, xim3, xweights);
        locator->InsertUniquePoint(xim3,xm27id);
        outputPD->InterpolatePoint(pd, xm27id, cell->PointIds, xweights);
        
        // internal sub-midpoint sm28 (39)
        vtkIdType xm28id;
        double xim4[3];
        for (int i=0; i<3; i++)
        {
            x[i] = 0.5 * (pcenter[i] + pnode8[i]);
        }
        cell->EvaluateLocation(subId,x, xim4, xweights);
        locator->InsertUniquePoint(xim4,xm28id);
        outputPD->InterpolatePoint(pd, xm28id, cell->PointIds, xweights);
        
        // internal sub-midpoint sm24 (40)
        vtkIdType xm29id;
        double xim5[3];
        for (int i=0; i<3; i++)
        {
            x[i] = 0.5 * (pcenter[i] + pnode9[i]);
        }
        cell->EvaluateLocation(subId,x, xim5, xweights);
        locator->InsertUniquePoint(xim5,xm29id);
        outputPD->InterpolatePoint(pd, xm29id, cell->PointIds, xweights);
        
        

        // Now create tetrahedra
        // First, four tetra from each vertex

        
        pts[0] = p0;
        pts[1] = m0;
        pts[2] = m2;
        pts[3] = m3;
        pts[4] = xm0id;
        pts[5] = xm1id;
        pts[6] = xm2id;
        pts[7] = xm3id;
        pts[8] = xm4id;
        pts[9] = xm5id;        
        output->InsertNextCell(VTK_QUADRATIC_TETRA, 10, pts);


        pts[0] = m0;
        pts[1] = p1;
        pts[2] = m1;
        pts[3] = m4;
        pts[4] = xm12id;
        pts[5] = xm13id;
        pts[6] = xm14id;
        pts[7] = xm15id;
        pts[8] = xm16id;
        pts[9] = xm17id;        
        output->InsertNextCell(VTK_QUADRATIC_TETRA, 10, pts);

        pts[0] = m1;
        pts[1] = p2;
        pts[2] = m2;
        pts[3] = m5;
        pts[4] = xm19id;
        pts[5] = xm20id;
        pts[6] = xm18id;
        pts[7] = xm22id;
        pts[8] = xm23id;
        pts[9] = xm21id;        
        output->InsertNextCell(VTK_QUADRATIC_TETRA, 10, pts);

        pts[0] = m3;
        pts[1] = m4;
        pts[2] = m5;
        pts[3] = p3;
        pts[4] = xm6id;
        pts[5] = xm7id;
        pts[6] = xm8id;
        pts[7] = xm9id;
        pts[8] = xm10id;
        pts[9] = xm11id;        
        output->InsertNextCell(VTK_QUADRATIC_TETRA, 10, pts);

        
        // Now four tetra from cut-off tetra corners
       
        pts[0] = m0;
        pts[1] = center;
        pts[2] = m2;
        pts[3] = m3;
        pts[4] = xm24id;
        pts[5] = xm26id;
        pts[6] = xm1id;
        pts[7] = xm4id;
        pts[8] = xm27id;
        pts[9] = xm5id;        
        output->InsertNextCell(VTK_QUADRATIC_TETRA, 10, pts);
        
        pts[0] = m0;
        pts[1] = m1;
        pts[2] = center;
        pts[3] = m4;
        pts[4] = xm14id;
        pts[5] = xm25id;
        pts[6] = xm24id;
        pts[7] = xm15id;
        pts[8] = xm17id;
        pts[9] = xm28id;        
        output->InsertNextCell(VTK_QUADRATIC_TETRA, 10, pts);

        pts[0] = m1;
        pts[1] = m2;
        pts[2] = center;
        pts[3] = m5;
        pts[4] = xm18id;
        pts[5] = xm26id;
        pts[6] = xm25id;
        pts[7] = xm22id;
        pts[8] = xm21id;
        pts[9] = xm29id;        
        output->InsertNextCell(VTK_QUADRATIC_TETRA, 10, pts);


        pts[0] = m3;
        pts[1] = m5;
        pts[2] = m4;
        pts[3] = center;
        pts[4] = xm8id;
        pts[5] = xm7id;
        pts[6] = xm6id;
        pts[7] = xm27id;
        pts[8] = xm29id;
        pts[9] = xm28id;        
        output->InsertNextCell(VTK_QUADRATIC_TETRA, 10, pts);
   

        
        // Now four tetra from triangles on tetra faces
        pts[0] = m0;
        pts[1] = m1;
        pts[2] = m2;
        pts[3] = center;
        pts[4] = xm14id;
        pts[5] = xm18id;
        pts[6] = xm1id;
        pts[7] = xm24id;
        pts[8] = xm25id;
        pts[9] = xm26id;        
        output->InsertNextCell(VTK_QUADRATIC_TETRA, 10, pts);

        pts[0] = m0;
        pts[1] = m4;
        pts[2] = m3;
        pts[3] = center;
        pts[4] = xm15id;
        pts[5] = xm6id;
        pts[6] = xm4id;
        pts[7] = xm24id;
        pts[8] = xm28id;
        pts[9] = xm27id;        
        output->InsertNextCell(VTK_QUADRATIC_TETRA, 10, pts);
        
    
        pts[0] = m1;
        pts[1] = m5;
        pts[2] = m4;
        pts[3] = center;
        pts[4] = xm22id;
        pts[5] = xm7id;
        pts[6] = xm17id;
        pts[7] = xm25id;
        pts[8] = xm29id;
        pts[9] = xm28id;        
        output->InsertNextCell(VTK_QUADRATIC_TETRA, 10, pts);

        pts[0] = m2;
        pts[1] = m3;
        pts[2] = m5;
        pts[3] = center;
        pts[4] = xm5id;
        pts[5] = xm8id;
        pts[6] = xm21id;
        pts[7] = xm26id;
        pts[8] = xm27id;
        pts[9] = xm29id;        
        output->InsertNextCell(VTK_QUADRATIC_TETRA, 10, pts);

     
    
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
        pts[4] = m0;
        pts[5] = m1;
        pts[6] = m2;
        pts[7] = m3;
        pts[8] = m4;
        pts[9] = m5;
        output->InsertNextCell(VTK_QUADRATIC_TETRA, 10, pts);
    
    }    
    
   
     
  } //for all cells
 
  cell->Delete();
  

  vtkDebugMacro(<<"Subdivided " << numsubCells << " cells");
   

  locator->Delete();  
  newPts->Delete();
  output->Squeeze();
 
}

void vtkSubdivideSelectQuadraticTetra::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}


