/*=========================================================================

  Program:   FEMRI
  Module:    $RCSfile: vtkfemriKSpaceZeroPadder.cxx,v $
  Language:  C++
  Date:      $Date: 2007/03/19 13:31:25 $
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

#include "vtkfemriKSpaceZeroPadder.h"
#include "vtkImageData.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkObjectFactory.h"

vtkCxxRevisionMacro(vtkfemriKSpaceZeroPadder, "$Revision: 1.1.1.1 $");
vtkStandardNewMacro(vtkfemriKSpaceZeroPadder);

vtkfemriKSpaceZeroPadder::vtkfemriKSpaceZeroPadder() 
{
  this->PadSize[0] = this->PadSize[1] = this->PadSize[2] = 0;
}

int vtkfemriKSpaceZeroPadder::RequestInformation (
  vtkInformation * vtkNotUsed(request),
  vtkInformationVector** inputVector,
  vtkInformationVector *outputVector)
{
  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation* outInfo = outputVector->GetInformationObject(0);

  int inputWholeExtent[6];
  int wholeExtent[6];

  inInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),inputWholeExtent);

  wholeExtent[0] = inputWholeExtent[0];
  wholeExtent[1] = inputWholeExtent[1] + this->PadSize[0];
  wholeExtent[2] = inputWholeExtent[2];
  wholeExtent[3] = inputWholeExtent[3] + this->PadSize[1];
  wholeExtent[4] = inputWholeExtent[4];
  wholeExtent[5] = inputWholeExtent[5] + this->PadSize[2];

  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),wholeExtent,6);
  
  vtkDataObject::SetPointDataActiveScalarInfo(outInfo, VTK_DOUBLE, 2);
  return 1;
}

void vtkfemriKSpaceZeroPadder::SimpleExecute(vtkImageData* input, vtkImageData* output)
{
  vtkIdType i, j, k, q; 
  vtkIdType id; 
  int inputIjk[3];
  int outputIjk[3];
  int inputDimensions[3];
  int outputDimensions[3];
  int inputMids[3];
  int outputMids[3];
//  int quadrantFactors[4][2] = {{1,1}, {-1,1}, {1,-1}, {-1,-1}};
  int octantFactors[8][3] = {{1,1,1}, {-1,1,1}, {1,-1,1}, {-1,-1,1}, {1,1,-1}, {-1,1,-1}, {1,-1,-1}, {-1,-1,-1}};

  double re, im;
  vtkDataArray *outputScalars;

  outputScalars = output->GetPointData()->GetScalars();
  outputScalars->FillComponent(0,0.0);
  outputScalars->FillComponent(1,0.0);

  input->GetDimensions(inputDimensions);
  output->GetDimensions(outputDimensions);

  inputMids[0] = (inputDimensions[0] ) / 2;
  inputMids[1] = (inputDimensions[1] ) / 2;
  inputMids[2] = (inputDimensions[2] ) / 2;

  outputMids[0] = (outputDimensions[0] ) / 2;
  outputMids[1] = (outputDimensions[1] ) / 2;
  outputMids[2] = (outputDimensions[2] ) / 2;

  for (k=0; k<=inputMids[2]; k++)
    {
    for (j=0; j<=inputMids[1]; j++)
      {
      for (i=0; i<=inputMids[0]; i++)
        {
        for (q=0; q<8; q++)
          {
          inputIjk[0] = (octantFactors[q][0] * i + inputDimensions[0]) % inputDimensions[0];
          inputIjk[1] = (octantFactors[q][1] * j + inputDimensions[1]) % inputDimensions[1];
          inputIjk[2] = (octantFactors[q][2] * k + inputDimensions[2]) % inputDimensions[2];
          
          outputIjk[0] = (octantFactors[q][0] * i + outputDimensions[0]) % outputDimensions[0];
          outputIjk[1] = (octantFactors[q][1] * j + outputDimensions[1]) % outputDimensions[1];
          outputIjk[2] = (octantFactors[q][2] * k + outputDimensions[2]) % outputDimensions[2];
          
          if (((i == inputMids[0]) && (inputDimensions[0] % 2 == 0) && (octantFactors[q][0] == -1)) || 
              ((j == inputMids[1]) && (inputDimensions[1] % 2 == 0) && (octantFactors[q][1] == -1)) || 
              ((k == inputMids[2]) && (inputDimensions[2] % 2 == 0) && (octantFactors[q][2] == -1)))
            {
            continue;
            }

          id = output->ComputePointId(outputIjk);
          
          re = input->GetScalarComponentAsDouble(inputIjk[0],inputIjk[1],inputIjk[2],0);
          im = input->GetScalarComponentAsDouble(inputIjk[0],inputIjk[1],inputIjk[2],1);
          
          outputScalars->SetComponent(id,0,re);
          outputScalars->SetComponent(id,1,im);
          }
        }
      }
    }
}

void vtkfemriKSpaceZeroPadder::PrintSelf(ostream& os, vtkIndent indent)
{
  Superclass::PrintSelf(os,indent);
}
