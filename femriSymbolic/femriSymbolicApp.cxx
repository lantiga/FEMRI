#include "vtkXMLUnstructuredGridReader.h"
#include "vtkUnstructuredGrid.h"
#include "vtkfemriUnstructuredGridKSpaceAcquisition.h"
#include "vtkMetaImageWriter.h"

#include "vtkImageRFFT.h"
#include "vtkImageMagnitude.h"
#include "vtkImageChangeInformation.h"

int main(int argc, char* argv[])
{
  char basePath[512];
  if (argc == 1)
    {
    sprintf(basePath,".");
    }
  else if (argc == 2)
    {
    sprintf(basePath,"%s",argv[1]);
    }
  else
    {
    cerr << "Wrong number of parameters." << endl;
    exit(1);
    }

  char meshFileName[] = "cubeMesh.vtk";
  char meshFilePath[512];
  sprintf(meshFilePath,"%s/%s",basePath,meshFileName);

  vtkXMLUnstructuredGridReader *meshReader;
  meshReader = vtkXMLUnstructuredGridReader::New();
  meshReader->SetFileName(meshFilePath);
  meshReader->Update();

  vtkUnstructuredGrid *mesh = vtkUnstructuredGrid::New();
  mesh->ShallowCopy(meshReader->GetOutput());
  mesh->Update();

  meshReader->Delete();
  
  int kSpaceDimensionality = 3;
  int kSpaceDimensions[3] = {16, 16, 16};
  double FOV[3] = {2.0, 2.0, 2.0};
  double origin[3] = {-1.0, -1.0, -1.0};
  char functionString[] = "1";

//   int kSpaceDimensionality = 2;
//   int kSpaceDimensions[3] = {16, 16, 1};
//   double FOV[3] = {2.0, 2.0, 1.0};
//   double origin[3] = {-1.0, -1.0, 0.0};
//   char functionString[] = "1";

  vtkfemriUnstructuredGridKSpaceAcquisition* kSpaceAcquisition = vtkfemriUnstructuredGridKSpaceAcquisition::New();
  kSpaceAcquisition->SetInput(mesh);
  kSpaceAcquisition->SetFunctionString(functionString);
  kSpaceAcquisition->SetKSpaceDimensionality(kSpaceDimensionality);
  kSpaceAcquisition->SetKSpaceDimensions(kSpaceDimensions);
  kSpaceAcquisition->SetFOV(FOV);
  kSpaceAcquisition->SetOrigin(origin);
  kSpaceAcquisition->UseElementIntegralCacheOn();
  kSpaceAcquisition->Update();

  vtkImageRFFT* ifft = vtkImageRFFT::New();
  ifft->SetInput(kSpaceAcquisition->GetOutput());
  ifft->SetDimensionality(kSpaceDimensionality);
  ifft->Update();

  vtkImageMagnitude* ifftMagnitude = vtkImageMagnitude::New();
  ifftMagnitude->SetInput(ifft->GetOutput()); 
  ifftMagnitude->Update();

  double spacing[3];
  spacing[0] = FOV[0] / double(kSpaceDimensions[0]);
  spacing[1] = FOV[1] / double(kSpaceDimensions[1]);
  spacing[2] = FOV[2] / double(kSpaceDimensions[2]);

  vtkImageChangeInformation* imageInformation = vtkImageChangeInformation::New();
  imageInformation->SetInput(ifftMagnitude->GetOutput());
  imageInformation->SetOutputSpacing(spacing);
  imageInformation->SetOutputOrigin(origin);
  imageInformation->Update();

  char imageFileName[] = "cubeImage.mhd";
  char imageFilePath[512];
  sprintf(imageFilePath,"%s/%s",basePath,imageFileName);

  vtkMetaImageWriter* imageWriter = vtkMetaImageWriter::New();
  imageWriter->SetInput(imageInformation->GetOutput());
  imageWriter->SetFileName(imageFilePath);
  imageWriter->Write();

  kSpaceAcquisition->Delete();
  ifft->Delete();
  ifftMagnitude->Delete();
  imageInformation->Delete();
  imageWriter->Delete();
  mesh->Delete();

  return 0;
}
