#include "ginac.h"
using namespace GiNaC;
using namespace std;
#include "vtkfemriUnstructuredGridFourierIntegrator.h"
#include "ginacfemriFourierIntegrators.h"

#include "vtkXMLUnstructuredGridReader.h"
#include "vtkUnstructuredGrid.h"

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

  meshReader->UnRegisterAllOutputs();
  meshReader->Delete();

  idx i(symbol("i"),3);
  ex K = symbolic_matrix(3,1,"k");
  ex X = symbolic_matrix(3,1,"x");
  ex function_3d =  exp(-I*(indexed(K,i) * indexed(X,i)).simplify_indexed());

  vtkfemriUnstructuredGridFourierIntegrator* meshIntegrator = vtkfemriUnstructuredGridFourierIntegrator::New();
  meshIntegrator->SetMesh(mesh);
  meshIntegrator->SetFunction(function_3d);
  meshIntegrator->SetX(X);
  meshIntegrator->SetK(K);

  lst freqs = lst(3*Pi,0,0);
  double fourier_integral_re, fourier_integral_im;
  meshIntegrator->ComputeIntegral(lst(K[0]==freqs[0], K[1]==freqs[1], K[2]==freqs[2]), fourier_integral_re, fourier_integral_im);
  ex mesh_integral = fourier_integral_re + I*fourier_integral_im;
  cout<<"Mesh integral: "<<mesh_integral<<endl;

  meshIntegrator->Delete();
  mesh->Delete();

  symbol c("c");
  symbol x("x");
  ex analytic_function = exp(-I*c*x);
  symbol v0("v0");
  calc_std_line_fourier_integral do_line_integral;
  ex line_integral_symbolic = do_line_integral(analytic_function.subs(lst(c==freqs[0],x==v0-numeric(1)/2)),v0);
  ex line_integral = line_integral_symbolic.evalf();
  cout<<"Line integral: "<< line_integral_symbolic << " = " << line_integral << endl;
  ex analytic_integral = I / c * exp(-I*c*x);
  ex analytic_definite_integral = (analytic_integral.subs(lst(c==freqs[0],x==numeric(1)/2)) - analytic_integral.subs(lst(c==freqs[0],x==-numeric(1)/2))).evalf();
  cout<<"Analytic definite 1d integral: "<<analytic_definite_integral<<endl;

  numeric tolerance = 1E-12;
  if (!(abs(mesh_integral - line_integral) < tolerance && abs(mesh_integral - analytic_definite_integral) < tolerance))
    {
    cout<<"Test failed."<<endl;
    cout<<"mesh_integral - line_integral = "<<mesh_integral - line_integral<<endl;
    cout<<"mesh_integral - analytic_definite_integral = "<<mesh_integral - line_integral<<endl;
    return 1;
    }

  cout<<"Test passed (tolerance: "<<tolerance<<")."<<endl;
  return 0;
}
