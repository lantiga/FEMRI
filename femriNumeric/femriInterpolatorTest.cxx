#include "femriIntegrationDesignCurve.h"
#include "femriIntegrationDesignCurveInterpolator.h"
#include <vector>
#include <iostream>

using namespace std;

int main()
{
#if 0  
  const unsigned int MAX_MAP_QUAD_ORDERS = 20;   
  
  unsigned int qarr[MAX_MAP_QUAD_ORDERS] = {1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39};
  
  double ce4[MAX_MAP_QUAD_ORDERS] = {0.03, 0.21, 0.6, 1.08, 1.59, 2.16, 2.73, 3.33, 3.93, 4.53, 5.13, 5.73, 6.33, 6.93, 7.53, 8.13, 8.73, 9.33, 9.93, 10.53};
  double ce3[MAX_MAP_QUAD_ORDERS] = {0.06, 0.36, 0.87, 1.47, 2.07, 2.7, 3.33, 3.96, 4.62, 5.28, 5.94, 6.6, 7.26, 7.92, 8.58, 9.24, 9.90, 10.56, 11.22, 11.88}; 

  vector<double> cycElemE4;
  vector<double> cycElemE3;
  vector<unsigned int> QOrders;
  
  for(int i=0; i < MAX_MAP_QUAD_ORDERS; i++)
  {
    cycElemE4.push_back(ce4[i]);
    cycElemE3.push_back(ce3[i]);
    QOrders.push_back(qarr[i]);
  }
  
  femriIntegrationDesignCurveInterpolator interp = femriIntegrationDesignCurveInterpolator();
  
  interp.AddDesignCurve(1e-4, cycElemE4, QOrders);
  interp.AddDesignCurve(1e-3, cycElemE3, QOrders);
  
  femriIntegrationDesignCurve CEint = femriIntegrationDesignCurve(0.5e-3);

  interp.Interpolate(&CEint);
  
  cerr << "Cycles/Element (E=" << CEint.GetNumError() << "):" << endl;
  for (unsigned int i = 0; i < CEint.GetNumEntries(); i++)
  {
    cerr << " " << CEint.GetCycles(i);  
  }
  cerr << endl;
  
  cerr << "Quadrature order (E=" << CEint.GetNumError() << "):" << endl;
  for (unsigned int i = 0; i < CEint.GetNumEntries(); i++)
  {
    cerr << " " << CEint.GetQOrder(i);  
  }
  cerr << endl;
#endif  
  
  return 0;
}
