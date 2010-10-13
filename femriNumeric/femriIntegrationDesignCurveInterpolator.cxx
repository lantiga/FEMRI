/*=========================================================================

  Program:   FEMRI
  Module:    $RCSfile: femriIntegrationDesignCurveInterpolator.cxx,v $
  Language:  C++
  Date:      $Date: 2008/10/28 18:19:59 $
  Version:   $Revision: 1.4 $

  Copyright (c) Luca Antiga, David Steinman. All rights reserved.
  See LICENCE file for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm 
  for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

// .NAME femriIntegrationDesignCurveInterpolator - ..
// .SECTION Description
// ..
// .SECTION Thanks
// This class was developed by Paul Simedrea\n
// Imaging Research Labs - Robarts Research Institute \n
// email: simedrea@imaging.robarts.ca - homepage: http://www.imaging.robarts.ca/~simedrea

#include "femriIntegrationDesignCurveInterpolator.h"
#include <iostream>

using namespace std;

femriIntegrationDesignCurveInterpolator::femriIntegrationDesignCurveInterpolator()
{
  
}

femriIntegrationDesignCurveInterpolator::~femriIntegrationDesignCurveInterpolator()
{
  this->ResetDesignCurves(); 
}

void femriIntegrationDesignCurveInterpolator::AddDesignCurve(const double num_err, vector<double> &cyc_elems, vector<unsigned int> &qorders)
{
  this->curve_set.push_back(new femriIntegrationDesignCurve(num_err,cyc_elems,qorders));
}

void femriIntegrationDesignCurveInterpolator::ResetDesignCurves()
{
  for (unsigned int i=0; i<this->curve_set.size(); i++)
    {
    delete this->curve_set.at(i);
    }
  this->curve_set.resize(0);
}

//void femriIntegrationDesignCurveInterpolator::AddDesignCurve(femriIntegrationDesignCurve *curve)
//{
//  vector<femriIntegrationDesignCurve *>::iterator iter;
  
  // would be best to compare error and insert
//  if (curve != NULL)
//    this->curve_set.push_back(curve);
//}

void femriIntegrationDesignCurveInterpolator::Interpolate(femriIntegrationDesignCurve *new_curve)
{
  unsigned int q1, q2, qx;  // quadrature orders
  double c1, c2, cx;      // cycles/element
  double E1, E2, Ex;
  femriIntegrationDesignCurve *curr_curve = NULL;
  femriIntegrationDesignCurve *next_curve = NULL;
  
  // iterate through each curve in set, find the one closest to argument error value
  for (unsigned int vindex = 0; vindex < this->curve_set.size(); vindex++)
  {
    curr_curve = this->curve_set.at(vindex);
    
    if (vindex == this->curve_set.size() - 1)
    {
      cerr << "femriIntegrationDesignCurveInterpolator: ERROR nearest interpolant not found." << endl;
    }  
    else if (*new_curve > *curr_curve) 
    {  
      next_curve = this->curve_set.at(vindex + 1);
      
      for (unsigned int i=0; i < curr_curve->GetNumEntries(); i++)
      {
        q1 = curr_curve->GetQOrder(i);
        c1 = curr_curve->GetCycles(i);
        c2 = next_curve->GetCycles(i);
        E1 = curr_curve->GetNumError();
        E2 = next_curve->GetNumError();
        Ex = new_curve->GetNumError();
        
        // interpolate cycles/element
        cx = ((E2 - Ex)/(E2 - E1))*(c2 - c1) + c1;
        // quadrature order remains constant
        qx = q1;
        // add points to curve argument
        new_curve->AddCurveEntry(cx, qx);
      }    
      break;
    }
  } 
}
