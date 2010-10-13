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

// .NAME vtkfemriUnstructuredGridMRSlicer - ..
// .SECTION Description
// ..
// .SECTION Thanks
// This class was developed by Paul Simedrea\n
// Imaging Research Labs - Robarts Research Institute \n
// email: simedrea@imaging.robarts.ca - homepage: http://www.imaging.robarts.ca/~simedrea

#ifndef FEMRIINTEGRATIONDESIGNCURVEINTERPOLATOR_H_
#define FEMRIINTEGRATIONDESIGNCURVEINTERPOLATOR_H_

#include "femriIntegrationDesignCurve.h"
#include <vector>

class femriIntegrationDesignCurveInterpolator
{
  public:
    /**
     * Add a design curve to the set.
     */
//    void AddDesignCurve(femriIntegrationDesignCurve *curve);
    void AddDesignCurve(const double num_err, vector<double> &cyc_elems, vector<unsigned int> &qorders);
   
    void ResetDesignCurves();
 
    /**
     * Interpolate and return the design curve for the given error value 
     */
    void Interpolate(femriIntegrationDesignCurve *new_curve);  
    
    /**
     *   Default constructor
     */
    femriIntegrationDesignCurveInterpolator();
    
    /**
     *  Default destructor
     */
    ~femriIntegrationDesignCurveInterpolator();
      
  protected:
    
    
  private:
    vector < femriIntegrationDesignCurve *> curve_set;
  
};

#endif /*FEMRIINTEGRATIONDESIGNCURVEINTERPOLATOR_H_*/
