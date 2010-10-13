/*=========================================================================

  Program:   femri design curve
  Module:    $RCSfile: femriIntegrationDesignCurve.h,v $
  Language:  C++
  Date:      $Date: 2008/10/28 18:19:59 $
  Version:   $Revision: 1.3 $

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

#ifndef FEMRIDESIGNCURVE_H_
#define FEMRIDESIGNCURVE_H_

#include <vector>

/**
 *  Class containing an integration design curve
 *   @author Paul Simedrea
 * 
 */

using namespace std;

class femriIntegrationDesignCurve 
{
  public:
    
    /**
     * Constructor
     */
    femriIntegrationDesignCurve(
      const double num_err, vector<double> &cyc_elems, vector<unsigned int> &qorders);
  
    /**
     * Overloaded constructor, requiring only the numerical error
     */
    femriIntegrationDesignCurve(const double num_err);
        
    /**
     * Set the quadrature order list
     */
    void SetQOrders(const vector<unsigned int> &qorders);
    
    /**
     * Add a new curve entry (must be sequential)
     */
    void AddCurveEntry(const double cyc_elems, const unsigned int qorder);     
    
    /**
     * Set the cycles/element list to the curve
     */
    void SetCycElems(vector<double> &cyc_elems);
    
    /**
     * Set the numerical error value of the curve
     */
    void SetNumError(const double num_err);
    
    /**
     * Get the numerical error value of the curve
     */
    const double GetNumError();
    
    /**
     * Get the number of discrete entries in the curve
     */
    const unsigned int GetNumEntries();
    
    /**
     * LT comparison operator
     */  
    bool operator<(const femriIntegrationDesignCurve &cmp_curve);
    
    /**
     * LTE comparison operator
     */  
    bool operator<=(const femriIntegrationDesignCurve &cmp_curve);
      
    /**
     * GTE comparison operator
     */  
    bool operator>(const femriIntegrationDesignCurve &cmp_curve);
    
    /**
     * GTE comparison operator
     */   
    bool operator>=(const femriIntegrationDesignCurve &cmp_curve);

    /**
     * EQ comparison operator
     */  
    bool operator==(const femriIntegrationDesignCurve &cmp_curve);
    
    /**
     * Get the cycles/elem at specified point in curve 
     */
    const double GetCycles(const unsigned int index);
    
    /**
     * Get the cycles/elem at specified point in curve 
     */
    const unsigned int GetQOrder(const unsigned int index);
    
    
    femriIntegrationDesignCurve();
    ~femriIntegrationDesignCurve();
    
  protected:
    /**
     * Add a quadrature order to the curve
     */
    void AddQOrder(const unsigned int new_qorder);
    
    /**
     * Add a cycles/element entry to the curve
     */
    void AddCycElem(const double new_cyc_elem);
    

    
  private:
  
    vector<unsigned int> qorders;  //!< holds quadrature orders
    vector<double> cyc_elems;    //!< holds cycles/element values
    
    unsigned int num_entries;    //!< number of discrete values in curve
    double num_err_val;        //!< numerical error value for curve  
    
};  
  

#endif /*FEMRIDESIGNCURVE_H_*/
