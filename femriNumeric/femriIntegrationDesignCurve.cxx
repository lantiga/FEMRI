/*=========================================================================

  Program:   FEMRI
  Module:    $RCSfile: femriIntegrationDesignCurve.cxx,v $
  Language:  C++
  Date:      $Date: 2008/10/28 16:46:59 $
  Version:   $Revision: 1.2 $

  Copyright (c) Luca Antiga, David Steinman. All rights reserved.
  See LICENCE file for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm 
  for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

// .NAME femriIntegrationDesignCurve - ..
// .SECTION Description
// ..
// .SECTION Thanks
// This class was developed by Paul Simedrea\n
// Imaging Research Labs - Robarts Research Institute \n
// email: simedrea@imaging.robarts.ca - homepage: http://www.imaging.robarts.ca/~simedrea

#include "femriIntegrationDesignCurve.h"

femriIntegrationDesignCurve::femriIntegrationDesignCurve()
{
  this->num_err_val = 1;  // give it a default value of 1;
}

femriIntegrationDesignCurve::~femriIntegrationDesignCurve()
{
  
}

femriIntegrationDesignCurve::femriIntegrationDesignCurve(
  const double num_err, vector<double> &cyc_elems, vector<unsigned int> &qorders)
{
  this->num_err_val = num_err;
  
  // must have same number of elems
  // should throw exception or something if not, ideally
  if (cyc_elems.size() == qorders.size())
  {
    this->cyc_elems = cyc_elems;
    this->qorders = qorders;
  }
  
}

femriIntegrationDesignCurve::femriIntegrationDesignCurve(const double num_err)
{  
  this->num_err_val = num_err;
}

const double femriIntegrationDesignCurve::GetNumError()
{
  return this->num_err_val;
}

void femriIntegrationDesignCurve::AddCurveEntry(const double cyc_elem, const unsigned int qorder)     
{
  this->AddQOrder(qorder);
  this->AddCycElem(cyc_elem);
}

void femriIntegrationDesignCurve::AddCycElem(const double cyc_elem)
{
  this->cyc_elems.push_back(cyc_elem);
}

void femriIntegrationDesignCurve::AddQOrder(const unsigned int qorder)
{
  this->qorders.push_back(qorder);
}

bool femriIntegrationDesignCurve::operator<(const femriIntegrationDesignCurve &cmp_curve)
{
  if (this->num_err_val < cmp_curve.num_err_val)
    return true;
  else
    return false;  
}

bool femriIntegrationDesignCurve::operator>(const femriIntegrationDesignCurve &cmp_curve)
{
  if (this->num_err_val > cmp_curve.num_err_val)
    return true;
  else
    return false;  
}

bool femriIntegrationDesignCurve::operator<=(const femriIntegrationDesignCurve &cmp_curve)
{
  if (this->num_err_val <= cmp_curve.num_err_val)
    return true;
  else
    return false;  
}

bool femriIntegrationDesignCurve::operator>=(const femriIntegrationDesignCurve &cmp_curve)
{
  if (this->num_err_val >= cmp_curve.num_err_val)
    return true;
  else
    return false;  
}

bool femriIntegrationDesignCurve::operator==(const femriIntegrationDesignCurve &cmp_curve)
{
  if (this->num_err_val == cmp_curve.num_err_val)
    return true;
  else
    return false;  
}


const double femriIntegrationDesignCurve::GetCycles(const unsigned int index)
{
  return this->cyc_elems.at(index);  
}    

const unsigned int femriIntegrationDesignCurve::GetQOrder(const unsigned int index)
{
  return this->qorders.at(index);
}

const unsigned int femriIntegrationDesignCurve::GetNumEntries()
{
  return this->qorders.size();
}
