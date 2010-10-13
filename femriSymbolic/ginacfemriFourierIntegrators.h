/*=========================================================================

  Program:   feMRI
  Module:    $RCSfile: ginacfemriFourierIntegrators.h,v $
  Language:  C++
  Date:      $Date: 2007/03/19 13:31:26 $
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

// .NAME ginacfeMRIFourierIntegrators - ..
// .SECTION Description
// ..
// .SECTION Thanks
// This class was developed by Luca Antiga, PhD \n
// Imaging Research Labs - Robarts Research Institute \n
// Bioengineering Dept - Mario Negri Institute for Pharmacological Research \n 
// email: lantiga@imaging.robarts.ca - homepage: http://www.imaging.robarts.ca/~lantiga

#include "ginac.h"

using namespace GiNaC;
using namespace std;

struct calc_fourier_integral
{
  ex operator()(const ex& input_function, const symbol& v)
    {
      ex expanded_input_function = input_function.expand();

      if (input_function == 0)
        {
        return input_function;
        }

      if (is_a<add>(expanded_input_function))
        {
        ex function_sum = ex_to<add>(expanded_input_function);
        ex function_sum_integral;
        for (size_t i=0; i<function_sum.nops(); i++)
          {
          function_sum_integral += operator()(function_sum.op(i), v);
          }
        return function_sum_integral;
        }

      ex exp_term, polynomial_term;
      exp_term = 1;
      polynomial_term = 1;

      if (is_a<mul>(input_function))
        {
        ex function_mul = ex_to<mul>(input_function);
        for (size_t i=0; i<function_mul.nops(); i++)
          {
          split_into_poly_times_exp(function_mul.op(i),v,polynomial_term,exp_term);
          }
        }
      else
        {
        split_into_poly_times_exp(input_function,v,polynomial_term,exp_term);
        }

      if (polynomial_term.degree(v) < 0)
        {
        throw invalid_argument("calc_exp_integral: invalid input function.");
        }

      ex integral;

      if (exp_term != 1)
        {
        ex exp_integral = exp_term * exp_term / exp_term.diff(v);
        integral = polynomial_term * exp_integral;
        if (polynomial_term.degree(v) > 0)
          {
          integral += - operator()(polynomial_term.diff(v) * exp_integral, v);
          }
        }
      else
        {
        integral = numeric(1) / (polynomial_term.degree(v)+1) * polynomial_term * v;
        }

      return integral;
    }

  ex operator()(const ex& input_function, const symbol& v, const ex& va, const ex& vb)
    {
      ex function_indef_integral = operator()(input_function,v);
      return function_indef_integral.subs(v==vb,subs_options::no_pattern) - function_indef_integral.subs(v==va,subs_options::no_pattern);
    }

  void split_into_poly_times_exp(const ex& input_function, const symbol& v, ex& polynomial_term, ex& exp_term, bool initialize = false)
    {
      if (initialize)
        {
        polynomial_term = 1;
        exp_term = 1;
        }
      if ((input_function.info(info_flags::polynomial)) || 
          (input_function.info(info_flags::rational_function)) || 
          (input_function.info(info_flags::numeric)) ||
          (input_function.diff(v)==0))
        {
        polynomial_term *= input_function;
        }
      else if (input_function.match(exp(wild(0))))
        {
        exp_term *= input_function;
        }
      else
        {
        throw invalid_argument("baus! calc_exp_integral: invalid input function.");
        }
    }
};

struct calc_std_hexahedron_fourier_integral : public calc_fourier_integral
{
  ex operator()(const ex& input_function, const symbol& v0, const symbol& v1, const symbol& v2)
    {
      return calc_fourier_integral::operator()(calc_fourier_integral::operator()(calc_fourier_integral::operator()(input_function,v2,0,1),v1,0,1),v0,0,1);
    }
};

struct calc_std_wedge_fourier_integral : public calc_fourier_integral
{
  ex operator()(const ex& input_function, const symbol& v0, const symbol& v1, const symbol& v2)
    {
      return calc_fourier_integral::operator()(calc_fourier_integral::operator()(calc_fourier_integral::operator()(input_function,v2,0,1),v1,0,1-v0),v0,0,1);
    }
};

struct calc_std_pyramid_fourier_integral : public calc_fourier_integral
{
  ex operator()(const ex& input_function, const symbol& v0, const symbol& v1, const symbol& v2)
    {
      return calc_fourier_integral::operator()(calc_fourier_integral::operator()(calc_fourier_integral::operator()(input_function,v2,0,1-v0-v1),v1,0,1),v0,0,1);
    }
};

struct calc_std_tetrahedron_fourier_integral : public calc_fourier_integral
{
  ex operator()(const ex& input_function, const symbol& v0, const symbol& v1, const symbol& v2)
    {
      return calc_fourier_integral::operator()(calc_fourier_integral::operator()(calc_fourier_integral::operator()(input_function,v2,0,1-v0-v1),v1,0,1-v0),v0,0,1);
    }
};

struct calc_std_quadrilateral_fourier_integral : public calc_fourier_integral
{
  ex operator()(const ex& input_function, const symbol& v0, const symbol& v1)
    {
      return calc_fourier_integral::operator()(calc_fourier_integral::operator()(input_function,v1,0,1),v0,0,1);
    }
};

struct calc_std_triangle_fourier_integral : public calc_fourier_integral
{
  ex operator()(const ex& input_function, const symbol& v0, const symbol& v1)
    {
      return calc_fourier_integral::operator()(calc_fourier_integral::operator()(input_function,v1,0,1-v0),v0,0,1);
    }
};

struct calc_std_line_fourier_integral : public calc_fourier_integral
{
  ex operator()(const ex& input_function, const symbol& v0)
    {
      return calc_fourier_integral::operator()(input_function,v0,0,1);
    }
};
