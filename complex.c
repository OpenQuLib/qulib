/* complex.c: Complex number functions

   Copyright 2003 Bjoern Butscher, Hendrik Weimer

   This file is part of libquantum

   libquantum is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published
   by the Free Software Foundation; either version 3 of the License,
   or (at your option) any later version.

   libquantum is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with libquantum; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
   MA 02110-1301, USA

*/

#include <math.h>
#include <complex.h>

#include "complex.h"
#include "config.h"

/* Return the complex conjugate of a complex number */

/*COMPLEX_FLOAT
quantum_conj(COMPLEX_FLOAT a)
{
  REAL_FLOAT r, i;

  r = quantum_real(a);
  i = quantum_imag(a);

  return r - IMAGINARY * i;
  }*/

/* Calculate the square of a complex number (i.e. the probability) */

double
quantum_prob(COMPLEX_FLOAT a)
{
  return quantum_prob_inline(a);
}

/* Calculate e^(i * phi) */

COMPLEX_FLOAT quantum_cexp(REAL_FLOAT phi)
{
  return cos(phi) + IMAGINARY * sin(phi);
}
