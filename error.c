/* error.c: Error handling

   Copyright 2005 Bjoern Butscher, Hendrik Weimer

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

#include <stdlib.h>
#include <stdio.h>

#include "error.h"

void *
quantum_error_handler(void *f(int))
{
  static void *errfunc = 0;
  
  if(f)
    errfunc = f;

  return errfunc;
}

const char *
quantum_strerr(int errno)
{
  switch(errno)
    {
    case QUANTUM_SUCCESS:
      return "success";
    case QUANTUM_FAILURE:
      return "failure";
    case QUANTUM_ENOMEM:
      return "malloc failed";
    case QUANTUM_EMLARGE:
      return "matrix too large";
    case QUANTUM_EMSIZE:
      return "wrong matrix size";
    case QUANTUM_EHASHFULL:
      return "hash table full";
    case QUANTUM_EHERMITIAN:
      return "matrix not Hermitian";
    case QUANTUM_ENOCONVERGE:
      return "method failed to converge";
    case QUANTUM_ENOLAPACK:
      return "LAPACK support not compiled in";
    case QUANTUM_ELAPACKARG:
      return "wrong arguments supplied to LAPACK";
    case QUANTUM_ELAPACKCONV:
      return "LAPACK failed to converge";
    case QUANTUM_EMCMATRIX:
      return "single-column matrix expected";
    case QUANTUM_EOPCODE:
      return "unknown opcode";
    default:
      return "unknown error code";
    }
}

void
quantum_error(int errno)
{
  void (*p)(int);
  p = quantum_error_handler(0);
  
  if(p)
    p(errno);
  else
    {
      fflush(stdout);
      fprintf(stderr, "ERROR: %s\n", quantum_strerr(errno));
      fflush(stderr);

      abort();
    }
} 
