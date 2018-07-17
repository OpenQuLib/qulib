/* qureg.h: Declarations for qureg.c and inline hashing functions

   Copyright 2003-2013 Bjoern Butscher, Hendrik Weimer

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

#ifndef __QUREG_H

#define __QUREG_H

#include <sys/types.h>

#include "config.h"
#include "matrix.h"
#include "error.h"

/* The quantum register */

struct quantum_qft_reg_struct
{
	int width;
	COMPLEX_FLOAT *amplitude;
};
typedef struct quantum_qft_reg_struct quantum_qft_reg;

extern quantum_qft_reg quantum_new_qft_reg(MAX_UNSIGNED initval, int width);
extern quantum_qft_reg quantum_new_qft_reg_mpi(MAX_UNSIGNED initval, int width);
extern void quantum_qft_print_reg(quantum_qft_reg *reg);
extern void quantum_qft_delete_qureg(quantum_qft_reg *reg);


#endif
