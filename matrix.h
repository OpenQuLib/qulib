/* matrix.h: Declarations for matrix.c

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

#ifndef __MATRIX_H
#define __MATRIX_H

#include "config.h"

/* A ROWS x COLS matrix with complex elements */

struct quantum_matrix_struct {
  int rows;
  int cols;
  COMPLEX_FLOAT *t;
};
typedef struct quantum_matrix_struct quantum_matrix;

#define M(m,x,y) m.t[(x)+(y)*m.cols]

unsigned long quantum_memman(long change);

quantum_matrix quantum_new_matrix(int cols, int rows);
void quantum_delete_matrix(quantum_matrix *m);
void quantum_print_matrix(quantum_matrix m);
quantum_matrix quantum_mmult(quantum_matrix A, quantum_matrix B);
void quantum_adjoint(quantum_matrix *m);

#endif
