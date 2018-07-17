/* quantum.h: Header file for libquantum

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

#ifndef __QUANTUM_H

#define __QUANTUM_H

#define COMPLEX_FLOAT float _Complex
#define MAX_UNSIGNED unsigned long long

#define quantum_density_operation(function, rho, ...) \
do{ \
  int quantum_int; \
  for(quantum_int=0; quantum_int < rho.num; quantum_int++) \
    function(__VA_ARGS__, &rho.reg[quantum_int]); \
} while(0)

/* A ROWS x COLS matrix with complex elements */

struct quantum_matrix_struct {
  int rows;
  int cols;
  COMPLEX_FLOAT *t;
};

typedef struct quantum_matrix_struct quantum_matrix;

struct quantum_qft_reg_struct
{
	int width;
	COMPLEX_FLOAT *amplitude;
};
typedef struct quantum_qft_reg_struct quantum_qft_reg;

enum {
  QUANTUM_SOLVER_LANCZOS,
  QUANTUM_SOLVER_LANCZOS_MODIFIED,
  QUANTUM_SOLVER_IMAGINARY_TIME
};



extern quantum_matrix quantum_new_matrix(int cols, int rows);
extern void quantum_delete_matrix(quantum_matrix *m);
extern quantum_matrix quantum_mmult(quantum_matrix A, quantum_matrix B);

extern double quantum_prob(COMPLEX_FLOAT a);

extern const char * quantum_get_version();

extern void *quantum_error_handler(void *f(int));
extern const char *quantum_strerr(int errno);
extern void quantum_error(int errno);

extern quantum_qft_reg quantum_new_qft_reg(MAX_UNSIGNED initval, int width);
extern quantum_qft_reg quantum_new_qft_reg_mpi(MAX_UNSIGNED initval, int width);

extern void quantum_qft_print_reg(quantum_qft_reg *reg);

extern int quantum_qft_get_state(MAX_UNSIGNED num,quantum_qft_reg *reg);

extern void quantum_qft_gate1(int target, quantum_matrix m, quantum_qft_reg *reg);

extern void quantum_qft_gate1_mpi(int target, quantum_matrix m, quantum_qft_reg *reg);

extern void quantum_qft_hadamard(int target, quantum_qft_reg *reg);

extern void quantum_qft_hadamard_mpi(int target, quantum_qft_reg *reg);

extern void quantum_qft_cond_phase(int control, int target, quantum_qft_reg *reg);

extern void quantum_qft_cond_phase_mpi(int control, int target, quantum_qft_reg *reg);

extern void quantum_qft_qft(int width, quantum_qft_reg *reg);

extern void quantum_qft_qft_mpi(int width, quantum_qft_reg *reg);

extern void quantum_qft_delete_qureg(quantum_qft_reg *reg);

#endif
