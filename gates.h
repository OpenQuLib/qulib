/* gates.h: Declarations for gates.c

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

#ifndef __GATES_H

#define __GATES_H

#include "matrix.h"
#include "qureg.h"

extern void quantum_qft_gate1(int target, quantum_matrix m, quantum_qft_reg *reg);

extern void quantum_qft_gate1_mpi(int target, quantum_matrix m, quantum_qft_reg *reg);

extern void quantum_qft_hadamard(int target, quantum_qft_reg *reg);

extern void quantum_qft_hadamard_mpi(int target, quantum_qft_reg *reg);

extern void quantum_qft_cond_phase(int control, int target, quantum_qft_reg *reg);

extern void quantum_qft_cond_phase_mpi(int control, int target, quantum_qft_reg *reg);

#endif
