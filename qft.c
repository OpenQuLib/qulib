/* qft.c: Quantum Fourier Transform
   
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
#include <mpi.h>
#include <math.h>
#include "gates.h"
#include "qureg.h"
#include "defs.h"

void
quantum_qft_qft(int width, quantum_qft_reg *reg)
{
	int i, j;

	for(i=width-1; i>=0; i--) {
		for(j=width-1; j>i; j--) {
			quantum_qft_cond_phase(j, i, reg);
		}
		quantum_qft_hadamard(i, reg);
	}
}

void
quantum_qft_qft_mpi(int width, quantum_qft_reg *reg)
{
	int i, j, id, numprocs;

	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	for(i=width-1; i>=0; i--) {
		for(j=width-1; j>i; j--) {
			quantum_qft_cond_phase_mpi(j, i, reg);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		quantum_qft_hadamard_mpi(i, reg);
		MPI_Barrier(MPI_COMM_WORLD);
	}
}