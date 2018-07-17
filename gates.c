/* gates.c: Basic gates for quantum register manipulation

   Copyright 2003, 2004 Bjoern Butscher, Hendrik Weimer

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
#include <math.h>
#include <stdarg.h>
#include <mpi.h>

#include "matrix.h"
#include "defs.h"
#include "complex.h"
#include "qureg.h"
#include "error.h"

int
quantum_qft_get_state(MAX_UNSIGNED num,quantum_qft_reg *reg){
	if(num > (1 << reg->width))
		return -2;
	if(quantum_prob_inline(reg->amplitude[num])==0)
		return -1;
	else
		return 1;
}

void
quantum_qft_gate1(int target, quantum_matrix m, quantum_qft_reg *reg)
{
	int j;
	COMPLEX_FLOAT t, tnot=0;
	//float limit;
	char *done;
	MAX_UNSIGNED i, ii, iset;

	if((m.cols != 2) || (m.rows != 2)) {
		quantum_error(QUANTUM_EMSIZE);
	}

	done = calloc(1 << reg->width, sizeof(char));
	//limit = (1.0 / ((MAX_UNSIGNED) 1 << reg->width)) * epsilon;

	for(i = 0; i < (1 << reg->width); i++) {
		if(quantum_prob_inline(reg->amplitude[i])==0) {
			continue;
		}
		if(!done[i]) {
			/* determine if the target of the basis state is set */
			iset = i & ((MAX_UNSIGNED)1 << target);
			tnot = 0;
			ii = i ^ ((MAX_UNSIGNED)1 << target);
			j = quantum_qft_get_state(ii, reg);
			//printf("i = %lli, ii = %lli, j = %i\n", i, ii, j);
			if(j >= 0) {
				tnot = reg->amplitude[ii];
			}

			t = reg->amplitude[i];
			if(iset) {
				reg->amplitude[i] = m.t[2] * tnot + m.t[3] * t;
			} else {
				reg->amplitude[i] = m.t[0] * t + m.t[1] * tnot;
				
			}
			if(j >= 0) {
				if(iset)
					reg->amplitude[ii] = m.t[0] * tnot + m.t[1] * t;
				else
					reg->amplitude[ii] = m.t[2] * t + m.t[3] * tnot;
			}
			else {  
				if((m.t[1] == 0) && (iset)) {
					break;
				}
				if((m.t[2] == 0) && !(iset)) {
					break; 
				}
				if(iset) {
					reg->amplitude[ii] = m.t[1] * t;
				}
				else {
					reg->amplitude[ii] = m.t[2] * t;
				}
			}
			done[ii] = 1;
		}
	}
	free(done);
	done = NULL;
}

void
quantum_qft_gate1_mpi(int target, quantum_matrix m, quantum_qft_reg *reg){
	int id, numprocs, l, j_state;
	COMPLEX_FLOAT t, tnot;
	char *done, doneTmp;
	MAX_UNSIGNED length, baseAddress, i, ii, iset, jj, procs, j;

	if((m.cols != 2) || (m.rows != 2)){
		quantum_error(QUANTUM_EMSIZE);
	}
	length = 1 << reg->width;
	done = calloc(length, sizeof(char));
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Barrier(MPI_COMM_WORLD);
	
	for(l = 0 ; l < numprocs ; l++) {
		for(i = 0; i < length; i++) {
			if(id == l){
				doneTmp = done[i];
			}
			MPI_Barrier(MPI_COMM_WORLD);		
			MPI_Bcast(&doneTmp, 1, MPI_CHAR, l, MPI_COMM_WORLD);
			//printf("%d, %d\n", flag, doneTmp);
			if (doneTmp == 0) 
			{
				tnot = 0;
				if(id == l) {
					t = reg->amplitude[i];
					baseAddress = id * length;
					ii = baseAddress + i;
					iset = ii & ((MAX_UNSIGNED)1 << target);
					jj = ii ^ ((MAX_UNSIGNED)1 << target);
				}
				MPI_Barrier(MPI_COMM_WORLD);
				MPI_Bcast(&t, 1, MPI_UNSIGNED_LONG, l, MPI_COMM_WORLD);
				MPI_Bcast(&iset, 1, MPI_UNSIGNED_LONG, l, MPI_COMM_WORLD);	
				MPI_Bcast(&jj, 1, MPI_UNSIGNED_LONG, l, MPI_COMM_WORLD);	
				
				procs = jj / length;
				if (id == procs) {
					j = jj % length;
					j_state = quantum_qft_get_state(j, reg);
					if(j_state >= 0) {
						tnot = reg->amplitude[j];
						if(iset) {
							reg->amplitude[j] = m.t[0] * tnot + m.t[1] * t;
						}
						else{
							reg->amplitude[j] = m.t[2] * t + m.t[3] * tnot;
						}
					}
					else {  
						if((m.t[1] == 0) && (iset)) {
							break;
						}
						if((m.t[2] == 0) && !(iset)) {
							break; 
						}
						if(iset) {
							reg->amplitude[j] = m.t[1] * t;
						}
						else{
							reg->amplitude[j] = m.t[2] * t;
						}
					}
					done[j] = 1;	
				}
				MPI_Barrier(MPI_COMM_WORLD);
				MPI_Bcast(&tnot, 1, MPI_UNSIGNED_LONG, procs, MPI_COMM_WORLD);

				if (id == l) {
					if(iset) {
						reg->amplitude[i] = m.t[2] * tnot + m.t[3] * t;
					}
					else {
						reg->amplitude[i] = m.t[0] * t + m.t[1] * tnot;
					}
					done[i] = 1;
				}
				MPI_Barrier(MPI_COMM_WORLD);
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		//MPI_Barrier(MPI_COMM_WORLD);
	}
	
	free(done);
	done = NULL;
}

void
quantum_qft_hadamard(int target, quantum_qft_reg *reg)
{
	quantum_matrix m;

	m = quantum_new_matrix(2, 2);
	m.t[0] = sqrt(1.0/2);
	m.t[1] = sqrt(1.0/2);
	m.t[2] = sqrt(1.0/2);
	m.t[3] = -sqrt(1.0/2);

	quantum_qft_gate1(target, m, reg);
	quantum_delete_matrix(&m);
}

void
quantum_qft_hadamard_mpi(int target, quantum_qft_reg *reg)
{
	quantum_matrix m;

	m = quantum_new_matrix(2, 2);

	m.t[0] = sqrt(1.0/2);
	m.t[1] = sqrt(1.0/2);
	m.t[2] = sqrt(1.0/2);
	m.t[3] = -sqrt(1.0/2);

	quantum_qft_gate1_mpi(target, m, reg);
	quantum_delete_matrix(&m);
}

void
quantum_qft_cond_phase(int control, int target, quantum_qft_reg *reg) {
	MAX_UNSIGNED i;
	COMPLEX_FLOAT z;
	float limit;

	limit = (1.0 / ((MAX_UNSIGNED) 1 << reg->width)) * epsilon;
	z = quantum_cexp(pi / ((MAX_UNSIGNED) 1 << (control - target)));

	for(i=0; i < (1 << reg->width); i++) {
		if(quantum_prob_inline(reg->amplitude[i]) < limit) {
			continue;
		}
		if( i & ((MAX_UNSIGNED) 1 << control)) {
			if( i & ((MAX_UNSIGNED) 1 << target)) {	
				reg->amplitude[i] *= z;
			}
		}
	}
}

void
quantum_qft_cond_phase_mpi(int control, int target, quantum_qft_reg *reg)
{
	MAX_UNSIGNED i, length;
	COMPLEX_FLOAT z;
	int id;
	float limit;
	 
	limit = (1.0 / ((MAX_UNSIGNED) 1 << reg->width)) * epsilon;
	z = quantum_cexp(pi / ((MAX_UNSIGNED) 1 << (control - target)));
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	length = 1 << reg->width;
   
	for(i=0; i < length; i++) {
		if(quantum_prob_inline(reg->amplitude[i]) < limit) {
			continue;
		}
		if((i + id * length) & ((MAX_UNSIGNED) 1 << control)) {
			if((i + id * length) & ((MAX_UNSIGNED) 1 << target)) {
				reg->amplitude[i] *= z;
			}
		}
	}
}