/* qureg.c: Quantum register management

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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "matrix.h"
#include "qureg.h"
#include "config.h"
#include "complex_q.h"
#include "error.h"

/* Create a new quantum register from scratch */
quantum_qft_reg
quantum_new_qft_reg(int width, int id, int r, float k)
{
	long long i;
	quantum_qft_reg  qr;
	qr.width = width;

	qr.amplitude = calloc(1 << width, sizeof(COMPLEX_FLOAT));

	i = (r - (id*(1<<width))%r)%r;
	for(; i<(1<<width); i+=r){
		qr.amplitude[i] = sqrt(1.0/k);
	}

	return qr;
}


void
quantum_qft_print_reg(quantum_qft_reg *reg,int id){
	int j;
	MAX_UNSIGNED i;
	for(i=0; i<(1<<reg->width); i++)
	{
		if(quantum_prob_inline(reg->amplitude[i])<0.0001) {
			continue;
		}
		printf("%d: % f %+fi|%lli> (%e) (|", id,quantum_real(reg->amplitude[i]), 
			quantum_imag(reg->amplitude[i]), i+id*(1<<reg->width), 
			quantum_prob_inline(reg->amplitude[i]));

		printf(">)\n");
	}
	printf("\n");
}

void
quantum_qft_delete_qureg(quantum_qft_reg *reg)
{
	if(reg->amplitude){
		free(reg->amplitude);
		reg->amplitude = NULL;
	}
}
