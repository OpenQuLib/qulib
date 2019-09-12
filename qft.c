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
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <mpi.h>
#include <athread.h>

#include "matrix.h"
#include "defs.h"
#include "complex_q.h"
#include "error.h"
#include "gates.h"
#include "qureg.h"

#define IMAGINARY _Complex_I

/* Perform a QFT on a quantum register. This is done by application of
   conditional phase shifts and hadamard gates. At the end, the
   position of the bits is reversed. */
float global_limit;
COMPLEX_FLOAT global_z;

int quantum_status = 0;
float quantum_lambda = 0;

void quantum_qft_gate1_wq(int target, quantum_matrix m, quantum_qft_reg *reg) 
{
	int j_state;
	COMPLEX_FLOAT t, tnot=0;
	char *done;
	MAX_UNSIGNED i, j, iset,k,step;
	MAX_UNSIGNED ii,jj;
	float limit;
	int myid,targetid;
	int numprocs;
	int count;
	int index;
	MPI_Status stat;
	count = ((MAX_UNSIGNED) 1 << reg->width) > 128 ? 128: ((MAX_UNSIGNED) 1 << reg->width);
	COMPLEX_FLOAT *recvamplitude;
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);    
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
	limit = (1.0 / ((MAX_UNSIGNED) 1 << reg->width)) * epsilon;
	targetid = ((MAX_UNSIGNED)1 << (target - reg->width));
		
	if(targetid != 0)
	{
		recvamplitude = calloc(count, sizeof(COMPLEX_FLOAT));
		if((myid / targetid) % 2 == 0)
			targetid = myid + targetid;
		else
			targetid = myid - targetid;
		for(index = 0 ; index < (1 << reg->width) ; index = index + count)
		{
			//if(myid == 0)printf("index:%d\ttargetid:%d\n",index,targetid);
			
			//	if(myid == 0) printf("%d:send message\n",target);
			MPI_Send(&reg->amplitude[index], count*2, MPI_FLOAT, targetid, myid , MPI_COMM_WORLD);
			//	if(myid == 0) printf("%d:recv message\n",target);
			MPI_Recv(&recvamplitude[0], count*2, MPI_FLOAT, targetid, targetid, MPI_COMM_WORLD, &stat);

			for( i = index ; i < index + count ; i ++)
			{
				if(quantum_prob_inline(reg->amplitude[i])<=limit && quantum_prob_inline(recvamplitude[i-index])<=limit)
					continue;
				else if(quantum_prob_inline(reg->amplitude[i])>limit && quantum_prob_inline(recvamplitude[i-index])<=limit)
				{

					iset = (MAX_UNSIGNED) (i +(MAX_UNSIGNED) myid * (1 << reg->width)) & ((MAX_UNSIGNED)1 << target);
					tnot = 0;
					t = reg->amplitude[i];
					if(iset)
						reg->amplitude[i] = m.t[2] * tnot + m.t[3] * t;
					else
						reg->amplitude[i] = m.t[0] * t + m.t[1] * tnot;
				}
				else if(quantum_prob_inline(reg->amplitude[i])<=limit && quantum_prob_inline(recvamplitude[i-index])>limit)
				{
					iset =  ((MAX_UNSIGNED)i + (MAX_UNSIGNED)targetid * (1 << reg->width)) & ((MAX_UNSIGNED)1 << target);
					t = recvamplitude[i-index];
					if(iset)
						reg->amplitude[i] = m.t[1] * t;
					else
						reg->amplitude[i] = m.t[2] * t;
				}
				else
				{
					if(myid < targetid)
					{
						iset =  ((MAX_UNSIGNED)i + (MAX_UNSIGNED)myid * (1 << reg->width)) & ((MAX_UNSIGNED)1 << target);
						tnot = recvamplitude[i-index];
						t = reg->amplitude[i];
						if(iset)
							reg->amplitude[i] = m.t[2] * tnot + m.t[3] * t;
						else
							reg->amplitude[i] = m.t[0] * t + m.t[1] * tnot;
					}
					else
					{
						iset =  (MAX_UNSIGNED)(i + (MAX_UNSIGNED)targetid * (1 << reg->width)) & ((MAX_UNSIGNED)1 << target);
						tnot = reg->amplitude[i];
						t = recvamplitude[i-index];
						if(iset) 
							reg->amplitude[i] = m.t[0] * tnot + m.t[1] * t;
						else 
							reg->amplitude[i] = m.t[2] * t + m.t[3] * tnot;

					}
				}
			}
		}
		free(recvamplitude);
		recvamplitude = NULL;
	}
	else
	{

		step = (MAX_UNSIGNED)1 << (target+1);
		for(k = 0; k < (1 << reg->width); k = k + step)
		{
			for(i = k ;  i < k + step / 2 ; i ++)
			{
				j = i + step/2;
				if(quantum_prob_inline(reg->amplitude[i])<=limit && quantum_prob_inline(reg->amplitude[j])<=limit)
					continue;
				ii = i;
				jj = j;
				if(quantum_prob_inline(reg->amplitude[i]) <= limit)
				{
					ii = j;
					jj = i;
				}
				iset = (MAX_UNSIGNED)(ii + (MAX_UNSIGNED)myid * (1 << reg->width)) & ((MAX_UNSIGNED)1 << target);
				tnot = 0;
				if(quantum_prob_inline(reg->amplitude[jj]) > limit)
				{
					tnot = reg->amplitude[jj];
				}
				t = reg->amplitude[ii];


				if(iset)
					reg->amplitude[ii] = m.t[2] * tnot + m.t[3] * t;
				else
					reg->amplitude[ii] = m.t[0] * t + m.t[1] * tnot;

				if(quantum_prob_inline(reg->amplitude[jj]) > limit)
				{   
					if(iset)
						reg->amplitude[jj] = m.t[0] * tnot + m.t[1] * t;
					else 
						reg->amplitude[jj] = m.t[2] * t + m.t[3] * tnot;
				}   
				else
				{   
					if(iset)
						reg->amplitude[jj] = m.t[1] * t;
					else
						reg->amplitude[jj] = m.t[2] * t;
				}
			}
		}
	}

}

extern SLAVE_FUN(func)();
quantum_qft_reg *global_reg;
extern quantum_qft_reg qr;

void quantum_qft_qft(int width, quantum_qft_reg *reg)
{
	int i, j;
	COMPLEX_FLOAT z;

	MAX_UNSIGNED k;
	quantum_matrix m;
	int myid,numprocs;
	global_reg = reg;
	global_reg->amplitude = reg->amplitude;
	global_reg->width = reg->width;
	global_reg->g_width = width;	
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);    
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);  
	m = quantum_new_matrix(2, 2);
	m.t[0] = sqrt(1.0/2);
	m.t[1] = sqrt(1.0/2);
	m.t[2] = sqrt(1.0/2);
	m.t[3] = -sqrt(1.0/2);

	athread_init();

	global_limit = (1.0 / ((MAX_UNSIGNED) 1 << width)) * epsilon;
	for(i = width - 1; i >= 0; i--) 
	{
		for(j = width - 1; j > i; j--) 
		{
/*
			global_reg->master_i = i;
			global_reg->master_j = j;
			global_reg->master_id = myid;
			*/
			qr.master_i = i;
			qr.master_j = j;
			qr.master_id = myid;
			
			global_z = quantum_cexp(pi / ((MAX_UNSIGNED) 1 << (j - i)));
			athread_spawn(func,0);
			athread_join();
			fflush(NULL);
/*
			z = quantum_cexp(pi / ((MAX_UNSIGNED) 1 << (j - i)));

			for(k = 0; k < (1 << reg->width); k ++) 
			{
				if( (MAX_UNSIGNED)(k + (MAX_UNSIGNED)myid * ((MAX_UNSIGNED)1 << reg->width))& ((MAX_UNSIGNED) 1 << j)) 
				{
					if( (MAX_UNSIGNED)(k + (MAX_UNSIGNED)myid * ((MAX_UNSIGNED)1 << reg->width)) & ((MAX_UNSIGNED) 1 << i)) 
					{
						reg->amplitude[k] *= z;
					}
				}
			}
*/			
		}
		quantum_qft_gate1_wq(i, m, reg);
	}
	quantum_delete_matrix(&m);
}

void quantum_mpi_cnot(int control, int target, quantum_qft_reg *reg)
{
	MAX_UNSIGNED i,j,k,step,index;
	int myid,targetid;
	int numprocs;
	int count;

	MPI_Status stat;
	COMPLEX_FLOAT temp;
	COMPLEX_FLOAT *recvamplitude;
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);    
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
	count = ((MAX_UNSIGNED) 1 << reg->width) > 256 ? 256: ((MAX_UNSIGNED) 1 << reg->width);
	targetid = ((MAX_UNSIGNED)1 << (target - reg->width));
	if(targetid != 0)
	{
		recvamplitude = calloc(count, sizeof(COMPLEX_FLOAT));
		if((myid / targetid) % 2 == 0)
			targetid = myid + targetid;
		else
			targetid = myid - targetid;

		for(index = 0 ; index < (1 << reg->width) ; index = index + count)
		{
			MPI_Send(&reg->amplitude[index], count*2, MPI_FLOAT, targetid, myid , MPI_COMM_WORLD);
			MPI_Recv(&recvamplitude[0], count*2, MPI_FLOAT, targetid, targetid, MPI_COMM_WORLD, &stat);
			for( i = 0; i < count ; i ++)
			{
			    if((MAX_UNSIGNED)(i + index + (MAX_UNSIGNED)myid * ((MAX_UNSIGNED)1 << reg->width))& ((MAX_UNSIGNED) 1 << control))
                            {   
                                reg->amplitude[i+index] = recvamplitude[i];
                            } 	
			}
		}
		free(recvamplitude);
		recvamplitude = NULL;
	}
	else
	{
		step = (MAX_UNSIGNED)1 << (target+1);
		for(k = 0; k < (1 << reg->width); k = k + step)
		{
			for(i = k ;  i < k + step / 2 ; i ++)
			{
				if((MAX_UNSIGNED)(i + (MAX_UNSIGNED)myid * ((MAX_UNSIGNED)1 << reg->width))& ((MAX_UNSIGNED) 1 << control))
				{
					j = i + step/2;
					temp = reg->amplitude[i];
					reg->amplitude[i] = reg->amplitude[j];
					reg->amplitude[j] = temp;
				}
			}
		}
	}
}

void quantum_mpi_sigma_x(int target, quantum_qft_reg *reg)
{
	MAX_UNSIGNED i,j,k,step,index;
	int myid,targetid;
	int numprocs;
	int count;

	MPI_Status stat;
	COMPLEX_FLOAT temp;
	COMPLEX_FLOAT *recvamplitude;
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
	count = ((MAX_UNSIGNED) 1 << reg->width) > 256 ? 256: ((MAX_UNSIGNED) 1 << reg->width);
	targetid = ((MAX_UNSIGNED)1 << (target - reg->width));
	if(targetid != 0)
	{
		recvamplitude = calloc(count, sizeof(COMPLEX_FLOAT));
		if((myid / targetid) % 2 == 0)
		        targetid = myid + targetid;
		else
		        targetid = myid - targetid;

		for(index = 0 ; index < (1 << reg->width) ; index = index + count)
		{
			MPI_Send(&reg->amplitude[index], count*2, MPI_FLOAT, targetid, myid , MPI_COMM_WORLD);
			MPI_Recv(&recvamplitude[0], count*2, MPI_FLOAT, targetid, targetid, MPI_COMM_WORLD, &stat);
			for( i = 0; i < count ; i ++)
			{
				reg->amplitude[i+index] = recvamplitude[i];
			}
		}
		free(recvamplitude);
		recvamplitude = NULL;
	}
	else
	{
		step = (MAX_UNSIGNED)1 << (target+1);
		for(k = 0; k < (1 << reg->width); k = k + step)
		{
			for(i = k ;  i < k + step / 2 ; i ++)
			{

				j = i + step/2;
				temp = reg->amplitude[i];
				reg->amplitude[i] = reg->amplitude[j];
				reg->amplitude[j] = temp; 
			}
		}
	}
}

void quantum_mpi_toffoli(int control1, int control2, int target, quantum_qft_reg *reg)
{
	MAX_UNSIGNED i,j,k,step,index;
	int myid,targetid;
	int numprocs;
	int count;

	MPI_Status stat;
	COMPLEX_FLOAT temp;
	COMPLEX_FLOAT *recvamplitude;
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
	count = ((MAX_UNSIGNED) 1 << reg->width) > 256 ? 256: ((MAX_UNSIGNED) 1 << reg->width);
	targetid = ((MAX_UNSIGNED)1 << (target - reg->width));
	if(targetid != 0)
	{
		recvamplitude = calloc(count, sizeof(COMPLEX_FLOAT));
		if((myid / targetid) % 2 == 0)
			targetid = myid + targetid;
		else
			targetid = myid - targetid;

		for(index = 0 ; index < (1 << reg->width) ; index = index + count)
		{
			MPI_Send(&reg->amplitude[index], count*2, MPI_FLOAT, targetid, myid , MPI_COMM_WORLD);
			MPI_Recv(&recvamplitude[0], count*2, MPI_FLOAT, targetid, targetid, MPI_COMM_WORLD, &stat);
			for( i = 0; i < count ; i ++)
			{
				if((MAX_UNSIGNED)(i + index + (MAX_UNSIGNED)myid * ((MAX_UNSIGNED)1 << reg->width))& ((MAX_UNSIGNED) 1 << control1))
				{
					if((MAX_UNSIGNED)(i + index + (MAX_UNSIGNED)myid * ((MAX_UNSIGNED)1 << reg->width))& ((MAX_UNSIGNED) 1 << control2))
					{
						reg->amplitude[i+index] = recvamplitude[i];
					}
				}
			}
		}
		free(recvamplitude);
		recvamplitude = NULL;
	}
	else
	{
		step = (MAX_UNSIGNED)1 << (target+1);
		for(k = 0; k < (1 << reg->width); k = k + step)
		{
			for(i = k ;  i < k + step / 2 ; i ++)
			{
				if((MAX_UNSIGNED)(i + index + (MAX_UNSIGNED)myid * ((MAX_UNSIGNED)1 << reg->width))& ((MAX_UNSIGNED) 1 << control1))
				{
					if((MAX_UNSIGNED)(i + index + (MAX_UNSIGNED)myid * ((MAX_UNSIGNED)1 << reg->width))& ((MAX_UNSIGNED) 1 << control2))
					{ 
						j = i + step/2;
						temp = reg->amplitude[i];
						reg->amplitude[i] = reg->amplitude[j];
						reg->amplitude[j] = temp;
					}
				}
			}
		}
	}
}

void quantum_mpi_sigma_y(int target, quantum_qft_reg *reg)
{
        MAX_UNSIGNED i,j,k,step,index;
        int myid,targetid;
        int numprocs;
        int count;

        MPI_Status stat;
        COMPLEX_FLOAT temp;
        COMPLEX_FLOAT *recvamplitude;
        MPI_Comm_rank(MPI_COMM_WORLD,&myid);
        MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
        count = ((MAX_UNSIGNED) 1 << reg->width) > 256 ? 256: ((MAX_UNSIGNED) 1 << reg->width);
        targetid = ((MAX_UNSIGNED)1 << (target - reg->width));
        if(targetid != 0)
        {   
                recvamplitude = calloc(count, sizeof(COMPLEX_FLOAT));
                if((myid / targetid) % 2 == 0)
                        targetid = myid + targetid;
                else
                        targetid = myid - targetid;

                for(index = 0 ; index < (1 << reg->width) ; index = index + count)
                {
                        MPI_Send(&reg->amplitude[index], count*2, MPI_FLOAT, targetid, myid , MPI_COMM_WORLD);
                        MPI_Recv(&recvamplitude[0], count*2, MPI_FLOAT, targetid, targetid, MPI_COMM_WORLD, &stat);
                        for( i = 0; i < count ; i ++)
                        {
                            reg->amplitude[i+index] = recvamplitude[i];
                        }
                }
                free(recvamplitude);
                recvamplitude = NULL;
        }
        else
        {
                step = (MAX_UNSIGNED)1 << (target+1);
                for(k = 0; k < (1 << reg->width); k = k + step)
                {
                        for(i = k ;  i < k + step / 2 ; i ++)
                        {

                             j = i + step/2;
                             temp = reg->amplitude[i];
                             reg->amplitude[i] = reg->amplitude[j];
                             reg->amplitude[j] = temp;
                        }
                }
        }
		for(i=0; i<(1 << reg->width); i++)
		{
			if((MAX_UNSIGNED)(i + (MAX_UNSIGNED)myid * ((MAX_UNSIGNED)1 << reg->width))& ((MAX_UNSIGNED) 1 << target))
					reg->amplitude[i] = reg->amplitude[i] * IMAGINARY;
			else
					reg->amplitude[i] = reg->amplitude[i] * (-IMAGINARY);
		}	
}

