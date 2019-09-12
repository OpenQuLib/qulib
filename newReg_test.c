#include <complex.h>
#include <stdio.h>
#include <mpi.h>
#include <sys/time.h>

#include "complex_q.h"
#include "qureg.h"
#include "qft.h"

int ln_2(int a)
{
	int result = 0;
	while(a != 0)
	{
		result ++;
		a = a >> 1;
	}
	result --;
	return result;
};

quantum_qft_reg  qr;

int main(int argc,char *argv[]){

	MPI_Init(&argc,&argv);
	int width;
	int myid,numprocs;
	int i,j,step;
	struct  timeval start_t;
	struct  timeval stop_t;
	float diff;
	width = atoi(argv[1]);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);		    
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
		
	qr = quantum_new_qft_reg((width - ln_2(numprocs)), myid, 368, 2850.0);
	qr.g_width = width - ln_2(numprocs);
//	quantum_qft_print_reg(&qr, myid);

	MPI_Barrier(MPI_COMM_WORLD);
	step = 1;
	MPI_Barrier(MPI_COMM_WORLD);
	gettimeofday(&start_t,NULL);
	quantum_qft_qft(width, &qr);
	gettimeofday(&stop_t,NULL);
	MPI_Barrier(MPI_COMM_WORLD);	
	for(i=0; i<width/2; i++)
	{
		quantum_mpi_cnot(i, width-i-1, &qr);
		quantum_mpi_cnot(width-i-1, i, &qr);
		quantum_mpi_cnot(i, width-i-1, &qr);
	}
	quantum_qft_print_reg(&qr, myid);
//	quantum_mpi_cnot(width-2,1, &qr);
//	quantum_mpi_sigma_x(0,&qr);
//	quantum_mpi_sigma_y(0,&qr);
//	quantum_mpi_toffoli(1, 0, width-1,&qr);

	MPI_Barrier(MPI_COMM_WORLD);
	if(myid == 0)
		printf("start_t.tv_sec: %d, \nstart_t.tv_usec: %d, \nstop_t.tv_sec: %d, \nstop_t.tv_usec: %d.\n",start_t.tv_sec,start_t.tv_usec,stop_t.tv_sec,stop_t.tv_usec);
	diff = (1000000 * (stop_t.tv_sec-start_t.tv_sec)+ stop_t.tv_usec-start_t.tv_usec)/1000000.0;
	if(myid ==0 ) 
		printf("%d Use Time:%fs\n", width,diff);
	quantum_qft_delete_qureg(&qr);
	MPI_Finalize();
	
	return 0;
}  
