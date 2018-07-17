#include <mpi.h>
#include <stdio.h>
#include "../quantum.h"

int bits(MAX_UNSIGNED a)
{
	int i=0;
	while(a != 0)
	{
		a = a >>1;
		i ++;
	}
	return i-1;
}

int main(int argc, char *argv[]) {
	int id, numprocs, width, myWidth, i, j;
	quantum_qft_reg qr;

	width = 4;
	
	MPI_Init(&argc, &argv);	
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	
	myWidth = width - bits(numprocs);
	qr = quantum_new_qft_reg_mpi(7, myWidth);
	MPI_Barrier(MPI_COMM_WORLD);
	
	if (id == 0) {
		printf("Before:\n");
	}
	MPI_Barrier(MPI_COMM_WORLD);

	for(i = 0; i < numprocs; i++) { 
		if(id == i) {
			quantum_qft_print_reg(&qr);
		}
		MPI_Barrier(MPI_COMM_WORLD);
  	}

  	// do some operation to qr
  	quantum_qft_qft_mpi(width, &qr);
  	//quantum_qft_cond_phase_mpi(1, 0, &qr);
	MPI_Barrier(MPI_COMM_WORLD);

  	if (id == 0) {
		printf("After:\n");
	}
	MPI_Barrier(MPI_COMM_WORLD);

	for(i = 0; i < numprocs; i++) { 
		if(id == i)	{
			quantum_qft_print_reg(&qr);
		}
		MPI_Barrier(MPI_COMM_WORLD);
  	}

	MPI_Finalize();
	return 0;
}
