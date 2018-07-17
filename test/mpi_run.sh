#!/bin/sh
cd ..
./mk_mpi.sh > /dev/null 
cd - > /dev/null
if [ $# -gt 0 ] ; then
	mpicc -c $1.c -o $1.o
	mpicc $1.c ../libOpenQu_MPI.a -lm -o a.out
else
	mpicc -c mpi_test.c 
	mpicc mpi_test.o ../libOpenQu_MPI.a -lm -o a.out
fi
