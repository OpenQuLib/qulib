#!/bin/sh
for i in `ls *.c`
do
	mpicc $i -c 
done
ar rcs libOpenQu_MPI.a *.o
rm *.o
