#!/bin/sh
rm *.a
rm *.o
for i in `ls *.c`
do
	if [ $i != "slave_qft.c" ]
	then
		echo "mpicc $i"	
		mpicc  $i -c 
	fi
done
sw5cc -slave -msimd -c slave_qft.c
echo "after sw5cc"
ar rcs libOpenQu.a *.o
#rm *.o
mpicc  -c newReg_test.c
#mpicc  newReg_test.o libOpenQu.a -lm -o a.out
mpicc  newReg_test.o libOpenQu.a slave_qft.o -lm -o a.out




