#include <stdio.h>
#include "../quantum.h"
int main() {
        quantum_qft_reg qr;
        int target, width;
        int initval;
        initval = 2;
        width = 4;
        target = 1;
        qr = quantum_new_qft_reg(initval, width);
	printf("B:\n");
        quantum_qft_print_reg(&qr);
        quantum_qft_hadamard(target, &qr);
        printf("A:\n");
	quantum_qft_print_reg(&qr);
	quantum_qft_delete_qureg(&qr);
	return 0;
}   
