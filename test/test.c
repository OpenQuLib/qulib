#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "../quantum.h"

int main(int argc, char **argv) {
	quantum_qft_reg qr;
	int width, initval;
	int x = 0;
	width = 3;
	initval = 5;
	while(1){
		printf("please input the width and init value of QuReg\n");
		printf("width: ");
		scanf("%d", &width);
		printf("the init value: ");
		scanf("%d", &initval);
		if(pow(2, width) <= initval) {
			printf("initval %d is bigger than 2^%d!\nplease input again!\n", initval, width);
			continue;
		}

		qr=quantum_new_qft_reg(initval, width);
		printf("Before H:\n");
		quantum_qft_print_reg(&qr);
		quantum_qft_qft(width, &qr);
		printf("After H:\n"); 
		quantum_qft_print_reg(&qr);
		quantum_qft_delete_qureg(&qr);
	}
	return 0;
}
