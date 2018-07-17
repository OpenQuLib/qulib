#include <stdio.h>
#include "../quantum.h"
int main(){
	quantum_qft_reg qr;
	int width, initval;
  int x = 0;
  width = 3;
  initval = 5;
 qr = quantum_new_qft_reg(initval, width);
  printf("Before QFT:\n");
  quantum_qft_print_reg(&qr);
quantum_qft_hadamard(1, &qr);
  printf("After H:\n");
  quantum_qft_print_reg(&qr);
//  printf("the size is %d\n", qr.size);
  quantum_qft_delete_qureg(&qr);

return 0;
}

