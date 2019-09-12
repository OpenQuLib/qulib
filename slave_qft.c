/* name: slave_qft.c
 * function: do matrix multiplication
 *           with 160kB local memory, only do 2^14 qubits in slave core,2^14 * 8 = 2^17 = 128kB
 */
#include <slave.h>

//#include <crts.h>
#include "complex_q.h"
#include "config.h"
#include "qureg.h"

#define SLAVE_WIDTH (12)
#define SLAVE_NUM (6)
#define MAX_MASTER_WIDTH (SLAVE_WIDTH + SLAVE_NUM)
#define MASTER_LENGTH (1<<MAX_MASTER_WIDTH)
#define SLAVE_NUMS (1<<SLAVE_NUM)
#define SLAVE_LENGTH (1<<SLAVE_WIDTH)
#define HEAD_LEN (20)

__thread_local volatile unsigned long get_reply, put_reply;
__thread_local int slave_id;
__thread_local COMPLEX_FLOAT z;
__thread_local float limit;

typedef struct{
	int width;
	int g_width;
	int master_i;
	int master_j;
	int master_id;
	COMPLEX_FLOAT amplitude[SLAVE_LENGTH];
}Reg_Struct;

__thread_local Reg_Struct slave_reg;

extern COMPLEX_FLOAT global_z;
extern float global_limit;
extern quantum_qft_reg qr;

void func()
{

	int width_outside;
	unsigned long long i,k,slave_index;
	
	//printf("slave_func starting... \n");

	slave_id = athread_get_id(-1);

	z = global_z;
	limit = global_limit;
	/*
	slave_reg.g_width = global_reg->g_width;
	slave_reg.width = global_reg->width;
	slave_reg.master_i = global_reg->master_i;
	slave_reg.master_j = global_reg->master_j;
	slave_reg.master_id = global_reg->master_id;
	*/
	slave_reg.g_width = qr.g_width;
	slave_reg.width = qr.width;
	slave_reg.master_i = qr.master_i;
	slave_reg.master_j = qr.master_j;
	slave_reg.master_id = qr.master_id;
	
//	if(global_reg->master_id == 0 && slave_id == 0)
//		printf("slave_id:%d;  limit:%d;  slave_reg.width:%d;  slave_reg.master_i:%d;  slave_reg.master_j:%d;  slave_reg.master_id:%d\n", \
		slave_id, limit, slave_reg.width, slave_reg.master_i, slave_reg.master_j, slave_reg.master_id);


	width_outside = slave_reg.width-MAX_MASTER_WIDTH;    //if width is 20, only execute once in 64 slave cores
	if(width_outside<0)
		return;

	
//	printf("3slave_func starting... slave_id:%d\n", slave_id);
	slave_index = (unsigned long long)(slave_id*SLAVE_LENGTH + slave_reg.master_id*(1<<slave_reg.width));
//	if(global_reg->master_id == 0 && slave_id == 0) printf("%d  %d\n", slave_index, width_outside);
	for(i=0; i<(1<<width_outside); i++)
	{

//		if(global_reg->master_id == 0 && slave_id == 0) 
//			printf("global_reg:%x\n", (global_reg->amplitude + i*MASTER_LENGTH*8 + slave_id*SLAVE_LENGTH*8));
//		printf("4slave_func starting... width_outside:%d;  i:%d\n",width_outside, i);
/*		if(qr.master_id == 0 && slave_id == 3) {
			printf("0,global_reg:%e\n", quantum_prob_inline((qr.amplitude + i*MASTER_LENGTH + slave_id*SLAVE_LENGTH)[0]));
			printf("0,slave_reg:%e\n", quantum_prob_inline(slave_reg.amplitude[0]));
		}
		*/
		get_reply = 0;
		athread_get(PE_MODE, qr.amplitude + i*MASTER_LENGTH + slave_id*SLAVE_LENGTH, \
			slave_reg.amplitude, SLAVE_LENGTH*8, &get_reply, 0, 0, 0);

		while(get_reply!=1);
/*		if(qr.master_id == 0 && slave_id == 3) {
			printf("0-1,global_reg:%e\n", quantum_prob_inline((qr.amplitude + i*MASTER_LENGTH + slave_id*SLAVE_LENGTH)[0]));
			printf("0-1,slave_reg:%e\n", quantum_prob_inline(slave_reg.amplitude[0]));
			printf("%lld,i:%lld,slave_index:%lld\n", i*MASTER_LENGTH + slave_index,i,slave_index);
		}
		*/
		//printf("5slave_func starting... \n");
		for(k = 0; k < SLAVE_LENGTH; k ++)
		{
			if(quantum_prob_inline(slave_reg.amplitude[k]) <= limit)
			{
				continue;
			}
			if( (k + i*MASTER_LENGTH + slave_index) & ((MAX_UNSIGNED)1 << slave_reg.master_j))
			{
				if( (k + i*MASTER_LENGTH + slave_index) & ((MAX_UNSIGNED)1 << slave_reg.master_i))
				{
					slave_reg.amplitude[k] *= z;
				}
			}
		}
/*		if(qr.master_id == 0 && slave_id == 3) {
			printf("1,global_reg:%e\n", quantum_prob_inline((qr.amplitude + i*MASTER_LENGTH + slave_id*SLAVE_LENGTH)[0]));
			printf("1,slave_reg:%e\n", quantum_prob_inline(slave_reg.amplitude[0]));
		}
*/
		put_reply = 0;
//		if(global_reg->master_id == 0 && slave_id == 0 && i == 31) 
//			printf("slave_reg:%f  %f\n", slave_reg.amplitude[0]);
		athread_put(PE_MODE,slave_reg.amplitude,qr.amplitude+ i*MASTER_LENGTH + slave_id*SLAVE_LENGTH, \
			SLAVE_LENGTH*8, &put_reply, 0, 0);

		while(put_reply!=1);
/*		if(qr.master_id == 0 && slave_id == 3) {
			printf("1-1,global_reg:%e\n", quantum_prob_inline((qr.amplitude + i*MASTER_LENGTH + slave_id*SLAVE_LENGTH)[0]));
			printf("1-1,slave_reg:%e\n", quantum_prob_inline(slave_reg.amplitude[0]));
		}
*/
	}

//	printf("6slave_func starting... \n");

}


