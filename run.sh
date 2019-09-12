bsub -I -b -q q_sw_expr -n 4 -cgsp 64 -share_size 6000 -host_stack 512 ./a.out 20 > semi_qft_20_8.log
bsub -I -b -q q_sw_expr -n 1 -cgsp 64 -share_size 6000 -host_stack 512 ./a.out 20 > semi_qft_20_1.log
