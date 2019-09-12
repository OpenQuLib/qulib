/* qft.h: Declarations for qft.c
   
   Copyright 2003 Bjoern Butscher, Hendrik Weimer

   This file is part of libquantum

   libquantum is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published
   by the Free Software Foundation; either version 3 of the License,
   or (at your option) any later version.

   libquantum is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with libquantum; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
   MA 02110-1301, USA

*/

#ifndef __QFT_H
#define __QFT_H

void quantum_qft_qft(int width, quantum_qft_reg *reg);
void quantum_mpi_cnot(int control, int target, quantum_qft_reg *reg);
void quantum_mpi_sigma_x(int target, quantum_qft_reg *reg);
void quantum_mpi_toffoli(int control1, int control2, int target, quantum_qft_reg *reg);
void quantum_mpi_sigma_y(int target, quantum_qft_reg *reg);

#endif
