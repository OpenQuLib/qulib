/* error.h: Declarations for error.c

   Copyright 2005 Bjoern Butscher, Hendrik Weimer

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

#ifndef __ERROR_H

#define __ERROR_H

enum {
  QUANTUM_SUCCESS      = 0,
  QUANTUM_FAILURE      = 1,
  QUANTUM_ENOMEM       = 2,
  QUANTUM_EMLARGE      = 3,
  QUANTUM_EMSIZE       = 4,
  QUANTUM_EHASHFULL    = 5,
  QUANTUM_EHERMITIAN   = 6, 
  QUANTUM_ENOCONVERGE  = 7,
  QUANTUM_ENOSOLVER    = 8,
  QUANTUM_ENOLAPACK    = 32768, /* LAPACK errors start at 32768 */
  QUANTUM_ELAPACKARG   = 32769,
  QUANTUM_ELAPACKCONV  = 32770,
  QUANTUM_EMCMATRIX    = 65536, /* internal errors start at 65536 */
  QUANTUM_EOPCODE      = 65537
};

extern void *quantum_error_handler(void *f(int));
extern const char *quantum_strerr(int errno);
extern void quantum_error(int errno);

#endif
