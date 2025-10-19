#ifndef MP_PRINT_H
#define MP_PRINT_H

#include "mp_int.h"

/* Multiplication algorithm selection based on mp_int size */

int mp_print_dec(mp_int *x);
int mp_print_bin(mp_int *x);
int mp_print_hex(mp_int *x);

#endif /* MP_PRINT_H */
