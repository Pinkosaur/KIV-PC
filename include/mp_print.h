#ifndef MP_PRINT_H
#define MP_PRINT_H

#include "mp_int.h"

/*
 * mp_print_dec
 * ------------
 * Print mp_int in decimal (base 10) to stdout.
 */
int mp_print_dec(mp_int *x);

/*
 * mp_print_bin
 * ------------
 * Print mp_int in binary (base 2) with explicit "0b" prefix.
 * Output uses minimal two's complement width logic.
 */
int mp_print_bin(mp_int *x);

/*
 * mp_print_hex
 * ------------
 * Print mp_int in hexadecimal (base 16) with explicit "0x" prefix.
 * Output uses minimal two's complement nibble width logic.
 */
int mp_print_hex(mp_int *x);

#endif /* MP_PRINT_H */