#ifndef MP_PRINT_H
#define MP_PRINT_H

#include "mp_int.h"

/*
 * mp_print_dec
 *   Print mp_int in decimal (base 10) to stdout.
 *   Uses repeated division by 1e9 to produce chunked output
 *   Returns SUCCESS on success, FAILURE on error.
 */
int mp_print_dec(mp_int *x);

/*
 * mp_print_bin
 *   Print mp_int in binary with an explicit sign bit.
 *   - zero       => "0b0"
 *   - positive   => "0b" followed by (bitlen(|x|)+1) bits where top bit = 0 (explicit sign bit)
 *   - negative   => minimal two's-complement representation with top bit = 1
 *
 *   Returns SUCCESS or FAILURE.
 */
int mp_print_bin(mp_int *x);

/*
 * mp_print_hex
 *   Print mp_int in hexadecimal with nibble-wise explicit sign semantics:
 *   - zero       => "0x0"
 *   - positive   => minimal hex digits, but if top hex digit >= 8, prepend a '0' nibble
 *                  so the sign nibble is 0.
 *   - negative   => minimal two's-complement nibble representation (at most one leading 'f' as sign)
 *
 *   Returns SUCCESS or FAILURE.
 */
int mp_print_hex(mp_int *x);

#endif /* MP_PRINT_H */
