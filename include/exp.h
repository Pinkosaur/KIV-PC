#ifndef EXP_H
#define EXP_H

#include "mp_int.h"

/*
 * mp_pow
 * ------
 * Computes r = a ^ b (integer exponentiation).
 *
 * Semantics:
 * - 0^0 = 1
 * - 0^negative = Error (Division by zero)
 * - negative exponent: Returns integer reciprocal (0 unless result is +/- 1).
 *
 * Returns SUCCESS or FAILURE.
 */
int mp_pow(mp_int *r, const mp_int *a, const mp_int *b);

#endif /* EXP_H */