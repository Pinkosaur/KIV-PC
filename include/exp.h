#ifndef EXP_H
#define EXP_H

#include "mp_int.h"

/*
 * mp_pow
 *   Compute r = a ^ b using integer exponentiation semantics.
 *
 *   Behavior:
 *     - If a == 0 and b == 0, defines 0^0 = 1.
 *     - If a == 0 and b < 0, error (division by zero).
 *     - If b < 0: computes integer reciprocal semantics:
 *         r = 1 / (a^|b|) using integer division rules (result may be 0 or ±1).
 *     - Exponentiation performed by binary exponentiation (square-and-multiply).
 *
 *   All arguments and result are arbitrary-precision mp_int. Function
 *   returns SUCCESS on success, FAILURE on error.
 */
int mp_pow(mp_int *r, const mp_int *a, const mp_int *b);

#endif /* EXP_H */
