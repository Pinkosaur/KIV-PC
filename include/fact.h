#ifndef FACT_H
#define FACT_H

#include "mp_int.h"

/*
 * FactEntry
 *   Associates an unsigned integer n with a decimal string containing n!.
 *   The fact_table[] array is expected to be terminated with {0, NULL}.
 */
typedef struct {
    unsigned int n;
    const char *fact_str;
} FactEntry;

/* External table of some precomputed factorials (provided by fact.c) */
extern const FactEntry fact_table[];

/*
 * mp_fact
 *   Compute factorial r = a! where 'a' is a non-negative integer mp_int.
 *
 *   - Accepts arbitrarily large 'a' as mp_int (no conversion to native integers).
 *   - Uses precomputed factorials for a small set of values and then multiplies
 *     upward from the precomputed base using repeated multiplication.
 *   - Fast-path: when the running multiplier fits in a single limb, uses
 *     mp_mul_small for improved speed.
 *
 *   Returns SUCCESS on success, FAILURE on error (including negative input).
 */
int mp_fact(mp_int *r, const mp_int *a);

#endif /* FACT_H */
