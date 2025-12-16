#ifndef FACT_H
#define FACT_H

#include "mp_int.h"

/* Compute factorial: r = a! .
   a must be an integer (non-negative). If a is too large to run the
   prime-swing algorithm (it doesn't fit in unsigned long on the host),
   the implementation falls back to a slower mp_int-based method.
   Returns SUCCESS/FAILURE.
*/
int mp_fact(mp_int *r, const mp_int *a);

#endif /* FACT_H */
