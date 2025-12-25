#ifndef FACT_H
#define FACT_H

#include "mp_int.h"

/* Structure for precomputed factorial table entries */
typedef struct {
    unsigned int n;
    const char *fact_str;
} fact_entry;

/* External table provided by fact.c */
extern const fact_entry fact_table[];

/*
 * mp_fact
 * -------
 * Computes r = a! (factorial).
 *
 * Optimization:
 * Uses the precomputed table to jump to the nearest factorial <= a,
 * then performs incremental multiplication up to a.
 *
 * Returns SUCCESS or FAILURE (e.g., if input is negative).
 */
int mp_fact(mp_int *r, const mp_int *a);

#endif /* FACT_H */