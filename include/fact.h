#ifndef FACT_H
#define FACT_H

#include "mp_int.h"

typedef struct {
    unsigned int n;
    const char *fact_str;
} FactEntry;

extern const FactEntry fact_table[];

int mp_fact(mp_int *r, const mp_int *a);

#endif /* FACT_H */
