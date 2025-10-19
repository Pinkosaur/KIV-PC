#ifndef PRECOMPUTED_FACT_H
#define PRECOMPUTED_FACT_H

typedef struct {
    unsigned int n;
    const char *fact_str;
} FactEntry;

extern const FactEntry fact_table[];

#endif /* PRECOMPUTED_FACT_H */
