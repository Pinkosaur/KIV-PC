/*
 * fact.c -- Factorial implementation.
 *
 * mp_fact computes a! for non-negative integer a (mp_int). The implementation:
 * - Uses a precomputed table for small factorials to jump-start calculation.
 * - Multiplies upward from the largest precomputed base using repeated multiplication.
 * - Uses fast-path mp_mul_small when the multiplier fits in a single limb.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mp_int.h"
#include "fact.h"
#include "error.h"

/* Precomputed factorial table */
const fact_entry fact_table[] = {
    {0, "1"}, {1, "1"}, {2, "2"}, {5, "120"}, {10, "3628800"},
    {50, "30414093201713378043612608166064768844377641568960512000000000000"},
    {75, "24809140811395398091946477116594033660926243886570122837795894512655842677572867409443815424000000000000000000"},
    {100, "93326215443944152681699238856266700490715968264381621468592963895217599993229915608941463976156518286253697920827223758251185210916864000000000000000000000000"},
    {200, "78865786736479050355236321393218506229513597768717326329474253324435944996340334292030428401198462390417721213891963883025764279024263710506192662495282993111346285727076331723739698894392244562145166424025403329186413122742829485327752424240757390324032125740557956866022603190417032406235170085879617892222278962370389737472000000000000000000000000000000000000000000000000"},
    {256, "8578177753428426541190822716812326251577815202794856198596556503772694525531475893774402913604514084503758853423365843061571968346936964753222892884974260256796373325633687864426752076267945601879688679715211433077020775266464514647091873261008328763257028189807736717814541702505230186084953190681382574810702528175594594769870346657127381392862052347568082188607012036110831520935019474371091017269682628616062636624350228409441914084246159360000000000000000000000000000000000000000"},
    {0, NULL}
};

/* Helper: create mp_int representing small unsigned value v */
static int mp_set_uint(mp_int *dst, unsigned int v)
{
    if (!dst) return FAILURE;
    mp_free(dst);
    mp_init(dst);
    if (v == 0) {
        dst->sign = 0;
        dst->length = 0;
        return SUCCESS;
    }
    if (mp_reserve(dst, 1) != SUCCESS) return FAILURE;
    dst->digits[0] = v;
    dst->length = 1;
    dst->sign = +1;
    return SUCCESS;
}

int mp_fact(mp_int *r, const mp_int *a)
{
    mp_int i, tmp, acc, n_mp;
    size_t best_idx, k;
    int found, rc = FAILURE;
    unsigned int mul;

    if (!r || !a) return FAILURE;

    /* Negative input -> error */
    if (a->sign < 0) {
        calc_error_set();
        printf("Input of factorial must not be negative!\n");
        return FAILURE;
    }

    /* 0! and 1! */
    if (a->sign == 0 || (a->length == 1 && a->digits[0] <= 1)) {
        mp_free(r); mp_init(r);
        return mp_set_uint(r, 1U);
    }

    /* Find largest table entry with n <= a */
    best_idx = 0; found = 0; mp_init(&n_mp);
    for (k = 0; fact_table[k].fact_str != NULL; ++k) {
        mp_free(&n_mp); mp_init(&n_mp);
        if (mp_set_uint(&n_mp, fact_table[k].n) != SUCCESS) {
            mp_free(&n_mp); return FAILURE;
        }
        if (mp_cmp_abs(a, &n_mp) >= 0) { best_idx = k; found = 1; }
    }
    mp_free(&n_mp);

    mp_init(&acc);
    if (found) {
        if (mp_from_str_dec(&acc, fact_table[best_idx].fact_str) != SUCCESS) {
            mp_free(&acc); return FAILURE;
        }
    } else {
        if (mp_set_uint(&acc, 1U) != SUCCESS) { mp_free(&acc); return FAILURE; }
        best_idx = 0;
    }

    mp_init(&i);
    if (mp_set_uint(&i, fact_table[best_idx].n + 1U) != SUCCESS) {
        mp_free(&acc); mp_free(&i); return FAILURE;
    }

    mp_init(&tmp);

    /* Loop multiplication: acc *= i, while i <= a */
    while (mp_cmp_abs(&i, a) <= 0) {
        if (mp_fits_uint(&i)) {
            /* fast path */
            mul = i.digits[0];
            if (mp_mul_small(&acc, mul) != SUCCESS) goto fail;
        } else {
            /* general multiply */
            if (mp_mul(&tmp, &acc, &i) != SUCCESS) goto fail;
            if (mp_copy(&acc, &tmp) != SUCCESS) goto fail;
        }
        if (mp_inc(&i) != SUCCESS) goto fail;
    }

    mp_free(r); mp_init(r);
    if (mp_copy(r, &acc) != SUCCESS) goto fail;
    rc = SUCCESS;

fail:
    mp_free(&i); mp_free(&tmp); mp_free(&acc);
    return rc;
}