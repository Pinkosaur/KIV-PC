/*
 * fact.c -- factorial implementation using precomputed table + incremental multiplication
 *
 * mp_fact computes a! for non-negative integer a (mp_int). The implementation:
 *  - finds the largest precomputed table entry n! with n <= a
 *  - sets acc = n! (from table) and i = n+1
 *  - multiplies acc *= i, incrementing i until i > a
 *
 * This avoids forcing 'a' into a native integer type; all comparisons and increments
 * are done using mp_int operations. The function uses a fast path (mp_mul_small)
 * when the multiplier i fits in a single limb.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mp_int.h"
#include "fact.h"
#include "error.h"

/* Precomputed factorial table */
const fact_entry fact_table[] = {
    {0, "1"},
    {1, "1"},
    {2, "2"},
    {5, "120"},
    {10, "3628800"},
    {50, "30414093201713378043612608166064768844377641568960512000000000000"},
    {75, "24809140811395398091946477116594033660926243886570122837795894512655842677572867409443815424000000000000000000"},
    {100, "93326215443944152681699238856266700490715968264381621468592963895217599993229915608941463976156518286253697920827223758251185210916864000000000000000000000000"},
    {200, "788657867364790503552363213932185062295135977687173263294742533244359449963403342920304284011984623904177212138919638830257642790242637105061926624952829931113462857270763317237396988943922445621451664240254033291864131227428294853277524242407573903240321257405579568660226031904170324062351700858796178922222789623703897374720000000000000000000000000000000000000000000000000"},
    {256, "857817775342842654119082271681232625157781520279485619859655650377269452553147589377440291360451408450375885342336584306157196834693696475322289288497426025679637332563368786442675207626794560187968867971521143307702077526646451464709187326100832876325702818980773671781454170250523018608495319068138257481070252817559459476987034665712738139286205234756808218860701203611083152093501947437109101726968262861606263662435022840944191408424615936000000000000000000000000000000000000000000000000000000000000000"},
    {0, NULL}
};

/* Helper: create mp_int representing small unsigned value v
   Returns SUCCESS or FAILURE. */
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

/*
 * mp_fact
 *   Compute factorial r = a!
 *
 *   This implementation does not require 'a' to fit in any native integer type.
 */
int mp_fact(mp_int *r, const mp_int *a)
{
    mp_int i;
    mp_int tmp;
    mp_int acc;
    size_t best_idx;
    int found;
    size_t k;
    mp_int n_mp;
    int rc = FAILURE;
    unsigned int mul;

    if (!r || !a) return FAILURE;

    /* negative -> error */
    if (a->sign < 0) {
        calc_error_set();
        printf("Input of factorial must not be negative!\n");
        return FAILURE;
    }

    /* handle 0! and 1! quickly */
    if (a->sign == 0) {
        mp_free(r);
        mp_init(r);
        if (mp_set_uint(r, 1U) != SUCCESS) return FAILURE;
        return SUCCESS;
    }
    if (a->length == 1 && a->digits[0] == 1U) {
        mp_free(r);
        mp_init(r);
        if (mp_set_uint(r, 1U) != SUCCESS) return FAILURE;
        return SUCCESS;
    }

    /* Find largest table entry with n <= a (compare as mp_int) */
    best_idx = 0;
    found = 0;
    mp_init(&n_mp);

    for (k = 0; fact_table[k].fact_str != NULL; ++k) {
        /* represent the integer n as mp_int in n_mp */
        mp_free(&n_mp);
        mp_init(&n_mp);
        if (mp_set_uint(&n_mp, fact_table[k].n) != SUCCESS) {
            mp_free(&n_mp);
            return FAILURE;
        }
        /* if a >= n_mp then table entry applies */
        if (mp_cmp_abs(a, &n_mp) >= 0) {
            best_idx = k;
            found = 1;
        }
    }
    mp_free(&n_mp);

    /* acc = precomputed factorial (or 1 if none found) */
    mp_init(&acc);
    if (found) {
        if (mp_from_str_dec(&acc, fact_table[best_idx].fact_str) != SUCCESS) {
            mp_free(&acc);
            return FAILURE;
        }
    } else {
        if (mp_set_uint(&acc, 1U) != SUCCESS) { mp_free(&acc); return FAILURE; }
        best_idx = 0;
    }

    /* i = fact_table[best_idx].n + 1 */
    mp_init(&i);
    if (mp_set_uint(&i, fact_table[best_idx].n + 1U) != SUCCESS) {
        mp_free(&acc);
        mp_free(&i);
        return FAILURE;
    }

    mp_init(&tmp);

    /* Multiply acc by i while i <= a */
    while (mp_cmp_abs(&i, a) <= 0) {
        if (mp_fits_uint(&i)) {
            /* fast path: multiplier fits in single limb */
            mul = i.digits[0];
            if (mp_mul_small(&acc, mul) != SUCCESS) goto fail;
        } else {
            /* general multiply */
            if (mp_mul(&tmp, &acc, &i) != SUCCESS) goto fail;
            mp_free(&acc);
            mp_init(&acc);
            if (mp_copy(&acc, &tmp) != SUCCESS) goto fail;
        }

        if (mp_inc(&i) != SUCCESS) goto fail;
    }

    /* move acc into r */
    mp_free(r);
    mp_init(r);
    if (mp_copy(r, &acc) != SUCCESS) goto fail;

    rc = SUCCESS;

    /* cleanup */
    mp_free(&i);
    mp_free(&tmp);
    mp_free(&acc);
    return rc;

fail:
    mp_free(&i);
    mp_free(&tmp);
    mp_free(&acc);
    return FAILURE;
}
