#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mp_int.h"
#include "fact.h"
#include "precomputed_fact.h"

/* ----------------- Factorial using precomputed table + productRange ----------------- */

const FactEntry fact_table[] = {
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
    /* Beyond 256 the strings are too long for ANSI C */
    {0, NULL} /* sentinel */
};

/* Helper: create mp_int representing small unsigned value */
static int mp_set_uint(mp_int *dst, unsigned int v) {
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

/* productRange(low, high, out)
   computes product of integers in [low..high] into out (mp_int).
   returns SUCCESS/FAILURE. ANSI C90-compliant. */
static int productRange(unsigned int low, unsigned int high, mp_int *out) {
    unsigned int mid;
    mp_int left, right, tmp;
    int status;

    if (!out) return FAILURE;

    mp_init(&left); mp_init(&right); mp_init(&tmp);

    if (low > high) {
        /* empty product = 1 */
        status = mp_set_uint(out, 1);
        mp_free(&left); mp_free(&right); mp_free(&tmp);
        return status;
    }

    if (low == high) {
        status = mp_set_uint(out, low);
        mp_free(&left); mp_free(&right); mp_free(&tmp);
        return status;
    }

    if (high - low == 1) {
        /* out = low * high */
        if (mp_set_uint(&left, low) != SUCCESS) { mp_free(&left); mp_free(&right); mp_free(&tmp); return FAILURE; }
        if (mp_set_uint(&right, high) != SUCCESS) { mp_free(&left); mp_free(&right); mp_free(&tmp); return FAILURE; }
        if (mp_mul(out, &left, &right) != SUCCESS) { mp_free(&left); mp_free(&right); mp_free(&tmp); return FAILURE; }
        mp_free(&left); mp_free(&right); mp_free(&tmp);
        return SUCCESS;
    }

    /* small range: iterative multiply using mp_mul_small for speed */
    if (high - low < 10) {
        if (mp_set_uint(out, low) != SUCCESS) { mp_free(&left); mp_free(&right); mp_free(&tmp); return FAILURE; }
        {
            unsigned int i;
            for (i = low + 1; i <= high; ++i) {
                if (mp_mul_small(out, i) != SUCCESS) {
                    mp_free(&left); mp_free(&right); mp_free(&tmp);
                    return FAILURE;
                }
            }
        }
        mp_free(&left); mp_free(&right); mp_free(&tmp);
        return SUCCESS;
    }

    /* divide and conquer */
    mid = (low + high) / 2;

    if (productRange(low, mid, &left) != SUCCESS) {
        mp_free(&left); mp_free(&right); mp_free(&tmp);
        return FAILURE;
    }
    if (productRange(mid + 1, high, &right) != SUCCESS) {
        mp_free(&left); mp_free(&right); mp_free(&tmp);
        return FAILURE;
    }

    status = mp_mul(&tmp, &left, &right);
    if (status != SUCCESS) {
        mp_free(&left); mp_free(&right); mp_free(&tmp);
        return FAILURE;
    }

    /* move tmp to out */
    mp_free(out);
    mp_init(out);
    if (mp_copy(out, &tmp) != SUCCESS) {
        mp_free(&left); mp_free(&right); mp_free(&tmp);
        return FAILURE;
    }

    mp_free(&left); mp_free(&right); mp_free(&tmp);
    return SUCCESS;
}

/* Convert mp_int (non-negative) to uint32_t if possible.
   Returns SUCCESS and writes value into *out if convertible, otherwise FAILURE. */
static int mp_to_uint32(const mp_int *x, unsigned int *out) {
    if (!x || !out) return FAILURE;
    if (x->sign < 0) return FAILURE;
    if (x->length == 0) { *out = 0U; return SUCCESS; }
    if (x->length > 1) {
        /* if any higher limb non-zero -> overflow */
        if (x->length >= 2) {
            size_t i;
            for (i = 1; i < x->length; ++i) {
                if (x->digits[i] != 0U) return FAILURE;
            }
        }
    }
    /* now lower limb holds value (may be > UINT_MAX but fits in unsigned int on typical platforms) */
    /* Be strict: if value > UINT32_MAX reject */
    if (x->digits[0] > 0xFFFFFFFFU) return FAILURE; /* technically always false for 32-bit unsigned int limb */
    *out = x->digits[0];
    return SUCCESS;
}

/* mp_fact: compute factorial of integer 'a' and store in r.
   Expects 'a' to be non-negative and reasonably small (fits into uint32_t).
   Uses precomputed factorials and productRange for remaining range.
*/
int mp_fact(mp_int *r, const mp_int *a) {
    unsigned int n;
    size_t best_idx;
    int found;
    mp_int prod;
    mp_int pre;
    mp_int tmp;

    if (!r || !a) return FAILURE;

    /* negative input -> error */
    if (a->sign < 0) {
        fprintf(stderr, "error: factorial of negative number\n");
        return FAILURE;
    }

    /* convert to uint32 */
    if (mp_to_uint32(a, &n) != SUCCESS) {
        fprintf(stderr, "error: factorial argument too large or not integer\n");
        return FAILURE;
    }

    /* trivial cases */
    if (n < 2U) {
        /* r = 1 */
        mp_free(r);
        mp_init(r);
        if (mp_set_uint(r, 1U) != SUCCESS) return FAILURE;
        return SUCCESS;
    }

    /* find largest precomputed <= n */
    best_idx = 0;
    found = 0;
    {
        size_t i;
        for (i = 0; fact_table[i].fact_str != NULL; ++i) {
            if (fact_table[i].n <= n) {
                best_idx = i;
                found = 1;
            }
        }
    }

    /* if found, set pre = precomputed[best_idx], else pre = 1 */
    mp_init(&pre);
    if (found) {
        if (mp_from_str_dec(&pre, fact_table[best_idx].fact_str) != SUCCESS) {
            mp_free(&pre);
            return FAILURE;
        }
    } else {
        if (mp_set_uint(&pre, 1U) != SUCCESS) { mp_free(&pre); return FAILURE; }
        best_idx = 0;
    }

    /* If precomputed base < n, compute productRange(base+1, n) and multiply */
    mp_init(&prod);
    mp_init(&tmp);

    if (fact_table[best_idx].n < n) {
        unsigned int low = fact_table[best_idx].n + 1U;
        unsigned int high = n;
        if (productRange(low, high, &prod) != SUCCESS) {
            mp_free(&pre); mp_free(&prod); mp_free(&tmp);
            return FAILURE;
        }
        /* tmp = pre * prod */
        if (mp_mul(&tmp, &pre, &prod) != SUCCESS) {
            mp_free(&pre); mp_free(&prod); mp_free(&tmp);
            return FAILURE;
        }
        /* move tmp to r */
        mp_free(r);
        mp_init(r);
        if (mp_copy(r, &tmp) != SUCCESS) {
            mp_free(&pre); mp_free(&prod); mp_free(&tmp);
            return FAILURE;
        }
    } else {
        /* pre already equals n! */
        mp_free(r);
        mp_init(r);
        if (mp_copy(r, &pre) != SUCCESS) {
            mp_free(&pre); mp_free(&prod); mp_free(&tmp);
            return FAILURE;
        }
    }

    /* cleanup */
    mp_free(&pre);
    mp_free(&prod);
    mp_free(&tmp);
    return SUCCESS;
}