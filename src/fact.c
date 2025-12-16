/*
  Prime-Swing factorial implementation (ANSI C90)

  Notes:
   - Uses an internal simple sieve of Eratosthenes on unsigned long n.
   - If the argument `a` does not fit in unsigned long, falls back to
     a slower mp_int iterative factorial (correct but much slower).
   - Uses existing mp_int API: mp_init/mp_free/mp_copy/mp_mul/mp_mul_small/mp_add_small/mp_reserve/mp_from_str_dec/mp_cmp_abs.
   - No <stdint.h> used. Works in both 32-bit and 64-bit hosts; the code
     detects whether unsigned long has enough bits to hold the requested n.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h> /* for ULONG_MAX */
#include "mp_int.h"
#include "fact.h"

/* Precomputed factorial table (keeps small factorials for speed) */
typedef struct {
    unsigned int n;
    const char *fact_str;
} FactEntry;

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
    {0, NULL}
};

/* -------------------------- small helpers --------------------------- */

/* Set mp_int to unsigned int v (works even if v == 0) */
static int mp_set_uint(mp_int *dst, unsigned int v) {
    if (!dst) return FAILURE;
    mp_free(dst);
    mp_init(dst);
    if (v == 0U) {
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

/* Set mp_int from an unsigned long value (portable).
   Produces the minimal number of 32-bit limbs required.
*/
static int mp_set_ulong(mp_int *dst, unsigned long v) {
    unsigned long tmpv;
    unsigned int limb;
    size_t needed = 0;
    size_t i;

    if (!dst) return FAILURE;
    mp_free(dst);
    mp_init(dst);

    /* quick path */
    if (v == 0UL) {
        dst->sign = 0;
        dst->length = 0;
        return SUCCESS;
    }

    /* Count limbs needed: each limb is 32 bits */
    tmpv = v;
#if ULONG_MAX <= 0xFFFFFFFFUL
    /* unsigned long is 32-bit or less: only one limb possible */
    needed = 1;
#else
    /* unsigned long is wider (likely 64-bit): count how many 32-bit chunks required */
    tmpv = v;
    while (tmpv != 0UL) {
        needed++;
        tmpv = (tmpv >> 32);
    }
#endif

    if (mp_reserve(dst, needed) != SUCCESS) return FAILURE;
    /* fill limbs (little-endian) */
#if ULONG_MAX <= 0xFFFFFFFFUL
    /* only one limb */
    dst->digits[0] = (unsigned int)(v & 0xFFFFFFFFUL);
    dst->length = 1;
    dst->sign = +1;
#else
    tmpv = v;
    for (i = 0; i < needed; ++i) {
        limb = (unsigned int)(tmpv & 0xFFFFFFFFUL);
        dst->digits[i] = limb;
        tmpv = (tmpv >> 32);
    }
    dst->length = needed;
    dst->sign = +1;
#endif
    return SUCCESS;
}

/* Increment mp_int by 1 using mp_add_small */
static int mp_inc(mp_int *x) {
    if (!x) return FAILURE;
    return mp_add_small(x, 1U);
}

/* Test if mp_int fits into unsigned long; if so, write value to *out and return SUCCESS.
   Conservative: computes bitlength and compares against host unsigned long width.
*/
static int mp_to_ulong_safe(const mp_int *x, unsigned long *out) {
    size_t bits = 0;
    unsigned int top;
    size_t top_idx;
    unsigned int tb;
    unsigned int i;
    unsigned long res = 0UL;
    size_t limb_count;
    size_t ulong_bits = sizeof(unsigned long) * 8;

    if (!x || !out) return FAILURE;
    if (x->sign < 0) return FAILURE;
    if (x->length == 0) { *out = 0UL; return SUCCESS; }

    /* compute bitlength */
    top_idx = x->length - 1;
    top = x->digits[top_idx];
    tb = 0;
    while (top) { top >>= 1; ++tb; }
    bits = (size_t)top_idx * 32 + (size_t)tb;

    if (bits > ulong_bits) return FAILURE;

    /* Reconstruct value: go from high limb to low */
#if ULONG_MAX <= 0xFFFFFFFFUL
    /* unsigned long is 32-bit => only one limb permitted */
    if (x->length > 1) return FAILURE;
    res = (unsigned long)(x->digits[0]);
    *out = res;
    return SUCCESS;
#else
    /* unsigned long is at least 64-bit; can support multiple 32-bit limbs */
    res = 0UL;
    for (i = (unsigned int)x->length; i-- > 0;) {
        /* shift res left by 32 bits and add limb */
        res = (res << 32) | (unsigned long)x->digits[i];
    }
    *out = res;
    return SUCCESS;
#endif
}

/* exponent of prime p in n! : E(n, p) = sum_{k>=1} floor(n / p^k) */
static unsigned long prime_factor_exponent(unsigned long n, unsigned long p) {
    unsigned long cnt = 0UL;
    unsigned long q = p;
    if (p == 0UL) return 0UL;
    while (n / p > 0UL) {
        n /= p;
        cnt += n;
    }
    return cnt;
}

/* Fast power: compute base^exp into 'res' (res is initialized by this func).
   base is unsigned long and may be larger than a single limb; we use mp_set_ulong.
   Implements exponentiation by squaring.
*/
static int mp_pow_ulong(mp_int *res, unsigned long base, unsigned long exp) {
    mp_int b, tmp;
    int status = FAILURE;

    mp_init(&b); mp_init(&tmp);
    if (mp_set_ulong(&b, base) != SUCCESS) goto cleanup;

    if (mp_set_uint(res, 1U) != SUCCESS) goto cleanup;

    while (exp > 0UL) {
        if (exp & 1UL) {
            if (mp_mul(&tmp, res, &b) != SUCCESS) goto cleanup;
            mp_free(res);
            mp_init(res);
            if (mp_copy(res, &tmp) != SUCCESS) goto cleanup;
            mp_free(&tmp);
            mp_init(&tmp);
        }
        exp >>= 1;
        if (exp) {
            if (mp_mul(&tmp, &b, &b) != SUCCESS) goto cleanup;
            mp_free(&b);
            mp_init(&b);
            if (mp_copy(&b, &tmp) != SUCCESS) goto cleanup;
            mp_free(&tmp);
            mp_init(&tmp);
        }
    }

    status = SUCCESS;

cleanup:
    mp_free(&b);
    mp_free(&tmp);
    return status;
}

/* Compute swing(n) as product of primes^e where e = E(n,p) - 2*E(n/2,p).
   n is unsigned long (must fit), primes[] is an array of primes <= n, prime_count count.
   Result is returned in 'swing' (initialized by caller).
*/
static int compute_swing_from_primes(unsigned long n,
                                    const unsigned long *primes, size_t prime_count,
                                    mp_int *swing)
{
    size_t i;
    mp_int tmp_pow;
    mp_int acc;
    int status = FAILURE;

    mp_init(&tmp_pow);
    mp_init(&acc);
    if (mp_set_uint(&acc, 1U) != SUCCESS) goto cleanup;

    for (i = 0; i < prime_count; ++i) {
        unsigned long p = primes[i];
        unsigned long e_n = 0UL, e_half = 0UL, e;
        unsigned long nn;

        /* compute exponent E(n, p) and E(n/2, p) */
        nn = n;
        while (nn) { nn /= p; e_n += nn; }

        nn = n / 2UL;
        while (nn) { nn /= p; e_half += nn; }

        if (e_n <= 2UL * e_half) continue;
        e = e_n - 2UL * e_half;
        if (e == 0UL) continue;

        /* compute p^e */
        if (mp_pow_ulong(&tmp_pow, p, e) != SUCCESS) goto cleanup;

        /* acc *= tmp_pow */
        if (mp_mul(&acc, &acc, &tmp_pow) != SUCCESS) goto cleanup;
        mp_free(&tmp_pow); mp_init(&tmp_pow);
    }

    /* move acc into swing */
    mp_free(swing);
    mp_init(swing);
    if (mp_copy(swing, &acc) != SUCCESS) goto cleanup;

    status = SUCCESS;

cleanup:
    mp_free(&tmp_pow);
    mp_free(&acc);
    return status;
}

/* Simple sieve producing primes <= n.
   Returns dynamically allocated array of unsigned long primes (caller must free).
   On success returns pointer and sets *out_count. On failure returns NULL.
*/
static unsigned long *sieve_primes(unsigned long n, size_t *out_count) {
    unsigned char *is_comp = NULL;
    unsigned long i, j;
    unsigned long limit;
    unsigned long *primes = NULL;
    size_t cap = 0, cnt = 0;

    if (n < 2UL) {
        *out_count = 0;
        return NULL;
    }

    is_comp = (unsigned char*) malloc((size_t)(n + 1UL));
    if (!is_comp) return NULL;
    memset(is_comp, 0, (size_t)(n + 1UL));

    /* 0 and 1 are not prime */
    is_comp[0] = 1; is_comp[1] = 1;

    /* sieve */
    limit = (unsigned long) ( (unsigned long) ( (unsigned long) ( (unsigned long)1) << (sizeof(unsigned long)*4) ) ); /* unused; placeholder */
    /* compute integer sqrt(n) in unsigned long without math.h */
    limit = 0UL;
    while ((limit+1UL) * (limit+1UL) <= n) ++limit;

    for (i = 2UL; i <= limit; ++i) {
        if (!is_comp[i]) {
            for (j = i * i; j <= n; j += i) {
                is_comp[j] = 1;
            }
        }
    }

    /* collect primes */
    cap = 1024;
    primes = (unsigned long*) malloc(cap * sizeof(unsigned long));
    if (!primes) { free(is_comp); return NULL; }
    cnt = 0;
    for (i = 2UL; i <= n; ++i) {
        if (!is_comp[i]) {
            if (cnt >= cap) {
                unsigned long *newp;
                size_t newcap = cap * 2;
                newp = (unsigned long*) realloc(primes, newcap * sizeof(unsigned long));
                if (!newp) { free(primes); free(is_comp); return NULL; }
                primes = newp;
                cap = newcap;
            }
            primes[cnt++] = i;
        }
    }

    free(is_comp);
    *out_count = cnt;
    return primes;
}

/* Recursive Prime-Swing factorial for n (unsigned long).
   Computes r = n! .
   Uses precomputed table for small n.
*/
static int mp_fact_swing_rec(unsigned long n, mp_int *r) {
    int status = FAILURE;
    mp_int half_fact;
    mp_int half_sq;
    mp_int swing;
    unsigned long *primes = NULL;
    size_t prime_count = 0;

    /* Base cases: small n handled by table */
    {
        size_t k;
        for (k = 0; fact_table[k].fact_str != NULL; ++k) {
            if (fact_table[k].n == (unsigned int) n) {
                /* exact match */
                mp_free(r);
                mp_init(r);
                if (mp_from_str_dec(r, fact_table[k].fact_str) != SUCCESS) return FAILURE;
                return SUCCESS;
            }
            /* if table entry greater than n, break - but table not necessarily sorted strictly; we simply continue */
        }
    }

    if (n < 2UL) {
        mp_free(r);
        mp_init(r);
        if (mp_set_uint(r, 1U) != SUCCESS) return FAILURE;
        return SUCCESS;
    }

    /* Recursively compute (n/2)! */
    mp_init(&half_fact);
    if (mp_fact_swing_rec(n / 2UL, &half_fact) != SUCCESS) {
        mp_free(&half_fact); return FAILURE;
    }

    /* half_sq = half_fact * half_fact */
    mp_init(&half_sq);
    if (mp_mul(&half_sq, &half_fact, &half_fact) != SUCCESS) {
        mp_free(&half_fact); mp_free(&half_sq); return FAILURE;
    }

    /* compute primes up to n (required for swing) */
    primes = sieve_primes(n, &prime_count);
    if (!primes && prime_count != 0) {
        mp_free(&half_fact); mp_free(&half_sq); return FAILURE;
    }

    /* compute swing(n) */
    mp_init(&swing);
    if (compute_swing_from_primes(n, primes, prime_count, &swing) != SUCCESS) {
        mp_free(&half_fact); mp_free(&half_sq); mp_free(&swing);
        if (primes) free(primes);
        return FAILURE;
    }

    /* r = half_sq * swing */
    mp_free(r);
    mp_init(r);
    if (mp_mul(r, &half_sq, &swing) != SUCCESS) {
        mp_free(&half_fact); mp_free(&half_sq); mp_free(&swing);
        if (primes) free(primes);
        return FAILURE;
    }

    /* cleanup */
    mp_free(&half_fact);
    mp_free(&half_sq);
    mp_free(&swing);
    if (primes) free(primes);

    return SUCCESS;
}

/* Slow fallback factorial: iterate i = start..a and multiply (works for giant mp_int arguments).
   This is a correct but slow O(n) approach; used when a does not fit in unsigned long.
*/
static int mp_fact_slow(mp_int *r, const mp_int *a) {
    mp_int i, prod, tmp;
    int status = FAILURE;

    if (!r || !a) return FAILURE;

    mp_init(&i); mp_init(&prod); mp_init(&tmp);

    /* prod = 1 */
    if (mp_set_uint(&prod, 1U) != SUCCESS) goto cleanup;

    /* i = 1 */
    if (mp_set_uint(&i, 1U) != SUCCESS) goto cleanup;

    while (mp_cmp_abs(&i, a) <= 0) {
        /* multiply prod *= i */
        if (mp_fits_uint(&i)) {
            if (mp_mul_small(&prod, i.digits[0]) != SUCCESS) goto cleanup;
        } else {
            if (mp_mul(&tmp, &prod, &i) != SUCCESS) goto cleanup;
            mp_free(&prod);
            mp_init(&prod);
            if (mp_copy(&prod, &tmp) != SUCCESS) goto cleanup;
        }

        if (mp_inc(&i) != SUCCESS) goto cleanup;
    }

    mp_free(r);
    mp_init(r);
    if (mp_copy(r, &prod) != SUCCESS) goto cleanup;

    status = SUCCESS;

cleanup:
    mp_free(&i); mp_free(&prod); mp_free(&tmp);
    return status;
}

/* Public wrapper: convert argument to unsigned long if possible and use Prime-Swing,
   otherwise fall back to slow method.
*/
int mp_fact(mp_int *r, const mp_int *a) {
    unsigned long n_ul;
    int ok;

    if (!r || !a) return FAILURE;
    if (a->sign < 0) {
        fprintf(stderr, "error: factorial of negative number\n");
        return FAILURE;
    }

    ok = mp_to_ulong_safe(a, &n_ul);
    if (!ok) {
        /* fallback slow method (a too large to sieve) */
        return mp_fact_slow(r, a);
    }

    /* Use Prime-Swing recursion for n_ul */
    return mp_fact_swing_rec(n_ul, r);
}
