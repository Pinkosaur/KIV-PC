#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mp_int.h"

#define NAIVE_THRESHOLD 16


/* ---------------------------------- Memory management helpers -------------------------------------- */

/* Initializes mp_int x to empty state */
int mp_init(mp_int *x) {
    if (!x) return FAILURE;
    x->sign = 0;
    x->length = 0;
    x->digits = NULL;
    x->capacity = 0;
    return SUCCESS;
}

/* Frees memory used by the digits and reinitializes the mp_int */
int mp_free(mp_int *x) {
    if (!x) return SUCCESS;
    free(x->digits);
    /* reinitialize */
    x->digits = NULL;
    x->length = 0;
    x->capacity = 0;
    x->sign = 0;
    return SUCCESS;
}

/* Reallocates memory for the new desired size of a mp_int */
int mp_reserve(mp_int *x, size_t capacity) {
    unsigned int *new_digits;
    if (!x) return FAILURE;
    if (capacity <= x->capacity)
        return SUCCESS;
    /* allocate at least 1 to avoid realloc(NULL,0) ambiguity */
    if (capacity == 0) capacity = 1;
    new_digits = (unsigned int*) realloc(x->digits, capacity * sizeof(unsigned int));
    if (!new_digits)
        return FAILURE; /* out of memory */
    x->digits = new_digits;
    x->capacity = capacity;
    return SUCCESS;
}

/* Copy mp_int src -> dst (dst must be initialized with mp_init) */
int mp_copy(mp_int *dst, const mp_int *src) {
    if (!dst || !src) return FAILURE;
    if (src->length == 0) {
        dst->length = 0;
        dst->sign = 0;
        return SUCCESS;
    }
    if (mp_reserve(dst, src->length) != SUCCESS) return FAILURE;
    memcpy(dst->digits, src->digits, src->length * sizeof(unsigned int));
    dst->length = src->length;
    dst->sign = src->sign;
    return SUCCESS;
}


/* ----------------------------------------- Arithmetic logic ------------------------------------------- */


/* Adds a mp_int and a regular unsigned int */
int mp_add_small(mp_int *x, unsigned int a)
{
    unsigned long carry;
    size_t i;
    unsigned long sum;

    if (!x) return FAILURE;
    carry = a;
    i = 0;
    while (carry && i < x->length) {
        sum = (unsigned long)x->digits[i] + carry;
        x->digits[i] = (unsigned int)(sum & 0xFFFFFFFFUL);
        carry = sum >> 32;
        i++;
    }
    if (carry) {
        if (mp_reserve(x, x->length + 1) != SUCCESS)
            return FAILURE;
        x->digits[x->length++] = (unsigned int)carry;
    }
    if (x->sign == 0 && x->length > 0)
        x->sign = +1;
    return SUCCESS;
}

/* Multiplies a mp_int by a regular unsigned int */
int mp_mul_small(mp_int *x, unsigned int m)
{
    unsigned long carry;
    size_t i;
    unsigned long prod;

    if (!x) return FAILURE;
    if (x->sign == 0)
        return SUCCESS; /* 0 * m = 0 */

    carry = 0;

    for (i = 0; i < x->length; ++i) {
        prod = (unsigned long)x->digits[i] * m + carry;
        x->digits[i] = (unsigned int)(prod & 0xFFFFFFFFUL);
        carry = prod >> 32;
    }
    if (carry) {
        if (mp_reserve(x, x->length + 1) != SUCCESS)
            return FAILURE;
        x->digits[x->length++] = (unsigned int)carry;
    }
    return SUCCESS;
}

/* Divides a mp_int by a regular unsigned int */
int mp_div_small(mp_int *result, const mp_int *a, unsigned int divisor, unsigned int *remainder) {
    unsigned long r;
    size_t i, len;
    unsigned long cur;

    if (!result || !a || divisor == 0) return FAILURE;
    if (a->sign == 0) {
        result->sign = 0;
        result->length = 0;
        if (remainder) *remainder = 0;
        return SUCCESS;
    }

    if (mp_reserve(result, a->length) != SUCCESS) return FAILURE;
    r = 0;
    i = a->length;
    while (i > 0) {
        i--;
        cur = (r << 32) | (unsigned long)a->digits[i];
        result->digits[i] = (unsigned int)(cur / divisor);
        r = cur % divisor;
    }

    len = a->length;
    while (len > 0 && result->digits[len - 1] == 0)
        len--;
    result->length = len;
    result->sign = (len == 0) ? 0 : a->sign;
    if (remainder) *remainder = (unsigned int)r;
    return SUCCESS;
}

/* Compare absolute values of two mp_ints */
/* Returns 1 if a is greater, -1 if b is greater, otherwise 0 */
int mp_cmp_abs(const mp_int *a, const mp_int *b) {
    size_t i;
    if (!a || !b) return 0;
    if (a->length != b->length)
        return (a->length > b->length) ? 1 : -1;

    i = a->length;
    while (i > 0) {
        i--;
        if (a->digits[i] > b->digits[i]) return 1;
        if (a->digits[i] < b->digits[i]) return -1;
    }
    return 0;
}

/* Unsigned addition of mp_int (absolute) */
int mp_add_abs(mp_int *result, const mp_int *a, const mp_int *b) {
    unsigned long carry, av, bv, sum;
    size_t i, n;

    if (!result || !a || !b) return FAILURE;
    n = (a->length > b->length) ? a->length : b->length;
    if (mp_reserve(result, n + 1) != SUCCESS) return FAILURE;

    carry = 0;

    for (i = 0; i < n; ++i) {
        av = (i < a->length) ? a->digits[i] : 0;
        bv = (i < b->length) ? b->digits[i] : 0;
        sum = av + bv + carry;
        result->digits[i] = (unsigned int)(sum & 0xFFFFFFFFu);
        carry = sum >> 32;
    }
    if (carry) {
        result->digits[n++] = (unsigned int)carry;
    }
    result->length = n;
    result->sign = (n == 0) ? 0 : +1;
    return SUCCESS;
}

/* Unsigned subtraction: result = |a| - |b|  (assumes |a| >= |b|) */
int mp_sub_abs(mp_int *result, const mp_int *a, const mp_int *b) {
    size_t n;
    unsigned long av, bv, diff, borrow;
    size_t i;

    if (!result || !a || !b) return FAILURE;
    n = a->length;
    if (mp_reserve(result, n) != SUCCESS) return FAILURE;

    borrow = 0;
    for (i = 0; i < n; ++i) {
        av = (unsigned long)a->digits[i];
        bv = (i < b->length) ? (unsigned long)b->digits[i] : 0;
        diff = av - bv - borrow;
        if (av < bv + borrow)
            borrow = 1;
        else
            borrow = 0;
        result->digits[i] = (unsigned int)(diff & 0xFFFFFFFFUL);
    }

    /* Remove leading zero limbs */
    while (n > 0 && result->digits[n - 1] == 0)
        n--;
    result->length = n;
    result->sign = (n == 0) ? 0 : +1;
    return SUCCESS;
}


/* Copy-based mp_add (signed) */
int mp_add(mp_int *result, const mp_int *a, const mp_int *b) {
    if (!result || !a || !b) return FAILURE;
    if (a->sign == 0) return mp_copy(result, b);
    if (b->sign == 0) return mp_copy(result, a);

    if (a->sign == b->sign) {
        mp_add_abs(result, a, b);
        result->sign = a->sign;
    } else {
        int cmp = mp_cmp_abs(a, b);
        if (cmp == 0) {
            /* a == -b */
            result->sign = 0;
            result->length = 0;
        } else if (cmp > 0) {
            mp_sub_abs(result, a, b);
            result->sign = a->sign;
        } else {
            mp_sub_abs(result, b, a);
            result->sign = b->sign;
        }
    }
    return SUCCESS;
}

/* Subtracts mp_int b from mp_int a */
int mp_sub(mp_int *result, const mp_int *a, const mp_int *b) {
    mp_int tmp;

    if (!result || !a || !b) return FAILURE;
    mp_init(&tmp);
    if (mp_copy(&tmp, b) != SUCCESS) { mp_free(&tmp); return FAILURE; }
    tmp.sign = -tmp.sign;
    mp_add(result, a, &tmp);
    mp_free(&tmp);
    return SUCCESS;
}

/* Naively multiplies two mp_ints, suitable for smaller numbers*/
int mp_mul_naive(mp_int *result, const mp_int *a, const mp_int *b) {
    size_t n, m, i, j, len;
    unsigned long carry, av, idx;

    if (!result || !a || !b) return FAILURE;
    if (a->sign == 0 || b->sign == 0) {
        result->sign = 0;
        result->length = 0;
        return SUCCESS;
    }

    n = a->length;
    m = b->length;
    if (mp_reserve(result, n + m) != SUCCESS) return FAILURE;
    memset(result->digits, 0, (n + m) * sizeof(unsigned int));

    for (i = 0; i < n; ++i) {
        carry = 0;
        av = a->digits[i];
        for (j = 0; j < m; ++j) {
            idx = (unsigned long)result->digits[i + j]
                               + av * (unsigned long)b->digits[j] + carry;
            result->digits[i + j] = (unsigned int)(idx & 0xFFFFFFFFUL);
            carry = idx >> 32;
        }
        if (carry)
            result->digits[i + m] += (unsigned int)carry;
    }

    len = n + m;
    while (len > 0 && result->digits[len - 1] == 0)
        len--;
    result->length = len;
    result->sign = a->sign * b->sign;
    return SUCCESS;
}

/* Logical left shift by k limbs (base 2^(32*k)) */
int mp_shift_left_words(mp_int *x, size_t k) {
    if (!x) return FAILURE;
    if (x->sign == 0 || k == 0) return SUCCESS;
    if (mp_reserve(x, x->length + k) != SUCCESS)
        return FAILURE;
    memmove(x->digits + k, x->digits, x->length * sizeof(unsigned int));
    memset(x->digits, 0, k * sizeof(unsigned int));
    x->length += k;
    return SUCCESS;
}

/* Helper: right-shift by k limbs (words). Destructive. */
int mp_shift_right_words(mp_int *x, size_t k) {
    size_t i;
    if (!x) return FAILURE;
    if (x->length == 0 || k == 0) return SUCCESS;
    if (k >= x->length) {
        /* becomes zero */
        x->length = 0;
        x->sign = 0;
        return SUCCESS;
    }
    /* shift down */
    for (i = 0; i + k < x->length; ++i) {
        x->digits[i] = x->digits[i + k];
    }
    /* clear high limbs */
    for (; i < x->length; ++i) x->digits[i] = 0;
    x->length -= k;
    /* trim leading zeros */
    while (x->length > 0 && x->digits[x->length - 1] == 0) x->length--;
    if (x->length == 0) x->sign = 0;
    return SUCCESS;
}

static int mp_is_zero(const mp_int *x) {
    if (!x) return 1;
    return (x->length == 0 || x->sign == 0) ? 1 : 0;
}

void mp_split(const mp_int *src, mp_int *low, mp_int *high, size_t m) {
    size_t i;

    if (!src || !low || !high) return;

    mp_reserve(low, m);
    mp_reserve(high, (src->length > m) ? src->length - m : 0);

    low->length = (src->length < m) ? src->length : m;
    for (i = 0; i < low->length; ++i)
        low->digits[i] = src->digits[i];
    low->sign = (low->length == 0) ? 0 : +1;

    if (src->length > m) {
        high->length = src->length - m;
        for (i = 0; i < high->length; ++i)
            high->digits[i] = src->digits[i + m];
        high->sign = +1;
    } else {
        high->length = 0;
        high->sign = 0;
    }
}

/* Karatsuba multiply (recursive) */
int mp_karatsuba_mul(mp_int *result, const mp_int *a, const mp_int *b)
{
    size_t n, m;
    mp_int x0, x1, y0, y1, z0, z1, z2, tmp1, tmp2, res_tmp;

    if (!result || !a || !b) return FAILURE;

    /* Base case: small numbers → naive multiply */
    if (a->length <= 16 || b->length <= 16) {
        return mp_mul_naive(result, a, b);
    }

    n = (a->length > b->length) ? a->length : b->length;
    m = n / 2;

    mp_init(&x0); mp_init(&x1);
    mp_init(&y0); mp_init(&y1);
    mp_init(&z0); mp_init(&z1); mp_init(&z2);
    mp_init(&tmp1); mp_init(&tmp2);
    mp_init(&res_tmp);

    mp_split(a, &x0, &x1, m);
    mp_split(b, &y0, &y1, m);

    /* z0 = x0 * y0 */
    mp_karatsuba_mul(&z0, &x0, &y0);

    /* z2 = x1 * y1 */
    mp_karatsuba_mul(&z2, &x1, &y1);

    /* tmp1 = x0 + x1, tmp2 = y0 + y1 */
    mp_add_abs(&tmp1, &x0, &x1);
    mp_add_abs(&tmp2, &y0, &y1);

    /* z1 = (x0+x1)*(y0+y1) - z2 - z0 */
    mp_karatsuba_mul(&z1, &tmp1, &tmp2);
    mp_sub_abs(&z1, &z1, &z2);
    mp_sub_abs(&z1, &z1, &z0);

    /* combine result = z0 + (z1 << (m words)) + (z2 << (2*m words)) */
    mp_init(&res_tmp);
    mp_add_abs(&res_tmp, &z0, &res_tmp);
    mp_shift_left_words(&z1, m);
    mp_add_abs(&res_tmp, &res_tmp, &z1);
    mp_shift_left_words(&z2, 2 * m);
    mp_add_abs(&res_tmp, &res_tmp, &z2);

    /* move res_tmp to result */
    mp_copy(result, &res_tmp);
    result->sign = a->sign * b->sign;

    /* free temps */
    mp_free(&x0); mp_free(&x1);
    mp_free(&y0); mp_free(&y1);
    mp_free(&z0); mp_free(&z1); mp_free(&z2);
    mp_free(&tmp1); mp_free(&tmp2);
    mp_free(&res_tmp);

    return SUCCESS;
}

/* Multiplies two mp_ints (hybrid) */
int mp_mul(mp_int *result, const mp_int *a, const mp_int *b) {
    if (!result || !a || !b) return FAILURE;
    if (a->length < NAIVE_THRESHOLD || b->length < NAIVE_THRESHOLD)
        return mp_mul_naive(result, a, b);
    else
        return mp_karatsuba_mul(result, a, b);
}

/* Divides mp_int b by mp_int a */
/* Note: full long division is not implemented; only small divisor supported via mp_div_small */
int mp_div(mp_int *result, const mp_int *a, const mp_int *b) {
    if (!result || !a || !b) return FAILURE;
    /* If divisor is a single limb, use mp_div_small */
    if (b->length == 1) {
        unsigned int rem;
        return mp_div_small(result, a, b->digits[0], &rem);
    }
    return FAILURE; /* not implemented */
}

/* mp_mod: r = a % b
   Remainder has same sign as 'a' (like C's %). Uses only ANSI C90 features. */
int mp_mod(mp_int *r, const mp_int *a, const mp_int *b)
{
    /* variables declared at top per C90 */
    size_t shift_words;
    int cmp;
    int status = FAILURE;

    mp_int rem;
    mp_int tmp;
    mp_int candidate;
    mp_int doubled;
    mp_int newrem;
    unsigned int rem_small;

    if (!r || !a || !b) return FAILURE;
    /* divisor zero */
    if (b->length == 0 || b->sign == 0) return FAILURE;

    /* dividend zero => remainder 0 */
    if (a->length == 0 || a->sign == 0) {
        mp_free(r);
        mp_init(r);
        r->sign = 0;
        r->length = 0;
        return SUCCESS;
    }

    /* fast path: single-limb divisor */
    if (b->length == 1) {
        mp_int q;
        mp_init(&q);
        if (mp_div_small(&q, a, b->digits[0], &rem_small) != SUCCESS) {
            mp_free(&q);
            return FAILURE;
        }
        mp_free(&q);

        mp_free(r);
        mp_init(r);
        if (rem_small == 0) {
            r->sign = 0; r->length = 0;
        } else {
            if (mp_reserve(r, 1) != SUCCESS) return FAILURE;
            r->digits[0] = rem_small;
            r->length = 1;
            r->sign = (a->sign < 0) ? -1 : +1;
        }
        return SUCCESS;
    }

    /* General multi-limb algorithm (word-shift + repeated-doubling subtraction) */

    /* initialize temporaries */
    mp_init(&rem);
    mp_init(&tmp);
    mp_init(&candidate);
    mp_init(&doubled);
    mp_init(&newrem);

    /* rem = |a| */
    if (mp_copy(&rem, a) != SUCCESS) goto cleanup;
    rem.sign = +1;

    /* tmp = |b| */
    if (mp_copy(&tmp, b) != SUCCESS) goto cleanup;
    tmp.sign = +1;

    /* if |rem| < |tmp| -> remainder is a */
    cmp = mp_cmp_abs(&rem, &tmp);
    if (cmp < 0) {
        mp_free(r);
        mp_init(r);
        if (mp_copy(r, &rem) != SUCCESS) goto cleanup;
        r->sign = (a->sign < 0) ? -1 : +1;
        status = SUCCESS;
        goto cleanup;
    }

    /* shift tmp left so its top limb lines up with rem top limb */
    shift_words = rem.length - tmp.length;
    if (shift_words > 0) {
        if (mp_shift_left_words(&tmp, shift_words) != SUCCESS) goto cleanup;
    }

    /* iterate from current shift down to 0 */
    for ( ; ; ) {
        /* while rem >= tmp, subtract the largest multiple of tmp we can get by doubling */
        while (mp_cmp_abs(&rem, &tmp) >= 0) {
            /* candidate = tmp */
            mp_free(&candidate);
            mp_init(&candidate);
            if (mp_copy(&candidate, &tmp) != SUCCESS) goto cleanup;

            /* double candidate repeatedly while doubled <= rem */
            while (1) {
                mp_free(&doubled);
                mp_init(&doubled);
                if (mp_add_abs(&doubled, &candidate, &candidate) != SUCCESS) goto cleanup;
                /* if doubled > rem break */
                if (mp_cmp_abs(&doubled, &rem) > 0) {
                    mp_free(&doubled);
                    break;
                }
                /* accept doubled as new candidate */
                mp_free(&candidate);
                mp_init(&candidate);
                if (mp_copy(&candidate, &doubled) != SUCCESS) goto cleanup;
                /* loop and try doubling again */
            }

            /* rem = rem - candidate  (candidate <= rem guaranteed) */
            mp_free(&newrem);
            mp_init(&newrem);
            if (mp_sub_abs(&newrem, &rem, &candidate) != SUCCESS) goto cleanup;
            /* move newrem -> rem */
            mp_free(&rem);
            mp_init(&rem);
            if (mp_copy(&rem, &newrem) != SUCCESS) goto cleanup;
            /* free candidate, will be recreated on next iteration if needed */
            mp_free(&candidate);
        }

        /* if we've shifted down to zero words stop */
        if (shift_words == 0) break;

        /* shift tmp right by one word and continue */
        if (mp_shift_right_words(&tmp, 1) != SUCCESS) goto cleanup;
        shift_words--;
    }

    /* rem is the absolute remainder; attach sign of a */
    mp_free(r);
    mp_init(r);
    if (mp_copy(r, &rem) != SUCCESS) goto cleanup;
    if (r->length == 0) r->sign = 0;
    else r->sign = (a->sign < 0) ? -1 : +1;

    status = SUCCESS;

cleanup:
    mp_free(&rem);
    mp_free(&tmp);
    mp_free(&candidate);
    mp_free(&doubled);
    mp_free(&newrem);
    return status;
}