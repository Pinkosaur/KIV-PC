/*
 * mp_int.c -- Multiple-Precision Integer Implementation
 *
 * This module implements arbitrary-precision arithmetic using a limb-based
 * representation. It supports signed integers with dynamic memory management.
 *
 * Algorithms:
 * - Addition/Subtraction: Standard carry/borrow propagation (Sign-Magnitude).
 * - Multiplication: Hybrid approach using Naive O(N^2) for small operands
 * and Karatsuba O(N^1.58) for larger operands.
 * - Division: Shift-and-subtract (restoring division) with memory optimization.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "mp_int.h"
#include "error.h"

/* Threshold (in limbs) for switching between naive multiplication and Karatsuba.
 * A value of 128 (approx 2048 bits on 16-bit systems) balances the overhead
 * of recursion against the asymptotic speedup.
 */
#define NAIVE_THRESHOLD 128

/* --- Helper Functions ----------------------------------------------------- */

/*
 * mp_normalize
 * Removes leading zero limbs from the number.
 * If the result is zero, ensures the sign is set to 0.
 */
static void mp_normalize(mp_int *x)
{
    size_t len;
    if (!x) return;
    len = x->length;
    while (len > 0 && x->digits[len - 1] == (mp_limb_t)0) len--;
    x->length = len;
    if (x->length == 0) x->sign = 0;
}

/* --- Memory Management ---------------------------------------------------- */

int mp_init(mp_int *x)
{
    if (!x) return FAILURE;
    x->sign = 0;
    x->length = 0;
    x->digits = NULL;
    x->capacity = 0;
    return SUCCESS;
}

int mp_free(mp_int *x)
{
    if (!x) return SUCCESS;
    if (x->digits) free(x->digits);
    x->digits = NULL;
    x->length = 0;
    x->capacity = 0;
    x->sign = 0;
    return SUCCESS;
}

int mp_reserve(mp_int *x, size_t capacity)
{
    mp_limb_t *new_digits;
    if (!x) return FAILURE;
    if (capacity <= x->capacity) return SUCCESS;
    
    if (capacity == 0) capacity = 1;
    
    new_digits = (mp_limb_t *) realloc(x->digits, capacity * sizeof(mp_limb_t));
    if (!new_digits) return FAILURE;
    
    x->digits = new_digits;
    x->capacity = capacity;
    return SUCCESS;
}

int mp_copy(mp_int *dst, const mp_int *src)
{
    if (!dst || !src) return FAILURE;
    if (src->length == 0) {
        dst->length = 0;
        dst->sign = 0;
        return SUCCESS;
    }
    if (mp_reserve(dst, src->length) != SUCCESS) return FAILURE;
    memcpy(dst->digits, src->digits, src->length * sizeof(mp_limb_t));
    dst->length = src->length;
    dst->sign = src->sign;
    return SUCCESS;
}

/* --- Arithmetic Primitives (Small Constants) ------------------------------ */

int mp_add_small(mp_int *x, unsigned int a)
{
    mp_double_t carry;
    size_t i;
    
    if (!x) return FAILURE;
    
    carry = (mp_double_t)a;
    i = 0;
    while (carry && i < x->length) {
        mp_double_t sum = (mp_double_t)x->digits[i] + carry;
        x->digits[i] = (mp_limb_t)(sum & (mp_double_t)MP_LIMB_MASK);
        carry = sum >> MP_LIMB_BITS;
        i++;
    }
    
    if (carry) {
        if (mp_reserve(x, x->length + 1) != SUCCESS) return FAILURE;
        x->digits[x->length++] = (mp_limb_t)carry;
    }
    
    if (x->sign == 0 && x->length > 0) x->sign = +1;
    return SUCCESS;
}

int mp_mul_small(mp_int *x, unsigned int m)
{
    mp_double_t carry;
    size_t i;
    
    if (!x) return FAILURE;
    if (x->sign == 0) return SUCCESS;

    carry = 0;
    for (i = 0; i < x->length; ++i) {
        mp_double_t prod = (mp_double_t)x->digits[i] * (mp_double_t)m + carry;
        x->digits[i] = (mp_limb_t)(prod & (mp_double_t)MP_LIMB_MASK);
        carry = prod >> MP_LIMB_BITS;
    }
    
    if (carry) {
        if (mp_reserve(x, x->length + 1) != SUCCESS) return FAILURE;
        x->digits[x->length++] = (mp_limb_t)carry;
    }
    return SUCCESS;
}

int mp_div_small(mp_int *result, const mp_int *a, unsigned int divisor, unsigned int *remainder)
{
    mp_double_t r;
    size_t i, len;
    mp_double_t cur;

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
        cur = (r << MP_LIMB_BITS) | (mp_double_t)a->digits[i];
        result->digits[i] = (mp_limb_t)(cur / (mp_double_t)divisor);
        r = cur % (mp_double_t)divisor;
    }

    len = a->length;
    while (len > 0 && result->digits[len - 1] == (mp_limb_t)0) len--;
    result->length = len;
    result->sign = (len == 0) ? 0 : a->sign;
    
    if (remainder) *remainder = (unsigned int)r;
    return SUCCESS;
}

int mp_inc(mp_int *x)
{
    return mp_add_small(x, 1U);
}

int mp_fits_uint(const mp_int *x)
{
    if (!x) return 0;
    if (x->sign < 0) return 0;
    if (x->length == 0) return 1;
    if (x->length > 1) return 0;
    return 1;
}

/* --- Unsigned (Absolute Value) Comparison & Arithmetic -------------------- */

int mp_cmp_abs(const mp_int *a, const mp_int *b)
{
    size_t i;
    if (!a || !b) return 0;
    if (a->length != b->length) return (a->length > b->length) ? 1 : -1;
    
    i = a->length;
    while (i > 0) {
        i--;
        if (a->digits[i] > b->digits[i]) return 1;
        if (a->digits[i] < b->digits[i]) return -1;
    }
    return 0;
}

int mp_add_abs(mp_int *result, const mp_int *a, const mp_int *b)
{
    mp_double_t carry;
    mp_double_t av, bv, sum;
    size_t i, n;

    if (!result || !a || !b) return FAILURE;
    
    n = (a->length > b->length) ? a->length : b->length;
    if (mp_reserve(result, n + 1) != SUCCESS) return FAILURE;

    carry = 0;
    for (i = 0; i < n; ++i) {
        av = (i < a->length) ? (mp_double_t)a->digits[i] : (mp_double_t)0;
        bv = (i < b->length) ? (mp_double_t)b->digits[i] : (mp_double_t)0;
        sum = av + bv + carry;
        result->digits[i] = (mp_limb_t)(sum & (mp_double_t)MP_LIMB_MASK);
        carry = sum >> MP_LIMB_BITS;
    }
    
    if (carry) {
        result->digits[n++] = (mp_limb_t)carry;
    }
    
    result->length = n;
    result->sign = (n == 0) ? 0 : +1;
    return SUCCESS;
}

int mp_sub_abs(mp_int *result, const mp_int *a, const mp_int *b)
{
    size_t n, i;
    mp_double_t av, bv, diff;
    mp_double_t borrow;

    if (!result || !a || !b) return FAILURE;
    
    n = a->length;
    if (mp_reserve(result, n) != SUCCESS) return FAILURE;

    borrow = 0;
    for (i = 0; i < n; ++i) {
        av = (mp_double_t)a->digits[i];
        bv = (i < b->length) ? (mp_double_t)b->digits[i] : (mp_double_t)0;
        
        diff = av - bv - borrow;
        if (av < bv + borrow) borrow = 1;
        else borrow = 0;
        
        result->digits[i] = (mp_limb_t)(diff & (mp_double_t)MP_LIMB_MASK);
    }

    while (n > 0 && result->digits[n - 1] == (mp_limb_t)0) n--;
    result->length = n;
    result->sign = (n == 0) ? 0 : +1;
    return SUCCESS;
}

/* --- Signed Arithmetic ---------------------------------------------------- */

int mp_add(mp_int *result, const mp_int *a, const mp_int *b)
{
    int cmp;
    if (!result || !a || !b) return FAILURE;
    
    if (a->sign == 0) return mp_copy(result, b);
    if (b->sign == 0) return mp_copy(result, a);

    if (a->sign == b->sign) {
        if (mp_add_abs(result, a, b) != SUCCESS) return FAILURE;
        result->sign = a->sign;
    } else {
        cmp = mp_cmp_abs(a, b);
        if (cmp == 0) {
            result->sign = 0;
            result->length = 0;
        } else if (cmp > 0) {
            if (mp_sub_abs(result, a, b) != SUCCESS) return FAILURE;
            result->sign = a->sign;
        } else {
            if (mp_sub_abs(result, b, a) != SUCCESS) return FAILURE;
            result->sign = b->sign;
        }
    }
    return SUCCESS;
}

int mp_sub(mp_int *result, const mp_int *a, const mp_int *b)
{
    mp_int tmp;
    if (!result || !a || !b) return FAILURE;
    
    mp_init(&tmp);
    if (mp_copy(&tmp, b) != SUCCESS) { mp_free(&tmp); return FAILURE; }
    
    tmp.sign = -tmp.sign;
    if (mp_add(result, a, &tmp) != SUCCESS) { mp_free(&tmp); return FAILURE; }
    
    mp_free(&tmp);
    return SUCCESS;
}

/* --- Multiplication ------------------------------------------------------- */

int mp_mul_naive(mp_int *result, const mp_int *a, const mp_int *b)
{
    size_t n, m, i, j;
    mp_double_t carry;
    
    if (!result || !a || !b) return FAILURE;
    
    if (a->sign == 0 || b->sign == 0) {
        result->sign = 0;
        result->length = 0;
        return SUCCESS;
    }

    n = a->length;
    m = b->length;
    if (mp_reserve(result, n + m) != SUCCESS) return FAILURE;
    memset(result->digits, 0, (n + m) * sizeof(mp_limb_t));

    for (i = 0; i < n; ++i) {
        carry = 0;
        for (j = 0; j < m; ++j) {
            mp_double_t acc = (mp_double_t)result->digits[i + j]
                              + (mp_double_t)a->digits[i] * (mp_double_t)b->digits[j]
                              + carry;
            result->digits[i + j] = (mp_limb_t)(acc & (mp_double_t)MP_LIMB_MASK);
            carry = acc >> MP_LIMB_BITS;
        }
        if (carry) {
            size_t pos = i + m;
            while (carry) {
                mp_double_t sum = (mp_double_t)((pos < result->capacity) ? result->digits[pos] : 0) + carry;
                if (pos >= result->capacity) {
                    if (mp_reserve(result, pos + 1) != SUCCESS) return FAILURE;
                }
                result->digits[pos] = (mp_limb_t)(sum & (mp_double_t)MP_LIMB_MASK);
                carry = sum >> MP_LIMB_BITS;
                pos++;
            }
        }
    }

    {
        size_t len = n + m;
        while (len > 0 && result->digits[len - 1] == (mp_limb_t)0) len--;
        result->length = len;
        result->sign = (len == 0) ? 0 : a->sign * b->sign;
    }
    return SUCCESS;
}

int mp_shift_left_words(mp_int *x, size_t k)
{
    if (!x) return FAILURE;
    if (x->sign == 0 || k == 0) return SUCCESS;
    
    if (mp_reserve(x, x->length + k) != SUCCESS) return FAILURE;
    
    memmove(x->digits + k, x->digits, x->length * sizeof(mp_limb_t));
    memset(x->digits, 0, k * sizeof(mp_limb_t));
    x->length += k;
    return SUCCESS;
}

int mp_shift_right_words(mp_int *x, size_t k)
{
    size_t i;
    if (!x) return FAILURE;
    if (x->length == 0 || k == 0) return SUCCESS;
    
    if (k >= x->length) {
        x->length = 0;
        x->sign = 0;
        return SUCCESS;
    }
    
    for (i = 0; i + k < x->length; ++i) x->digits[i] = x->digits[i + k];
    for (; i < x->length; ++i) x->digits[i] = (mp_limb_t)0;
    
    x->length -= k;
    mp_normalize(x);
    return SUCCESS;
}

void mp_split(const mp_int *src, mp_int *low, mp_int *high, size_t m)
{
    size_t i;
    if (!src || !low || !high) return;

    mp_reserve(low, m);
    mp_reserve(high, (src->length > m) ? src->length - m : 0);

    low->length = (src->length < m) ? src->length : m;
    for (i = 0; i < low->length; ++i) low->digits[i] = src->digits[i];
    low->sign = (low->length == 0) ? 0 : +1;

    if (src->length > m) {
        high->length = src->length - m;
        for (i = 0; i < high->length; ++i) high->digits[i] = src->digits[i + m];
        high->sign = +1;
    } else {
        high->length = 0;
        high->sign = 0;
    }
}

int mp_karatsuba_mul(mp_int *result, const mp_int *a, const mp_int *b)
{
    size_t n, m;
    mp_int x0, x1, y0, y1, z0, z1, z2, t1, t2, res_tmp, tmp;

    if (!result || !a || !b) return FAILURE;

    if (a->length <= NAIVE_THRESHOLD || b->length <= NAIVE_THRESHOLD) {
        return mp_mul_naive(result, a, b);
    }

    n = (a->length > b->length) ? a->length : b->length;
    m = n / 2u;

    mp_init(&x0); mp_init(&x1); mp_init(&y0); mp_init(&y1);
    mp_init(&z0); mp_init(&z1); mp_init(&z2); mp_init(&t1); mp_init(&t2);
    mp_init(&res_tmp); mp_init(&tmp);

    mp_split(a, &x0, &x1, m);
    mp_split(b, &y0, &y1, m);

    if (mp_karatsuba_mul(&z0, &x0, &y0) != SUCCESS) goto cleanup;
    if (mp_karatsuba_mul(&z2, &x1, &y1) != SUCCESS) goto cleanup;

    if (mp_add_abs(&t1, &x0, &x1) != SUCCESS) goto cleanup;
    if (mp_add_abs(&t2, &y0, &y1) != SUCCESS) goto cleanup;

    if (mp_karatsuba_mul(&z1, &t1, &t2) != SUCCESS) goto cleanup;
    if (mp_sub_abs(&z1, &z1, &z2) != SUCCESS) goto cleanup;
    if (mp_sub_abs(&z1, &z1, &z0) != SUCCESS) goto cleanup;

    if (mp_copy(&res_tmp, &z0) != SUCCESS) goto cleanup;

    if (mp_shift_left_words(&z1, m) != SUCCESS) goto cleanup;
    if (mp_add_abs(&tmp, &res_tmp, &z1) != SUCCESS) goto cleanup;
    if (mp_copy(&res_tmp, &tmp) != SUCCESS) goto cleanup;

    if (mp_shift_left_words(&z2, 2u * m) != SUCCESS) goto cleanup;
    if (mp_add_abs(&tmp, &res_tmp, &z2) != SUCCESS) goto cleanup;
    if (mp_copy(&res_tmp, &tmp) != SUCCESS) goto cleanup;

    if (mp_copy(result, &res_tmp) != SUCCESS) goto cleanup;
    result->sign = a->sign * b->sign;

    mp_free(&x0); mp_free(&x1); mp_free(&y0); mp_free(&y1);
    mp_free(&z0); mp_free(&z1); mp_free(&z2); mp_free(&t1); mp_free(&t2);
    mp_free(&res_tmp); mp_free(&tmp);
    return SUCCESS;

cleanup:
    mp_free(&x0); mp_free(&x1); mp_free(&y0); mp_free(&y1);
    mp_free(&z0); mp_free(&z1); mp_free(&z2); mp_free(&t1); mp_free(&t2);
    mp_free(&res_tmp); mp_free(&tmp);
    return FAILURE;
}

int mp_mul(mp_int *result, const mp_int *a, const mp_int *b)
{
    if (!result || !a || !b) return FAILURE;
    if (a->length < NAIVE_THRESHOLD || b->length < NAIVE_THRESHOLD)
        return mp_mul_naive(result, a, b);
    else
        return mp_karatsuba_mul(result, a, b);
}

/* --- Division & Modulus --------------------------------------------------- */

int mp_div(mp_int *result, const mp_int *a, const mp_int *b)
{
    mp_int rem, tmp, candidate, doubled, newrem, quotient, mult, newq;
    size_t shift_words;
    int status = FAILURE;
    size_t i;

    if (!result || !a || !b) return FAILURE;
    
    if (b->length == 0 || b->sign == 0) {
        calc_error_set();
        printf("Division by zero!\n");
        return FAILURE;
    }
    
    if (a->length == 0 || a->sign == 0) {
        mp_free(result);
        mp_init(result);
        result->sign = 0;
        result->length = 0;
        return SUCCESS;
    }

    if (b->length == 1) {
        mp_int q;
        unsigned int rem_small;
        mp_init(&q);
        if (mp_div_small(&q, a, (unsigned int)b->digits[0], &rem_small) != SUCCESS) {
            mp_free(&q);
            return FAILURE;
        }
        if (q.length == 0) q.sign = 0;
        else q.sign = (a->sign * b->sign);
        
        mp_free(result);
        mp_init(result);
        if (mp_copy(result, &q) != SUCCESS) { mp_free(&q); return FAILURE; }
        mp_free(&q);
        return SUCCESS;
    }

    mp_init(&rem); mp_init(&tmp); mp_init(&candidate); mp_init(&doubled);
    mp_init(&newrem); mp_init(&quotient); mp_init(&mult); mp_init(&newq);

    if (mp_copy(&rem, a) != SUCCESS) goto cleanup;
    rem.sign = +1;
    if (mp_copy(&tmp, b) != SUCCESS) goto cleanup;
    tmp.sign = +1;

    if (mp_cmp_abs(&rem, &tmp) < 0) {
        mp_free(result);
        mp_init(result);
        result->sign = 0;
        result->length = 0;
        status = SUCCESS;
        goto cleanup;
    }

    shift_words = 0;
    if (rem.length > tmp.length) shift_words = rem.length - tmp.length;
    if (shift_words > 0) {
        if (mp_shift_left_words(&tmp, shift_words) != SUCCESS) goto cleanup;
    }

    quotient.length = 0;
    quotient.sign = 0;

    for (;;) {
        /* Main division loop */
        while (mp_cmp_abs(&rem, &tmp) >= 0) {
            if (mp_copy(&candidate, &tmp) != SUCCESS) goto cleanup;

            if (mp_reserve(&mult, shift_words + 1) != SUCCESS) goto cleanup;
            for (i = 0; i < shift_words; ++i) mult.digits[i] = (mp_limb_t)0;
            mult.digits[shift_words] = (mp_limb_t)1;
            mult.length = shift_words + 1;
            mult.sign = +1;

            /* Doubling step to find largest multiple */
            for (;;) {
                if (mp_add_abs(&doubled, &candidate, &candidate) != SUCCESS) goto cleanup;
                if (mp_cmp_abs(&doubled, &rem) > 0) {
                    break;
                }
                if (mp_copy(&candidate, &doubled) != SUCCESS) goto cleanup;

                if (mp_add_abs(&newq, &mult, &mult) != SUCCESS) goto cleanup;
                if (mp_copy(&mult, &newq) != SUCCESS) goto cleanup;
            }

            if (mp_sub_abs(&newrem, &rem, &candidate) != SUCCESS) goto cleanup;
            if (mp_copy(&rem, &newrem) != SUCCESS) goto cleanup;

            if (mp_add_abs(&newq, &quotient, &mult) != SUCCESS) goto cleanup;
            if (mp_copy(&quotient, &newq) != SUCCESS) goto cleanup;
        }

        if (shift_words == 0) break;

        if (mp_shift_right_words(&tmp, 1) != SUCCESS) goto cleanup;
        shift_words--;
    }

    mp_free(result);
    mp_init(result);
    if (mp_copy(result, &quotient) != SUCCESS) goto cleanup;
    if (result->length == 0) result->sign = 0;
    else result->sign = (a->sign * b->sign);

    status = SUCCESS;

cleanup:
    mp_free(&rem); mp_free(&tmp); mp_free(&candidate);
    mp_free(&doubled); mp_free(&newrem); mp_free(&quotient);
    mp_free(&mult); mp_free(&newq);
    return status;
}

int mp_mod(mp_int *r, const mp_int *a, const mp_int *b)
{
    mp_int rem, tmp, candidate, doubled, newrem;
    size_t shift_words;
    int cmp;
    int status = FAILURE;

    if (!r || !a || !b) return FAILURE;
    
    if (b->length == 0 || b->sign == 0) {
        calc_error_set();
        printf("Division by zero!\n");
        return FAILURE;
    }

    if (a->length == 0 || a->sign == 0) {
        mp_free(r);
        mp_init(r);
        r->sign = 0;
        r->length = 0;
        return SUCCESS;
    }

    if (b->length == 1) {
        mp_int q;
        unsigned int rem_small;
        mp_init(&q);
        if (mp_div_small(&q, a, (unsigned int)b->digits[0], &rem_small) != SUCCESS) {
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
            r->digits[0] = (mp_limb_t)rem_small;
            r->length = 1;
            r->sign = (a->sign < 0) ? -1 : +1;
        }
        return SUCCESS;
    }

    mp_init(&rem); mp_init(&tmp); mp_init(&candidate); mp_init(&doubled); mp_init(&newrem);

    if (mp_copy(&rem, a) != SUCCESS) goto cleanup;
    rem.sign = +1;

    if (mp_copy(&tmp, b) != SUCCESS) goto cleanup;
    tmp.sign = +1;

    cmp = mp_cmp_abs(&rem, &tmp);
    if (cmp < 0) {
        mp_free(r);
        mp_init(r);
        if (mp_copy(r, &rem) != SUCCESS) goto cleanup;
        r->sign = (a->sign < 0) ? -1 : +1;
        status = SUCCESS;
        goto cleanup;
    }

    shift_words = rem.length - tmp.length;
    if (shift_words > 0) if (mp_shift_left_words(&tmp, shift_words) != SUCCESS) goto cleanup;

    for (;;) {
        while (mp_cmp_abs(&rem, &tmp) >= 0) {
            if (mp_copy(&candidate, &tmp) != SUCCESS) goto cleanup;

            while (1) {
                if (mp_add_abs(&doubled, &candidate, &candidate) != SUCCESS) goto cleanup;
                if (mp_cmp_abs(&doubled, &rem) > 0) { break; }
                if (mp_copy(&candidate, &doubled) != SUCCESS) goto cleanup;
            }

            if (mp_sub_abs(&newrem, &rem, &candidate) != SUCCESS) goto cleanup;
            if (mp_copy(&rem, &newrem) != SUCCESS) goto cleanup;
        }

        if (shift_words == 0) break;
        if (mp_shift_right_words(&tmp, 1) != SUCCESS) goto cleanup;
        shift_words--;
    }

    mp_free(r);
    mp_init(r);
    if (mp_copy(r, &rem) != SUCCESS) goto cleanup;
    if (r->length == 0) r->sign = 0;
    else r->sign = (a->sign < 0) ? -1 : +1;

    status = SUCCESS;

cleanup:
    mp_free(&rem); mp_free(&tmp); mp_free(&candidate);
    mp_free(&doubled); mp_free(&newrem);
    return status;
}

/* --- String Conversions --------------------------------------------------- */

int mp_from_str_dec(mp_int *x, const char *str)
{
    const char *p;
    int sign;
    unsigned int d;

    if (!x || !str) return FAILURE;
    mp_free(x);
    mp_init(x);

    p = str;
    while (isspace((unsigned char)*p)) p++;

    sign = +1;
    if (*p == '-') { sign = -1; p++; }
    else if (*p == '+') { p++; }

    while (*p == '0') p++;

    if (*p == '\0') {
        x->sign = 0;
        x->length = 0;
        return SUCCESS;
    }

    if (!isdigit((unsigned char)*p)) {
        return FAILURE;
    }

    if (mp_reserve(x, 1) != SUCCESS) return FAILURE;
    x->digits[0] = 0;
    x->length = 1;
    x->sign = +1;

    for (; *p; ++p) {
        if (!isdigit((unsigned char)*p)) break;
        d = (unsigned int)(*p - '0');
        if (mp_mul_small(x, 10U) != SUCCESS) return FAILURE;
        if (mp_add_small(x, d) != SUCCESS) return FAILURE;
    }

    while (isspace((unsigned char)*p)) p++;
    if (*p != '\0') return FAILURE;

    if (x->length == 1 && x->digits[0] == 0)
        x->sign = 0;
    else
        x->sign = sign;
    return SUCCESS;
}

int mp_from_str_bin(mp_int *x, const char *str) {
    size_t i;
    const char *p;
    const char *digits;
    size_t digits_len;
    unsigned int bit;
    int explicit_sign = +1;
    mp_int pow2;
    mp_int tmp;

    if (!x || !str) return FAILURE;
    mp_free(x);
    mp_init(x);

    p = str;
    while (isspace((unsigned char)*p)) p++;

    if (*p == '-') { explicit_sign = -1; p++; }
    else if (*p == '+') { p++; }

    if (p[0] == '0' && (p[1] == 'b' || p[1] == 'B')) p += 2;

    digits = p;
    digits_len = 0;
    while (p[digits_len] == '0' || p[digits_len] == '1') digits_len++;

    if (digits_len == 0) return FAILURE;

    if (mp_reserve(x, 1) != SUCCESS) return FAILURE;
    x->digits[0] = 0;
    x->length = 1;
    x->sign = 1;

    for (i = 0; i < digits_len; ++i) {
        bit = (unsigned int)(digits[i] - '0');
        if (mp_mul_small(x, 2U) != SUCCESS) return FAILURE;
        if (mp_add_small(x, bit) != SUCCESS) return FAILURE;
    }

    p = digits + digits_len;
    while (isspace((unsigned char)*p)) p++;
    if (*p != '\0') return FAILURE;

    if (digits[0] == '1') {
        mp_init(&pow2);
        mp_init(&tmp);

        if (mp_reserve(&pow2, 1) != SUCCESS) { mp_free(&pow2); mp_free(&tmp); return FAILURE; }
        pow2.digits[0] = 1;
        pow2.length = 1;
        pow2.sign = 1;
        for (i = 0; i < digits_len; ++i) {
            if (mp_mul_small(&pow2, 2U) != SUCCESS) { mp_free(&pow2); mp_free(&tmp); return FAILURE; }
        }

        if (mp_sub(&tmp, x, &pow2) != SUCCESS) { mp_free(&pow2); mp_free(&tmp); return FAILURE; }

        mp_free(x);
        mp_init(x);
        if (mp_copy(x, &tmp) != SUCCESS) { mp_free(&pow2); mp_free(&tmp); return FAILURE; }

        mp_free(&pow2);
        mp_free(&tmp);
    }

    if (explicit_sign == -1) x->sign = -x->sign;

    /* Normalize zero sign */
    if (x->length == 1 && x->digits[0] == 0) x->sign = 0;

    return SUCCESS;
}

int mp_from_str_hex(mp_int *x, const char *str) {
    size_t i;
    const char *p;
    const char *digits;
    size_t digits_len;
    int val;
    int explicit_sign = +1;
    int first_nibble;

    if (!x || !str) return FAILURE;
    mp_free(x);
    mp_init(x);

    p = str;
    while (isspace((unsigned char)*p)) p++;

    if (*p == '-') { explicit_sign = -1; p++; }
    else if (*p == '+') { p++; }

    if (p[0] == '0' && (p[1] == 'x' || p[1] == 'X')) p += 2;

    digits = p;
    digits_len = 0;
    while (isxdigit((unsigned char)p[digits_len])) digits_len++;

    if (digits_len == 0) return FAILURE;

    p = digits + digits_len;
    while (isspace((unsigned char)*p)) p++;
    if (*p != '\0') return FAILURE;

    if (mp_reserve(x, 1) != SUCCESS) return FAILURE;
    x->digits[0] = 0;
    x->length = 1;
    x->sign = 1;

    for (i = 0; i < digits_len; ++i) {
        char c = digits[i];
        if (isdigit((unsigned char)c)) val = c - '0';
        else if (isupper((unsigned char)c)) val = c - 'A' + 10;
        else val = c - 'a' + 10;

        if (mp_mul_small(x, 16U) != SUCCESS) return FAILURE;
        if (mp_add_small(x, (unsigned int)val) != SUCCESS) return FAILURE;
    }

    first_nibble = 0;
    if (isxdigit((unsigned char)digits[0])) {
        char c = digits[0];
        if (isdigit((unsigned char)c)) first_nibble = c - '0';
        else if (isupper((unsigned char)c)) first_nibble = c - 'A' + 10;
        else first_nibble = c - 'a' + 10;
    }

    if (first_nibble & 0x8) {
        mp_int pow2;
        mp_int tmp;
        mp_init(&pow2);
        mp_init(&tmp);

        if (mp_reserve(&pow2, 1) != SUCCESS) { mp_free(&pow2); mp_free(&tmp); return FAILURE; }
        pow2.digits[0] = 1;
        pow2.length = 1;
        pow2.sign = 1;

        for (i = 0; i < (4 * digits_len); ++i) {
            if (mp_mul_small(&pow2, 2U) != SUCCESS) { mp_free(&pow2); mp_free(&tmp); return FAILURE; }
        }

        if (mp_sub(&tmp, x, &pow2) != SUCCESS) { mp_free(&pow2); mp_free(&tmp); return FAILURE; }

        mp_free(x);
        mp_init(x);
        if (mp_copy(x, &tmp) != SUCCESS) { mp_free(&pow2); mp_free(&tmp); return FAILURE; }

        mp_free(&pow2);
        mp_free(&tmp);
    }

    if (explicit_sign == -1) x->sign = -x->sign;

    /* Normalize zero sign */
    if (x->length == 1 && x->digits[0] == 0) x->sign = 0;

    return SUCCESS;
}