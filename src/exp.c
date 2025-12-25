/*
 * exp.c -- integer exponentiation for mp_int
 *
 * Implements binary exponentiation and handles negative exponents with
 * integer-division semantics (i.e., 1/(a^k) truncated to integer).
 *
 * ANSI C90 compatible, 32-bit safe.
 */

#include <stdio.h>
#include <stdlib.h>
#include "mp_int.h"
#include "exp.h"
#include "error.h"

/*
 * Compute r = a ^ b using integer exponentiation semantics.
 *
 * Implementation notes:
 *  - Works on temporaries to avoid altering the inputs.
 *  - Uses mp_div_small(..., 2, &rem) to obtain exponent parity each iteration.
 *  - For negative exponent, only integer reciprocal results are returned:
 *      if |a^|b|| == 1 => result is ±1, else 0
 *
 * Returns SUCCESS or FAILURE.
 */
int mp_pow(mp_int *r, const mp_int *a, const mp_int *b)
{
    mp_int base;
    mp_int exp;
    mp_int res;
    mp_int tmp;
    mp_int tmpb;
    mp_int q;
    unsigned int rem;
    unsigned int is_exp_odd;
    int negative_exp = 0;
    int status = FAILURE;

    if (!r || !a || !b) return FAILURE;

    /* Special cases for base == 0 */
    if (a->sign == 0) {
        if (b->sign == 0) {
            /* define 0^0 = 1 */
            mp_free(r);
            mp_init(r);
            if (mp_reserve(r, 1) != SUCCESS) return FAILURE;
            r->digits[0] = 1;
            r->length = 1;
            r->sign = 1;
            return SUCCESS;
        }
        /* 0^negative -> division by zero */
        if (b->sign < 0) {
            calc_error_set();
            printf("Division by zero!\n");
            return FAILURE;
        }
        /* 0^positive = 0 */
        mp_free(r);
        mp_init(r);
        r->sign = 0;
        r->length = 0;
        return SUCCESS;
    }

    /* exponent == 0 -> 1 */
    if (b->sign == 0) {
        mp_free(r);
        mp_init(r);
        if (mp_reserve(r, 1) != SUCCESS) return FAILURE;
        r->digits[0] = 1;
        r->length = 1;
        r->sign = 1;
        return SUCCESS;
    }

    /* prepare temporaries */
    mp_init(&base);
    mp_init(&exp);
    mp_init(&res);
    mp_init(&tmp);
    mp_init(&tmpb);
    mp_init(&q);

    /* make local copies; if mp_copy fails we'll go to cleanup */
    if (mp_copy(&base, a) != SUCCESS) goto cleanup;
    if (mp_copy(&exp, b) != SUCCESS) goto cleanup;

    /* remember if exponent was negative, work with |exp| inside loop */
    if (exp.sign < 0) {
        negative_exp = 1;
        exp.sign = +1;
    }

    /* res = 1 */
    if (mp_reserve(&res, 1) != SUCCESS) { status = FAILURE; goto cleanup; }
    res.digits[0] = 1;
    res.length = 1;
    res.sign = 1;

    /* Binary exponentiation loop: while exp > 0 */
    while (exp.sign != 0) {
        /* q = exp / 2, rem = exp % 2 */
        rem = 0;
        mp_free(&q);
        mp_init(&q);
        if (mp_div_small(&q, &exp, 2U, &rem) != SUCCESS) {
            status = FAILURE; goto cleanup;
        }

        /* if current exponent bit is 1: res = res * base */
        if (rem == 1U) {
            mp_free(&tmp);
            mp_init(&tmp);
            if (mp_mul(&tmp, &res, &base) != SUCCESS) { status = FAILURE; goto cleanup; }
            mp_free(&res);
            mp_init(&res);
            if (mp_copy(&res, &tmp) != SUCCESS) { status = FAILURE; goto cleanup; }
            mp_free(&tmp);
        }

        /* base = base * base */
        mp_free(&tmp);
        mp_init(&tmp);
        if (mp_mul(&tmp, &base, &base) != SUCCESS) { status = FAILURE; goto cleanup; }
        mp_free(&base);
        mp_init(&base);
        if (mp_copy(&base, &tmp) != SUCCESS) { status = FAILURE; goto cleanup; }
        mp_free(&tmp);

        /* exp = q (floor division) */
        mp_free(&exp);
        mp_init(&exp);
        if (mp_copy(&exp, &q) != SUCCESS) { status = FAILURE; goto cleanup; }
    }

    /* Determine parity of original exponent (for sign when base negative) */
    is_exp_odd = 0;
    mp_free(&tmpb);
    mp_init(&tmpb);
    if (mp_copy(&tmpb, b) != SUCCESS) { status = FAILURE; goto cleanup; }
    rem = 0;
    mp_free(&q);
    mp_init(&q);
    if (mp_div_small(&q, &tmpb, 2U, &rem) != SUCCESS) { status = FAILURE; goto cleanup; }
    if (rem == 1U) is_exp_odd = 1U;

    /* If base negative and exponent odd, ensure negative sign on result */
    if (a->sign < 0 && is_exp_odd) {
        res.sign = -1;
    }

    /* If exponent was negative, return integer reciprocal semantics:
       - if |res| == 1 => result is ±1
       - else => integer division gives 0
    */
    if (negative_exp) {
        if (res.length == 0 || res.sign == 0) {
            fprintf(stderr, "error: division by zero in negative exponent\n");
            status = FAILURE;
            goto cleanup;
        }

        if (res.length == 1 && res.digits[0] == 1U) {
            /* result is ±1 */
            mp_free(r);
            mp_init(r);
            if (mp_reserve(r, 1) != SUCCESS) { status = FAILURE; goto cleanup; }
            r->digits[0] = 1U;
            r->length = 1;
            r->sign = (res.sign < 0) ? -1 : +1;
            status = SUCCESS;
            goto cleanup;
        } else {
            /* integer reciprocal truncates to zero */
            mp_free(r);
            mp_init(r);
            r->sign = 0;
            r->length = 0;
            status = SUCCESS;
            goto cleanup;
        }
    }

    /* positive exponent: copy res to r */
    mp_free(r);
    mp_init(r);
    if (mp_copy(r, &res) != SUCCESS) { status = FAILURE; goto cleanup; }
    status = SUCCESS;

cleanup:
    /* free temporaries */
    mp_free(&base);
    mp_free(&exp);
    mp_free(&res);
    mp_free(&tmp);
    mp_free(&tmpb);
    mp_free(&q);
    return status;
}
