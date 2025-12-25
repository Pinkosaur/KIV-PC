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
            /* 0^0 is defined as 1 here */
            mp_free(r);
            mp_init(r);
            if (mp_reserve(r, 1) != SUCCESS) return FAILURE;
            r->digits[0] = 1;
            r->length = 1;
            r->sign = 1;
            return SUCCESS;
        }
        if (b->sign < 0) {
            calc_error_set();
            printf("Division by zero!\n");
            return FAILURE;
        }
        mp_free(r);
        mp_init(r);
        r->sign = 0;
        r->length = 0;
        return SUCCESS;
    }

    /* Optimization: If base is 1 or -1, result is always 1 or -1 */
    if (a->length == 1 && a->digits[0] == 1) {
        mp_free(r);
        mp_init(r);
        if (mp_reserve(r, 1) != SUCCESS) return FAILURE;
        r->digits[0] = 1;
        r->length = 1;
        
        if (a->sign > 0) {
            r->sign = 1;
        } else {
            /* -1 ^ exp: check if exp is odd */
            mp_init(&tmpb); mp_init(&q);
            if (mp_copy(&tmpb, b) != SUCCESS) { mp_free(&tmpb); mp_free(&q); return FAILURE; }
            if (mp_div_small(&q, &tmpb, 2U, &rem) != SUCCESS) { mp_free(&tmpb); mp_free(&q); return FAILURE; }
            r->sign = (rem == 1) ? -1 : 1;
            mp_free(&tmpb); mp_free(&q);
        }
        return SUCCESS;
    }

    mp_init(&base);
    mp_init(&exp);
    mp_init(&res);
    mp_init(&tmp);
    mp_init(&tmpb);
    mp_init(&q);

    if (mp_copy(&base, a) != SUCCESS) goto cleanup;
    if (mp_copy(&exp, b) != SUCCESS) goto cleanup;

    if (exp.sign < 0) {
        negative_exp = 1;
        exp.sign = +1;
    }

    if (mp_reserve(&res, 1) != SUCCESS) { status = FAILURE; goto cleanup; }
    res.digits[0] = 1;
    res.length = 1;
    res.sign = 1;

    /* Binary exponentiation loop */
    while (exp.sign != 0) {
        rem = 0;
        /* q = exp / 2, rem = exp % 2 */
        if (mp_div_small(&q, &exp, 2U, &rem) != SUCCESS) {
            status = FAILURE; goto cleanup;
        }

        if (rem == 1U) {
            /* res = res * base */
            if (mp_mul(&tmp, &res, &base) != SUCCESS) { status = FAILURE; goto cleanup; }
            if (mp_copy(&res, &tmp) != SUCCESS) { status = FAILURE; goto cleanup; }
        }

        /* base = base * base (if we need to continue) */
        if (q.sign != 0) {
            if (mp_mul(&tmp, &base, &base) != SUCCESS) { status = FAILURE; goto cleanup; }
            if (mp_copy(&base, &tmp) != SUCCESS) { status = FAILURE; goto cleanup; }
        }

        if (mp_copy(&exp, &q) != SUCCESS) { status = FAILURE; goto cleanup; }
    }

    /* Handle sign for base < 0 */
    is_exp_odd = 0;
    if (mp_copy(&tmpb, b) != SUCCESS) { status = FAILURE; goto cleanup; }
    if (mp_div_small(&q, &tmpb, 2U, &rem) != SUCCESS) { status = FAILURE; goto cleanup; }
    if (rem == 1U) is_exp_odd = 1U;

    if (a->sign < 0 && is_exp_odd) {
        res.sign = -1;
    }

    if (negative_exp) {
        /* Reciprocal logic: integer division 1 / result */
        if (res.length == 0 || res.sign == 0) {
            calc_error_set();
            printf("Division by zero!\n");
            status = FAILURE;
            goto cleanup;
        }
        if (res.length == 1 && res.digits[0] == 1U) {
            mp_free(r); mp_init(r);
            if (mp_reserve(r, 1) != SUCCESS) { status = FAILURE; goto cleanup; }
            r->digits[0] = 1U; r->length = 1;
            r->sign = (res.sign < 0) ? -1 : +1;
        } else {
            mp_free(r); mp_init(r);
            r->sign = 0; r->length = 0;
        }
        status = SUCCESS;
    } else {
        mp_free(r);
        mp_init(r);
        if (mp_copy(r, &res) != SUCCESS) { status = FAILURE; goto cleanup; }
        status = SUCCESS;
    }

cleanup:
    mp_free(&base); mp_free(&exp); mp_free(&res);
    mp_free(&tmp); mp_free(&tmpb); mp_free(&q);
    return status;
}