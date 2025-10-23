#include <stdio.h>
#include <stdlib.h>
#include "mp_int.h"
#include "exp.h"

int mp_pow(mp_int *r, const mp_int *a, const mp_int *b)
{
    /* declarations first (ANSI C90) */
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

    /* Special-cases: base == 0 */
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
        /* if exponent negative -> 1/0 -> error */
        if (b->sign < 0) {
            fprintf(stderr, "error: division by zero (0 raised to negative exponent)\n");
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

    mp_copy(&base, a);
    mp_copy(&exp, b);

    /* remember if exponent was negative and work with |exp| for the loop */
    if (exp.sign < 0) {
        negative_exp = 1;
        exp.sign = +1;
    }

    /* res = 1 */
    if (mp_reserve(&res, 1) != SUCCESS) { status = FAILURE; goto cleanup; }
    res.digits[0] = 1;
    res.length = 1;
    res.sign = 1;

    /* main loop: binary exponentiation while exp > 0 */
    while (exp.sign != 0) {
        rem = 0;
        mp_free(&q);
        mp_init(&q);
        if (mp_div_small(&q, &exp, 2U, &rem) != SUCCESS) {
            status = FAILURE; goto cleanup;
        }

        /* if exponent bit is 1: res = res * base */
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

        /* exp = q (the quotient) */
        mp_free(&exp);
        mp_init(&exp);
        if (mp_copy(&exp, &q) != SUCCESS) { status = FAILURE; goto cleanup; }
    }

    /* Determine parity of original exponent (needed only for sign when base negative) */
    is_exp_odd = 0;
    mp_free(&tmpb);
    mp_init(&tmpb);
    if (mp_copy(&tmpb, b) != SUCCESS) { status = FAILURE; goto cleanup; }
    rem = 0;
    mp_free(&q);
    mp_init(&q);
    if (mp_div_small(&q, &tmpb, 2U, &rem) != SUCCESS) { status = FAILURE; goto cleanup; }
    if (rem == 1U) is_exp_odd = 1U;

    /* mp_mul already set the sign of 'res' correctly in multiplications.
       However in some codepaths it's safe to ensure sign for odd negative-base exponents: */
    if (a->sign < 0 && is_exp_odd) {
        /* ensure result sign is negative */
        res.sign = -1;
    }

    /* Handle negative exponent by integer division: r = 1 / res (integer division) */
    if (negative_exp) {
        /* division by zero check */
        if (res.length == 0 || res.sign == 0) {
            fprintf(stderr, "error: division by zero in negative exponent\n");
            status = FAILURE;
            goto cleanup;
        }

        /* if |res| == 1, 1 / res is +/-1; else integer division gives 0 */
        if (res.length == 1 && res.digits[0] == 1U) {
            mp_free(r);
            mp_init(r);
            if (mp_reserve(r, 1) != SUCCESS) { status = FAILURE; goto cleanup; }
            r->digits[0] = 1U;
            r->length = 1;
            r->sign = (res.sign < 0) ? -1 : +1;
            status = SUCCESS;
            goto cleanup;
        } else {
            /* 1 / |res| where |res| > 1 -> integer result 0 (C truncates toward zero) */
            mp_free(r);
            mp_init(r);
            r->sign = 0;
            r->length = 0;
            status = SUCCESS;
            goto cleanup;
        }
    }

    /* positive exponent: move res into r */
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