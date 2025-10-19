#include <stdio.h>
#include <stdlib.h>
#include "mp_int.h"
#include "exp.h"

int mp_pow(mp_int *r, const mp_int *a, const mp_int *b)
{
    mp_int base, exp, res, tmp, tmpb, qq, q;
    unsigned int is_exp_odd, rem;

    if (!r || !a || !b) return FAILURE;

    /* Negative exponent not supported */
    if (b->sign < 0) {
        fprintf(stderr, "error: negative exponent\n");
        return FAILURE;
    }

    /* Handle 0^0 and 0^positive cases */
    if (a->sign == 0) {
        if (b->sign == 0) {
            /* define 0^0 = 1 */
            mp_free(r);
            mp_init(r);
            mp_reserve(r, 1);
            r->digits[0] = 1;
            r->length = 1;
            r->sign = 1;
            return SUCCESS;
        } else {
            /* 0^n = 0 */
            mp_free(r);
            mp_init(r);
            r->sign = 0;
            r->length = 0;
            return SUCCESS;
        }
    }

    /* Handle exponent == 0 => 1 */
    if (b->sign == 0) {
        mp_free(r);
        mp_init(r);
        mp_reserve(r, 1);
        r->digits[0] = 1;
        r->length = 1;
        r->sign = 1;
        return SUCCESS;
    }

    mp_init(&base);
    mp_init(&exp);
    mp_init(&res);
    mp_init(&tmp);

    mp_copy(&base, a);
    mp_copy(&exp, b);

    /* res = 1 */
    mp_reserve(&res, 1);
    res.digits[0] = 1;
    res.length = 1;
    res.sign = 1;

    /* main loop: while exp > 0 */
    while (exp.sign != 0) {
        rem = 0;
        mp_init(&q);

        /* q = exp / 2, rem = exp % 2 */
        if (mp_div_small(&q, &exp, 2U, &rem) != SUCCESS) {
            mp_free(&q);
            mp_free(&base); mp_free(&exp); mp_free(&res); mp_free(&tmp);
            return FAILURE;
        }

        /* if exponent is odd, multiply result by base */
        if (rem == 1U) {
            mp_init(&tmp);
            if (mp_mul(&tmp, &res, &base) != SUCCESS) {
                mp_free(&tmp); mp_free(&q);
                mp_free(&base); mp_free(&exp); mp_free(&res);
                return FAILURE;
            }
            mp_free(&res);
            mp_init(&res);
            mp_copy(&res, &tmp);
            mp_free(&tmp);
        }

        /* square base */
        mp_init(&tmp);
        if (mp_mul(&tmp, &base, &base) != SUCCESS) {
            mp_free(&tmp); mp_free(&q);
            mp_free(&base); mp_free(&exp); mp_free(&res);
            return FAILURE;
        }
        mp_free(&base);
        mp_init(&base);
        mp_copy(&base, &tmp);
        mp_free(&tmp);

        /* move quotient back into exp */
        mp_free(&exp);
        mp_init(&exp);
        mp_copy(&exp, &q);
        mp_free(&q);
    }

    /* Final sign: only relevant if base was negative and exponent is odd */
    is_exp_odd = 0;
    {
        mp_init(&tmpb);
        mp_copy(&tmpb, b);
        rem = 0;
        mp_init(&qq);
        mp_div_small(&qq, &tmpb, 2U, &rem);
        if (rem == 1) is_exp_odd = 1;
        mp_free(&qq);
        mp_free(&tmpb);
    }

    if (a->sign < 0 && is_exp_odd)
        res.sign = -1;

    /* move res into r */
    mp_free(r);
    mp_init(r);
    mp_copy(r, &res);

    mp_free(&base);
    mp_free(&exp);
    mp_free(&res);
    mp_free(&tmp);
    return SUCCESS;
}