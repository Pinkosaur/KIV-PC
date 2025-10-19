#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mp_int.h"
#include "mp_print.h"

/* Helper: copy absolute value of x into dst (dst initialized by caller) */
static int mp_abs_copy(mp_int *dst, const mp_int *x) {
    if (!dst || !x) return FAILURE;
    if (mp_copy(dst, x) != SUCCESS) return FAILURE;
    if (dst->length == 0) { dst->sign = 0; return SUCCESS; }
    dst->sign = +1;
    return SUCCESS;
}

/* Helper: create pow2 = 2^k */
static int mp_make_pow2(mp_int *pow2, size_t k) {
    if (!pow2) return FAILURE;
    mp_free(pow2);
    mp_init(pow2);
    if (mp_reserve(pow2, 1) != SUCCESS) return FAILURE;
    pow2->digits[0] = 1;
    pow2->length = 1;
    pow2->sign = 1;
    /* multiply by 2 k times */
    while (k--) {
        if (mp_mul_small(pow2, 2U) != SUCCESS) return FAILURE;
    }
    return SUCCESS;
}

/* Prints the mp_int @x in decimal format (unchanged behavior) */
int mp_print_dec(mp_int *x) {
    mp_int tmp, q;
    size_t approx;
    char *buf;
    size_t idx;
    unsigned int rem;

    if (!x) return FAILURE;
    if (x->sign == 0) {
        printf("0");
        return SUCCESS;
    }

    mp_init(&tmp); mp_init(&q);
    if (mp_copy(&tmp, x) != SUCCESS) { mp_free(&tmp); mp_free(&q); return FAILURE; }
    tmp.sign = +1;

    approx = tmp.length * 10 + 4;
    buf = (char*) malloc(approx);
    if (!buf) { mp_free(&tmp); mp_free(&q); return FAILURE; }
    idx = 0;

    while (tmp.sign != 0) {
        mp_init(&q);
        if (mp_div_small(&q, &tmp, 10U, &rem) != SUCCESS) {
            free(buf); mp_free(&tmp); mp_free(&q);
            return FAILURE;
        }
        buf[idx++] = (char)('0' + (rem % 10));
        mp_free(&tmp);
        mp_init(&tmp);
        mp_copy(&tmp, &q);
        mp_free(&q);
    }

    if (x->sign < 0) putchar('-');
    if (idx == 0) putchar('0');
    else while (idx > 0) putchar(buf[--idx]);

    free(buf);
    mp_free(&tmp);
    return SUCCESS;
}

/* Prints the mp_int @x in binary (two's complement if negative), minimal width */
int mp_print_bin(mp_int *x) {
    mp_int tmp, q;
    size_t approx;
    char *buf;
    size_t idx;
    unsigned int rem;

    if (!x) return FAILURE;
    if (x->sign == 0) { printf("0b0"); return SUCCESS; }

    if (x->sign > 0) {
        /* Positive: print minimal binary (no leading zeros) */
        mp_init(&tmp); mp_init(&q);
        mp_copy(&tmp, x);
        tmp.sign = +1;

        approx = tmp.length * 32 + 4;
        buf = (char*) malloc(approx);
        if (!buf) { mp_free(&tmp); mp_free(&q); return FAILURE; }
        idx = 0;

        while (tmp.sign != 0) {
            mp_init(&q);
            if (mp_div_small(&q, &tmp, 2U, &rem) != SUCCESS) {
                free(buf); mp_free(&tmp); mp_free(&q);
                return FAILURE;
            }
            buf[idx++] = (char)('0' + (rem & 1));
            mp_free(&tmp);
            mp_init(&tmp);
            mp_copy(&tmp, &q);
            mp_free(&q);
        }

        printf("0b");
        if (idx == 0) putchar('0');
        else while (idx > 0) putchar(buf[--idx]);

        free(buf);
        mp_free(&tmp);
        return SUCCESS;
    } else {
        /* Negative: choose minimal w such that x >= -2^(w-1) (i.e. fits in w-bit two's complement).
           Then compute val = 2^w + x and print exactly w bits. */
        mp_int absx;
        size_t w;
        mp_int pow;
        int cmp;
        mp_int poww, val;
        mp_int q2;
        mp_init(&absx);
        if (mp_abs_copy(&absx, x) != SUCCESS) { mp_free(&absx); return FAILURE; }

        /* find minimal w (start from 1) */
        w = 1;
        for (;;) {
            mp_init(&pow);
            if (mp_make_pow2(&pow, w - 1) != SUCCESS) { mp_free(&absx); mp_free(&pow); return FAILURE; }
            /* condition: abs(x) <= 2^(w-1) */
            cmp = mp_cmp_abs(&absx, &pow);
            mp_free(&pow);
            if (cmp <= 0) break;
            w++;
            /* safety: cap w to a reasonable max based on absx length */
            if (w > absx.length * 32 + 64) break;
        }

        /* compute val = 2^w + x */
        mp_init(&poww); mp_init(&val);
        if (mp_make_pow2(&poww, w) != SUCCESS) { mp_free(&absx); mp_free(&poww); mp_free(&val); return FAILURE; }
        /* val = poww + x  (x is negative, addition yields the two's complement representation) */
        if (mp_add(&val, &poww, x) != SUCCESS) { mp_free(&absx); mp_free(&poww); mp_free(&val); return FAILURE; }

        /* extract bits of val */
        approx = w + 4;
        buf = (char*) malloc(approx);
        if (!buf) { mp_free(&absx); mp_free(&poww); mp_free(&val); return FAILURE; }
        idx = 0;
        mp_init(&tmp);
        mp_copy(&tmp, &val);
        tmp.sign = +1;
        while (tmp.sign != 0) {
            mp_init(&q2);
            if (mp_div_small(&q2, &tmp, 2U, &rem) != SUCCESS) {
                free(buf); mp_free(&tmp); mp_free(&q2); mp_free(&absx); mp_free(&poww); mp_free(&val);
                return FAILURE;
            }
            buf[idx++] = (char)('0' + (rem & 1));
            mp_free(&tmp);
            mp_init(&tmp);
            mp_copy(&tmp, &q2);
            mp_free(&q2);
        }
        /* pad to exactly w bits */
        while (idx < w) buf[idx++] = '1'; /* two's complement negative should pad with ones */

        printf("0b");
        while (idx > 0) putchar(buf[--idx]);

        free(buf);
        mp_free(&tmp);
        mp_free(&poww);
        mp_free(&val);
        mp_free(&absx);
        return SUCCESS;
    }
}

/* Prints the mp_int @x in hexadecimal (two's complement if negative), minimal nibble width */
int mp_print_hex(mp_int *x) {
    mp_int tmp;
    size_t approx;
    char *buf;
    size_t idx;
    unsigned int rem;
    mp_int q2;

    if (!x) return FAILURE;
    if (x->sign == 0) { printf("0x0"); return SUCCESS; }

    if (x->sign > 0) {
        /* Positive: straightforward minimal hex */
        mp_init(&tmp);
        mp_copy(&tmp, x);
        tmp.sign = +1;

        approx = tmp.length * 8 + 4;
        buf = (char*) malloc(approx);
        if (!buf) { mp_free(&tmp); return FAILURE; }
        idx = 0;
        while (tmp.sign != 0) {
            mp_init(&q2);
            if (mp_div_small(&q2, &tmp, 16U, &rem) != SUCCESS) {
                free(buf); mp_free(&tmp); mp_free(&q2);
                return FAILURE;
            }
            buf[idx++] = (rem < 10) ? (char)('0' + rem) : (char)('a' + rem - 10);
            mp_free(&tmp);
            mp_init(&tmp);
            mp_copy(&tmp, &q2);
            mp_free(&q2);
        }

        printf("0x");
        if (idx == 0) putchar('0');
        else while (idx > 0) putchar(buf[--idx]);

        free(buf);
        mp_free(&tmp);
        return SUCCESS;
    } else {
        /* Negative: pick minimal number of nibbles (4-bit groups) so that x fits in that many bits */
        mp_int absx;
        size_t w;
        mp_int pow;
        int cmp;
        size_t nibbles, bits;
        mp_int powb, val;
        mp_init(&absx);
        if (mp_abs_copy(&absx, x) != SUCCESS) { mp_free(&absx); return FAILURE; }

        /* find minimal w bits as before */
        w = 1;
        for (;;) {
            mp_init(&pow);
            if (mp_make_pow2(&pow, w - 1) != SUCCESS) { mp_free(&absx); mp_free(&pow); return FAILURE; }
            cmp = mp_cmp_abs(&absx, &pow);
            mp_free(&pow);
            if (cmp <= 0) break;
            w++;
            if (w > absx.length * 32 + 64) break;
        }
        /* convert w to nibble count */
        nibbles = (w + 3) / 4;
        bits = nibbles * 4;

        /* val = 2^bits + x */
        mp_init(&powb); mp_init(&val);
        if (mp_make_pow2(&powb, bits) != SUCCESS) { mp_free(&absx); mp_free(&powb); mp_free(&val); return FAILURE; }
        if (mp_add(&val, &powb, x) != SUCCESS) { mp_free(&absx); mp_free(&powb); mp_free(&val); return FAILURE; }

        /* Extract hex digits of val and pad to nibbles */
        approx = nibbles + 4;
        buf = (char*) malloc(approx);
        if (!buf) { mp_free(&absx); mp_free(&powb); mp_free(&val); return FAILURE; }
        idx = 0;
        mp_init(&tmp);
        mp_copy(&tmp, &val);
        tmp.sign = +1;
        while (tmp.sign != 0) {
            mp_init(&q2);
            if (mp_div_small(&q2, &tmp, 16U, &rem) != SUCCESS) {
                free(buf); mp_free(&tmp); mp_free(&q2); mp_free(&absx); mp_free(&powb); mp_free(&val);
                return FAILURE;
            }
            buf[idx++] = (rem < 10) ? (char)('0' + rem) : (char)('a' + rem - 10);
            mp_free(&tmp);
            mp_init(&tmp);
            mp_copy(&tmp, &q2);
            mp_free(&q2);
        }
        while (idx < nibbles) buf[idx++] = 'f'; /* pad leading Fs for negative */
        printf("0x");
        while (idx > 0) putchar(buf[--idx]);

        free(buf);
        mp_free(&tmp);
        mp_free(&powb);
        mp_free(&val);
        mp_free(&absx);
        return SUCCESS;
    }
}
