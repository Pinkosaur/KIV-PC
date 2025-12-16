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

int mp_print_dec(mp_int *x)
{
    mp_int tmp, q;
    unsigned int rem;
    char **chunks;
    size_t num_chunks, chunk_cap;
    const unsigned int BASE = 1000000000U; /* 1e9 */
    size_t i;

    if (!x) return FAILURE;
    if (x->sign == 0) {
        printf("0");
        return SUCCESS;
    }

    mp_init(&tmp);
    mp_init(&q);
    if (mp_copy(&tmp, x) != SUCCESS) goto cleanup;
    tmp.sign = +1;

    /* allocate initial chunk array */
    chunk_cap = 64;
    chunks = (char**)malloc(chunk_cap * sizeof(char*));
    if (!chunks) goto cleanup;
    num_chunks = 0;

    /* repeatedly divide by 1e9 instead of 10 */
    while (tmp.sign != 0) {
        if (num_chunks >= chunk_cap) {
            char **newchunks;
            chunk_cap *= 2;
            newchunks = (char**)realloc(chunks, chunk_cap * sizeof(char*));
            if (!newchunks) { free(chunks); goto cleanup; }
            chunks = newchunks;
        }

        /* q = tmp / BASE, rem = tmp % BASE */
        if (mp_div_small(&q, &tmp, BASE, &rem) != SUCCESS) {
            free(chunks);
            goto cleanup;
        }

        /* store remainder as string chunk */
        chunks[num_chunks] = (char*)malloc(10);
        if (!chunks[num_chunks]) { free(chunks); goto cleanup; }
        sprintf(chunks[num_chunks], "%09u", rem); /* pad with zeros */
        num_chunks++;

        mp_free(&tmp);
        mp_init(&tmp);
        mp_copy(&tmp, &q);
    }

    /* print sign */
    if (x->sign < 0) putchar('-');

    /* print most significant chunk without padding zeros */
    if (num_chunks > 0) {
        unsigned long ms_val = strtoul(chunks[num_chunks - 1], NULL, 10);
        printf("%lu", ms_val);
    }

    /* print remaining chunks padded */
    for (i = num_chunks - 1; i > 0; --i)
        printf("%s", chunks[i - 1]);

    /* cleanup */
    for (i = 0; i < num_chunks; ++i)
        free(chunks[i]);
    free(chunks);
    mp_free(&tmp);
    mp_free(&q);
    return SUCCESS;

cleanup:
    mp_free(&tmp);
    mp_free(&q);
    return FAILURE;
}


int mp_print_bin(mp_int *x)
{
    /* New implementation prints a single sign bit plus minimal magnitude bits.
       For positive x: prints w = bitlen(x) + 1 bits, leading bit = 0.
       For negative x: prints w = bitlen(|x|) + 1 bits as two's-complement (leading bit = 1).
    */

    size_t bitlen, w, bit_index;
    unsigned int bit;
    mp_int absx, pow2, val;
    size_t top_idx;
    unsigned int top;
    unsigned int tb;

    if (!x) return FAILURE;
    if (x->sign == 0) {
        printf("0b0");
        return SUCCESS;
    }

    printf("0b");

    /* Helper: compute bit length of absolute value (number of significant bits)
       zero -> 0, otherwise highest_bit_index + 1. */
    if (x->sign > 0) {
        /* positive: compute bitlen directly from x */
        if (x->length == 0) {
            bitlen = 0;
        } else {
            top_idx = x->length - 1;
            top = x->digits[top_idx];
            tb = 0;
            while (top) { top >>= 1; ++tb; }
            bitlen = ((size_t)top_idx) * 32 + (size_t)tb;
        }

        /* w = bitlen + 1 (one sign bit '0') */
        w = bitlen + 1;

        /* Print exactly w bits, MSB first */
        for (bit_index = w; bit_index-- > 0; ) {
            size_t limb = bit_index / 32;
            unsigned int pos = (unsigned int)(bit_index % 32);
            unsigned int limbval = (limb < x->length) ? x->digits[limb] : 0U;
            bit = (limbval >> pos) & 1U;
            putchar('0' + (int)bit);
        }

    } else {
        /* negative: compute bitlen of absolute value and print w = bitlen+1 bits
           of (2^w + x) which yields the two's complement representation. */
        mp_init(&absx); mp_init(&pow2); mp_init(&val);

        mp_abs_copy(&absx, x);

        if (absx.length == 0) {
            bitlen = 0;
        } else {
            top_idx = absx.length - 1;
            top = absx.digits[top_idx];
            tb = 0;
            while (top) { top >>= 1; ++tb; }
            bitlen = ((size_t)top_idx) * 32 + (size_t)tb;
        }

        w = bitlen + 1;

        /* pow2 = 2^w */
        if (mp_make_pow2(&pow2, w) != SUCCESS) {
            mp_free(&absx);
            mp_free(&pow2);
            mp_free(&val);
            return FAILURE;
        }

        /* val = pow2 + x  (x is negative, so this yields two's-complement magnitude) */
        if (mp_add(&val, &pow2, x) != SUCCESS) {
            mp_free(&absx);
            mp_free(&pow2);
            mp_free(&val);
            return FAILURE;
        }

        /* print exactly w bits of val (MSB first) */
        for (bit_index = w; bit_index-- > 0; ) {
            size_t limb = bit_index / 32;
            unsigned int pos = (unsigned int)(bit_index % 32);
            unsigned int limbval = (limb < val.length) ? val.digits[limb] : 0U;
            bit = (limbval >> pos) & 1U;
            putchar('0' + (int)bit);
        }

        mp_free(&absx);
        mp_free(&pow2);
        mp_free(&val);
    }

    return SUCCESS;
}



int mp_print_hex(mp_int *x)
{
    static const char hex_arr[17] = "0123456789abcdef";
    size_t bitlen;
    size_t mag_nibbles;   /* nibbles required to represent magnitude */
    size_t total_nibbles; /* final number of nibbles to print */
    size_t nib_index;
    mp_int absx, powb, mpval;
    size_t top_idx;
    unsigned int top;
    unsigned int tb;

    if (!x) return FAILURE;
    if (x->sign == 0) {
        printf("0x0");
        return SUCCESS;
    }

    printf("0x");

    /* compute bit length of absolute value (number of significant bits) */
    if (x->sign > 0) {
        if (x->length == 0) {
            bitlen = 0;
        } else {
            top_idx = x->length - 1;
            top = x->digits[top_idx];
            tb = 0;
            while (top) { top >>= 1; ++tb; }
            bitlen = ((size_t)top_idx) * 32 + (size_t)tb;
        }

        /* minimal nibbles to represent magnitude (ceil(bitlen / 4)) */
        if (bitlen == 0) mag_nibbles = 1;
        else mag_nibbles = (bitlen + 3) / 4;

        /* If most-significant nibble would have its high bit set (>=8),
           add one leading 0 nibble so sign nibble = 0. */
        {
            /* find value of the current top nibble (without extra leading nibble) */
            size_t top_nib_index = mag_nibbles - 1;
            size_t limb = top_nib_index / 8;        /* 8 nibbles per 32-bit limb */
            unsigned int nib_shift = (unsigned int)((top_nib_index % 8) * 4);
            unsigned int limbval = (limb < x->length) ? x->digits[limb] : 0U;
            unsigned int top_nibble = (limbval >> nib_shift) & 0xF;
            if (top_nibble >= 8U) {
                total_nibbles = mag_nibbles + 1; /* add a leading 0 nibble */
            } else {
                total_nibbles = mag_nibbles;
            }
        }

        /* Print exactly total_nibbles hex digits (MSB first).
           We index from most-significant nibble (total_nibbles-1) down to 0. */
        for (nib_index = total_nibbles; nib_index-- > 0; ) {
            size_t limb = nib_index / 8;
            unsigned int nib_shift = (unsigned int)((nib_index % 8) * 4);
            unsigned int limbval = (limb < x->length) ? x->digits[limb] : 0U;
            unsigned int digit = (limbval >> nib_shift) & 0xF;
            putchar(hex_arr[digit]);
        }
    } else {
        /* Negative: compute minimal two's-complement nibble width and print exactly that many nibbles:
           wbits = bitlen(|x|) + 1  => nibbles = ceil(wbits / 4)
           mpval = 2^(4*nibbles) + x  (x negative)  => print nibbles hex digits of mpval
        */
        mp_init(&absx);
        mp_init(&powb);
        mp_init(&mpval);

        mp_abs_copy(&absx, x);

        if (absx.length == 0) {
            bitlen = 0;
        } else {
            top_idx = absx.length - 1;
            top = absx.digits[top_idx];
            tb = 0;
            while (top) { top >>= 1; ++tb; }
            bitlen = ((size_t)top_idx) * 32 + (size_t)tb;
        }

        /* wbits = bitlen + 1  => minimal nibbles to hold sign bit */
        {
            size_t wbits = bitlen + 1;
            size_t nibbles = (wbits + 3) / 4; /* ceil(wbits/4) */
            if (nibbles == 0) nibbles = 1;
            total_nibbles = nibbles;
        }

        /* powb = 2^(4 * total_nibbles) */
        if (mp_make_pow2(&powb, total_nibbles * 4) != SUCCESS) {
            mp_free(&absx);
            mp_free(&powb);
            mp_free(&mpval);
            return FAILURE;
        }

        /* mpval = powb + x  (note x is negative) */
        if (mp_add(&mpval, &powb, x) != SUCCESS) {
            mp_free(&absx);
            mp_free(&powb);
            mp_free(&mpval);
            return FAILURE;
        }

        /* Print exactly total_nibbles hex digits from mpval (MSB first) */
        for (nib_index = total_nibbles; nib_index-- > 0; ) {
            size_t limb = nib_index / 8;
            unsigned int nib_shift = (unsigned int)((nib_index % 8) * 4);
            unsigned int limbval = (limb < mpval.length) ? mpval.digits[limb] : 0U;
            unsigned int digit = (limbval >> nib_shift) & 0xF;
            putchar(hex_arr[digit]);
        }

        mp_free(&absx);
        mp_free(&powb);
        mp_free(&mpval);
    }

    return SUCCESS;
}
