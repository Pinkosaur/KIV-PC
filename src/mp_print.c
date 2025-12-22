/*
 * mp_print.c -- printing helpers for mp_int
 *
 * Uses mp_limb_t and MP_LIMB_BITS from mp_int.h.
 *
 * The printing functions compute bit / nibble widths from MP_LIMB_BITS
 * so they work correctly regardless of chosen limb size (16-bit or 32-bit).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mp_int.h"
#include "mp_print.h"

/* Local helper: copy absolute value of x into dst (dst must be initialized) */
static int mp_abs_copy(mp_int *dst, const mp_int *x)
{
    if (!dst || !x) return FAILURE;
    if (mp_copy(dst, x) != SUCCESS) return FAILURE;
    if (dst->length == 0) { dst->sign = 0; return SUCCESS; }
    dst->sign = +1;
    return SUCCESS;
}

/* Local helper: build pow2 = 2^k (k bits). */
static int mp_make_pow2(mp_int *pow2, size_t k)
{
    if (!pow2) return FAILURE;
    mp_free(pow2);
    mp_init(pow2);

    /* allocate at least one limb so mp_mul_small will work */
    if (mp_reserve(pow2, 1) != SUCCESS) return FAILURE;
    pow2->digits[0] = (mp_limb_t)1;
    pow2->length = 1;
    pow2->sign = 1;

    /* multiply by two k times */
    while (k--) {
        if (mp_mul_small(pow2, 2U) != SUCCESS) return FAILURE;
    }
    return SUCCESS;
}

/* ----------------------- Decimal printing --------------------------------
 * Strategy:
 *   Repeatedly divide the absolute value by BASE = 1e9 and store remainders
 *   as zero-padded chunks. Then print the most-significant chunk without
 *   padding and the rest padded. This reduces the number of division calls.
 */
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
    tmp.sign = +1; /* operate on magnitude */

    /* allocate initial chunk array */
    chunk_cap = 64;
    chunks = (char**) malloc(chunk_cap * sizeof(char*));
    if (!chunks) goto cleanup;
    num_chunks = 0;

    /* repeatedly divide by BASE storing remainder chunks */
    while (tmp.sign != 0) {
        if (num_chunks >= chunk_cap) {
            char **newchunks;
            chunk_cap *= 2;
            newchunks = (char**) realloc(chunks, chunk_cap * sizeof(char*));
            if (!newchunks) { free(chunks); goto cleanup; }
            chunks = newchunks;
        }

        if (mp_div_small(&q, &tmp, BASE, &rem) != SUCCESS) {
            free(chunks);
            goto cleanup;
        }

        chunks[num_chunks] = (char*) malloc(10); /* 9 digits + null */
        if (!chunks[num_chunks]) { free(chunks); goto cleanup; }
        sprintf(chunks[num_chunks], "%09u", rem);
        num_chunks++;

        mp_free(&tmp);
        mp_init(&tmp);
        if (mp_copy(&tmp, &q) != SUCCESS) { free(chunks); goto cleanup; }
    }

    /* print sign */
    if (x->sign < 0) putchar('-');

    /* print most-significant chunk without leading zeros */
    if (num_chunks > 0) {
        unsigned long ms_val = strtoul(chunks[num_chunks - 1], NULL, 10);
        printf("%lu", ms_val);
    }

    /* print remaining chunks padded */
    for (i = num_chunks - 1; i > 0; --i) {
        printf("%s", chunks[i - 1]);
    }

    /* regular cleanup */
    for (i = 0; i < num_chunks; ++i) free(chunks[i]);
    free(chunks);
    mp_free(&tmp);
    mp_free(&q);
    return SUCCESS;

    /* error cleanup */
cleanup:
    if(chunks) free(chunks);
    mp_free(&tmp);
    mp_free(&q);
    return FAILURE;
}

/* ----------------------- Binary printing ---------------------------------
 * For positive x: print '0b' followed by w = bitlen(x)+1 bits, top bit is 0
 *   (explicit sign bit). This ensures a positive number always starts with 0.
 * For negative x: compute minimal w = bitlen(|x|)+1 and print (2^w + x) in w bits.
 *
 * The implementation uses MP_LIMB_BITS (limb width) so it works for 16/32-bit limbs.
 */
int mp_print_bin(mp_int *x)
{
    size_t bitlen, w, bit_index;
    unsigned int bit;
    mp_int absx, pow2, val;
    size_t top_idx;
    mp_limb_t top;
    unsigned int tb;

    if (!x) return FAILURE;
    if (x->sign == 0) {
        printf("0b0");
        return SUCCESS;
    }

    printf("0b");

    /* compute bit length of absolute value */
    if (x->sign > 0) {
        if (x->length == 0) {
            bitlen = 0;
        } else {
            top_idx = x->length - 1;
            top = x->digits[top_idx];
            tb = 0;
            while (top) { top >>= 1; ++tb; }
            bitlen = ((size_t)top_idx) * (size_t)MP_LIMB_BITS + (size_t)tb;
        }

        /* add explicit sign bit */
        w = bitlen + 1;

        /* print exactly w bits MSB-first */
        for (bit_index = w; bit_index-- > 0; ) {
            size_t limb = bit_index / (size_t)MP_LIMB_BITS;
            unsigned int pos = (unsigned int)(bit_index % (size_t)MP_LIMB_BITS);
            mp_limb_t limbval = (limb < x->length) ? x->digits[limb] : (mp_limb_t)0;
            bit = (unsigned int)((limbval >> pos) & (mp_limb_t)1);
            putchar('0' + (int)bit);
        }

    } else {
        /* negative: compute bitlen(|x|), choose w = bitlen+1 and print two's complement */
        mp_init(&absx); mp_init(&pow2); mp_init(&val);

        mp_abs_copy(&absx, x);

        if (absx.length == 0) {
            bitlen = 0;
        } else {
            top_idx = absx.length - 1;
            top = absx.digits[top_idx];
            tb = 0;
            while (top) { top >>= 1; ++tb; }
            bitlen = ((size_t)top_idx) * (size_t)MP_LIMB_BITS + (size_t)tb;
        }

        w = bitlen + 1;

        /* pow2 = 2^w */
        if (mp_make_pow2(&pow2, w) != SUCCESS) {
            mp_free(&absx); mp_free(&pow2); mp_free(&val);
            return FAILURE;
        }

        /* val = pow2 + x  (x negative -> yields two's complement magnitude) */
        if (mp_add(&val, &pow2, x) != SUCCESS) {
            mp_free(&absx); mp_free(&pow2); mp_free(&val);
            return FAILURE;
        }

        /* print exactly w bits of val */
        for (bit_index = w; bit_index-- > 0; ) {
            size_t limb = bit_index / (size_t)MP_LIMB_BITS;
            unsigned int pos = (unsigned int)(bit_index % (size_t)MP_LIMB_BITS);
            mp_limb_t limbval = (limb < val.length) ? val.digits[limb] : (mp_limb_t)0;
            bit = (unsigned int)((limbval >> pos) & (mp_limb_t)1);
            putchar('0' + (int)bit);
        }

        mp_free(&absx); mp_free(&pow2); mp_free(&val);
    }

    return SUCCESS;
}

/* ----------------------- Hexadecimal printing -----------------------------
 * Hex printing is nibble-aware. We compute minimal nibble width required to
 * represent the signed value with an explicit sign nibble (for positives),
 * or minimal two's-complement nibble-width (for negatives).
 *
 * The implementation uses MP_LIMB_BITS to derive the number of nibbles per limb.
 */
int mp_print_hex(mp_int *x)
{
    static const char HEX_CHARS[17] = "0123456789abcdef";
    size_t bitlen;
    size_t mag_nibbles;
    size_t total_nibbles;
    size_t nib_index;
    mp_int absx, powb, mpval;
    size_t top_idx;
    mp_limb_t top;
    unsigned int tb;
    unsigned int nibbles_per_limb;

    if (!x) return FAILURE;
    if (x->sign == 0) {
        printf("0x0");
        return SUCCESS;
    }

    printf("0x");

    /* nibbles per limb = MP_LIMB_BITS / 4 (MP_LIMB_BITS is multiple of 4 in our setups) */
    nibbles_per_limb = (unsigned int)(MP_LIMB_BITS / 4u);

    /* positive values: minimal magnitude nibbles, but if top nibble >= 8 add leading 0 nibble */
    if (x->sign > 0) {
        if (x->length == 0) {
            bitlen = 0;
        } else {
            top_idx = x->length - 1;
            top = x->digits[top_idx];
            tb = 0;
            while (top) { top >>= 1; ++tb; }
            bitlen = ((size_t)top_idx) * (size_t)MP_LIMB_BITS + (size_t)tb;
        }

        /* minimal nibbles to represent magnitude */
        if (bitlen == 0) mag_nibbles = 1;
        else mag_nibbles = (bitlen + 3) / 4;

        /* determine whether top nibble's high bit would be 1 -> prepend 0 nibble */
        {
            size_t top_nib_index = mag_nibbles - 1;
            size_t limb = top_nib_index / (size_t)nibbles_per_limb;
            unsigned int nib_shift = (unsigned int)((top_nib_index % nibbles_per_limb) * 4u);
            mp_limb_t limbval = (limb < x->length) ? x->digits[limb] : (mp_limb_t)0;
            unsigned int top_nibble = (unsigned int)((limbval >> nib_shift) & (mp_limb_t)0xFu);
            if (top_nibble >= 8u) total_nibbles = mag_nibbles + 1;
            else total_nibbles = mag_nibbles;
        }

        /* print exactly total_nibbles hex digits (MSB first) */
        for (nib_index = total_nibbles; nib_index-- > 0; ) {
            size_t limb = nib_index / (size_t)nibbles_per_limb;
            unsigned int nib_shift = (unsigned int)((nib_index % nibbles_per_limb) * 4u);
            mp_limb_t limbval = (limb < x->length) ? x->digits[limb] : (mp_limb_t)0;
            unsigned int digit = (unsigned int)((limbval >> nib_shift) & (mp_limb_t)0xFu);
            putchar(HEX_CHARS[digit]);
        }

    } else {
        /* Negative: minimal two's-complement nibbles */
        mp_init(&absx); mp_init(&powb); mp_init(&mpval);

        mp_abs_copy(&absx, x);

        if (absx.length == 0) {
            bitlen = 0;
        } else {
            top_idx = absx.length - 1;
            top = absx.digits[top_idx];
            tb = 0;
            while (top) { top >>= 1; ++tb; }
            bitlen = ((size_t)top_idx) * (size_t)MP_LIMB_BITS + (size_t)tb;
        }

        /* wbits = bitlen + 1; compute minimal nibbles = ceil(wbits/4) */
        {
            size_t wbits = bitlen + 1;
            size_t nibbles = (wbits + 3) / 4;
            if (nibbles == 0) nibbles = 1;
            total_nibbles = nibbles;
        }

        /* powb = 2^(4 * total_nibbles) */
        if (mp_make_pow2(&powb, total_nibbles * 4) != SUCCESS) {
            mp_free(&absx); mp_free(&powb); mp_free(&mpval);
            return FAILURE;
        }

        /* mpval = powb + x  (x negative) */
        if (mp_add(&mpval, &powb, x) != SUCCESS) {
            mp_free(&absx); mp_free(&powb); mp_free(&mpval);
            return FAILURE;
        }

        /* print exactly total_nibbles digits from mpval */
        for (nib_index = total_nibbles; nib_index-- > 0; ) {
            size_t limb = nib_index / (size_t)nibbles_per_limb;
            unsigned int nib_shift = (unsigned int)((nib_index % nibbles_per_limb) * 4u);
            mp_limb_t limbval = (limb < mpval.length) ? mpval.digits[limb] : (mp_limb_t)0;
            unsigned int digit = (unsigned int)((limbval >> nib_shift) & (mp_limb_t)0xFu);
            putchar(HEX_CHARS[digit]);
        }

        mp_free(&absx); mp_free(&powb); mp_free(&mpval);
    }

    return SUCCESS;
}