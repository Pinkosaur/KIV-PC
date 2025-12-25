/*
 * mp_print.c -- Output formatting for multiple-precision integers.
 *
 * This module implements printing of mp_int structures in Decimal, Binary,
 * and Hexadecimal formats.
 *
 * Design:
 * - Decimal printing uses repeated division by a large power of 10 (BASE)
 * to produce output chunks in reverse order.
 * - Binary and Hexadecimal printing use bitwise operations and two's complement
 * arithmetic to generate precise, minimal representations for negative numbers.
 * - The code detects the limb size (16-bit vs 32-bit) to ensure safe arithmetic
 * and prevent overflows during printing operations on all platforms.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mp_int.h"
#include "mp_print.h"

/*
 * Helper: Copy absolute value of x into dst.
 * dst must be initialized before calling.
 */
static int mp_abs_copy(mp_int *dst, const mp_int *x)
{
    if (!dst || !x) return FAILURE;
    if (mp_copy(dst, x) != SUCCESS) return FAILURE;
    if (dst->length == 0) { dst->sign = 0; return SUCCESS; }
    dst->sign = +1;
    return SUCCESS;
}

/*
 * Helper: Build a power of 2 (pow2 = 2^k).
 * Used for calculating two's complement offsets.
 */
static int mp_make_pow2(mp_int *pow2, size_t k)
{
    if (!pow2) return FAILURE;
    mp_free(pow2);
    mp_init(pow2);

    /* Allocate at least one limb so mp_mul_small will work */
    if (mp_reserve(pow2, 1) != SUCCESS) return FAILURE;
    pow2->digits[0] = (mp_limb_t)1;
    pow2->length = 1;
    pow2->sign = 1;

    /* Multiply by two k times */
    while (k--) {
        if (mp_mul_small(pow2, 2U) != SUCCESS) return FAILURE;
    }
    return SUCCESS;
}

/*
 * Helper: Check if mp_int is a power of two (has exactly one bit set).
 * Assumes x is positive/unsigned. Returns 1 if true, 0 otherwise.
 */
static int is_power_of_two(const mp_int *x)
{
    size_t i;
    mp_limb_t top;
    
    if (x->length == 0) return 0;

    /* Check that all limbs except the most significant are zero */
    for (i = 0; i < x->length - 1; i++) {
        if (x->digits[i] != 0) return 0;
    }

    /* Check that the top limb is a power of 2 */
    top = x->digits[x->length - 1];
    return (top != 0) && ((top & (top - 1)) == 0);
}

/* ----------------------- Decimal printing -------------------------------- */

int mp_print_dec(mp_int *x)
{
    mp_int tmp, q;
    unsigned int rem;
    char **chunks = NULL;
    char **newchunks;
    size_t num_chunks, chunk_cap;
    size_t i;
    unsigned long ms_val;

    /* Define the printing base based on the limb size to prevent overflow.
     * On 16-bit limb systems (e.g. Windows), (rem << 16) must fit in 32-bit.
     */
    #if defined(MP_LIMB_BITS) && (MP_LIMB_BITS == 16)
        const unsigned int BASE = 10000U;       /* 10^4 */
        const char *FMT = "%04u";
    #else
        const unsigned int BASE = 1000000000U;  /* 10^9 */
        const char *FMT = "%09u";
    #endif

    if (!x) return FAILURE;
    if (x->sign == 0) {
        printf("0");
        return SUCCESS;
    }

    mp_init(&tmp);
    mp_init(&q);
    
    if (mp_copy(&tmp, x) != SUCCESS) goto cleanup;
    tmp.sign = +1; /* operate on magnitude */

    /* Allocate initial chunk array */
    chunk_cap = 64;
    chunks = (char**) malloc(chunk_cap * sizeof(char*));
    if (!chunks) goto cleanup;
    num_chunks = 0;

    /* Repeatedly divide by BASE and store remainder chunks */
    while (tmp.sign != 0) {
        if (num_chunks >= chunk_cap) {
            chunk_cap *= 2;
            newchunks = (char**) realloc(chunks, chunk_cap * sizeof(char*));
            if (!newchunks) { free(chunks); chunks = NULL; goto cleanup; }
            chunks = newchunks;
        }

        if (mp_div_small(&q, &tmp, BASE, &rem) != SUCCESS) {
            goto cleanup;
        }

        chunks[num_chunks] = (char*) malloc(16); 
        if (!chunks[num_chunks]) { goto cleanup; }
        sprintf(chunks[num_chunks], FMT, rem);
        num_chunks++;

        mp_free(&tmp);
        mp_init(&tmp);
        if (mp_copy(&tmp, &q) != SUCCESS) { goto cleanup; }
    }

    /* Print sign */
    if (x->sign < 0) putchar('-');

    /* Print most-significant chunk (no leading zeros) */
    if (num_chunks > 0) {
        ms_val = strtoul(chunks[num_chunks - 1], NULL, 10);
        printf("%lu", ms_val);
    }

    /* Print remaining chunks (padded with zeros) */
    for (i = num_chunks - 1; i > 0; --i) {
        printf("%s", chunks[i - 1]);
    }

    /* Cleanup */
    for (i = 0; i < num_chunks; ++i) free(chunks[i]);
    free(chunks);
    mp_free(&tmp);
    mp_free(&q);
    return SUCCESS;

cleanup:
    if (chunks) {
        for (i = 0; i < num_chunks; ++i) if(chunks[i]) free(chunks[i]);
        free(chunks);
    }
    mp_free(&tmp);
    mp_free(&q);
    return FAILURE;
}

/* ----------------------- Binary printing --------------------------------- */

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

    if (x->sign > 0) {
        /* Calculate bit length */
        if (x->length == 0) {
            bitlen = 0;
        } else {
            top_idx = x->length - 1;
            top = x->digits[top_idx];
            tb = 0;
            while (top) { top >>= 1; ++tb; }
            bitlen = ((size_t)top_idx) * (size_t)MP_LIMB_BITS + (size_t)tb;
        }

        /* Positive number: always prepend 0 to indicate sign */
        w = bitlen + 1;

        for (bit_index = w; bit_index-- > 0; ) {
            size_t limb = bit_index / (size_t)MP_LIMB_BITS;
            unsigned int pos = (unsigned int)(bit_index % (size_t)MP_LIMB_BITS);
            mp_limb_t limbval = (limb < x->length) ? x->digits[limb] : (mp_limb_t)0;
            bit = (unsigned int)((limbval >> pos) & (mp_limb_t)1);
            putchar('0' + (int)bit);
        }

    } else {
        /* Negative number: print as two's complement */
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

        /* Determine minimal bit width */
        if (is_power_of_two(&absx)) w = bitlen;
        else w = bitlen + 1;

        /* Compute 2^w */
        if (mp_make_pow2(&pow2, w) != SUCCESS) {
            mp_free(&absx); mp_free(&pow2); mp_free(&val);
            return FAILURE;
        }

        /* Value = 2^w - |x| = 2^w + x (since x is negative) */
        if (mp_add(&val, &pow2, x) != SUCCESS) {
            mp_free(&absx); mp_free(&pow2); mp_free(&val);
            return FAILURE;
        }

        /* Print bits */
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

/* ----------------------- Hexadecimal printing ----------------------------- */

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

    nibbles_per_limb = (unsigned int)(MP_LIMB_BITS / 4u);

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

        if (bitlen == 0) mag_nibbles = 1;
        else mag_nibbles = (bitlen + 3) / 4;

        /* Determine if a leading zero nibble is needed */
        {
            size_t top_nib_index = mag_nibbles - 1;
            size_t limb = top_nib_index / (size_t)nibbles_per_limb;
            unsigned int nib_shift = (unsigned int)((top_nib_index % nibbles_per_limb) * 4u);
            mp_limb_t limbval = (limb < x->length) ? x->digits[limb] : (mp_limb_t)0;
            unsigned int top_nibble = (unsigned int)((limbval >> nib_shift) & (mp_limb_t)0xFu);
            if (top_nibble >= 8u) total_nibbles = mag_nibbles + 1;
            else total_nibbles = mag_nibbles;
        }

        for (nib_index = total_nibbles; nib_index-- > 0; ) {
            size_t limb = nib_index / (size_t)nibbles_per_limb;
            unsigned int nib_shift = (unsigned int)((nib_index % nibbles_per_limb) * 4u);
            mp_limb_t limbval = (limb < x->length) ? x->digits[limb] : (mp_limb_t)0;
            unsigned int digit = (unsigned int)((limbval >> nib_shift) & (mp_limb_t)0xFu);
            putchar(HEX_CHARS[digit]);
        }

    } else {
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

        {
            size_t wbits;
            size_t nibbles;
            
            if (is_power_of_two(&absx)) wbits = bitlen;
            else wbits = bitlen + 1;

            nibbles = (wbits + 3) / 4;
            if (nibbles == 0) nibbles = 1;
            total_nibbles = nibbles;
        }

        if (mp_make_pow2(&powb, total_nibbles * 4) != SUCCESS) {
            mp_free(&absx); mp_free(&powb); mp_free(&mpval);
            return FAILURE;
        }

        if (mp_add(&mpval, &powb, x) != SUCCESS) {
            mp_free(&absx); mp_free(&powb); mp_free(&mpval);
            return FAILURE;
        }

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