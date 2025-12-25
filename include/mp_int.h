#ifndef MP_INT_H
#define MP_INT_H

/* * ANSI C90 portable multiple-precision integer header.
 *
 * This header handles platform-specific limb size selection:
 * - Systems with 64-bit long (e.g., Linux x64): Use 32-bit limbs.
 * - Systems with 32-bit long (e.g., Windows x64/x86): Use 16-bit limbs.
 */

#include <stddef.h> /* size_t */
#include <limits.h> /* ULONG_MAX */

/* Return codes */
#define SUCCESS 1
#define FAILURE 0

/* Output format selectors */
#define DEC 0
#define BIN 1
#define HEX 2

/* --- Limb Configuration -------------------------------------------------- */

#if ULONG_MAX > 0xFFFFFFFFUL
  /* Platform where unsigned long >= 64 bits (e.g. LP64) */
  typedef unsigned int mp_limb_t;    /* Limb: 32 bits */
  typedef unsigned long mp_double_t; /* Accumulator: 64 bits */
  #define MP_LIMB_BITS 32u
  #define MP_LIMB_MASK 0xFFFFFFFFUL
#else
  /* Typical 32-bit platform (or Windows x64 LLP64) */
  typedef unsigned short mp_limb_t;  /* Limb: 16 bits */
  typedef unsigned long mp_double_t; /* Accumulator: 32 bits */
  #define MP_LIMB_BITS 16u
  #define MP_LIMB_MASK 0xFFFFu
#endif

/* --- Type Definition ----------------------------------------------------- */

typedef struct {
    int sign;           /* 0 for zero, +1 or -1 for non-zero */
    mp_limb_t *digits;  /* Little-endian array of limbs */
    size_t length;      /* Number of used limbs */
    size_t capacity;    /* Allocated capacity in limbs */
} mp_int;

/* --- Memory Management --------------------------------------------------- */

/* Initialize an mp_int to zero. */
int mp_init(mp_int *x);

/* Free memory and reset mp_int to zero state. */
int mp_free(mp_int *x);

/* Ensure internal buffer has capacity for at least 'capacity' limbs. */
int mp_reserve(mp_int *x, size_t capacity);

/* Deep copy: dst = src. */
int mp_copy(mp_int *dst, const mp_int *src);

/* --- String Conversions -------------------------------------------------- */

int mp_from_str_dec(mp_int *x, const char *str);
int mp_from_str_bin(mp_int *x, const char *str);
int mp_from_str_hex(mp_int *x, const char *str);

/* --- Small Operand Arithmetic -------------------------------------------- */

int mp_add_small(mp_int *x, unsigned int a);
int mp_mul_small(mp_int *x, unsigned int m);
int mp_div_small(mp_int *result, const mp_int *a, unsigned int divisor, unsigned int *remainder);
int mp_inc(mp_int *x);

/* Returns 1 if x fits in a single limb/uint, 0 otherwise. */
int mp_fits_uint(const mp_int *x);

/* --- Absolute Value Arithmetic ------------------------------------------- */

int mp_cmp_abs(const mp_int *a, const mp_int *b);
int mp_add_abs(mp_int *result, const mp_int *a, const mp_int *b);
int mp_sub_abs(mp_int *result, const mp_int *a, const mp_int *b);

/* --- Signed Arithmetic --------------------------------------------------- */

int mp_add(mp_int *result, const mp_int *a, const mp_int *b);
int mp_sub(mp_int *result, const mp_int *a, const mp_int *b);
int mp_mul(mp_int *result, const mp_int *a, const mp_int *b);
int mp_div(mp_int *result, const mp_int *a, const mp_int *b);
int mp_mod(mp_int *r, const mp_int *a, const mp_int *b);

/* --- Bitwise Utilities --------------------------------------------------- */

int mp_shift_left_words(mp_int *x, size_t k);
int mp_shift_right_words(mp_int *x, size_t k);

#endif /* MP_INT_H */