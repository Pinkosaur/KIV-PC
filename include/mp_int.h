#ifndef MP_INT_H
#define MP_INT_H

/* ANSI C90 portable multiple-precision integer header.
 *
 * This header chooses an internal limb size at compile time:
 *  - If 'unsigned long' is wider than 32 bits (ULONG_MAX > 0xFFFFFFFF),
 *    we use 32-bit limbs (type: unsigned int) and unsigned long as the
 *    double-width accumulator.
 *  - Otherwise (typical 32-bit platforms), we use 16-bit limbs
 *    (type: unsigned short) and unsigned long as the accumulator.
 */

#include <stddef.h> /* size_t */
#include <limits.h> /* ULONG_MAX, USHRT_MAX */
#include <limits.h>

/* return codes */
#define SUCCESS 1
#define FAILURE 0

/* Output format selectors (used by print helpers elsewhere) */
#define DEC 0
#define BIN 1
#define HEX 2

/* --- Limb configuration -------------------------------------------------- */
/* Choose limb size based on ULONG_MAX. If unsigned long is wider than
   32 bits, we may safely use 32-bit limbs; otherwise use 16-bit limbs. */
#if ULONG_MAX > 0xFFFFFFFFUL
  /* Platform where unsigned long >= 64 bits (e.g. LP64) */
  typedef unsigned int mp_limb_t;    /* limb (32 bits) */
  typedef unsigned long mp_double_t; /* accumulator (>=64 bits) */
  #define MP_LIMB_BITS 32u
  #define MP_LIMB_MASK 0xFFFFFFFFUL
#else
  /* Typical 32-bit platform: use 16-bit limbs to ensure accumulator
     (unsigned long, usually 32-bit) holds limb*limb without overflow. */
  typedef unsigned short mp_limb_t;  /* limb (16 bits) */
  typedef unsigned long mp_double_t; /* accumulator (>=32 bits) */
  #define MP_LIMB_BITS 16u
  #define MP_LIMB_MASK 0xFFFFu
#endif

/* --- mp_int definition -------------------------------------------------- */

/* Multiple-precision integer type.
 *  - sign:  0 for zero, +1 or -1 for non-zero values
 *  - digits: little-endian array of limbs (least-significant limb first)
 *  - length: number of used limbs
 *  - capacity: allocated size (in limbs)
 */
typedef struct {
    int sign;
    mp_limb_t *digits;
    size_t length;
    size_t capacity;
} mp_int;

/* ---------------- Memory & low-level helpers ------------------------------
 *
 * All functions return SUCCESS (1) or FAILURE (0).
 */

/* Initialize mp_int object to zero/empty state. */
int mp_init(mp_int *x);

/* Free memory held by mp_int and reinitialize it to empty. */
int mp_free(mp_int *x);

/* Ensure capacity of at least 'capacity' limbs (no shrink). */
int mp_reserve(mp_int *x, size_t capacity);

/* Copy src into dst. dst must have been initialized with mp_init. */
int mp_copy(mp_int *dst, const mp_int *src);

/* ---------------- Conversions from strings ------------------------------- */

/* Parse decimal string into x (leading +/-, whitespace allowed). */
int mp_from_str_dec(mp_int *x, const char *str);

/* Parse binary string (optional sign, optional 0b/0B prefix). */
int mp_from_str_bin(mp_int *x, const char *str);

/* Parse hexadecimal string (optional sign, optional 0x/0X prefix). */
int mp_from_str_hex(mp_int *x, const char *str);

/* ---------------- Small-operand helpers ---------------------------------- */

/* Add small unsigned number 'a' to x (in-place): x += a. */
int mp_add_small(mp_int *x, unsigned int a);

/* Multiply x in-place by small unsigned m: x *= m. */
int mp_mul_small(mp_int *x, unsigned int m);

/* result = a / divisor (integer division). remainder optionally returned. */
int mp_div_small(mp_int *result, const mp_int *a, unsigned int divisor, unsigned int *remainder);

/* Increment x by 1 (convenience wrapper). */
int mp_inc(mp_int *x);

/* Return 1 if x is non-negative and fits into single limb (i.e. fits
   into an unsigned int as previously used). Otherwise 0. */
int mp_fits_uint(const mp_int *x);

/* ---------------- Absolute-value helpers (unsigned arithmetic) ---------- */

/* Compare absolute values: returns 1 if |a|>|b|, -1 if |a|<|b|, 0 if equal */
int mp_cmp_abs(const mp_int *a, const mp_int *b);

/* result = |a| + |b| (unsigned addition). */
int mp_add_abs(mp_int *result, const mp_int *a, const mp_int *b);

/* result = |a| - |b| (unsigned subtraction). Assumes |a| >= |b|. */
int mp_sub_abs(mp_int *result, const mp_int *a, const mp_int *b);

/* ---------------- Arithmetic (signed) ----------------------------------- */

/* result = a + b (signed). */
int mp_add(mp_int *result, const mp_int *a, const mp_int *b);

/* result = a - b (signed). */
int mp_sub(mp_int *result, const mp_int *a, const mp_int *b);

/* result = a * b (signed). Uses hybrid naive/Karatsuba. */
int mp_mul(mp_int *result, const mp_int *a, const mp_int *b);

/* result = a / b (integer division). */
int mp_div(mp_int *result, const mp_int *a, const mp_int *b);

/* r = a % b (remainder, sign follows a). */
int mp_mod(mp_int *r, const mp_int *a, const mp_int *b);

/* ---------------- Shifts / utility -------------------------------------- */

/* Shift left by k limbs (words): multiplies by (base)^(k) where base = 2^(MP_LIMB_BITS) */
int mp_shift_left_words(mp_int *x, size_t k);

/* Destructive shift right by k limbs (words). */
int mp_shift_right_words(mp_int *x, size_t k);

#endif /* MP_INT_H */
