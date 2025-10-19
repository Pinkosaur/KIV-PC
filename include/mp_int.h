#ifndef MP_INT_H
#define MP_INT_H

#include <stddef.h>

/* return codes */
#define SUCCESS 1
#define FAILURE 0

typedef struct {
    int sign;
    unsigned int *digits;
    size_t length;
    size_t capacity;
} mp_int;

/* memory & low-level helpers */
int mp_init(mp_int *x);
int mp_free(mp_int *x);
int mp_reserve(mp_int *x, size_t capacity);
int mp_copy(mp_int *dst, const mp_int *src);

/* conversion from strings */
int mp_from_str_dec(mp_int *x, const char *str);
int mp_from_str_bin(mp_int *x, const char *str);
int mp_from_str_hex(mp_int *x, const char *str);

/* small-operand helpers */
int mp_add_small(mp_int *x, unsigned int a);
int mp_mul_small(mp_int *x, unsigned int m);
int mp_div_small(mp_int *result, const mp_int *a, unsigned int divisor, unsigned int *remainder);

/* absolute ops used internally / optionally public */
int mp_cmp_abs(const mp_int *a, const mp_int *b);
int mp_add_abs(mp_int *result, const mp_int *a, const mp_int *b);
int mp_sub_abs(mp_int *result, const mp_int *a, const mp_int *b);

/* arithmetic ops */
int mp_add(mp_int *result, const mp_int *a, const mp_int *b);
int mp_sub(mp_int *result, const mp_int *a, const mp_int *b);
int mp_mul(mp_int *result, const mp_int *a, const mp_int *b);
int mp_div(mp_int *result, const mp_int *a, const mp_int *b);
int mp_mod(mp_int *r, const mp_int *a, const mp_int *b);

/* shifts / utility */
int mp_shift_left_words(mp_int *x, size_t k);
int mp_shift_right_words(mp_int *x, size_t k);

#endif /* MP_INT_H */
