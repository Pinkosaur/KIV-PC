#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

/* Output format selectors */
#define DEC 0
#define BIN 1
#define HEX 2

/* function returns */
#define SUCCESS 1
#define FAILURE 0

/* Multiplication algorithm selection based on mp_int size */
#define NAIVE_THRESHOLD 16

/* Multiple precision int representation */
typedef struct {
    int sign;
    unsigned int *digits;
    size_t length, capacity; /* number of limbs used (current, max) */
} mp_int;

/* Forward declarations (to avoid implicit-declaration warnings) */
int parse_operand_to_mp(mp_int *dst, const char *s);
int mp_from_str_bin(mp_int *x, const char *str);
int mp_from_str_hex(mp_int *x, const char *str);
int mp_print_dec(mp_int *x);
int mp_print_bin(mp_int *x);
int mp_print_hex(mp_int *x);

/* ---------------------------------- Memory management helpers -------------------------------------- */

/* Initializes mp_int x to empty state */
int mp_init(mp_int *x) {
    if (!x) return FAILURE;
    x->sign = 0;
    x->length = 0;
    x->digits = NULL;
    x->capacity = 0;
    return SUCCESS;
}

/* Frees memory used by the digits and reinitializes the mp_int */
int mp_free(mp_int *x) {
    if (!x) return SUCCESS;
    free(x->digits);
    /* reinitialize */
    x->digits = NULL;
    x->length = 0;
    x->capacity = 0;
    x->sign = 0;
    return SUCCESS;
}

/* Reallocates memory for the new desired size of a mp_int */
int mp_reserve(mp_int *x, size_t capacity) {
    unsigned int *new_digits;
    if (!x) return FAILURE;
    if (capacity <= x->capacity)
        return SUCCESS;
    /* allocate at least 1 to avoid realloc(NULL,0) ambiguity */
    if (capacity == 0) capacity = 1;
    new_digits = (unsigned int*) realloc(x->digits, capacity * sizeof(unsigned int));
    if (!new_digits)
        return FAILURE; /* out of memory */
    x->digits = new_digits;
    x->capacity = capacity;
    return SUCCESS;
}

/* Copy mp_int src -> dst (dst must be initialized with mp_init) */
int mp_copy(mp_int *dst, const mp_int *src) {
    if (!dst || !src) return FAILURE;
    if (src->length == 0) {
        dst->length = 0;
        dst->sign = 0;
        return SUCCESS;
    }
    if (mp_reserve(dst, src->length) != SUCCESS) return FAILURE;
    memcpy(dst->digits, src->digits, src->length * sizeof(unsigned int));
    dst->length = src->length;
    dst->sign = src->sign;
    return SUCCESS;
}

/* ----------------------------------------- Arithmetic logic ------------------------------------------- */

/* Compare absolute values of two mp_ints */
/* Returns 1 if a is greater, -1 if b is greater, otherwise 0 */
int mp_cmp_abs(const mp_int *a, const mp_int *b) {
    size_t i;
    if (!a || !b) return 0;
    if (a->length != b->length)
        return (a->length > b->length) ? 1 : -1;

    i = a->length;
    while (i > 0) {
        i--;
        if (a->digits[i] > b->digits[i]) return 1;
        if (a->digits[i] < b->digits[i]) return -1;
    }
    return 0;
}

/* Adds a mp_int and a regular unsigned int */
int mp_add_small(mp_int *x, unsigned int a)
{
    unsigned long carry;
    size_t i;
    unsigned long sum;

    if (!x) return FAILURE;
    carry = a;
    i = 0;
    while (carry && i < x->length) {
        sum = (unsigned long)x->digits[i] + carry;
        x->digits[i] = (unsigned int)(sum & 0xFFFFFFFFUL);
        carry = sum >> 32;
        i++;
    }
    if (carry) {
        if (mp_reserve(x, x->length + 1) != SUCCESS)
            return FAILURE;
        x->digits[x->length++] = (unsigned int)carry;
    }
    if (x->sign == 0 && x->length > 0)
        x->sign = +1;
    return SUCCESS;
}

/* Multiplies a mp_int by a regular unsigned int */
int mp_mul_small(mp_int *x, unsigned int m)
{
    unsigned long carry;
    size_t i;
    unsigned long prod;

    if (!x) return FAILURE;
    if (x->sign == 0)
        return SUCCESS; /* 0 * m = 0 */

    carry = 0;

    for (i = 0; i < x->length; ++i) {
        prod = (unsigned long)x->digits[i] * m + carry;
        x->digits[i] = (unsigned int)(prod & 0xFFFFFFFFUL);
        carry = prod >> 32;
    }
    if (carry) {
        if (mp_reserve(x, x->length + 1) != SUCCESS)
            return FAILURE;
        x->digits[x->length++] = (unsigned int)carry;
    }
    return SUCCESS;
}

/* Divides a mp_int by a regular unsigned int */
int mp_div_small(mp_int *result, const mp_int *a, unsigned int divisor, unsigned int *remainder) {
    unsigned long r;
    size_t i, len;
    unsigned long cur;

    if (!result || !a || divisor == 0) return FAILURE;
    if (a->sign == 0) {
        result->sign = 0;
        result->length = 0;
        if (remainder) *remainder = 0;
        return SUCCESS;
    }

    if (mp_reserve(result, a->length) != SUCCESS) return FAILURE;
    r = 0;
    i = a->length;
    while (i > 0) {
        i--;
        cur = (r << 32) | (unsigned long)a->digits[i];
        result->digits[i] = (unsigned int)(cur / divisor);
        r = cur % divisor;
    }

    len = a->length;
    while (len > 0 && result->digits[len - 1] == 0)
        len--;
    result->length = len;
    result->sign = (len == 0) ? 0 : a->sign;
    if (remainder) *remainder = (unsigned int)r;
    return SUCCESS;
}

/* Unsigned addition of mp_int (absolute) */
int mp_add_abs(mp_int *result, const mp_int *a, const mp_int *b) {
    unsigned long carry, av, bv, sum;
    size_t i, n;

    if (!result || !a || !b) return FAILURE;
    n = (a->length > b->length) ? a->length : b->length;
    if (mp_reserve(result, n + 1) != SUCCESS) return FAILURE;

    carry = 0;

    for (i = 0; i < n; ++i) {
        av = (i < a->length) ? a->digits[i] : 0;
        bv = (i < b->length) ? b->digits[i] : 0;
        sum = av + bv + carry;
        result->digits[i] = (unsigned int)(sum & 0xFFFFFFFFu);
        carry = sum >> 32;
    }
    if (carry) {
        result->digits[n++] = (unsigned int)carry;
    }
    result->length = n;
    result->sign = (n == 0) ? 0 : +1;
    return SUCCESS;
}

/* Unsigned subtraction: result = |a| - |b|  (assumes |a| >= |b|) */
int mp_sub_abs(mp_int *result, const mp_int *a, const mp_int *b) {
    size_t n;
    unsigned long av, bv, diff, borrow;
    size_t i;

    if (!result || !a || !b) return FAILURE;
    n = a->length;
    if (mp_reserve(result, n) != SUCCESS) return FAILURE;

    borrow = 0;
    for (i = 0; i < n; ++i) {
        av = (unsigned long)a->digits[i];
        bv = (i < b->length) ? (unsigned long)b->digits[i] : 0;
        diff = av - bv - borrow;
        if (av < bv + borrow)
            borrow = 1;
        else
            borrow = 0;
        result->digits[i] = (unsigned int)(diff & 0xFFFFFFFFUL);
    }

    /* Remove leading zero limbs */
    while (n > 0 && result->digits[n - 1] == 0)
        n--;
    result->length = n;
    result->sign = (n == 0) ? 0 : +1;
    return SUCCESS;
}

/* Copy-based mp_add (signed) */
int mp_add(mp_int *result, const mp_int *a, const mp_int *b) {
    if (!result || !a || !b) return FAILURE;
    if (a->sign == 0) return mp_copy(result, b);
    if (b->sign == 0) return mp_copy(result, a);

    if (a->sign == b->sign) {
        mp_add_abs(result, a, b);
        result->sign = a->sign;
    } else {
        int cmp = mp_cmp_abs(a, b);
        if (cmp == 0) {
            /* a == -b */
            result->sign = 0;
            result->length = 0;
        } else if (cmp > 0) {
            mp_sub_abs(result, a, b);
            result->sign = a->sign;
        } else {
            mp_sub_abs(result, b, a);
            result->sign = b->sign;
        }
    }
    return SUCCESS;
}

/* Subtracts mp_int b from mp_int a */
int mp_sub(mp_int *result, const mp_int *a, const mp_int *b) {
    mp_int tmp;

    if (!result || !a || !b) return FAILURE;
    mp_init(&tmp);
    if (mp_copy(&tmp, b) != SUCCESS) { mp_free(&tmp); return FAILURE; }
    tmp.sign = -tmp.sign;
    mp_add(result, a, &tmp);
    mp_free(&tmp);
    return SUCCESS;
}

/* Naively multiplies two mp_ints, suitable for smaller numbers*/
int mp_mul_naive(mp_int *result, const mp_int *a, const mp_int *b) {
    size_t n, m, i, j, len;
    unsigned long carry, av, idx;

    if (!result || !a || !b) return FAILURE;
    if (a->sign == 0 || b->sign == 0) {
        result->sign = 0;
        result->length = 0;
        return SUCCESS;
    }

    n = a->length;
    m = b->length;
    if (mp_reserve(result, n + m) != SUCCESS) return FAILURE;
    memset(result->digits, 0, (n + m) * sizeof(unsigned int));

    for (i = 0; i < n; ++i) {
        carry = 0;
        av = a->digits[i];
        for (j = 0; j < m; ++j) {
            idx = (unsigned long)result->digits[i + j]
                               + av * (unsigned long)b->digits[j] + carry;
            result->digits[i + j] = (unsigned int)(idx & 0xFFFFFFFFUL);
            carry = idx >> 32;
        }
        if (carry)
            result->digits[i + m] += (unsigned int)carry;
    }

    len = n + m;
    while (len > 0 && result->digits[len - 1] == 0)
        len--;
    result->length = len;
    result->sign = a->sign * b->sign;
    return SUCCESS;
}

/* Logical left shift by k limbs (base 2^(32*k)) */
int mp_shift_left_words(mp_int *x, size_t k) {
    if (!x) return FAILURE;
    if (x->sign == 0 || k == 0) return SUCCESS;
    if (mp_reserve(x, x->length + k) != SUCCESS)
        return FAILURE;
    memmove(x->digits + k, x->digits, x->length * sizeof(unsigned int));
    memset(x->digits, 0, k * sizeof(unsigned int));
    x->length += k;
    return SUCCESS;
}

/* Helper: right-shift by k limbs (words). Destructive. */
int mp_shift_right_words(mp_int *x, size_t k) {
    size_t i;
    if (!x) return FAILURE;
    if (x->length == 0 || k == 0) return SUCCESS;
    if (k >= x->length) {
        /* becomes zero */
        x->length = 0;
        x->sign = 0;
        return SUCCESS;
    }
    /* shift down */
    for (i = 0; i + k < x->length; ++i) {
        x->digits[i] = x->digits[i + k];
    }
    /* clear high limbs */
    for (; i < x->length; ++i) x->digits[i] = 0;
    x->length -= k;
    /* trim leading zeros */
    while (x->length > 0 && x->digits[x->length - 1] == 0) x->length--;
    if (x->length == 0) x->sign = 0;
    return SUCCESS;
}

static int mp_is_zero(const mp_int *x) {
    if (!x) return 1;
    return (x->length == 0 || x->sign == 0) ? 1 : 0;
}

void mp_split(const mp_int *src, mp_int *low, mp_int *high, size_t m) {
    size_t i;

    if (!src || !low || !high) return;

    mp_reserve(low, m);
    mp_reserve(high, (src->length > m) ? src->length - m : 0);

    low->length = (src->length < m) ? src->length : m;
    for (i = 0; i < low->length; ++i)
        low->digits[i] = src->digits[i];
    low->sign = (low->length == 0) ? 0 : +1;

    if (src->length > m) {
        high->length = src->length - m;
        for (i = 0; i < high->length; ++i)
            high->digits[i] = src->digits[i + m];
        high->sign = +1;
    } else {
        high->length = 0;
        high->sign = 0;
    }
}

/* Karatsuba multiply (recursive) */
int mp_karatsuba_mul(mp_int *result, const mp_int *a, const mp_int *b)
{
    size_t n, m;
    mp_int x0, x1, y0, y1, z0, z1, z2, tmp1, tmp2, res_tmp;

    if (!result || !a || !b) return FAILURE;

    /* Base case: small numbers → naive multiply */
    if (a->length <= 16 || b->length <= 16) {
        return mp_mul_naive(result, a, b);
    }

    n = (a->length > b->length) ? a->length : b->length;
    m = n / 2;

    mp_init(&x0); mp_init(&x1);
    mp_init(&y0); mp_init(&y1);
    mp_init(&z0); mp_init(&z1); mp_init(&z2);
    mp_init(&tmp1); mp_init(&tmp2);
    mp_init(&res_tmp);

    mp_split(a, &x0, &x1, m);
    mp_split(b, &y0, &y1, m);

    /* z0 = x0 * y0 */
    mp_karatsuba_mul(&z0, &x0, &y0);

    /* z2 = x1 * y1 */
    mp_karatsuba_mul(&z2, &x1, &y1);

    /* tmp1 = x0 + x1, tmp2 = y0 + y1 */
    mp_add_abs(&tmp1, &x0, &x1);
    mp_add_abs(&tmp2, &y0, &y1);

    /* z1 = (x0+x1)*(y0+y1) - z2 - z0 */
    mp_karatsuba_mul(&z1, &tmp1, &tmp2);
    mp_sub_abs(&z1, &z1, &z2);
    mp_sub_abs(&z1, &z1, &z0);

    /* combine result = z0 + (z1 << (m words)) + (z2 << (2*m words)) */
    mp_init(&res_tmp);
    mp_add_abs(&res_tmp, &z0, &res_tmp);
    mp_shift_left_words(&z1, m);
    mp_add_abs(&res_tmp, &res_tmp, &z1);
    mp_shift_left_words(&z2, 2 * m);
    mp_add_abs(&res_tmp, &res_tmp, &z2);

    /* move res_tmp to result */
    mp_copy(result, &res_tmp);
    result->sign = a->sign * b->sign;

    /* free temps */
    mp_free(&x0); mp_free(&x1);
    mp_free(&y0); mp_free(&y1);
    mp_free(&z0); mp_free(&z1); mp_free(&z2);
    mp_free(&tmp1); mp_free(&tmp2);
    mp_free(&res_tmp);

    return SUCCESS;
}

/* Multiplies two mp_ints (hybrid) */
int mp_mul(mp_int *result, const mp_int *a, const mp_int *b) {
    if (!result || !a || !b) return FAILURE;
    if (a->length < NAIVE_THRESHOLD || b->length < NAIVE_THRESHOLD)
        return mp_mul_naive(result, a, b);
    else
        return mp_karatsuba_mul(result, a, b);
}

/* Divides mp_int b by mp_int a */
/* Note: full long division is not implemented; only small divisor supported via mp_div_small */
int mp_div(mp_int *result, const mp_int *a, const mp_int *b) {
    if (!result || !a || !b) return FAILURE;
    /* If divisor is a single limb, use mp_div_small */
    if (b->length == 1) {
        unsigned int rem;
        return mp_div_small(result, a, b->digits[0], &rem);
    }
    return FAILURE; /* not implemented */
}

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


/* ----------------- Factorial using precomputed table + productRange ----------------- */

/* Precomputed table of factorials (n -> n! as decimal string).
   You can extend this table with more entries if desired. */
typedef struct {
    unsigned int n;
    const char *fact_str;
} FactEntry;

static FactEntry fact_table[] = {
    {0, "1"},
    {1, "1"},
    {2, "2"},
    {5, "120"},
    {10, "3628800"},
    {50, "30414093201713378043612608166064768844377641568960512000000000000"},
    {75, "24809140811395398091946477116594033660926243886570122837795894512655842677572867409443815424000000000000000000"},
    {100, "93326215443944152681699238856266700490715968264381621468592963895217599993229915608941463976156518286253697920827223758251185210916864000000000000000000000000"},
    {200, "788657867364790503552363213932185062295135977687173263294742533244359449963403342920304284011984623904177212138919638830257642790242637105061926624952829931113462857270763317237396988943922445621451664240254033291864131227428294853277524242407573903240321257405579568660226031904170324062351700858796178922222789623703897374720000000000000000000000000000000000000000000000000"},
    {0, NULL} /* sentinel */
};

/* Helper: create mp_int representing small unsigned value */
static int mp_set_uint(mp_int *dst, unsigned int v) {
    if (!dst) return FAILURE;
    mp_free(dst);
    mp_init(dst);
    if (v == 0) {
        dst->sign = 0;
        dst->length = 0;
        return SUCCESS;
    }
    if (mp_reserve(dst, 1) != SUCCESS) return FAILURE;
    dst->digits[0] = v;
    dst->length = 1;
    dst->sign = +1;
    return SUCCESS;
}

/* productRange(low, high, out)
   computes product of integers in [low..high] into out (mp_int).
   returns SUCCESS/FAILURE. ANSI C90-compliant. */
static int productRange(unsigned int low, unsigned int high, mp_int *out) {
    unsigned int mid;
    mp_int left, right, tmp;
    int status;

    if (!out) return FAILURE;

    mp_init(&left); mp_init(&right); mp_init(&tmp);

    if (low > high) {
        /* empty product = 1 */
        status = mp_set_uint(out, 1);
        mp_free(&left); mp_free(&right); mp_free(&tmp);
        return status;
    }

    if (low == high) {
        status = mp_set_uint(out, low);
        mp_free(&left); mp_free(&right); mp_free(&tmp);
        return status;
    }

    if (high - low == 1) {
        /* out = low * high */
        if (mp_set_uint(&left, low) != SUCCESS) { mp_free(&left); mp_free(&right); mp_free(&tmp); return FAILURE; }
        if (mp_set_uint(&right, high) != SUCCESS) { mp_free(&left); mp_free(&right); mp_free(&tmp); return FAILURE; }
        if (mp_mul(out, &left, &right) != SUCCESS) { mp_free(&left); mp_free(&right); mp_free(&tmp); return FAILURE; }
        mp_free(&left); mp_free(&right); mp_free(&tmp);
        return SUCCESS;
    }

    /* small range: iterative multiply using mp_mul_small for speed */
    if (high - low < 10) {
        if (mp_set_uint(out, low) != SUCCESS) { mp_free(&left); mp_free(&right); mp_free(&tmp); return FAILURE; }
        {
            unsigned int i;
            for (i = low + 1; i <= high; ++i) {
                if (mp_mul_small(out, i) != SUCCESS) {
                    mp_free(&left); mp_free(&right); mp_free(&tmp);
                    return FAILURE;
                }
            }
        }
        mp_free(&left); mp_free(&right); mp_free(&tmp);
        return SUCCESS;
    }

    /* divide and conquer */
    mid = (low + high) / 2;

    if (productRange(low, mid, &left) != SUCCESS) {
        mp_free(&left); mp_free(&right); mp_free(&tmp);
        return FAILURE;
    }
    if (productRange(mid + 1, high, &right) != SUCCESS) {
        mp_free(&left); mp_free(&right); mp_free(&tmp);
        return FAILURE;
    }

    status = mp_mul(&tmp, &left, &right);
    if (status != SUCCESS) {
        mp_free(&left); mp_free(&right); mp_free(&tmp);
        return FAILURE;
    }

    /* move tmp to out */
    mp_free(out);
    mp_init(out);
    if (mp_copy(out, &tmp) != SUCCESS) {
        mp_free(&left); mp_free(&right); mp_free(&tmp);
        return FAILURE;
    }

    mp_free(&left); mp_free(&right); mp_free(&tmp);
    return SUCCESS;
}

/* Convert mp_int (non-negative) to uint32_t if possible.
   Returns SUCCESS and writes value into *out if convertible, otherwise FAILURE. */
static int mp_to_uint32(const mp_int *x, unsigned int *out) {
    if (!x || !out) return FAILURE;
    if (x->sign < 0) return FAILURE;
    if (x->length == 0) { *out = 0U; return SUCCESS; }
    if (x->length > 1) {
        /* if any higher limb non-zero -> overflow */
        if (x->length >= 2) {
            size_t i;
            for (i = 1; i < x->length; ++i) {
                if (x->digits[i] != 0U) return FAILURE;
            }
        }
    }
    /* now lower limb holds value (may be > UINT_MAX but fits in unsigned int on typical platforms) */
    /* Be strict: if value > UINT32_MAX reject */
    if (x->digits[0] > 0xFFFFFFFFU) return FAILURE; /* technically always false for 32-bit unsigned int limb */
    *out = x->digits[0];
    return SUCCESS;
}


int mp_from_str_dec(mp_int *x, const char *str);
/* mp_fact: compute factorial of integer 'a' and store in r.
   Expects 'a' to be non-negative and reasonably small (fits into uint32_t).
   Uses precomputed factorials and productRange for remaining range.
*/
int mp_fact(mp_int *r, const mp_int *a) {
    unsigned int n;
    size_t best_idx;
    int found;
    mp_int prod;
    mp_int pre;
    mp_int tmp;

    if (!r || !a) return FAILURE;

    /* negative input -> error */
    if (a->sign < 0) {
        fprintf(stderr, "error: factorial of negative number\n");
        return FAILURE;
    }

    /* convert to uint32 */
    if (mp_to_uint32(a, &n) != SUCCESS) {
        fprintf(stderr, "error: factorial argument too large or not integer\n");
        return FAILURE;
    }

    /* trivial cases */
    if (n < 2U) {
        /* r = 1 */
        mp_free(r);
        mp_init(r);
        if (mp_set_uint(r, 1U) != SUCCESS) return FAILURE;
        return SUCCESS;
    }

    /* find largest precomputed <= n */
    best_idx = 0;
    found = 0;
    {
        size_t i;
        for (i = 0; fact_table[i].fact_str != NULL; ++i) {
            if (fact_table[i].n <= n) {
                best_idx = i;
                found = 1;
            }
        }
    }

    /* if found, set pre = precomputed[best_idx], else pre = 1 */
    mp_init(&pre);
    if (found) {
        if (mp_from_str_dec(&pre, fact_table[best_idx].fact_str) != SUCCESS) {
            mp_free(&pre);
            return FAILURE;
        }
    } else {
        if (mp_set_uint(&pre, 1U) != SUCCESS) { mp_free(&pre); return FAILURE; }
        best_idx = 0;
    }

    /* If precomputed base < n, compute productRange(base+1, n) and multiply */
    mp_init(&prod);
    mp_init(&tmp);

    if (fact_table[best_idx].n < n) {
        unsigned int low = fact_table[best_idx].n + 1U;
        unsigned int high = n;
        if (productRange(low, high, &prod) != SUCCESS) {
            mp_free(&pre); mp_free(&prod); mp_free(&tmp);
            return FAILURE;
        }
        /* tmp = pre * prod */
        if (mp_mul(&tmp, &pre, &prod) != SUCCESS) {
            mp_free(&pre); mp_free(&prod); mp_free(&tmp);
            return FAILURE;
        }
        /* move tmp to r */
        mp_free(r);
        mp_init(r);
        if (mp_copy(r, &tmp) != SUCCESS) {
            mp_free(&pre); mp_free(&prod); mp_free(&tmp);
            return FAILURE;
        }
    } else {
        /* pre already equals n! */
        mp_free(r);
        mp_init(r);
        if (mp_copy(r, &pre) != SUCCESS) {
            mp_free(&pre); mp_free(&prod); mp_free(&tmp);
            return FAILURE;
        }
    }

    /* cleanup */
    mp_free(&pre);
    mp_free(&prod);
    mp_free(&tmp);
    return SUCCESS;
}


/* mp_mod: r = a % b
   Remainder has same sign as 'a' (like C's %). Uses only ANSI C90 features. */
int mp_mod(mp_int *r, const mp_int *a, const mp_int *b)
{
    /* variables declared at top per C90 */
    size_t shift_words;
    int cmp;
    int status = FAILURE;

    mp_int rem;
    mp_int tmp;
    mp_int candidate;
    mp_int doubled;
    mp_int newrem;
    unsigned int rem_small;

    if (!r || !a || !b) return FAILURE;
    /* divisor zero */
    if (b->length == 0 || b->sign == 0) return FAILURE;

    /* dividend zero => remainder 0 */
    if (a->length == 0 || a->sign == 0) {
        mp_free(r);
        mp_init(r);
        r->sign = 0;
        r->length = 0;
        return SUCCESS;
    }

    /* fast path: single-limb divisor */
    if (b->length == 1) {
        mp_int q;
        mp_init(&q);
        if (mp_div_small(&q, a, b->digits[0], &rem_small) != SUCCESS) {
            mp_free(&q);
            return FAILURE;
        }
        mp_free(&q);

        mp_free(r);
        mp_init(r);
        if (rem_small == 0) {
            r->sign = 0; r->length = 0;
        } else {
            if (mp_reserve(r, 1) != SUCCESS) return FAILURE;
            r->digits[0] = rem_small;
            r->length = 1;
            r->sign = (a->sign < 0) ? -1 : +1;
        }
        return SUCCESS;
    }

    /* General multi-limb algorithm (word-shift + repeated-doubling subtraction) */

    /* initialize temporaries */
    mp_init(&rem);
    mp_init(&tmp);
    mp_init(&candidate);
    mp_init(&doubled);
    mp_init(&newrem);

    /* rem = |a| */
    if (mp_copy(&rem, a) != SUCCESS) goto cleanup;
    rem.sign = +1;

    /* tmp = |b| */
    if (mp_copy(&tmp, b) != SUCCESS) goto cleanup;
    tmp.sign = +1;

    /* if |rem| < |tmp| -> remainder is a */
    cmp = mp_cmp_abs(&rem, &tmp);
    if (cmp < 0) {
        mp_free(r);
        mp_init(r);
        if (mp_copy(r, &rem) != SUCCESS) goto cleanup;
        r->sign = (a->sign < 0) ? -1 : +1;
        status = SUCCESS;
        goto cleanup;
    }

    /* shift tmp left so its top limb lines up with rem top limb */
    shift_words = rem.length - tmp.length;
    if (shift_words > 0) {
        if (mp_shift_left_words(&tmp, shift_words) != SUCCESS) goto cleanup;
    }

    /* iterate from current shift down to 0 */
    for ( ; ; ) {
        /* while rem >= tmp, subtract the largest multiple of tmp we can get by doubling */
        while (mp_cmp_abs(&rem, &tmp) >= 0) {
            /* candidate = tmp */
            mp_free(&candidate);
            mp_init(&candidate);
            if (mp_copy(&candidate, &tmp) != SUCCESS) goto cleanup;

            /* double candidate repeatedly while doubled <= rem */
            while (1) {
                mp_free(&doubled);
                mp_init(&doubled);
                if (mp_add_abs(&doubled, &candidate, &candidate) != SUCCESS) goto cleanup;
                /* if doubled > rem break */
                if (mp_cmp_abs(&doubled, &rem) > 0) {
                    mp_free(&doubled);
                    break;
                }
                /* accept doubled as new candidate */
                mp_free(&candidate);
                mp_init(&candidate);
                if (mp_copy(&candidate, &doubled) != SUCCESS) goto cleanup;
                /* loop and try doubling again */
            }

            /* rem = rem - candidate  (candidate <= rem guaranteed) */
            mp_free(&newrem);
            mp_init(&newrem);
            if (mp_sub_abs(&newrem, &rem, &candidate) != SUCCESS) goto cleanup;
            /* move newrem -> rem */
            mp_free(&rem);
            mp_init(&rem);
            if (mp_copy(&rem, &newrem) != SUCCESS) goto cleanup;
            /* free candidate, will be recreated on next iteration if needed */
            mp_free(&candidate);
        }

        /* if we've shifted down to zero words stop */
        if (shift_words == 0) break;

        /* shift tmp right by one word and continue */
        if (mp_shift_right_words(&tmp, 1) != SUCCESS) goto cleanup;
        shift_words--;
    }

    /* rem is the absolute remainder; attach sign of a */
    mp_free(r);
    mp_init(r);
    if (mp_copy(r, &rem) != SUCCESS) goto cleanup;
    if (r->length == 0) r->sign = 0;
    else r->sign = (a->sign < 0) ? -1 : +1;

    status = SUCCESS;

cleanup:
    mp_free(&rem);
    mp_free(&tmp);
    mp_free(&candidate);
    mp_free(&doubled);
    mp_free(&newrem);
    return status;
}

/* ----------------------------------------------- Parsing -------------------------------------------------- */

/* Operator precedence and associativity */
typedef struct {
    char op;
    int precedence;
    int right_assoc;
    int unary;
} OpInfo;

static OpInfo op_table[] = {
    {'!', 5, 1, 1},  /* factorial, unary postfix */
    {'^', 4, 1, 0},  /* exponentiation */
    {'*', 3, 0, 0},
    {'/', 3, 0, 0},
    {'%', 3, 0, 0},
    {'+', 2, 0, 0},
    {'-', 2, 0, 0},
    {'(', 1, 0, 0},
    {0, 0, 0, 0}
};

static int get_precedence(char op) {
    int i;
    for (i = 0; op_table[i].op; ++i)
        if (op_table[i].op == op) return op_table[i].precedence;
    return 0;
}

static int is_right_assoc(char op) {
    int i;
    for (i = 0; op_table[i].op; ++i)
        if (op_table[i].op == op) return op_table[i].right_assoc;
    return 0;
}

static int is_unary(char op) {
    int i;
    for (i = 0; op_table[i].op; ++i)
        if (op_table[i].op == op) return op_table[i].unary;
    return 0;
}

static int is_operator(char c) {
    return strchr("+-*/%^!", c) != NULL;
}

/* ---------- Tokenizer ---------- */
typedef struct {
    char tokens[256][64];
    int count;
} TokenList;

/*
 Tokenizer rules:
  - Produces tokens that are either:
      * an operator single-char: "+", "-", "*", "/", "%", "^", "!"
      * "(" or ")"
      * a full operand string (may include a leading + or - sign, and hex/bin prefixes)
  - A leading + / - is considered part of a number when it appears:
      * at the very start of the expression, or
      * immediately after '(' or another operator token.
*/
static void tokenize(const char *expr, TokenList *out) {
    const char *p;
    char buf[64];
    int i;
    int prev_was_operand; /* 0 = start/after-operator/'(', 1 = after-operand or ')' */

    out->count = 0;
    prev_was_operand = 0; /* at start, treat +/ - as sign if they appear */

    p = expr;
    while (*p) {
        /* skip whitespace */
        while (isspace((unsigned char)*p)) p++;
        if (*p == '\0') break;

        /* safety: don't overflow tokens */
        if (out->count >= 256) break;

        /* If we see + or - we must decide: operator or sign of number */
        if ((*p == '+' || *p == '-') && !prev_was_operand) {
            /* treat + / - as sign for a following number token if next chars form a number */
            const char *q = p + 1;
            /* if next is whitespace then it's ambiguous; treat as operator */
            if (*q != '\0' && !isspace((unsigned char)*q)) {
                /* gather a number token starting with sign */
                i = 0;
                if (i < 63) buf[i++] = *p; /* leading sign */
                p++; /* consume sign */
                while (*p && !isspace((unsigned char)*p) && !is_operator(*p) && *p != '(' && *p != ')') {
                    if (i < 63) buf[i++] = *p;
                    p++;
                }
                buf[i] = '\0';
                /* store token */
                strncpy(out->tokens[out->count], buf, 64);
                out->tokens[out->count][63] = '\0';
                out->count++;
                prev_was_operand = 1;
                continue;
            }
            /* else fall through to treat + / - as operator */
        }

        /* If it's a single-char operator or parentheses, emit operator token */
        if (is_operator(*p) || *p == '(' || *p == ')') {
            out->tokens[out->count][0] = *p;
            out->tokens[out->count][1] = '\0';
            out->count++;
            /* update prev flag: '(' is not an operand, ')' is */
            if (*p == '(') prev_was_operand = 0;
            else if (*p == ')') prev_was_operand = 1;
            else prev_was_operand = 0; /* operator */
            p++;
            continue;
        }

        /* Otherwise it's the start of a plain operand token (no leading sign) */
        i = 0;
        while (*p && !isspace((unsigned char)*p) && !is_operator(*p) && *p != '(' && *p != ')') {
            if (i < 63) buf[i++] = *p;
            p++;
        }
        buf[i] = '\0';
        strncpy(out->tokens[out->count], buf, 64);
        out->tokens[out->count][63] = '\0';
        out->count++;
        prev_was_operand = 1;
    }
}

/* ---------- Shunting Yard: infix -> postfix ---------- */
static void to_postfix(TokenList *infix, TokenList *postfix) {
    char stack[256];
    int i;
    int sp = 0;
    postfix->count = 0;

    for (i = 0; i < infix->count; ++i) {
        char *tok = infix->tokens[i];
        /* operator token (single char) */
        if (tok[0] != '\0' && tok[1] == '\0' && is_operator(tok[0])) {
            char op = tok[0];
            if (is_unary(op)) {
                /* postfix unary operators: they behave like operators with precedence */
                while (sp > 0 && get_precedence(stack[sp-1]) >= get_precedence(op)) {
                    postfix->tokens[postfix->count][0] = stack[--sp];
                    postfix->tokens[postfix->count][1] = '\0';
                    postfix->count++;
                }
                stack[sp++] = op;
            } else {
                /* binary operator: pop according to precedence/associativity */
                while (sp > 0) {
                    char top = stack[sp-1];
                    if (is_operator(top) &&
                        ((is_right_assoc(op) && get_precedence(op) < get_precedence(top)) ||
                         (!is_right_assoc(op) && get_precedence(op) <= get_precedence(top)))) {
                        postfix->tokens[postfix->count][0] = stack[--sp];
                        postfix->tokens[postfix->count][1] = '\0';
                        postfix->count++;
                    } else break;
                }
                stack[sp++] = op;
            }
        } else if (tok[0] == '(' && tok[1] == '\0') {
            stack[sp++] = '(';
        } else if (tok[0] == ')' && tok[1] == '\0') {
            while (sp > 0 && stack[sp-1] != '(') {
                postfix->tokens[postfix->count][0] = stack[--sp];
                postfix->tokens[postfix->count][1] = '\0';
                postfix->count++;
            }
            if (sp > 0 && stack[sp-1] == '(') sp--; /* pop '(' */
        } else {
            /* operand (possibly multi-char) */
            strncpy(postfix->tokens[postfix->count], tok, 64);
            postfix->tokens[postfix->count][63] = '\0';
            postfix->count++;
        }
    }
    while (sp > 0) {
        postfix->tokens[postfix->count][0] = stack[--sp];
        postfix->tokens[postfix->count][1] = '\0';
        postfix->count++;
    }
}

/* ---------- Postfix evaluator ---------- */
static int eval_postfix(TokenList *postfix, mp_int *result) {
    mp_int stack[128];
    int i;
    int sp = 0;

    for (i = 0; i < postfix->count; ++i) {
        char *tok = postfix->tokens[i];
        /* operator token (single char) */
        if (tok[0] != '\0' && tok[1] == '\0' && is_operator(tok[0])) {
            char op = tok[0];
            if (is_unary(op)) {
                /* postfix unary: pop one operand, apply, push result */
                if (sp < 1) return FAILURE;
                {
                    mp_int a;
                    mp_int out;
                    mp_init(&a); mp_init(&out);
                    /* pop */
                    mp_copy(&a, &stack[--sp]);
                    mp_free(&stack[sp]);
                    /* evaluate */
                    if (op == '!') {
                        if (mp_fact(&out, &a) != SUCCESS) {
                            mp_free(&a); mp_free(&out); return FAILURE;
                        }
                    } else {
                        /* unknown unary - should not happen */
                        mp_free(&a); mp_free(&out); return FAILURE;
                    }
                    mp_free(&a);
                    /* push out */
                    mp_init(&stack[sp]);
                    mp_copy(&stack[sp], &out);
                    mp_free(&out);
                    sp++;
                }
            } else {
                /* binary operator */
                if (sp < 2) return FAILURE;
                {
                    mp_int a;
                    mp_int b;
                    mp_int out;
                    mp_init(&a); mp_init(&b); mp_init(&out);

                    /* pop b then a */
                    mp_copy(&b, &stack[--sp]);
                    mp_free(&stack[sp]);
                    mp_copy(&a, &stack[--sp]);
                    mp_free(&stack[sp]);

                    /* compute out = a op b */
                    if (op == '+') {
                        if (mp_add(&out, &a, &b) != SUCCESS) { mp_free(&a); mp_free(&b); mp_free(&out); return FAILURE; }
                    } else if (op == '-') {
                        if (mp_sub(&out, &a, &b) != SUCCESS) { mp_free(&a); mp_free(&b); mp_free(&out); return FAILURE; }
                    } else if (op == '*') {
                        if (mp_mul(&out, &a, &b) != SUCCESS) { mp_free(&a); mp_free(&b); mp_free(&out); return FAILURE; }
                    } else if (op == '/') {
                        if (mp_div(&out, &a, &b) != SUCCESS) { mp_free(&a); mp_free(&b); mp_free(&out); return FAILURE; }
                    } else if (op == '%') {
                        if (mp_mod(&out, &a, &b) != SUCCESS) { mp_free(&a); mp_free(&b); mp_free(&out); return FAILURE; }
                    } else if (op == '^') {
                        if (mp_pow(&out, &a, &b) != SUCCESS) { mp_free(&a); mp_free(&b); mp_free(&out); return FAILURE; }
                    } else {
                        mp_free(&a); mp_free(&b); mp_free(&out); return FAILURE;
                    }

                    mp_free(&a);
                    mp_free(&b);

                    /* push result */
                    mp_init(&stack[sp]);
                    mp_copy(&stack[sp], &out);
                    mp_free(&out);
                    sp++;
                }
            }
        } else {
            /* operand: parse and push */
            if (sp >= 128) return FAILURE;
            mp_init(&stack[sp]);
            if (parse_operand_to_mp(&stack[sp], tok) != SUCCESS) {
                /* cleanup partially built stack */
                int j;
                for (j = 0; j <= sp; ++j) mp_free(&stack[j]);
                return FAILURE;
            }
            sp++;
        }
    }

    if (sp != 1) {
        int j;
        for (j = 0; j < sp; ++j) mp_free(&stack[j]);
        return FAILURE;
    }

    mp_copy(result, &stack[0]);
    mp_free(&stack[0]);
    return SUCCESS;
}

/* ---------- Decimal ---------- */
int mp_from_str_dec(mp_int *x, const char *str)
{
    const char *p;
    int sign;
    unsigned int d;

    if (!x || !str) return FAILURE;
    mp_free(x);
    mp_init(x);

    p = str;
    while (isspace((unsigned char)*p)) p++;

    sign = +1;
    if (*p == '-') { sign = -1; p++; }
    else if (*p == '+') { p++; }

    while (*p == '0') p++;

    if (!isdigit((unsigned char)*p)) {
        x->sign = 0;
        x->length = 0;
        return SUCCESS;
    }

    if (mp_reserve(x, 1) != SUCCESS) return FAILURE;
    x->digits[0] = 0;
    x->length = 1;
    x->sign = +1;

    for (; *p; ++p) {
        if (!isdigit((unsigned char)*p)) break;
        d = (unsigned int)(*p - '0');
        if (mp_mul_small(x, 10U) != SUCCESS) return FAILURE;
        if (mp_add_small(x, d) != SUCCESS) return FAILURE;
    }

    if (x->length == 1 && x->digits[0] == 0)
        x->sign = 0;
    else
        x->sign = sign;
    return SUCCESS;
}

/* ---------- Binary (with implicit sign bit handling) ---------- */
int mp_from_str_bin(mp_int *x, const char *str) {
    size_t i;
    const char *p;
    const char *start;
    unsigned int bit;

    if (!x || !str) return FAILURE;
    mp_free(x);
    mp_init(x);

    p = str;
    while (isspace((unsigned char)*p)) p++;

    {
        int explicit_sign = +1;
        if (*p == '-') { explicit_sign = -1; p++; }
        else if (*p == '+') { p++; }

        if (p[0] == '0' && (p[1] == 'b' || p[1] == 'B')) p += 2;

        start = p;
        while (*p == '0') p++;

        if (*p != '0' && *p != '1' && *p != '\0') {
            x->sign = 0; x->length = 0;
            return SUCCESS;
        }

        if (mp_reserve(x, 1) != SUCCESS) return FAILURE;
        x->digits[0] = 0;
        x->length = 1;

        for (; *p; ++p) {
            if (*p != '0' && *p != '1') break;
            bit = (unsigned int)(*p - '0');
            if (mp_mul_small(x, 2U) != SUCCESS) return FAILURE;
            if (mp_add_small(x, bit) != SUCCESS) return FAILURE;
        }

        /* Two's complement detection */
        {
            size_t bitwidth = strlen(start);
            if (start[0] == '1') {
                mp_int pow2;
                mp_int tmp;
                mp_init(&pow2);
                if (mp_reserve(&pow2, 1) != SUCCESS) { mp_free(&pow2); return FAILURE; }
                pow2.sign = 1; pow2.digits[0] = 1; pow2.length = 1;

                for (i = 0; i < bitwidth; ++i)
                    mp_mul_small(&pow2, 2U);
                mp_init(&tmp);
                mp_sub(&tmp, x, &pow2);
                mp_copy(x, &tmp);
                mp_free(&pow2);
                mp_free(&tmp);
            }
        }

        if (explicit_sign == -1)
            x->sign = -x->sign;
    }

    return SUCCESS;
}

/* ---------- Hexadecimal (with implicit sign bit handling) ---------- */
int mp_from_str_hex(mp_int *x, const char *str) {
    size_t i;
    size_t len;
    const char *start;
    int first_digit;
    size_t bitwidth;
    mp_int pow2;
    mp_int tmp;
    int val;
    const char *p;

    if (!x || !str) return FAILURE;
    mp_free(x);
    mp_init(x);

    p = str;
    while (isspace((unsigned char)*p)) p++;

    {
        int explicit_sign = +1;
        if (*p == '-') { explicit_sign = -1; p++; }
        else if (*p == '+') { p++; }

        if (p[0] == '0' && (p[1] == 'x' || p[1] == 'X')) p += 2;

        start = p;
        while (*p == '0') p++;
        len = strlen(p);

        if (!isxdigit((unsigned char)*p)) {
            x->sign = 0; x->length = 0;
            return SUCCESS;
        }

        if (mp_reserve(x, 1) != SUCCESS) return FAILURE;
        x->digits[0] = 0;
        x->length = 1;

        for (; *p; ++p) {
            if (isdigit((unsigned char)*p)) val = *p - '0';
            else if (isupper((unsigned char)*p)) val = *p - 'A' + 10;
            else val = *p - 'a' + 10;
            if (!isxdigit((unsigned char)*p)) break;
            if (mp_mul_small(x, 16U) != SUCCESS) return FAILURE;
            if (mp_add_small(x, (unsigned int)val) != SUCCESS) return FAILURE;
        }

        /* Two's complement detection */
        first_digit = 0;
        if (isxdigit((unsigned char)start[0])) {
            if (isdigit((unsigned char)start[0])) first_digit = start[0] - '0';
            else if (isupper((unsigned char)start[0])) first_digit = start[0] - 'A' + 10;
            else first_digit = start[0] - 'a' + 10;
        }
        bitwidth = len * 4;
        if (first_digit & 0x8) {
            mp_init(&pow2);
            if (mp_reserve(&pow2, 1) != SUCCESS) { mp_free(&pow2); return FAILURE; }
            pow2.sign = 1; pow2.digits[0] = 1; pow2.length = 1;
            for (i = 0; i < bitwidth; ++i)
                mp_mul_small(&pow2, 2U);
            mp_init(&tmp);
            mp_sub(&tmp, x, &pow2);
            mp_copy(x, &tmp);
            mp_free(&pow2);
            mp_free(&tmp);
        }

        if (explicit_sign == -1)
            x->sign = -x->sign;
    }

    return SUCCESS;
}


/* ---------------------------------------------------- Printing ---------------------------------------------------- */

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

/* Helper: copy absolute value of x into dst (dst initialized by caller) */
static int mp_abs_copy(mp_int *dst, const mp_int *x) {
    if (!dst || !x) return FAILURE;
    if (mp_copy(dst, x) != SUCCESS) return FAILURE;
    if (dst->length == 0) { dst->sign = 0; return SUCCESS; }
    dst->sign = +1;
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

/* ------------------------------------------------- File processing ----------------------------------------------- */

/* parses the input file name, checks validity and processes each line as a decimal number */
int process_file(char str[]) {
    FILE *f;
    char line[1024];
    mp_int val;

    if (!str) return FAILURE;
    f = fopen(str, "r");
    if (!f) return FAILURE;
    mp_init(&val);

    while (fgets(line, sizeof(line), f)) {
        /* trim newline */
        size_t L = strlen(line);
        while (L > 0 && (line[L - 1] == '\n' || line[L - 1] == '\r')) { line[--L] = '\0'; }

        if (L == 0) continue;
        /* detect hex or bin prefix */
        if (line[0] == '0' && (line[1] == 'x' || line[1] == 'X')) {
            if (mp_from_str_hex(&val, line) == SUCCESS) {
                mp_print_hex(&val); putchar('\n');
            } else {
                printf("parse error: %s\n", line);
            }
        } else if (line[0] == '0' && (line[1] == 'b' || line[1] == 'B')) {
            if (mp_from_str_bin(&val, line) == SUCCESS) {
                mp_print_bin(&val); putchar('\n');
            } else {
                printf("parse error: %s\n", line);
            }
        } else {
            if (mp_from_str_dec(&val, line) == SUCCESS) {
                mp_print_dec(&val); putchar('\n');
            } else {
                printf("parse error: %s\n", line);
            }
        }
    }

    mp_free(&val);
    fclose(f);
    return SUCCESS;
}

/* -------------------- Expression parsing & evaluation helpers -------------------- */

/* Trim in-place leading and trailing whitespace */
static void trim_inplace(char *s) {
    char *p = s;
    size_t len;
    while (*p && isspace((unsigned char)*p)) p++;
    if (p != s) memmove(s, p, strlen(p) + 1);
    len = strlen(s);
    while (len > 0 && isspace((unsigned char)s[len - 1])) s[--len] = '\0';
}

/*
  split_expr:
    - input: a nul-terminated string (may contain leading/trailing spaces)
    - outputs: lhs (buffer), rhs (buffer), op (single char)
    - returns 1 if an operator was found and split succeeded, 0 if no binary operator detected
  Notes:
    - it accepts forms like "123+456", "  -2   *  0x10", "0b101 - -0b1"
    - does NOT implement parenthesis or operator precedence
*/
int split_expr(const char *input, char *lhs, size_t lhs_cap, char *op, char *rhs, size_t rhs_cap) {
    char buf[1024];
    size_t i, n;
    size_t left_len;
    size_t right_len;

    if (!input || !lhs || !rhs || !op) return 0;
    /* copy and trim */
    strncpy(buf, input, sizeof(buf) - 1);
    buf[sizeof(buf) - 1] = '\0';
    trim_inplace(buf);
    if (buf[0] == '\0') return 0;

    /* find first operator + - * / that is not the leading sign of the first operand */
    n = strlen(buf);
    i = 0;
    {
        int found = 0;
        for (i = 0; i < n; ++i) {
            char ch = buf[i];
            if (ch == '+' || ch == '-' || ch == '*' || ch == '/') {
                if (i == 0) {
                    /* This is sign of the first operand -> skip */
                    continue;
                }
                /* We accept this as operator */
                found = 1;
                break;
            }
        }
        if (!found) return 0;
    }

    /* split */
    left_len = i;
    right_len = (n > i + 1) ? n - (i + 1) : 0;
    /* copy left */
    if (left_len >= lhs_cap) left_len = lhs_cap - 1;
    strncpy(lhs, buf, left_len);
    lhs[left_len] = '\0';
    trim_inplace(lhs);
    /* copy operator */
    *op = buf[i];
    /* copy right */
    if (right_len >= rhs_cap) right_len = rhs_cap - 1;
    strncpy(rhs, buf + i + 1, right_len);
    rhs[right_len] = '\0';
    trim_inplace(rhs);

    /* If rhs is empty -> not a valid binary expression */
    if (rhs[0] == '\0') return 0;
    /* lhs must not be empty either */
    if (lhs[0] == '\0') return 0;

    return 1;
}

/* parse_operand_to_mp: parses a single operand string into dst (auto detects base) */
int parse_operand_to_mp(mp_int *dst, const char *s) {
    char buf[1024];
    if (!dst || !s) return FAILURE;
    /* trim into local buffer */
    strncpy(buf, s, sizeof(buf) - 1);
    buf[sizeof(buf) - 1] = '\0';
    trim_inplace(buf);
    if (buf[0] == '\0') return FAILURE;

    /* detect base and call appropriate parser */
    if (buf[0] == '0' && (buf[1] == 'x' || buf[1] == 'X')) {
        return mp_from_str_hex(dst, buf);
    } else if (buf[0] == '0' && (buf[1] == 'b' || buf[1] == 'B')) {
        return mp_from_str_bin(dst, buf);
    } else {
        /* default decimal (handles leading + / -) */
        return mp_from_str_dec(dst, buf);
    }
}

/* ---------------------------------------------------- Main ------------------------------------------------------- */

int main(int argc, char *argv[]) {
    /* Format selection */
    int (*mp_print[])(mp_int *) = { mp_print_dec, mp_print_bin, mp_print_hex };
    int print_format = DEC;

    mp_int a, b, c;
    char input_file[300];
    char buf[1024];
    TokenList infix, postfix;
    size_t L;
    char *p;

    if (argc > 2) {
        printf("Error: invalid input\nUsage: calc.exe [<input-file>]\n");
        return EXIT_FAILURE;
    }
    if (argc == 2) {
        strncpy(input_file, argv[1], sizeof(input_file) - 1);
        input_file[sizeof(input_file) - 1] = '\0';
        return process_file(input_file) != FAILURE ? EXIT_SUCCESS : EXIT_FAILURE;
    }

    mp_init(&a);
    mp_init(&b);
    mp_init(&c);

    while (1) {
        printf("> ");
        if (!fgets(buf, sizeof(buf), stdin)) break;

        /* trim newline */
        L = strlen(buf);
        while (L > 0 && (buf[L - 1] == '\n' || buf[L - 1] == '\r')) { buf[--L] = '\0'; }

        /* Trim leading spaces */
        p = buf;
        while (*p && isspace((unsigned char)*p)) p++;

        if (p[0] == '\0') continue; /* empty line -> prompt again */

        /* exit command */
        if (strcmp(p, "exit") == 0 || strcmp(p, "quit") == 0) break;

        if (strcmp(p, "bin") == 0) {
            print_format = BIN;
            printf("bin\n");
            continue;
        }
        if (strcmp(p, "hex") == 0) {
            print_format = HEX;
            printf("hex\n");
            continue;
        }
        if (strcmp(p, "dec") == 0) {
            print_format = DEC;
            printf("dec\n");
            continue;
        }
        if (strcmp(p, "out") == 0) {
            switch(print_format) {
                case 0: printf("dec\n"); break;
                case 1: printf("bin\n"); break;
                case 2: printf("hex\n"); break;
            }
            continue;
        }

        tokenize(p, &infix);
        to_postfix(&infix, &postfix);
        mp_free(&c);
        mp_init(&c);
        if (eval_postfix(&postfix, &c) != SUCCESS) {
            printf("parse or eval error\n");
            continue;
        }

        mp_print[print_format](&c);
        putchar('\n');


        /* Not an expression -> treat as a single literal to parse and echo */
        mp_free(&c);
        mp_init(&c);
        if (p[0] == '0' && (p[1] == 'x' || p[1] == 'X')) {
            if (mp_from_str_hex(&c, p) != SUCCESS) { 
                printf("parse error\n"); 
                continue; 
            }
        } else if (p[0] == '0' && (p[1] == 'b' || p[1] == 'B')) {
            if (mp_from_str_bin(&c, p) != SUCCESS) { 
                printf("parse error\n"); 
                continue; 
            }
        } else {
            if (mp_from_str_dec(&c, p) != SUCCESS) { 
                printf("parse error\n"); 
                continue; 
            }
        }
        putchar('\n');
    }

    mp_free(&a);
    mp_free(&b);
    mp_free(&c);
    return EXIT_SUCCESS;
}
