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

/* Multiplication algorithm selection based on MpInt size */
#define NAIVE_THRESHOLD 16

/* Multiple precision int representation */
typedef struct {
    int sign;
    unsigned int *digits;
    size_t length, capacity; /* number of limbs used (current, max) */
} MpInt;

/* Forward declarations (to avoid implicit-declaration warnings) */
int parse_operand_to_mp(MpInt *dst, const char *s);
int mp_from_str_bin(MpInt *x, const char *str);
int mp_from_str_hex(MpInt *x, const char *str);
int mp_print_dec(MpInt *x);
int mp_print_bin(MpInt *x);
int mp_print_hex(MpInt *x);

/* ---------------------------------- Memory management helpers -------------------------------------- */

/* Initializes MpInt x to empty state */
int mp_init(MpInt *x) {
    if (!x) return FAILURE;
    x->sign = 0;
    x->length = 0;
    x->digits = NULL;
    x->capacity = 0;
    return SUCCESS;
}

/* Frees memory used by the digits and reinitializes the MpInt */
int mp_free(MpInt *x) {
    if (!x) return SUCCESS;
    free(x->digits);
    /* reinitialize */
    x->digits = NULL;
    x->length = 0;
    x->capacity = 0;
    x->sign = 0;
    return SUCCESS;
}

/* Reallocates memory for the new desired size of a MpInt */
int mp_reserve(MpInt *x, size_t capacity) {
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

/* Copy MpInt src -> dst (dst must be initialized with mp_init) */
int mp_copy(MpInt *dst, const MpInt *src) {
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

/* Compare absolute values of two MpInts */
/* Returns 1 if a is greater, -1 if b is greater, otherwise 0 */
int mp_cmp_abs(const MpInt *a, const MpInt *b) {
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

/* Adds a MpInt and a regular unsigned int */
int mp_add_small(MpInt *x, unsigned int a)
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

/* Multiplies a MpInt by a regular unsigned int */
int mp_mul_small(MpInt *x, unsigned int m)
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

/* Divides a MpInt by a regular unsigned int */
int mp_div_small(MpInt *result, const MpInt *a, unsigned int divisor, unsigned int *remainder) {
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

/* Unsigned addition of MpInt (absolute) */
int mp_add_abs(MpInt *result, const MpInt *a, const MpInt *b) {
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
int mp_sub_abs(MpInt *result, const MpInt *a, const MpInt *b) {
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
int mp_add(MpInt *result, const MpInt *a, const MpInt *b) {
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

/* Subtracts MpInt b from MpInt a */
int mp_sub(MpInt *result, const MpInt *a, const MpInt *b) {
    MpInt tmp;

    if (!result || !a || !b) return FAILURE;
    mp_init(&tmp);
    if (mp_copy(&tmp, b) != SUCCESS) { mp_free(&tmp); return FAILURE; }
    tmp.sign = -tmp.sign;
    mp_add(result, a, &tmp);
    mp_free(&tmp);
    return SUCCESS;
}

/* Naively multiplies two MpInts, suitable for smaller numbers*/
int mp_mul_naive(MpInt *result, const MpInt *a, const MpInt *b) {
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
int mp_shift_left_words(MpInt *x, size_t k) {
    if (!x) return FAILURE;
    if (x->sign == 0 || k == 0) return SUCCESS;
    if (mp_reserve(x, x->length + k) != SUCCESS)
        return FAILURE;
    memmove(x->digits + k, x->digits, x->length * sizeof(unsigned int));
    memset(x->digits, 0, k * sizeof(unsigned int));
    x->length += k;
    return SUCCESS;
}

void mp_split(const MpInt *src, MpInt *low, MpInt *high, size_t m) {
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
int mp_karatsuba_mul(MpInt *result, const MpInt *a, const MpInt *b)
{
    size_t n, m;
    MpInt x0, x1, y0, y1, z0, z1, z2, tmp1, tmp2, res_tmp;

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

/* Multiplies two MpInts (hybrid) */
int mp_mul(MpInt *result, const MpInt *a, const MpInt *b) {
    if (!result || !a || !b) return FAILURE;
    if (a->length < NAIVE_THRESHOLD || b->length < NAIVE_THRESHOLD)
        return mp_mul_naive(result, a, b);
    else
        return mp_karatsuba_mul(result, a, b);
}

/* Divides MpInt b by MpInt a */
/* Note: full long division is not implemented; only small divisor supported via mp_div_small */
int mp_div(MpInt *result, const MpInt *a, const MpInt *b) {
    if (!result || !a || !b) return FAILURE;
    /* If divisor is a single limb, use mp_div_small */
    if (b->length == 1) {
        unsigned int rem;
        return mp_div_small(result, a, b->digits[0], &rem);
    }
    return FAILURE; /* not implemented */
}

int mp_pow(MpInt *r, const MpInt *a, const MpInt *b) {
    /* TODO: implement exponentiation */
    printf("[warning] exponentiation not yet implemented\n");
    mp_copy(r, a);
    return SUCCESS;
}

int mp_fact(MpInt *r, const MpInt *a) {
    /* TODO: implement factorial */
    printf("[warning] factorial not yet implemented\n");
    mp_copy(r, a);
    return SUCCESS;
}

int mp_mod(MpInt *r, const MpInt *a, const MpInt *b) {
    /* TODO: implement modulo */
    printf("[warning] modulo not yet implemented\n");
    mp_copy(r, a);
    return SUCCESS;
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

static void tokenize(const char *expr, TokenList *out) {
    const char *p = expr;
    int i;
    char buf[64];

    out->count = 0;

    while (*p) {
        while (isspace((unsigned char)*p)) p++;
        if (*p == '\0') break;

        if (is_operator(*p) || *p == '(' || *p == ')') {
            out->tokens[out->count][0] = *p;
            out->tokens[out->count][1] = '\0';
            out->count++;
            p++;
        } else {
            /* number or identifier */
            i = 0;
            while (*p && !isspace((unsigned char)*p)
                   && !is_operator(*p)
                   && *p != '(' && *p != ')') {
                if (i < 63) buf[i++] = *p;
                p++;
            }
            buf[i] = '\0';
            /* copy into token list */
            strncpy(out->tokens[out->count], buf, 64);
            out->tokens[out->count][63] = '\0';
            out->count++;
        }
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
        if (is_operator(tok[0]) && tok[1] == '\0') {
            char op = tok[0];
            if (is_unary(op)) {
                while (sp > 0 && get_precedence(stack[sp-1]) >= get_precedence(op)) {
                    postfix->tokens[postfix->count][0] = stack[--sp];
                    postfix->tokens[postfix->count][1] = '\0';
                    postfix->count++;
                }
                stack[sp++] = op;
            } else {
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
            if (sp > 0 && stack[sp-1] == '(') sp--;
        } else {
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
static int eval_postfix(TokenList *postfix, MpInt *result) {
    MpInt stack[128];
    int i;
    int sp = 0;

    for (i = 0; i < postfix->count; ++i) {
        char *tok = postfix->tokens[i];
        if (is_operator(tok[0]) && tok[1] == '\0') {
            char op = tok[0];
            if (is_unary(op)) {
                if (sp < 1) return FAILURE;
                mp_init(&stack[sp]); /* ensure slot exists */
                /* pop one operand into a */
                mp_copy(&stack[sp], &stack[sp-1]);
                mp_free(&stack[sp-1]);
                /* evaluate unary */
                if (op == '!') {
                    mp_fact(&stack[sp], &stack[sp]);
                }
                sp = sp + 0; /* stack size unchanged: replaced top-1 with result in place */
            } else {
                if (sp < 2) return FAILURE;
                /* pop b then a; evaluate into stack[sp-2] */
                mp_init(&stack[sp-2]); /* ensure slots valid */
                /* we already have copies in stack[sp-2], stack[sp-1] */
                /* perform operation and store into stack[sp-2] */
                if (op == '+') {
                    mp_add(&stack[sp-2], &stack[sp-2], &stack[sp-1]);
                } else if (op == '-') {
                    mp_sub(&stack[sp-2], &stack[sp-2], &stack[sp-1]);
                } else if (op == '*') {
                    mp_mul(&stack[sp-2], &stack[sp-2], &stack[sp-1]);
                } else if (op == '/') {
                    mp_div(&stack[sp-2], &stack[sp-2], &stack[sp-1]);
                } else if (op == '%') {
                    mp_mod(&stack[sp-2], &stack[sp-2], &stack[sp-1]);
                } else if (op == '^') {
                    mp_pow(&stack[sp-2], &stack[sp-2], &stack[sp-1]);
                }
                /* free the slot stack[sp-1] and reduce stack size by 1 */
                mp_free(&stack[sp-1]);
                sp--;
            }
        } else {
            /* operand: parse and push */
            mp_init(&stack[sp]);
            if (parse_operand_to_mp(&stack[sp], tok) != SUCCESS)
                return FAILURE;
            sp++;
        }
    }
    if (sp != 1) return FAILURE;
    mp_copy(result, &stack[0]);
    for (i = 0; i < sp; ++i) mp_free(&stack[i]);
    return SUCCESS;
}

/* ---------- Decimal ---------- */
int mp_from_str_dec(MpInt *x, const char *str)
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
int mp_from_str_bin(MpInt *x, const char *str) {
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
        /* size_t len = strlen(p); -- len not used; removed */

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
                MpInt pow2;
                MpInt tmp;
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
int mp_from_str_hex(MpInt *x, const char *str) {
    size_t i;
    size_t len;
    const char *start;
    int first_digit;
    size_t bitwidth;
    MpInt pow2;
    MpInt tmp;
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
static int mp_make_pow2(MpInt *pow2, size_t k) {
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
static int mp_abs_copy(MpInt *dst, const MpInt *x) {
    if (!dst || !x) return FAILURE;
    if (mp_copy(dst, x) != SUCCESS) return FAILURE;
    if (dst->length == 0) { dst->sign = 0; return SUCCESS; }
    dst->sign = +1;
    return SUCCESS;
}

/* Prints the MpInt @x in decimal format (unchanged behavior) */
int mp_print_dec(MpInt *x) {
    MpInt tmp, q;
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

/* Prints the MpInt @x in binary (two's complement if negative), minimal width */
int mp_print_bin(MpInt *x) {
    MpInt tmp, q;
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
        MpInt absx;
        size_t w;
        MpInt pow;
        int cmp;
        MpInt poww, val;
        MpInt q2;
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

/* Prints the MpInt @x in hexadecimal (two's complement if negative), minimal nibble width */
int mp_print_hex(MpInt *x) {
    MpInt tmp;
    size_t approx;
    char *buf;
    size_t idx;
    unsigned int rem;
    MpInt q2;

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
        MpInt absx;
        size_t w;
        MpInt pow;
        int cmp;
        size_t nibbles, bits;
        MpInt powb, val;
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
    MpInt val;

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
int parse_operand_to_mp(MpInt *dst, const char *s) {
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
    int (*mp_print[])(MpInt *) = { mp_print_dec, mp_print_bin, mp_print_hex };
    int print_format = DEC;

    MpInt a, b, c;
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
        printf(">>> ");
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
            continue;
        }
        if (strcmp(p, "hex") == 0) {
            print_format = HEX;
            continue;
        }
        if (strcmp(p, "dec") == 0) {
            print_format = DEC;
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
            if (mp_from_str_hex(&c, p) != SUCCESS) { printf("parse error\n"); continue; }
        } else if (p[0] == '0' && (p[1] == 'b' || p[1] == 'B')) {
            if (mp_from_str_bin(&c, p) != SUCCESS) { printf("parse error\n"); continue; }
        } else {
            if (mp_from_str_dec(&c, p) != SUCCESS) { printf("parse error\n"); continue; }
        }
        /*mp_print[print_format](&c);*/
        putchar('\n');
    }

    mp_free(&a);
    mp_free(&b);
    mp_free(&c);
    return EXIT_SUCCESS;
}
