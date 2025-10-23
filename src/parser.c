#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "parser.h"
#include "mp_int.h"
#include "mp_print.h"     /* if process_file prints results */
#include "exp.h"
#include "fact.h"


static int is_operator(char c) {
    return strchr("+-*/%^!", c) != NULL;
}

/* Trim in-place leading and trailing whitespace */
static void trim_inplace(char *s) {
    char *p = s;
    size_t len;
    while (*p && isspace((unsigned char)*p)) p++;
    if (p != s) memmove(s, p, strlen(p) + 1);
    len = strlen(s);
    while (len > 0 && isspace((unsigned char)s[len - 1])) s[--len] = '\0';
}

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
void tokenize(const char *expr, TokenList *out) {
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
void to_postfix(TokenList *infix, TokenList *postfix) {
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
int eval_postfix(TokenList *postfix, mp_int *result) {
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
/* parse_operand_to_mp: parses a single operand string into dst (auto detects base)
   Accepts operands that may start with an optional + or - sign, followed by
   decimal, 0b... binary, or 0x... hex. */
int parse_operand_to_mp(mp_int *dst, const char *s) {
    char buf[1024];
    const char *q;
    int status;

    if (!dst || !s) return FAILURE;

    /* trim into local buffer */
    strncpy(buf, s, sizeof(buf) - 1);
    buf[sizeof(buf) - 1] = '\0';
    trim_inplace(buf);
    if (buf[0] == '\0') return FAILURE;

    /* Look past optional leading sign to detect 0x/0b prefixes */
    q = buf;
    if (*q == '+' || *q == '-') q++;

    /* Now test the prefix at q */
    if (q[0] == '0' && (q[1] == 'x' || q[1] == 'X')) {
        status = mp_from_str_hex(dst, buf);
        return status;
    } else if (q[0] == '0' && (q[1] == 'b' || q[1] == 'B')) {
        status = mp_from_str_bin(dst, buf);
        return status;
    } else {
        /* default decimal (handles leading + / -) */
        status = mp_from_str_dec(dst, buf);
        return status;
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
    /* ANSI C90 declarations at top */
    size_t i;
    const char *p;
    const char *digits;
    size_t digits_len;
    unsigned int bit;
    int explicit_sign = +1;

    if (!x || !str) return FAILURE;
    mp_free(x);
    mp_init(x);

    /* skip leading whitespace */
    p = str;
    while (isspace((unsigned char)*p)) p++;

    /* optional explicit sign */
    if (*p == '-') { explicit_sign = -1; p++; }
    else if (*p == '+') { p++; }

    /* optional 0b/0B prefix */
    if (p[0] == '0' && (p[1] == 'b' || p[1] == 'B')) p += 2;

    /* collect contiguous binary digits */
    digits = p;
    digits_len = 0;
    while (p[digits_len] == '0' || p[digits_len] == '1') digits_len++;

    if (digits_len == 0) {
        /* no binary digits -> treat as zero/invalid */
        x->sign = 0;
        x->length = 0;
        return SUCCESS;
    }

    /* initialize x = 0 */
    if (mp_reserve(x, 1) != SUCCESS) return FAILURE;
    x->digits[0] = 0;
    x->length = 1;
    x->sign = 1; /* magnitude building */

    /* parse digits left-to-right */
    for (i = 0; i < digits_len; ++i) {
        bit = (unsigned int)(digits[i] - '0');
        if (mp_mul_small(x, 2U) != SUCCESS) return FAILURE;
        if (mp_add_small(x, bit) != SUCCESS) return FAILURE;
    }

    /* two's complement detection: if the highest provided bit is 1, interpret as negative value
       (value = magnitude - 2^digits_len). */
    if (digits[0] == '1') {
        mp_int pow2;
        mp_int tmp;

        mp_init(&pow2);
        mp_init(&tmp);

        /* pow2 = 2^digits_len */
        if (mp_reserve(&pow2, 1) != SUCCESS) { mp_free(&pow2); mp_free(&tmp); return FAILURE; }
        pow2.digits[0] = 1;
        pow2.length = 1;
        pow2.sign = 1;
        for (i = 0; i < digits_len; ++i) {
            if (mp_mul_small(&pow2, 2U) != SUCCESS) { mp_free(&pow2); mp_free(&tmp); return FAILURE; }
        }

        /* tmp = x - pow2  (may be negative) */
        if (mp_sub(&tmp, x, &pow2) != SUCCESS) { mp_free(&pow2); mp_free(&tmp); return FAILURE; }

        /* copy tmp back into x */
        mp_free(x);
        mp_init(x);
        if (mp_copy(x, &tmp) != SUCCESS) { mp_free(&pow2); mp_free(&tmp); return FAILURE; }

        mp_free(&pow2);
        mp_free(&tmp);
    }

    /* apply explicit sign: explicit '-' flips the computed sign (even if two's-complement produced negative) */
    if (explicit_sign == -1) {
        x->sign = -x->sign;
    }

    return SUCCESS;
}

/* ---------- Hexadecimal (with implicit sign bit handling) ---------- */
int mp_from_str_hex(mp_int *x, const char *str) {
    /* ANSI C90 declarations at top */
    size_t i;
    const char *p;
    const char *digits;
    size_t digits_len;
    int val;
    int explicit_sign = +1;
    int first_nibble;

    if (!x || !str) return FAILURE;
    mp_free(x);
    mp_init(x);

    /* skip leading whitespace */
    p = str;
    while (isspace((unsigned char)*p)) p++;

    /* optional explicit sign */
    if (*p == '-') { explicit_sign = -1; p++; }
    else if (*p == '+') { p++; }

    /* optional 0x/0X prefix */
    if (p[0] == '0' && (p[1] == 'x' || p[1] == 'X')) p += 2;

    /* collect contiguous hex digits */
    digits = p;
    digits_len = 0;
    while (isxdigit((unsigned char)p[digits_len])) digits_len++;

    if (digits_len == 0) {
        x->sign = 0;
        x->length = 0;
        return SUCCESS;
    }

    /* initialize x = 0 */
    if (mp_reserve(x, 1) != SUCCESS) return FAILURE;
    x->digits[0] = 0;
    x->length = 1;
    x->sign = 1;

    for (i = 0; i < digits_len; ++i) {
        char c = digits[i];
        if (isdigit((unsigned char)c)) val = c - '0';
        else if (isupper((unsigned char)c)) val = c - 'A' + 10;
        else val = c - 'a' + 10;

        if (mp_mul_small(x, 16U) != SUCCESS) return FAILURE;
        if (mp_add_small(x, (unsigned int)val) != SUCCESS) return FAILURE;
    }

    /* detect sign-bit from first nibble */
    first_nibble = 0;
    if (isxdigit((unsigned char)digits[0])) {
        char c = digits[0];
        if (isdigit((unsigned char)c)) first_nibble = c - '0';
        else if (isupper((unsigned char)c)) first_nibble = c - 'A' + 10;
        else first_nibble = c - 'a' + 10;
    }

    if (first_nibble & 0x8) {
        /* interpret as two's complement negative: x = x - 2^(4*digits_len) */
        mp_int pow2;
        mp_int tmp;
        mp_init(&pow2);
        mp_init(&tmp);

        if (mp_reserve(&pow2, 1) != SUCCESS) { mp_free(&pow2); mp_free(&tmp); return FAILURE; }
        pow2.digits[0] = 1;
        pow2.length = 1;
        pow2.sign = 1;

        /* pow2 = 2^(4 * digits_len) */
        for (i = 0; i < (4 * digits_len); ++i) {
            if (mp_mul_small(&pow2, 2U) != SUCCESS) { mp_free(&pow2); mp_free(&tmp); return FAILURE; }
        }

        if (mp_sub(&tmp, x, &pow2) != SUCCESS) { mp_free(&pow2); mp_free(&tmp); return FAILURE; }

        mp_free(x);
        mp_init(x);
        if (mp_copy(x, &tmp) != SUCCESS) { mp_free(&pow2); mp_free(&tmp); return FAILURE; }

        mp_free(&pow2);
        mp_free(&tmp);
    }

    /* apply explicit sign if present */
    if (explicit_sign == -1) x->sign = -x->sign;

    return SUCCESS;
}
