/*
 * parser.c -- tokenizer, shunting-yard, postfix evaluator, file processing
 *
 * Important behaviors explained:
 *  - Tokenizer: recognizes operands (optionally signed literals with 0x/0b prefixes),
 *    single-char operators (+ - * / % ^ !) and parentheses.
 *  - Unary rules:
 *      * A leading '+'/'-' that appears at start or immediately after '('/another operator
 *        is generally treated as a sign for a following numeric literal.
 *      * Exception: when a signed literal would be immediately followed by a postfix '!' (factorial),
 *        the tokenizer will *not* attach the sign to the literal; it rewrites e.g. "-49!" as:
 *           "0" "-" "49" "!"
 *        This ensures factorial has higher precedence than unary minus (so -49! == -(49!)).
 *      * For "-(...)" forms the tokenizer rewrites to "0" "-" "(" ... to allow parsing "-(...)".
 *
 *  - Shunting-yard: supports postfix unary '!' with highest precedence (per op_table).
 *  - Postfix evaluator: uses an mp_int stack (fixed depth 128). That limit is conservative but
 *    should be adjusted if you expect extremely deep expression trees.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "parser.h"
#include "mp_int.h"
#include "mp_print.h"
#include "exp.h"
#include "fact.h"
#include "error.h"

/* Helper: returns true if the character is one of the single-char operators we support. */
static int is_operator(char c) {
    return strchr("+-*/%^!", c) != NULL;
}

/* Trim in-place leading and trailing whitespace.
   Implemented carefully to use only C90 features and avoid temporary allocations. */
static void trim_inplace(char *s) {
    char *p = s;
    size_t len;
    while (*p && isspace((unsigned char)*p)) p++;
    if (p != s) memmove(s, p, strlen(p) + 1); /* shift left if needed */
    len = strlen(s);
    while (len > 0 && isspace((unsigned char)s[len - 1])) s[--len] = '\0';
}

/* Operator metadata table:
   - op: character
   - precedence: higher = evaluated before
   - right_assoc: 1 if operator is right-associative (e.g., exponentiation '^')
   - unary: 1 if operator is unary (we use this primarily for postfix factorial '!')
*/
static OpInfo op_table[] = {
    {'!', 5, 1, 1},  /* factorial, unary postfix (highest precedence) */
    {'^', 4, 1, 0},  /* exponentiation (right associative) */
    {'*', 3, 0, 0},
    {'/', 3, 0, 0},
    {'%', 3, 0, 0},
    {'+', 2, 0, 0},
    {'-', 2, 0, 0},
    {'(', 1, 0, 0},
    {0, 0, 0, 0}
};

/* Accessor helpers for op_table */
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
      * a single-character operator: "+", "-", "*", "/", "%", "^", "!"
      * "(" or ")"
      * an operand string (may include a leading + or - sign and binary/hex prefixes)
  - Leading + / - is considered part of a literal when it appears at expression start
    or immediately after '(' or another operator — *except* when the literal is followed
    immediately (possibly after spaces) by a postfix '!' token. In that case the sign
    must be treated as a separate binary operator so that factorial binds more tightly.
*/
void tokenize(const char *expr, TokenList *out) {
    const char *p;
    char buf[64];
    int i;
    int prev_was_operand; /* 0 = start/after-operator/'(', 1 = after-operand or ')' */

    out->count = 0;
    prev_was_operand = 0; /* treat leading +/ - as sign if they appear */

    p = expr;
    while (*p) {
        /* skip whitespace */
        while (isspace((unsigned char)*p)) p++;
        if (*p == '\0') break;

        /* safety: do not exceed token storage */
        if (out->count >= 256) break;

        /* Handle leading + / - that may be a sign or operator.
           Several cases handled:
             - "-("  => rewrite as "0" "-" "(" so unary minus before parenthesis works.
             - "-number" where 'number' not followed immediately by '!' => treat as signed literal.
             - "-number !" => treat sign as separate binary operator: "0" "-" "number" "!"
        */
        if ((*p == '+' || *p == '-') && !prev_was_operand) {
            const char *q = p + 1;
            /* skip spaces when checking for parentheses case */
            while (isspace((unsigned char)*q)) q++;
            if (*q == '(') {
                /* rewrite "-(" as tokens: "0" "-" "(" */
                strncpy(out->tokens[out->count], "0", 64);
                out->tokens[out->count][63] = '\0';
                out->count++;

                out->tokens[out->count][0] = *p;
                out->tokens[out->count][1] = '\0';
                out->count++;

                /* consume sign; '(' will be handled in next loop iteration */
                p++;
                continue;
            }

            /* If next char looks like a number token, check whether the signed literal
               would be immediately followed by a postfix '!' — if so, do not absorb sign
               into the literal (we want factorial to bind before unary minus). */
            q = p + 1;
            if (*q != '\0' && !isspace((unsigned char)*q)) {
                const char *r = q;
                /* scan number body */
                while (*r && !isspace((unsigned char)*r) &&
                    !is_operator(*r) && *r != '(' && *r != ')') { r++; }

                /* scan past spaces to see if '!' follows */
                {
                    const char *s = r;
                    while (isspace((unsigned char)*s)) s++;
                    if (*s == '!') {
                        /* rewrite "-49!" as "0" "-" "49" "!" */
                        strncpy(out->tokens[out->count], "0", 64);
                        out->tokens[out->count][63] = '\0';
                        out->count++;

                        out->tokens[out->count][0] = *p; /* '+' or '-' */
                        out->tokens[out->count][1] = '\0';
                        out->count++;

                        p++; /* consume only the sign; number will be tokenized in next pass */
                        continue;
                    }
                }

                /* Otherwise collect a signed literal token (e.g. "-123", "+0x1f") */
                i = 0;
                if (i < 63) buf[i++] = *p; /* leading sign */
                p++;
                while (*p && !isspace((unsigned char)*p) &&
                    !is_operator(*p) && *p != '(' && *p != ')') {
                    if (i < 63) buf[i++] = *p;
                    p++;
                }
                buf[i] = '\0';
                strncpy(out->tokens[out->count], buf, 64);
                out->tokens[out->count][63] = '\0';
                out->count++;
                prev_was_operand = 1;
                continue;
            }
            /* otherwise fallthrough to emit sign as operator token */
        }

        /* Single-char operator or parentheses */
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

        /* Otherwise parse a plain operand token (no leading sign) */
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

/*
 * to_postfix -- convert infix token list to postfix using the shunting-yard algorithm
 *
 * Notes:
 *  - Supports binary operators and a postfix unary operator '!' (factorial).
 *  - Uses operator precedence and associativity from op_table.
 *  - The stack capacity is 256 characters (sufficient for single-char operators and parentheses).
 */
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
                /* postfix unary operators behave like operators with precedence:
                   pop while top has >= precedence, then push unary operator. */
                while (sp > 0 && get_precedence(stack[sp-1]) >= get_precedence(op)) {
                    postfix->tokens[postfix->count][0] = stack[--sp];
                    postfix->tokens[postfix->count][1] = '\0';
                    postfix->count++;
                }
                stack[sp++] = op;
            } else {
                /* binary operator: pop according to precedence and associativity */
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
            /* pop until '(' */
            while (sp > 0 && stack[sp-1] != '(') {
                postfix->tokens[postfix->count][0] = stack[--sp];
                postfix->tokens[postfix->count][1] = '\0';
                postfix->count++;
            }
            if (sp > 0 && stack[sp-1] == '(') sp--; /* pop '(' */
        } else {
            /* operand (possibly multi-char) -> emit to postfix */
            strncpy(postfix->tokens[postfix->count], tok, 64);
            postfix->tokens[postfix->count][63] = '\0';
            postfix->count++;
        }
    }

    /* flush operator stack */
    while (sp > 0) {
        postfix->tokens[postfix->count][0] = stack[--sp];
        postfix->tokens[postfix->count][1] = '\0';
        postfix->count++;
    }
}

/*
 * eval_postfix -- evaluate a postfix token list and produce an mp_int result.
 *
 * Stack usage:
 *   - Uses mp_int stack[128]. Each mp_int is small (contains pointers) but the
 *     stack depth limit is 128 items; enlarge if you expect deep nesting.
 *
 * Operator semantics:
 *   - Binary: + - * / % ^  (power uses mp_pow)
 *   - Postfix unary: ! (factorial) uses mp_fact
 */
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

                    /* pop a from stack */
                    mp_copy(&a, &stack[--sp]);
                    mp_free(&stack[sp]);

                    /* apply unary operator */
                    if (op == '!') {
                        if (mp_fact(&out, &a) != SUCCESS) {
                            mp_free(&a); mp_free(&out);
                            return FAILURE;
                        }
                    } else {
                        mp_free(&a); mp_free(&out);
                        return FAILURE; /* unknown unary */
                    }
                    mp_free(&a);

                    /* push result */
                    mp_init(&stack[sp]);
                    mp_copy(&stack[sp], &out);
                    mp_free(&out);
                    sp++;
                }
            } else {
                /* binary operator: pop two operands (b then a), compute a op b */
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

                    /* compute */
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
                        mp_free(&a); mp_free(&b); mp_free(&out);
                        return FAILURE;
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
            /* operand: parse and push onto stack */
            if (sp >= 128) return FAILURE;
            mp_init(&stack[sp]);
            if (parse_operand_to_mp(&stack[sp], tok) != SUCCESS) {
                /* cleanup partially built stack */
                {
                    int j;
                    for (j = 0; j <= sp; ++j) mp_free(&stack[j]);
                }
                return FAILURE;
            }
            sp++;
        }
    }

    /* final result must be single value on stack */
    if (sp != 1) {
        {
            int j;
            for (j = 0; j < sp; ++j) mp_free(&stack[j]);
        }
        return FAILURE;
    }

    mp_copy(result, &stack[0]);
    mp_free(&stack[0]);
    return SUCCESS;
}

/*
  split_expr:
    - utility that splits a simple two-operand expression into lhs, op, rhs.
    - This helper does NOT support parentheses or operator precedence;
      it exists for simpler one-line parsing use-cases.
*/
int split_expr(const char *input, char *lhs, size_t lhs_cap, char *op, char *rhs, size_t rhs_cap) {
    char buf[1024];
    size_t i, n;
    size_t left_len;
    size_t right_len;

    if (!input || !lhs || !rhs || !op) return 0;

    /* copy and trim input into local buffer */
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
                    /* leading sign, skip */
                    continue;
                }
                found = 1;
                break;
            }
        }
        if (!found) return 0;
    }

    /* split into left, op, right; trim both sides */
    left_len = i;
    right_len = (n > i + 1) ? n - (i + 1) : 0;

    if (left_len >= lhs_cap) left_len = lhs_cap - 1;
    strncpy(lhs, buf, left_len);
    lhs[left_len] = '\0';
    trim_inplace(lhs);

    *op = buf[i];

    if (right_len >= rhs_cap) right_len = rhs_cap - 1;
    strncpy(rhs, buf + i + 1, right_len);
    rhs[right_len] = '\0';
    trim_inplace(rhs);

    if (rhs[0] == '\0') return 0;
    if (lhs[0] == '\0') return 0;

    return 1;
}

/* parse_operand_to_mp: parse one operand string into mp_int (auto-detect base)
   Accepts an optional leading + / - followed by decimal, 0b... binary, or 0x... hex.
*/
int parse_operand_to_mp(mp_int *dst, const char *s) {
    char buf[1024];
    const char *q;
    int status;

    if (!dst || !s) return FAILURE;

    /* copy and trim local */
    strncpy(buf, s, sizeof(buf) - 1);
    buf[sizeof(buf) - 1] = '\0';
    trim_inplace(buf);
    if (buf[0] == '\0') return FAILURE;

    /* skip optional sign when detecting prefix */
    q = buf;
    if (*q == '+' || *q == '-') q++;

    if (q[0] == '0' && (q[1] == 'x' || q[1] == 'X')) {
        status = mp_from_str_hex(dst, buf);
        return status;
    } else if (q[0] == '0' && (q[1] == 'b' || q[1] == 'B')) {
        status = mp_from_str_bin(dst, buf);
        return status;
    } else {
        status = mp_from_str_dec(dst, buf);
        return status;
    }
}

/* ------------------------------------------------- File processing -----------------------------------------------
   process_file:
     - Reads the file line by line
     - Echoes the line (prints "> <line>")
     - Supports file-local commands (bin/hex/dec/out)
     - Otherwise treats the line as an expression and evaluates it using the same
       tokenize -> to_postfix -> eval_postfix pipeline as the REPL.
*/
int process_file(char str[]) {
    FILE *f;
    char line[1024];
    TokenList infix, postfix;
    mp_int result;
    int print_format = DEC;
    int (*mp_print[])(mp_int *) = { mp_print_dec, mp_print_bin, mp_print_hex };

    if (!str) return FAILURE;
    f = fopen(str, "r");
    if (!f) {
        /* note: original code exits here; returning FAILURE is generally
           better for embedding/automated use. Currently mirrors original behavior. */
        printf("Invalid input file!\n");
        exit(EXIT_FAILURE);
    }

    mp_init(&result);

    while (fgets(line, sizeof(line), f)) {
        size_t L = strlen(line);
        /* Trim newline / CR */
        while (L > 0 && (line[L - 1] == '\n' || line[L - 1] == '\r')) { line[--L] = '\0'; }

        /* Trim leading/trailing whitespace */
        trim_inplace(line);

        /* Skip empty lines */
        if (line[0] == '\0') continue;

        /* Echo input line before evaluating (requested behavior) */
        printf("> %s\n", line);

        /* Mirror REPL file-local commands */
        if (strcmp(line, "bin") == 0) {
            print_format = BIN;
            printf("bin\n");
            continue;
        }
        if (strcmp(line, "hex") == 0) {
            print_format = HEX;
            printf("hex\n");
            continue;
        }
        if (strcmp(line, "dec") == 0) {
            print_format = DEC;
            printf("dec\n");
            continue;
        }
        if (strcmp(line, "out") == 0) {
            switch (print_format) {
                case DEC: printf("dec\n"); break;
                case BIN: printf("bin\n"); break;
                case HEX: printf("hex\n"); break;
                default: printf("dec\n"); break;
            }
            continue;
        }

        /* Evaluate expression line */
        tokenize(line, &infix);
        to_postfix(&infix, &postfix);

        mp_free(&result);
        mp_init(&result);

        if (eval_postfix(&postfix, &result) == SUCCESS) {
            mp_print[print_format](&result);
            putchar('\n');
        }

        /* reset result and continue */
        mp_free(&result);
        mp_init(&result);
    }

    mp_free(&result);
    fclose(f);
    return SUCCESS;
}
