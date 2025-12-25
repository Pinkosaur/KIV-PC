/*
 * parser.c -- Tokenizer, Shunting-yard parser, Postfix evaluator, and File processing.
 *
 * This module implements the core parsing logic for the arbitrary precision calculator.
 *
 * Logic Overview:
 * - Tokenizer: Scans the input string and identifies operators, parentheses, and operands.
 * It handles the ambiguity of the '-' character (binary subtraction vs unary negation)
 * by tracking the context. Unary minus is converted to a special internal operator '~'.
 *
 * - Shunting-yard: Converts the infix stream of tokens into Reverse Polish Notation (Postfix).
 * It respects operator precedence and associativity (including right-associative unary operators).
 *
 * - Evaluator: Processes the postfix queue using a stack of mp_int values to compute the result.
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

/* Helper: returns true if the character is one of the single-char operators we support.
 * Includes '~', the internal representation for unary negation.
 */
static int is_operator(char c) {
    return strchr("+-*/%^!~", c) != NULL;
}

/*
 * Trim in-place leading and trailing whitespace.
 */
static void trim_inplace(char *s) {
    char *p = s;
    size_t len;
    
    while (*p && isspace((unsigned char)*p)) p++;
    
    if (p != s) memmove(s, p, strlen(p) + 1); /* shift left if needed */
    
    len = strlen(s);
    while (len > 0 && isspace((unsigned char)s[len - 1])) s[--len] = '\0';
}

/* Operator metadata table.
 *
 * Fields:
 * - op:           The operator character.
 * - precedence:   Higher value means tighter binding.
 * - right_assoc:  1 if right-associative (e.g., ^, ~), 0 if left-associative.
 * - unary:        1 if the operator is unary (prefix or postfix).
 *
 * Notes on Precedence:
 * - '!' (Factorial, Postfix) has the highest precedence (6).
 * - '~' (Unary Minus, Prefix) has precedence 5.
 * - '^' (Exponentiation) has precedence 4.
 * This ensures that -x! is parsed as -(x!) rather than (-x)!.
 */
static OpInfo op_table[] = {
    {'!', 6, 0, 1},  /* factorial, unary postfix */
    {'~', 5, 1, 1},  /* unary minus, prefix (internal operator) */
    {'^', 4, 1, 0},  /* exponentiation */
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

/* Scans the expression string and populates the TokenList.
 *
 * Context-Sensitive Logic:
 * - Detects whether a '+' or '-' is a binary operator or a unary sign based on
 * the previous token (start of expression, or following an operator/parenthesis).
 * - A unary '-' is emitted as the special operator '~'.
 * - A unary '+' is ignored (no-op).
 * - Numeric literals are parsed as positive strings; their sign is determined
 * by the operators acting on them.
 */
void tokenize(const char *expr, TokenList *out) {
    const char *p;
    char buf[64];
    int i;
    int prev_was_operand; /* 0 = start/after-operator/'(', 1 = after-operand or ')' */

    out->count = 0;
    prev_was_operand = 0; 

    p = expr;
    while (*p) {
        /* skip whitespace */
        while (isspace((unsigned char)*p)) p++;
        if (*p == '\0') break;

        /* safety: do not exceed token storage */
        if (out->count >= 256) break;

        /* Handle unary signs context */
        if ((*p == '+' || *p == '-') && !prev_was_operand) {
            if (*p == '-') {
                /* Emit unary minus operator '~' */
                out->tokens[out->count][0] = '~';
                out->tokens[out->count][1] = '\0';
                out->count++;
                /* prev_was_operand remains 0 because '~' is an operator */
            }
            /* If '+', we simply skip it (unary plus is a no-op) */
            p++;
            continue;
        }

        /* Single-char operator or parentheses */
        if (is_operator(*p) || *p == '(' || *p == ')') {
            out->tokens[out->count][0] = *p;
            out->tokens[out->count][1] = '\0';
            out->count++;
            
            /* update prev flag */
            if (*p == ')') prev_was_operand = 1;
            else prev_was_operand = 0; /* '(' or operators */
            
            p++;
            continue;
        }

        /* Operand (numeric literal) */
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

/* Converts infix token list to postfix using the Shunting-yard algorithm.
 * Handles operator precedence and associativity (Right vs Left).
 */
void to_postfix(TokenList *infix, TokenList *postfix) {
    char stack[256];
    int i;
    int sp = 0;
    char *tok;
    char op;
    char top;
    int p_op, p_top;

    postfix->count = 0;

    for (i = 0; i < infix->count; ++i) {
        tok = infix->tokens[i];

        /* operator token (single char) */
        if (tok[0] != '\0' && tok[1] == '\0' && is_operator(tok[0])) {
            op = tok[0];
            
            /* Pop operators from stack based on precedence and associativity */
            while (sp > 0) {
                top = stack[sp-1];
                if (top == '(') break; /* Stop at open paren */

                p_op = get_precedence(op);
                p_top = get_precedence(top);

                /* Pop if:
                   - Top has higher precedence
                   - OR Top has equal precedence AND op is Left Associative */
                if (p_top > p_op || (p_top == p_op && !is_right_assoc(op))) {
                    postfix->tokens[postfix->count][0] = stack[--sp];
                    postfix->tokens[postfix->count][1] = '\0';
                    postfix->count++;
                } else {
                    break;
                }
            }
            stack[sp++] = op;
            
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
            /* operand */
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

/* Evaluates the postfix token list.
 * Supports binary operators (+, -, *, /, %, ^) and unary operators (!, ~).
 */
int eval_postfix(TokenList *postfix, mp_int *result) {
    mp_int stack[128], a, b, out, zero;
    int i;
    int sp = 0;
    char *tok;
    char op;
    int j;
    int op_status = SUCCESS;

    for (i = 0; i < postfix->count; ++i) {
        tok = postfix->tokens[i];

        /* operator token (single char) */
        if (tok[0] != '\0' && tok[1] == '\0' && is_operator(tok[0])) {
            op = tok[0];

            if (is_unary(op)) {
                /* postfix/prefix unary: pop one operand */
                if (sp < 1) {
                    /* Not enough operands */
                    for (j = 0; j < sp; ++j) mp_free(&stack[j]);
                    return FAILURE;
                }
                
                mp_init(&a); 
                mp_init(&out);

                mp_copy(&a, &stack[--sp]);
                mp_free(&stack[sp]);

                op_status = SUCCESS;
                if (op == '!') {
                    op_status = mp_fact(&out, &a);
                } else if (op == '~') {
                    /* Negation: 0 - a */
                    mp_init(&zero); /* 0 by default */
                    op_status = mp_sub(&out, &zero, &a);
                    mp_free(&zero);
                } else {
                    op_status = FAILURE; /* unknown unary */
                }

                mp_free(&a);

                if (op_status != SUCCESS) {
                    mp_free(&out);
                    /* Leak fix: free remaining stack items */
                    for (j = 0; j < sp; ++j) mp_free(&stack[j]);
                    return FAILURE;
                }

                mp_init(&stack[sp]);
                mp_copy(&stack[sp], &out);
                mp_free(&out);
                sp++;
            } else {
                /* binary operator: pop two operands */
                if (sp < 2) {
                    /* Not enough operands (and potentially stack[0] leaking) */
                    for (j = 0; j < sp; ++j) mp_free(&stack[j]);
                    return FAILURE;
                }
                                  
                mp_init(&a); 
                mp_init(&b); 
                mp_init(&out);

                /* pop b then a */
                mp_copy(&b, &stack[--sp]);
                mp_free(&stack[sp]);
                mp_copy(&a, &stack[--sp]);
                mp_free(&stack[sp]);

                /* compute */
                op_status = SUCCESS;
                if (op == '+')      op_status = mp_add(&out, &a, &b);
                else if (op == '-') op_status = mp_sub(&out, &a, &b);
                else if (op == '*') op_status = mp_mul(&out, &a, &b);
                else if (op == '/') op_status = mp_div(&out, &a, &b);
                else if (op == '%') op_status = mp_mod(&out, &a, &b);
                else if (op == '^') op_status = mp_pow(&out, &a, &b);
                else                op_status = FAILURE;

                mp_free(&a);
                mp_free(&b);

                if (op_status != SUCCESS) {
                    mp_free(&out);
                    /* Leak fix: free remaining stack items */
                    for (j = 0; j < sp; ++j) mp_free(&stack[j]);
                    return FAILURE;
                }

                mp_init(&stack[sp]);
                mp_copy(&stack[sp], &out);
                mp_free(&out);
                sp++;
            }
        } else {
            /* operand */
            if (sp >= 128) {
                /* Stack overflow - clean up */
                for (j = 0; j < sp; ++j) mp_free(&stack[j]);
                return FAILURE;
            }
            mp_init(&stack[sp]);
            if (parse_operand_to_mp(&stack[sp], tok) != SUCCESS) {
                /* cleanup partially built stack */
                for (j = 0; j <= sp; ++j) mp_free(&stack[j]);
                return FAILURE;
            }
            sp++;
        }
    }

    if (sp != 1) {
        for (j = 0; j < sp; ++j) mp_free(&stack[j]);
        return FAILURE;
    }

    mp_copy(result, &stack[0]);
    mp_free(&stack[0]);
    return SUCCESS;
}

/* Utility to split a simple two-operand expression (lhs op rhs) without parentheses.
 * Useful for simple internal tests or one-line parsing logic that guarantees no complex syntax.
 */
int split_expr(const char *input, char *lhs, size_t lhs_cap, char *op, char *rhs, size_t rhs_cap) {
    char buf[1024], ch;
    size_t i, n;
    size_t left_len;
    size_t right_len;
    int found = 0;

    if (!input || !lhs || !rhs || !op) return 0;

    strncpy(buf, input, sizeof(buf) - 1);
    buf[sizeof(buf) - 1] = '\0';
    trim_inplace(buf);
    if (buf[0] == '\0') return 0;

    n = strlen(buf);
    i = 0;
    
    for (i = 0; i < n; ++i) {
        ch = buf[i];
        /* * Naive split only looks for binary operators. 
         * It treats the first operator found (ignoring leading sign) as the split point.
         */
        if (ch == '+' || ch == '-' || ch == '*' || ch == '/') {
            if (i == 0) continue; 
            found = 1;
            break;
        }
    }
    if (!found) return 0;

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

/* Parses a string operand into an mp_int.
 * Automatically detects decimal, binary (0b), and hexadecimal (0x) formats.
 */
int parse_operand_to_mp(mp_int *dst, const char *s) {
    char buf[1024];
    const char *q;
    int status;

    if (!dst || !s) return FAILURE;

    strncpy(buf, s, sizeof(buf) - 1);
    buf[sizeof(buf) - 1] = '\0';
    trim_inplace(buf);
    if (buf[0] == '\0') return FAILURE;

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

/* Reads an input file line-by-line and processes expressions or commands.
 * Mirrors the behavior of the interactive REPL.
 */
int process_file(char str[]) {
    FILE *f;
    char line[1024];
    TokenList infix, postfix;
    mp_int result;
    size_t L;
    int print_format = DEC;
    int (*mp_print[])(mp_int *) = { mp_print_dec, mp_print_bin, mp_print_hex };

    if (!str) return FAILURE;
    f = fopen(str, "r");
    if (!f) {
        printf("Invalid input file!\n");
        exit(EXIT_FAILURE);
    }

    mp_init(&result);

    while (fgets(line, sizeof(line), f)) {
        L = strlen(line);
        while (L > 0 && (line[L - 1] == '\n' || line[L - 1] == '\r')) { line[--L] = '\0'; }
        trim_inplace(line);
        if (line[0] == '\0') continue;

        printf("> %s\n", line);

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

        tokenize(line, &infix);
        to_postfix(&infix, &postfix);

        mp_free(&result);
        mp_init(&result);

        if (eval_postfix(&postfix, &result) == SUCCESS) {
            mp_print[print_format](&result);
            putchar('\n');
        }

        mp_free(&result);
        mp_init(&result);
    }

    mp_free(&result);
    fclose(f);
    return SUCCESS;
}