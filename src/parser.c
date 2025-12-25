/*
 * parser.c -- Tokenizer, Shunting-yard parser, Postfix evaluator, and File processing.
 *
 * This module implements the core parsing logic for the arbitrary precision calculator.
 *
 * Logic Overview:
 * - Tokenizer: Scans the input string and identifies operators, parentheses, and operands.
 * It handles the ambiguity of the '-' character (binary subtraction vs unary negation)
 * by tracking the context. Unary minus is converted to a special internal operator '~'.
 * - Shunting-yard: Converts the infix stream of tokens into Reverse Polish Notation (Postfix).
 * - Evaluator: Processes the postfix queue using a stack of mp_int values.
 * - Validation: Pre-validates tokens to prevent execution of syntax errors.
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

static void trim_inplace(char *s) {
    char *p = s;
    size_t len;
    
    while (*p && isspace((unsigned char)*p)) p++;
    
    if (p != s) memmove(s, p, strlen(p) + 1);
    
    len = strlen(s);
    while (len > 0 && isspace((unsigned char)s[len - 1])) s[--len] = '\0';
}

/* Operator table with precedence and associativity */
static OpInfo op_table[] = {
    {'!', 6, 0, 1},
    {'~', 5, 1, 1},
    {'^', 4, 1, 0},
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

/* Internal helper: Check if string is a valid numeric literal */
static int is_valid_literal(const char *s) {
    const char *p = s;
    if (!s || !*s) return 0;
    
    /* Skip optional sign */
    if (*p == '+' || *p == '-') p++;
    
    /* Check Hex */
    if (p[0] == '0' && (p[1] == 'x' || p[1] == 'X')) {
        p += 2;
        if (*p == '\0') return 0;
        while (*p) {
            if (!isxdigit((unsigned char)*p)) return 0;
            p++;
        }
        return 1;
    }
    /* Check Binary */
    if (p[0] == '0' && (p[1] == 'b' || p[1] == 'B')) {
        p += 2;
        if (*p == '\0') return 0;
        while (*p) {
            if (*p != '0' && *p != '1') return 0;
            p++;
        }
        return 1;
    }
    /* Check Decimal */
    if (*p == '\0') return 0;
    while (*p) {
        if (!isdigit((unsigned char)*p)) return 0;
        p++;
    }
    return 1;
}

/* Public validator */
int validate_tokens(const TokenList *tokens) {
    int i;
    for (i = 0; i < tokens->count; ++i) {
        const char *t = tokens->tokens[i];
        /* Skip operators and parens */
        if (t[1] == '\0' && (is_operator(t[0]) || t[0] == '(' || t[0] == ')')) {
            continue;
        }
        /* Check operand syntax */
        if (!is_valid_literal(t)) {
            return FAILURE;
        }
    }
    return SUCCESS;
}

/* Tokenizer */
void tokenize(const char *expr, TokenList *out) {
    const char *p;
    char buf[64];
    int i;
    int prev_was_operand; 

    out->count = 0;
    prev_was_operand = 0; 

    p = expr;
    while (*p) {
        while (isspace((unsigned char)*p)) p++;
        if (*p == '\0') break;

        if (out->count >= 256) break;

        /* Handle + and - based on context (unary vs binary) */
        if ((*p == '+' || *p == '-') && !prev_was_operand) {
            if (*p == '-') {
                out->tokens[out->count][0] = '~';
                out->tokens[out->count][1] = '\0';
                out->count++;
            }
            /* Unary + is ignored */
            p++;
            continue;
        }

        /* Handle operators and parentheses */
        if (is_operator(*p) || *p == '(' || *p == ')') {
            out->tokens[out->count][0] = *p;
            out->tokens[out->count][1] = '\0';
            out->count++;
            
            if (*p == ')') {
                prev_was_operand = 1;
            } else if (*p == '!') {
                /* Postfix operator '!' acts as an operand for the next operator.
                   This ensures "10! - 5" is parsed as binary minus, not unary. */
                prev_was_operand = 1;
            } else {
                prev_was_operand = 0;
            }
            
            p++;
            continue;
        }

        /* Handle operands (literals) */
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

/* Shunting-yard */
void to_postfix(TokenList *infix, TokenList *postfix) {
    char stack[256];
    int i, sp = 0, p_op, p_top;
    char *tok, op, top;

    postfix->count = 0;

    for (i = 0; i < infix->count; ++i) {
        tok = infix->tokens[i];
        if (tok[0] != '\0' && tok[1] == '\0' && is_operator(tok[0])) {
            op = tok[0];
            while (sp > 0) {
                top = stack[sp-1];
                if (top == '(') break;
                p_op = get_precedence(op);
                p_top = get_precedence(top);
                if (p_top > p_op || (p_top == p_op && !is_right_assoc(op))) {
                    postfix->tokens[postfix->count][0] = stack[--sp];
                    postfix->tokens[postfix->count][1] = '\0';
                    postfix->count++;
                } else break;
            }
            stack[sp++] = op;
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

/* Postfix Evaluator */
int eval_postfix(TokenList *postfix, mp_int *result) {
    mp_int stack[128], a, b, out, zero;
    int i, sp = 0, j, op_status;
    char *tok, op;

    for (i = 0; i < postfix->count; ++i) {
        tok = postfix->tokens[i];
        if (tok[0] != '\0' && tok[1] == '\0' && is_operator(tok[0])) {
            op = tok[0];
            if (is_unary(op)) {
                if (sp < 1) { for (j=0; j<sp; ++j) mp_free(&stack[j]); return FAILURE; }
                mp_init(&a); mp_init(&out);
                mp_copy(&a, &stack[--sp]); mp_free(&stack[sp]);
                op_status = SUCCESS;
                if (op == '!') op_status = mp_fact(&out, &a);
                else if (op == '~') { mp_init(&zero); op_status = mp_sub(&out, &zero, &a); mp_free(&zero); }
                else op_status = FAILURE;
                mp_free(&a);
                if (op_status != SUCCESS) { mp_free(&out); for (j=0; j<sp; ++j) mp_free(&stack[j]); return FAILURE; }
                mp_init(&stack[sp]); mp_copy(&stack[sp], &out); mp_free(&out); sp++;
            } else {
                if (sp < 2) { for (j=0; j<sp; ++j) mp_free(&stack[j]); return FAILURE; }
                mp_init(&a); mp_init(&b); mp_init(&out);
                mp_copy(&b, &stack[--sp]); mp_free(&stack[sp]);
                mp_copy(&a, &stack[--sp]); mp_free(&stack[sp]);
                op_status = SUCCESS;
                if (op == '+') op_status = mp_add(&out, &a, &b);
                else if (op == '-') op_status = mp_sub(&out, &a, &b);
                else if (op == '*') op_status = mp_mul(&out, &a, &b);
                else if (op == '/') op_status = mp_div(&out, &a, &b);
                else if (op == '%') op_status = mp_mod(&out, &a, &b);
                else if (op == '^') op_status = mp_pow(&out, &a, &b);
                else op_status = FAILURE;
                mp_free(&a); mp_free(&b);
                if (op_status != SUCCESS) { mp_free(&out); for (j=0; j<sp; ++j) mp_free(&stack[j]); return FAILURE; }
                mp_init(&stack[sp]); mp_copy(&stack[sp], &out); mp_free(&out); sp++;
            }
        } else {
            if (sp >= 128) { for (j=0; j<sp; ++j) mp_free(&stack[j]); return FAILURE; }
            mp_init(&stack[sp]);
            if (parse_operand_to_mp(&stack[sp], tok) != SUCCESS) { for (j=0; j<=sp; ++j) mp_free(&stack[j]); return FAILURE; }
            sp++;
        }
    }
    if (sp != 1) { for (j=0; j<sp; ++j) mp_free(&stack[j]); return FAILURE; }
    mp_copy(result, &stack[0]); mp_free(&stack[0]);
    return SUCCESS;
}

int parse_operand_to_mp(mp_int *dst, const char *s) {
    char buf[1024]; const char *q; int status;
    if (!dst || !s) return FAILURE;
    strncpy(buf, s, sizeof(buf)-1); buf[sizeof(buf)-1]='\0'; trim_inplace(buf);
    if (buf[0] == '\0') return FAILURE;
    q = buf; if (*q == '+' || *q == '-') q++;
    if (q[0] == '0' && (q[1] == 'x' || q[1] == 'X')) status = mp_from_str_hex(dst, buf);
    else if (q[0] == '0' && (q[1] == 'b' || q[1] == 'B')) status = mp_from_str_bin(dst, buf);
    else status = mp_from_str_dec(dst, buf);
    return status;
}

/* File Processing */
int process_file(char str[]) {
    FILE *f; char line[1024]; TokenList infix, postfix; mp_int result;
    int print_format = DEC;
    int (*mp_print[])(mp_int *) = { mp_print_dec, mp_print_bin, mp_print_hex };

    if (!str) return FAILURE;
    f = fopen(str, "r");
    if (!f) { printf("Invalid input file!\n"); exit(EXIT_FAILURE); }
    mp_init(&result);

    while (fgets(line, sizeof(line), f)) {
        size_t L = strlen(line);
        calc_error_clear();
        while (L > 0 && (line[L-1] == '\n' || line[L-1] == '\r')) line[--L] = '\0';
        trim_inplace(line);
        if (line[0] == '\0') continue;
        printf("> %s\n", line);

        /* Handle Commands */
        if (!strcmp(line, "quit") || !strcmp(line, "exit")) {
            printf("%s\n", line); /* Echo the quit command to output */
            break; /* Stop processing the file */
        }
        if (!strcmp(line, "bin")) { print_format=BIN; printf("bin\n"); continue; }
        if (!strcmp(line, "hex")) { print_format=HEX; printf("hex\n"); continue; }
        if (!strcmp(line, "dec")) { print_format=DEC; printf("dec\n"); continue; }
        if (!strcmp(line, "out")) {
            if (print_format==DEC) printf("dec\n");
            else if (print_format==BIN) printf("bin\n");
            else printf("hex\n");
            continue;
        }

        tokenize(line, &infix);
        
        /* Validate syntax before evaluation */
        if (validate_tokens(&infix) != SUCCESS) {
            if (infix.count == 1 && isalpha((unsigned char)infix.tokens[0][0])) {
                printf("Invalid command \"%s\"!\n", infix.tokens[0]);
            } else {
                printf("Syntax error!\n");
            }
            continue;
        }

        to_postfix(&infix, &postfix);
        mp_free(&result); mp_init(&result);

        if (eval_postfix(&postfix, &result) == SUCCESS) {
            mp_print[print_format](&result); putchar('\n');
        } else {
            if (!calc_error_was_set()) printf("Syntax error!\n");
        }
        mp_free(&result); mp_init(&result);
    }
    mp_free(&result); fclose(f); return SUCCESS;
}