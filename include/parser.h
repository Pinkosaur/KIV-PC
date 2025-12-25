#ifndef PARSER_H
#define PARSER_H

#include "mp_int.h"

/*
 * TokenList
 *   A simple fixed-size token container.
 *   Each token is stored as a null-terminated string.
 *
 *   Limits:
 *    - max tokens per expression: 256
 *    - max token length: 63 chars + '\0'
 *
 *   This is sufficient for all supported expressions and ANSI C compatibility.
 */
typedef struct {
    char tokens[256][64];
    int count;
} TokenList;

/*
 * OpInfo
 *   Describes operator properties for the shunting-yard algorithm.
 *
 *   op           : operator character
 *   precedence   : higher value = higher precedence
 *   right_assoc  : nonzero if right-associative (e.g. '^')
 *   unary        : nonzero if unary (currently only postfix '!')
 */
typedef struct {
    char op;
    int precedence;
    int right_assoc;
    int unary;
} OpInfo;

/* Tokenize an infix expression into a flat token list */
void tokenize(const char *expr, TokenList *out);

/* Convert infix token list to postfix (Reverse Polish) form */
void to_postfix(TokenList *infix, TokenList *postfix);

/* Evaluate postfix expression and store result */
int eval_postfix(TokenList *postfix, mp_int *result);

/*
 * Parse a single operand token into mp_int.
 * Automatically detects:
 *   - decimal
 *   - binary (0b / 0B)
 *   - hexadecimal (0x / 0X)
 * Handles optional leading '+' or '-'.
 */
int parse_operand_to_mp(mp_int *dst, const char *s);

/*
 * Process an input file line-by-line.
 * Supports:
 *   - expressions
 *   - format commands: dec / bin / hex / out
 * Uses the same logic as the interactive REPL.
 */
int process_file(char str[]);

#endif /* PARSER_H */
