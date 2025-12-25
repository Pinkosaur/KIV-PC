#ifndef PARSER_H
#define PARSER_H

#include "mp_int.h"

/*
 * TokenList
 * ---------
 * Container for tokens extracted from an expression string.
 *
 * Limits:
 * - Max tokens per expression: 256
 * - Max length per token: 64 characters (including null terminator)
 */
typedef struct {
    char tokens[256][64];
    int count;
} TokenList;

/*
 * OpInfo
 * ------
 * Describes the properties of an operator for the Shunting-yard algorithm.
 *
 * Fields:
 * - op:           The operator character (e.g., '+', '*', '!', '~').
 * - precedence:   Higher values indicate higher precedence.
 * - right_assoc:  1 if the operator is right-associative (e.g., '^', '~'), 0 otherwise.
 * - unary:        1 if the operator is unary (e.g., '!', '~'), 0 if binary.
 */
typedef struct {
    char op;
    int precedence;
    int right_assoc;
    int unary;
} OpInfo;

/*
 * tokenize
 * --------
 * Splits an infix expression string into a list of tokens.
 * Handles whitespace, operators, parentheses, and numeric literals.
 * Detects unary minus context and emits '~' instead of '-'.
 */
void tokenize(const char *expr, TokenList *out);

/*
 * validate_tokens
 * ---------------
 * Checks the syntax of all tokens in the list BEFORE parsing begins.
 * Returns SUCCESS if all tokens are valid operators or strict numeric literals.
 * Returns FAILURE if any token is malformed (e.g., "0b102" or invalid chars).
 */
int validate_tokens(const TokenList *tokens);

/*
 * to_postfix
 * ----------
 * Converts an infix token list into Reverse Polish Notation (postfix)
 * using the Shunting-yard algorithm.
 */
void to_postfix(TokenList *infix, TokenList *postfix);

/*
 * eval_postfix
 * ------------
 * Evaluates a postfix token list using a stack of mp_int values.
 * Returns SUCCESS and populates 'result' on completion.
 * Returns FAILURE on division by zero, stack errors, or memory failure.
 */
int eval_postfix(TokenList *postfix, mp_int *result);

/*
 * parse_operand_to_mp
 * -------------------
 * Parses a single string token into a multiple-precision integer.
 * Automatically detects base:
 * - "0x..." -> Hexadecimal
 * - "0b..." -> Binary
 * - Otherwise -> Decimal
 */
int parse_operand_to_mp(mp_int *dst, const char *s);

/*
 * process_file
 * ------------
 * Reads a file line-by-line, treating each line as a command or expression.
 * Echoes input to stdout and prints the result.
 */
int process_file(char str[]);

#endif /* PARSER_H */