#ifndef PARSER_H
#define PARSER_H

#include "mp_int.h"

typedef struct {
    char tokens[256][64];
    int count;
} TokenList;

/* Operator precedence and associativity */
typedef struct {
    char op;
    int precedence;
    int right_assoc;
    int unary;
} OpInfo;

void tokenize(const char *expr, TokenList *out);
void to_postfix(TokenList *infix, TokenList *postfix);
int eval_postfix(TokenList *postfix, mp_int *result);

/* helper used by parser */
int parse_operand_to_mp(mp_int *dst, const char *s);

/* a simple file processing helper */
int process_file(char str[]);

#endif /* PARSER_H */
