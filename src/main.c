/*
 * main.c -- Calculator REPL Entry Point
 *
 * Handles argument processing, interactive REPL loop, and dispatching
 * to file processing or line evaluation.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "parser.h"
#include "mp_print.h"
#include "mp_int.h"
#include "error.h"

int main(int argc, char *argv[]) {
    int (*mp_print[])(mp_int *) = { mp_print_dec, mp_print_bin, mp_print_hex };
    int print_format = DEC;
    mp_int c;
    char input_file[300];
    char buf[1024];
    TokenList infix, postfix;
    size_t L;
    char *p;

    /* Handle file argument */
    if (argc > 2) {
        printf("Error: invalid input\nUsage: calc.exe [<input-file>]\n");
        return EXIT_FAILURE;
    }
    if (argc == 2) {
        strncpy(input_file, argv[1], sizeof(input_file) - 1);
        input_file[sizeof(input_file) - 1] = '\0';
        return process_file(input_file) != FAILURE ? EXIT_SUCCESS : EXIT_FAILURE;
    }

    mp_init(&c);

    /* Interactive REPL loop */
    while (1) {
        printf("> ");
        if (!fgets(buf, sizeof(buf), stdin)) break;

        calc_error_clear();

        L = strlen(buf);
        while (L > 0 && (buf[L - 1] == '\n' || buf[L - 1] == '\r')) { buf[--L] = '\0'; }

        p = buf;
        while (*p && isspace((unsigned char)*p)) p++;

        if (p[0] == '\0') continue;

        if (strcmp(p, "exit") == 0 || strcmp(p, "quit") == 0) break;

        if (strcmp(p, "bin") == 0) { print_format = BIN; printf("bin\n"); continue; }
        if (strcmp(p, "hex") == 0) { print_format = HEX; printf("hex\n"); continue; }
        if (strcmp(p, "dec") == 0) { print_format = DEC; printf("dec\n"); continue; }
        if (strcmp(p, "out") == 0) {
            if (print_format == DEC) printf("dec\n");
            else if (print_format == BIN) printf("bin\n");
            else printf("hex\n");
            continue;
        }

        tokenize(p, &infix);

        /* Validate syntax BEFORE evaluation to prevent crashes on invalid inputs */
        if (validate_tokens(&infix) != SUCCESS) {
            /* Check for invalid command */
            if (infix.count == 1 && isalpha((unsigned char)infix.tokens[0][0])) {
                printf("Invalid command \"%s\"!\n", infix.tokens[0]);
            } else {
                printf("Syntax error!\n");
            }
            continue;
        }

        to_postfix(&infix, &postfix);

        mp_free(&c);
        mp_init(&c);

        if (eval_postfix(&postfix, &c) == SUCCESS) {
            mp_print[print_format](&c);
            putchar('\n');
            mp_free(&c);
        } else {
            if (!calc_error_was_set()) {
                printf("Syntax error!\n");
            }
        }
    }

    mp_free(&c);
    return EXIT_SUCCESS;
}