/*
 * main.c -- calculator REPL / file-processing entry point
 *
 * Behavior:
 * - With no arguments: interactive REPL (commands: bin/hex/dec/out/exit/quit).
 * - With one argument: treat as filename and call process_file(filename).
 *
 * Important change:
 * - Do NOT attempt to parse the input again as a literal after a successful
 * expression evaluation. Only try literal parsing when expression parsing
 * / evaluation fails (and there was no runtime error).
 * - Literal parsing is strict: optional sign, optional 0x/0b prefix, then
 * only digits of the selected base until the end of the token.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "parser.h"
#include "mp_print.h"
#include "mp_int.h"
#include "error.h" /* Included to access calc_error functions */

/* Validate that 's' is a strict standalone numeric literal (no trailing garbage).
   Returns:
     10 for strict decimal literal,
      2 for strict binary literal (optionally "0b" prefix),
     16 for strict hex literal (optionally "0x" prefix),
      0 if not a valid strict literal.
   Accepts optional leading '+' or '-' (leading whitespace is trimmed before calling).
*/
static int is_strict_literal(const char *s) {
    const char *p;
    int saw_digit = 0;

    if (!s) return 0;
    p = s;

    /* skip leading spaces (caller generally trimmed, but be robust) */
    while (*p && isspace((unsigned char)*p)) p++;

    /* optional sign */
    if (*p == '+' || *p == '-') p++;

    /* detect 0x / 0b prefix (only if there is at least one more char) */
    if (p[0] == '0' && (p[1] == 'x' || p[1] == 'X')) {
        p += 2;
        /* require at least one hex digit */
        for ( ; *p; ++p) {
            if (!isxdigit((unsigned char)*p)) return 0;
            saw_digit = 1;
        }
        return saw_digit ? 16 : 0; /* hexadecimal or decimal */
    }
    if (p[0] == '0' && (p[1] == 'b' || p[1] == 'B')) {
        p += 2;
        for (; *p; ++p) {
            if (*p != '0' && *p != '1') return 0;
            saw_digit = 1;
        }
        return saw_digit ? 2 : 0; /* binary or decimal */
    }

    /* decimal: require at least one digit, and only digits until end */
    for (; *p; ++p) {
        if (!isdigit((unsigned char)*p)) return 0;
        saw_digit = 1;
    }
    return saw_digit ? 10 : 0;
}

int main(int argc, char *argv[]) {
    /* Format selection: function pointer table for printing functions.
       Each entry is int func(mp_int *). */
    int (*mp_print[])(mp_int *) = { mp_print_dec, mp_print_bin, mp_print_hex };
    int print_format = DEC;
    int base;
    /* Temporaries used by REPL and expression evaluation */
    mp_int a, b, c;
    /* Buffers sized for typical use; modest for 32-bit stack safety */
    char input_file[300];
    char buf[1024];
    TokenList infix, postfix;
    size_t L;
    char *p;

    /* Argument handling */
    if (argc > 2) {
        printf("Error: invalid input\nUsage: calc.exe [<input-file>]\n");
        return EXIT_FAILURE;
    }
    if (argc == 2) {
        /* copy filename safely and delegate to file processor */
        strncpy(input_file, argv[1], sizeof(input_file) - 1);
        input_file[sizeof(input_file) - 1] = '\0';
        return process_file(input_file) != FAILURE ? EXIT_SUCCESS : EXIT_FAILURE;
    }

    /* Initialize temporaries */
    mp_init(&a);
    mp_init(&b);
    mp_init(&c);

    /* REPL loop */
    while (1) {
        printf("> ");
        if (!fgets(buf, sizeof(buf), stdin)) break;

        /* Clear the global error flag before processing new input */
        calc_error_clear();

        /* Trim trailing newline / CR */
        L = strlen(buf);
        while (L > 0 && (buf[L - 1] == '\n' || buf[L - 1] == '\r')) { buf[--L] = '\0'; }

        /* Trim leading spaces before parsing commands/expressions */
        p = buf;
        while (*p && isspace((unsigned char)*p)) p++;

        if (p[0] == '\0') continue; /* empty line -> prompt again */

        /* Exit commands */
        if (strcmp(p, "exit") == 0 || strcmp(p, "quit") == 0) break;

        /* Change output format */
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
            /* Echo the currently selected output format */
            switch(print_format) {
                case DEC: printf("dec\n"); break;
                case BIN: printf("bin\n"); break;
                case HEX: printf("hex\n"); break;
                default: printf("dec\n"); break;
            }
            continue;
        }

        /*
         * Expression handling:
         * - Tokenize input
         * - Convert to postfix (shunting-yard)
         * - Evaluate postfix, print the evaluated result
         *
         * The parser/tokenizer handle unary minus, postfix factorial precedence, parentheses, etc.
         */

        tokenize(p, &infix);
        to_postfix(&infix, &postfix);

        /* ensure 'c' is empty and ready to receive the evaluation result */
        mp_free(&c);
        mp_init(&c);

        if (eval_postfix(&postfix, &c) == SUCCESS) {
            /* Expression evaluated successfully -> print and SKIP literal parsing */
            mp_print[print_format](&c);
            putchar('\n');
            mp_free(&c);
            continue;
        } else {
            /* Evaluation failed. 
             * Check if it was due to a runtime error (like Division by Zero) that was already reported.
             */
            if (calc_error_was_set()) {
                /* Error message already printed by lower-level function, so don't print "Syntax error!" */
                continue;
            }

            /* Try strict literal parsing (only if the entire input is a valid literal) */
            {
                base = is_strict_literal(p);
                if (base == 10) {
                    mp_free(&c);
                    mp_init(&c);
                    if (mp_from_str_dec(&c, p) == SUCCESS) {
                        mp_print[print_format](&c);
                        putchar('\n');
                        mp_free(&c);
                        continue;
                    } else {
                        printf("Syntax error!\n");
                        continue;
                    }
                } else if (base == 2) {
                    mp_free(&c);
                    mp_init(&c);
                    if (mp_from_str_bin(&c, p) == SUCCESS) {
                        mp_print[print_format](&c);
                        putchar('\n');
                        mp_free(&c);
                        continue;
                    } else {
                        printf("Syntax error!\n");
                        continue;
                    }
                } else if (base == 16) {
                    mp_free(&c);
                    mp_init(&c);
                    if (mp_from_str_hex(&c, p) == SUCCESS) {
                        mp_print[print_format](&c);
                        putchar('\n');
                        mp_free(&c);
                        continue;
                    } else {
                        printf("Syntax error!\n");
                        continue;
                    }
                } else {
                    /* not a strict literal and expression parsing failed -> syntax error */
                    printf("Syntax error!\n");
                    continue;
                }
            }
        }
    }

    /* Cleanup and exit */
    mp_free(&a);
    mp_free(&b);
    mp_free(&c);
    return EXIT_SUCCESS;
}