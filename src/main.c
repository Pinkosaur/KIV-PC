#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "parser.h"
#include "mp_print.h"
#include "mp_int.h"

int main(int argc, char *argv[]) {
    /* Format selection */
    int (*mp_print[])(mp_int *) = { mp_print_dec, mp_print_bin, mp_print_hex };
    int print_format = DEC;

    mp_int a, b, c;
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
        printf("> ");
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
            if (mp_from_str_hex(&c, p) != SUCCESS) { 
                printf("parse error\n"); 
                continue; 
            }
        } else if (p[0] == '0' && (p[1] == 'b' || p[1] == 'B')) {
            if (mp_from_str_bin(&c, p) != SUCCESS) { 
                printf("parse error\n"); 
                continue; 
            }
        } else {
            if (mp_from_str_dec(&c, p) != SUCCESS) { 
                printf("parse error\n"); 
                continue; 
            }
        }
        putchar('\n');
    }

    mp_free(&a);
    mp_free(&b);
    mp_free(&c);
    return EXIT_SUCCESS;
}
