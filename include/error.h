#ifndef ERROR_H
#define ERROR_H

/* * Global Runtime Error Flag
 * -------------------------
 * This flag indicates that a lower-level routine (arithmetic, factorial, etc.)
 * has already reported a specific error (e.g., "Division by zero!").
 * The main loop checks this flag to avoid printing a generic "Syntax error!"
 * on top of the specific error message.
 */

extern int calc_error_reported;

void calc_error_set(void);
void calc_error_clear(void);
int calc_error_was_set(void);

#endif /* ERROR_H */