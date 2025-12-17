#ifndef ERROR_H
#define ERROR_H

/* A global runtime-error flag used to indicate that a lower-level
   routine already printed an error message (division by zero, negative
   factorial etc.). This avoids duplicating "Syntax error!" after a runtime
   error has already been reported. */

extern int calc_error_reported;

void calc_error_set(void);
void calc_error_clear(void);
int calc_error_was_set(void);

#endif /* ERROR_H */
