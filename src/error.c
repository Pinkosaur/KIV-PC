#include "error.h"

int calc_error_reported = 0;

void calc_error_set(void)      { calc_error_reported = 1; }
void calc_error_clear(void)    { calc_error_reported = 0; }
int calc_error_was_set(void)   { return calc_error_reported; }