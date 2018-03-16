#include "structure.h"

void my_gsl_error(const char * reason, const char * file, int line, int gsl_errno);

// Integrate an even function f from x = -1 to 1 assuming f(-x) = f(x)
double integral(double (* f)(double x, void * params), Fn_Param *params, double *error, unsigned int *intervals);
