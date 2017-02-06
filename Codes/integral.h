#include "structure.h"

double integral(const double (* f)(double x, void * params), Fn_Param *params, double *error, unsigned int *intervals);
