#include "structure.h"

double integral(double (* f)(double x, void * params), Fn_Param *params, double *error, unsigned int *intervals);
