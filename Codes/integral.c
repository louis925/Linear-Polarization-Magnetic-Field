/*
GSL Integration Example: http://www.gnu.org/software/gsl/manual/html_node/Numerical-integration-examples.html
*/

#include <stdio.h>
#include <gsl/gsl_integration.h>
#include "global_value.h"
#include "integral.h"
#include "physics_coeff.h"

double integral(const double (* f)(double x, void * params), Fn_Param *params, double *error, unsigned int *intervals)
{ 
	// We use integral() to integral only the inclination angle (from 0 to pi) 
	
	static gsl_function F;
	double result;
	//w = gsl_integration_workspace_alloc (3000); //have already been allocated in main.c
	F.function = f;
	F.params = params;
	
	gsl_integration_qags (&F, 0, 1, EpsAbs, EpsRel, Gsl_Integ_Space, w, &result, error); //Integrate from 0 to 1
	
	*intervals = w->size;

	return result*2;
}
