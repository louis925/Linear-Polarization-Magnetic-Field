/*
GSL Integration Example: http://www.gnu.org/software/gsl/manual/html_node/Numerical-integration-examples.html
*/

#include <stdio.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include "global_value.h"
#include "integral.h"
#include "parameters.h"
#include "tools.h"

void my_gsl_error(const char * reason, const char * file, int line, int gsl_errno)
{
	printf("\nGSL Error %d: ", gsl_errno);
	gsl_stream_printf("", file, line, reason);
	fflush(stdout);
	//fprintf(stderr, "Default GSL error handler invoked.\n");
	//fflush(stderr);
	pause();
	//abort();
}

double integral(double (* f)(double x, void * params), Fn_Param *params, double *error, unsigned int *intervals)
{ 
	// We use integral() to integral only the inclination angle (from 0 to pi) 
	
	static gsl_function F;
	double result;
	int status = 0;
	//w = gsl_integration_workspace_alloc (3000); //have already been allocated in main.c
	F.function = f;
	F.params = params;
	
	// Modify the GSL error handler so that the program doesn't crash when it encounters numerical error.
	gsl_set_error_handler(&my_gsl_error);

	// Integrate from 0 to 1
#if GSL_INTEGRAL_QNG
	status = gsl_integration_qng(&F, 0, 1, EpsAbs, EpsRel, &result, error, intervals);	
#else
	status = gsl_integration_qags (&F, 0, 1, EpsAbs, EpsRel, Gsl_Integ_Space, w, &result, error); 
	//printf("%d ", w->size);
	*intervals = w->size;
#endif

	if (status == 21) {
		printf("\nIntegrand: ");
		for (int j = 0; j <= 20; j++) {
			printf("%.3e ", f(j/20., params));
		}
	}

	return result*2;
}
