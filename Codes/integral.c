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

// Integrate an even function f from x = -1 to 1 assuming f(-x) = f(x)
double integral(double (* f)(double x, void * params), Fn_Param *params, double *error, unsigned int *intervals) { 
	static gsl_function F;
	double result;
	int status = 0;
	//w = gsl_integration_workspace_alloc (3000); //have already been allocated in main.c
	F.function = f;
	F.params = params;
	
	// Modify the GSL error handler so that the program doesn't abort when encountering numerical error.
	//gsl_set_error_handler(&my_gsl_error);  // Already done in main.c


	// Integrate from x = 0 to 1
#if GSL_INTEGRAL_QNG
	status = gsl_integration_qng(&F, 0, 1, EpsAbs, EpsRel, &result, error, intervals);
#elif GSL_INTEGRAL_CQUAD
	status = gsl_integration_cquad(&F, 0, 1, EpsAbs, EpsRel, ws, &result, error, NULL);
	//printf("%d ", ws->size);
	*intervals = ws->size;
#else
	status = gsl_integration_qags (&F, 0, 1, EpsAbs, EpsRel, Gsl_Integ_Space, w, &result, error); 
	//printf("%d ", w->size);
	*intervals = w->size;
#endif

	if (status == GSL_ESING) {  // If apparent singularity detected
		printf("\nIntegrand: ");
		for (int j = 0; j <= 20; j++) {
			printf("%.3e ", f(j/20., params));
		}
		print_n(params->n);
	}

	return 2 * result;  // Double the result to include contribution from x = -1 to 0
}
