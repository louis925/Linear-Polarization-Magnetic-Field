#include "parameters.h"

#ifndef _FN_PARAM
typedef struct _fn_param
{
	double *n;
	double tau;
	int j; //J -> j, J-1 = j = J' //change at 2009.11.13
	double k0;
}Fn_Param;
#define _FN_PARAM
#endif

/*
#ifndef _PHYS_COEFF
typedef struct _phys_coeff
{
	double *n;  //nJM = n[J(J+1)/2+M], n[] = {n00,n10,n11,n20,n21,n22...}
	double *A;  //Einstein coefficients, AJJ' = A[J'] (J-1 = J')
	double *C;  //CJ->J' = C[(J-1)J/2+J'], C[] = {C10,C20,C21,C30,C31...}
	double *E;  //exp(-hv/kT)
	double *F;  //2h(vJJ')^3/c^2, FJJ' = F[J'] (J-1 = J')
	double *v;  //frequency from J to J'(vJJ'), vJJ' = v[J'] (J-1 = J')
	double tau; //Optical depth in this loop (tau = tau(0,1->0,PI/2))
}Phys_Coeff;
#define _PHYS_COEFF
#endif
*/