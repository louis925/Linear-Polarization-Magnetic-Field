#ifndef _GLOBAL_VALUE_H_
#include <gsl/gsl_integration.h>
extern double A[];            //Einstein coefficients, AJJ' = A[J'] (J-1 = J')
extern double C[];            //CJ->J' = C[(J-1)J/2+J'], C[] = {C10,C20,C21,C30,C31...}
extern double E[];            //exp(-hv/kT)
extern double F[];            //2h(vJJ')^3/c^2, FJJ' = F[J'] (J-1 = J')
extern double v[];            //frequency (GHz) from J to J'(vJJ'), vJJ' = v[J'] (J-1 = J'), GHz = 10^9 Hz
extern double Br_n[];         //normalized Cosmic Blackbody Radiation intensity in each level frequency
extern double energy_level[]; //Energy from ground state(J=0) for each level(J)
extern double a_matrix_i[];   //initial a_matrix[], use 1D array to simulate 2D matrix
extern double a_matrix[];     //use 1D array to simulate 2D matrix
extern unsigned int loop_count;        //count of loops
extern unsigned int interval_count;    //count for integral intervals number
extern double T;              //temperature of cloud
extern gsl_integration_workspace *w;
#define _GLOBAL_VALUE_H_
#endif
