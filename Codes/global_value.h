#ifndef _GLOBAL_VALUE_H_
#include <gsl/gsl_integration.h>
extern double A[];            //Einstein coefficients, AJJ' = A[J'] (J-1 = J')
extern double C[];            // Collisional excitation rates for J -> J', CJJ' = C[(J-1)J/2+J'], C[] = {C10,C20,C21,C30,C31...}
extern double E[];            //exp(-hv/kT)
extern double F[];            //2h(vJJ')^3/c^2, FJJ' = F[J'] (J-1 = J')
extern double v[];            //frequency (GHz) from J to J'(vJJ'), vJJ' = v[J'] (J-1 = J'), GHz = 10^9 Hz
extern double Br_n[];         //normalized Cosmic Blackbody Radiation intensity in each level frequency
extern double S_ext_n[];      // Normalized intensity from external source for J -> J'=J-1, S_ext_JJ' / FJJ' = S_ext_n[J']
							  // S_ext_n = (1 - exp(-TAU_ext)) / (exp(hv/kT_ext) - 1)
extern double energy_level[]; //Energy from ground state(J=0) for each level(J)
extern double a_matrix_i[];   //initial a_matrix[], use 1D array to simulate 2D matrix
extern double a_matrix[];     //use 1D array to simulate 2D matrix
extern unsigned long long int loop_count;        //count of loops
extern unsigned long long int interval_count;    //count for integral intervals number
extern double T;              //temperature of cloud
extern gsl_integration_workspace *w;
extern gsl_integration_cquad_workspace *ws;  // GSL integration workspace for CQUAD method
#define _GLOBAL_VALUE_H_
#endif
