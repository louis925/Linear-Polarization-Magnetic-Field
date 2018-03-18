#include "structure.h"

/*
In this code, q is defined as polarization
q = 0 : parallel
q = 1 : perpendicular
c.f. Goldreich and Kylafis (1981)

Notation: j = J', transition is J -> j
*/

// Index conversion functions
int indexN(int J, int M);   //(J,M ) to J(J+1)/2+M , for n[]
int indexJJ(int J, int Jp); //(J,J') to (J-1)J/2+J', for C[], E[]

// Einstein A coefficients for the transition (J, M) -> (J' = J-1, M' = M+dm)
double A_coeff(int J, int M, int dM);

// Normalized source function
// Normalized source function is defined as S_JJ' / F_JJ' where S_JJ' is the source function 
// for J -> J' = J - 1 and F_JJ' = 2hv_JJ'^3 / c^2.
// Given population n[], angle=cos(\theta), polarization q, higher level J, return the value 
// of normalized source function
double source_f_n(const double n[TOTAL_N], double angle, int q, int J);

// Normalized absorption coefficients
// Absorption coefficients divided by the normalized factor, phi(v-v_JJ') * c_light_speed^2 * 3 / 8pi,
// for the transition j+1 -> j with polarization q, viewing angle, angle = cos(theta), and population n[],
// where phi(v-v_JJ') is the line profile function.
// Note lower j = J' = J-1.
double k_f_n(const double n[TOTAL_N], double angle, int q, int j);

// Escape probability
// beta_f(tau=0) = 1; beta_f(tau=+inf) = 0; 
double beta_f(double tau);

// Normalized profile-averaged specific intensity (divided by FJJ')
double I_pa_n(const double n[TOTAL_N], double angle, int q, double tau_d, int j);

// Calculate stimulated emission coefficient factors R[][]
// Given n[], tau[0], write result to R[][]
// Note: R[][] are dimensionless.
void R_cal(const double n[TOTAL_N], double tau, double R[LEVEL_N-1][2]);

// Output integrand of R to file 
// Given n[], tau[0], number of sample points n_sample, write integrand of R to file output_filename
void output_integrand_R(const double n[TOTAL_N], double tau, int n_sample, const char* output_filename);

// Integrand of R_JJ' for dM = 0 (R[J'][0]) without the factor 3/2
// Return sin^2 * I_pa_n[0]
double i_f_r0m(double x, void * params);

// Integrand of R_JJ' for dM = 1 (R[J'][1]) without the factor 3/4
// Return cos^2 * I_pa_n[0] + I_pa_n[1]
double i_f_r1m(double x, void * params);

// Construct the optical depth array, tau[][]
// Given n[], TAU, calculate all tau[][].
void tau_array(double TAU, double tau[][2], const double n[TOTAL_N]);

// Normalized emerge specific intensity, (intensity divided by F = 2hv^3/c^2)
void I_emerge_n(const double n[], const double tau[][2], double I[][2]);

// Testing ====================================================================================

// Test source function for isotropic case
int test_source_f_n_iso();

#if LEVEL_N == 3
// Test source_f_n() for LEVEL_N = 3
int test_source_f_n_3();

// Test k_f_n() for LEVEL_N = 3
int test_k_f_n_3();
#endif

// Test beta_f()
int test_beta_f();

// Test R_cal() for isotropic case
int test_R_cal_iso();