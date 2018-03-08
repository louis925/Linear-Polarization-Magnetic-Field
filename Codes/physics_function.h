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

// Normalized source function (divided by FJJ')
double source_f_n(const double n[TOTAL_N], double angle, int q, int J); 

// Normalized absorption coefficients (divided by FJJ')
double k_f_n(const double n[TOTAL_N], double angle, int q, int j);

// Escape probability
double beta_f(double tau);

// Normalized profile-averaged specific intensity (divided by FJJ')
double I_pa_n(const double n[TOTAL_N], double angle, int q, double tau_d, int j);

// Calculate multi-level normalized rate functions R[][]
// Given n[], tau[0], write result to R[][]
void Rate_f_n_cal(const double n[TOTAL_N], double tau, double R[LEVEL_N-1][2]); 

double i_f_r0m(double x, void * params); //for R0, integral of I_pa_n[0]*sin^2
double i_f_r1m(double x, void * params); //for R1, integral of I_pa_n[0]*cos^2 + I_pa_n[1]

// Construct the optical depth array, tau[][]
// Given n[], TAU, calculate all tau[][].
void tau_array(double TAU, double tau[][2], const double n[TOTAL_N]);

// Normalized emerge specific intensity, (intensity divided by F = 2hv^3/c^2)
void I_emerge_n(const double n[], const double tau[][2], double I[][2]);
