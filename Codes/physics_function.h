#include "structure.h"

//In this code, q is defined as polarization
//q = 0 : parallel
//q = 1 : perpendicular
//c.f. Goldreich and Kylafis (1981)

//Notation: J-1 = j = J', transition is J -> j

//Index Conversion Function-----------------------
int indexN(int J, int M);   //(J,M ) to J(J+1)/2+M , for n[]
int indexJJ(int J, int Jp); //(J,J') to (J-1)J/2+J', for C[], E[]
//----------------------------------------------//

double A_coeff(int J, int M, int dM);//sublevel A coefficient for (J,M)->(j=J-1,M+dm)

double source_f_n(const double n[TOTAL_N], double angle, int q, int J); //normalized Source function , = source_function/FJJ'
double k_f_n(const double n[TOTAL_N], double angle, int q, int j); //normalized Absorption efficients
double beta_f(double tau); //escape probability function
double I_pa_n(const double n[TOTAL_N], double angle, int q, double tau_d, int j); //normalized profile-averaged specific intensity(divided by _F)

void Rate_f_n_cal(const double n[TOTAL_N], double tau, double R[LEVEL_N-1][2]); //calculation for multi-level normalized rate function
double i_f_r0m(double x, void * params); //for R0, integral of I_pa_n[0]*sin^2
double i_f_r1m(double x, void * params); //for R1, integral of I_pa_n[0]*cos^2 + I_pa_n[1]

void tau_array(double TAU, double tau[][2], const double n[TOTAL_N]); //calculate all tau[] from tau[0]
