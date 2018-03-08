/*
Made by Louis
NTHU
2009.04.19
*/

// Jp denotes J' = j
// dM denotes delta M

#include <stdio.h>
#include <math.h>
#include "parameters.h"
#include "global_value.h"
#include "integral.h"
#include "physics_function.h"

// Index conversion function for n[]. (J, M) to J(J+1)/2 + |M|.
int indexN(int J, int M) {
	if(M >= 0) {
		return (J*(J+1))/2 + M;
	}
	else {
		return (J*(J+1))/2 - M;
	}
}

// Index conversion function for C[], E[]. (J, J') to (J-1)J/2 + J'.
int indexJJ(int J, int Jp) {
	return (J*(J-1))/2 + Jp;
}

// Einstein A coefficients for the transition (J, M) -> (J' = J-1, M' = M+dm)
double A_coeff(int J, int M, int dM) {
	switch(dM)
	{
	case 0:
		return A[J-1]*(J*J-M*M)/(2*J-1)/J;
		break;
	case -1:
		return A[J-1]*(J+M)*(J-1+M)/2/(2*J-1)/J;
		break;
	case 1:
		return A[J-1]*(J-M)*(J-1-M)/2/(2*J-1)/J;
		break;
	default:
		return 0.0;
	}
}

// Normalized source function
// Given n[], angle, q, J, return the source function / F_JJ' for the transition J -> J' = J-1
double source_f_n(const double n[TOTAL_N], double angle, int q, int J) 
{
	double c, s;
	int M;
	double sum_1, sum_2; // Summation temps
	int A1_n, A0_n;      // Normalized Einstein A coefficients for q = 1 and 0
	int dA;              // Next changed in A

	switch(q)
	{
	case 1: // Perpendicular, there is no direction dependent.
		M = -J;
		sum_1 = 0; // A_JMJ'M' * n_JM
		sum_2 = 0; // A_JMJ'M' * n_J'M'
		A1_n = J*(2*J-1); // Ratio factor of Einstein A coefficient for different M
		dA = 2*J-1;

		while(M < (J-1)) // M: -J ~ J-2
		{
			// J, M -> J-1, M+1
			sum_1 += n[indexN(J,M)]*A1_n;
			sum_2 += n[indexN(J-1,M+1)]*A1_n;

			A1_n -= dA;
			dA--;
			M++;
		}
		return sum_1/(sum_2 - sum_1) / 2;
		break;
	case 0: // Parallel
		c = angle;        // angle = cos(angle between magnetic field and line of sight)
		c *= c;           // c = cos^2
		s = fabs(1-c);
				
		M = -J;
		sum_1 = 0;
		sum_2 = 0;
		A1_n = J*(2*J-1);
		dA = 2*J-1;

		while(M < J)
		{
			A0_n = J*J - M*M;
			sum_1 += n[indexN(J, M)] * (A1_n*c + A0_n*s);
			sum_2 += (n[indexN(J - 1, M + 1)] * A1_n*c + n[indexN(J - 1, M)] * A0_n*s);

			A1_n -= dA;
			dA--;
			M++;
		}
		return sum_1/(sum_2 - sum_1) / 2;
		break;
	default: // Error
		printf("q is out of range in source_f()!\n");
		return 0;
	}
}

// Normalized absorption efficients
double k_f_n(const double n[TOTAL_N], double angle, int q, int j)
{//J -> j, J-1 = j = J', change at 2009.11.13
	double c,s;
	int M;
	int J = j+1;
	int A1_n, A0_n;
	int dA;
	double sum;

	switch(q)
	{
	case 1: // Perpendicular
		M = -J;
		sum = 0;
		A1_n = J*(2*J-1);
		dA = 2*J-1;

		while(M < (J-1))
		{
			sum += (n[indexN(J-1,M+1)]-n[indexN(J,M)])*A1_n;
			
			A1_n -= dA;
			dA--;
			M++;
		}
		
		return sum / (v[j] * v[j]) / (2 * J - 1) / J * A[j]; // Before 2018
		//return 2 * sum / (v[j] * v[j]) / (2 * J - 1) / J * A[j]; // test, original eq in DW's paper, wrong(?)
		break;
	case 0: // Parallel
		c = angle;        // angle = angle between magnetic field and line of sight
		c *= c;           // c = cos^2
		s = fabs(1-c);
		
		M = -J;
		sum = 0;
		A1_n = J*(2*J-1);
		dA = 2*J-1;

		while(M < J)
		{
			A0_n = J*J - M*M;
			sum += ((n[indexN(J-1,M+1)]-n[indexN(J,M)])*A1_n*c + (n[indexN(J-1, M)]-n[indexN(J,M)])*A0_n*s);
			
			A1_n -= dA;
			dA--;
			M++;
		}
		return sum / (v[j] * v[j]) / (2 * J - 1) / J * A[j]; // modified at 2009.12.30
		//return 2 * sum / (v[j] * v[j]) / (2 * J - 1) / J * A[j]; // test on 2018.02.28
		break;
	default: // Error
		printf("q is out of range in k_f()!\n");
		return 0;
	}
}

// Escape probability function
double beta_f(double tau) { return (1 - exp(-1*tau)) / tau; }

// Normalized profile-averaged specific intensity (divided by F_JJ')
double I_pa_n(const double n[TOTAL_N], double angle, int q, double tau_d, int j) 
{ //from J -> J', J-1 = J' = j
	double beta; // Escape probability 
	beta = beta_f(tau_d);
	return source_f_n(n, angle, q, j+1) * (1 - beta) + Br_n[j]/2 * beta; 
}

// Calculate multi-level normalized rate functions R[][]
// Given n[], tau[0], write result to R[][]
void Rate_f_n_cal(const double n[TOTAL_N], double tau, double R[LEVEL_N-1][2]) {
	// R[j][dM] for the transition from J to j = J' = J-1 with |M - M'| = dM
	int j; // j = J' = J-1
	double error;
	unsigned int intervals;
	Fn_Param params;

	params.n = n;
	params.tau = tau; //the total optical depth tau[0] in this case(computation)
	params.k0 = k_f_n(n, cos(TAU_ANG), 0, 0); //[2010.01.22] OBS_ANG was replaced by TAU_ANG
		
	for (j = 0; j < (LEVEL_N - 1); j++) //j = 0 ~ LEVEL_N-2
	{
		params.j = j;
		R[j][0] = integral(i_f_r0m, &params, &error, &intervals) * 3.0 / 2.0;
		// R has already integrated over azimuth angle. Thus, it gains a 1/2 factor
		// The unit of R are the same as that in Deguchi and Watson's paper(1984)
		interval_count += (unsigned long long)intervals;
		
		R[j][1] = integral(i_f_r1m, &params, &error, &intervals) * 3.0 / 4.0;
		// R has already integrated over azimuth angle. Thus, it gains a 1/2 factor
		interval_count += (unsigned long long)intervals;
		
		//R[j][0] = max(R[j][0], 0.);  // Prevent negative intensity [2018.03.01]
		//R[j][1] = max(R[j][1], 0.);  // Prevent negative intensity [2018.03.01]		
	}
}

//[2010.01.22] OBS_ANG was replaced by TAU_ANG
double i_f_r0m(double x, void * params) //for R0, integral of I_pa_n[0]*sin^2
{
	static double *n;
	static double c,s;
	static double tau_d; //tau_d is the tau in this direction(angle) and polariation
	static int j; //J -> j, J-1 = j = J', change at 2009.11.13
	static double k0;
	n = ((Fn_Param *) params)->n;
	j = ((Fn_Param *) params)->j;
	tau_d = ((Fn_Param *) params)->tau;
	k0 =  ((Fn_Param *) params)->k0;
	c = x*x; //x = cos
	s = fabs(1-c);

#if TwoD
	tau_d = tau_d * (k_f_n(n,x,0,j)/k0) * sin(TAU_ANG)*sin(TAU_ANG)/ s;
#elif OneD
	tau_d = tau_d * (k_f_n(n,x,0,j)/k0) * (cos(TAU_ANG)*cos(TAU_ANG)/ c);
#elif Mix
	tau_d = tau_d * (k_f_n(n,x,0,j)/k0) * ((cos(TAU_ANG)*cos(TAU_ANG) + MixRatio*sin(TAU_ANG)*sin(TAU_ANG))/ (c + MixRatio*s));
#else //Isotropic
	tau_d = tau_d * (k_f_n(n,x,0,j)/k0);
#endif
	return s * I_pa_n(n,x,0,tau_d,j);
}

//[2010.01.22] OBS_ANG was replaced by TAU_ANG
double i_f_r1m(double x, void * params) //for R1, integral of I_pa_n[0]*cos^2+I_pa_n[1]
{
	static double *n;
	static double c,s;
	static double tau_d0, tau_d1; //tau_d is the tau in this direction(angle) and polariation
	static int j; //J -> j, J-1 = j = J', change at 2009.11.13
	static double k0;
	n = ((Fn_Param *) params)->n;
	j = ((Fn_Param *) params)->j;
	tau_d0 = ((Fn_Param *) params)->tau;
	tau_d1 = tau_d0;
	k0 =  ((Fn_Param *) params)->k0;
	c = x*x; //x = cos
	s = fabs(1-c);

#if TwoD //use s*s for 2D velocity field, or c*c for 1D velocity field
	tau_d0 = tau_d0 * (k_f_n(n,x,0,j)/k0) * sin(TAU_ANG)*sin(TAU_ANG)/ s;
	tau_d1 = tau_d1 * (k_f_n(n,x,1,j)/k0) * sin(TAU_ANG)*sin(TAU_ANG)/ s;
#elif OneD
	tau_d0 = tau_d0 * (k_f_n(n,x,0,j)/k0) * (cos(TAU_ANG)*cos(TAU_ANG)/ c);
	tau_d1 = tau_d1 * (k_f_n(n,x,1,j)/k0) * (cos(TAU_ANG)*cos(TAU_ANG)/ c);
#elif Mix
	tau_d0 = tau_d0 * (k_f_n(n,x,0,j)/k0) * ((cos(TAU_ANG)*cos(TAU_ANG) + MixRatio*sin(TAU_ANG)*sin(TAU_ANG))/ (c + MixRatio*s));
	tau_d1 = tau_d1 * (k_f_n(n,x,1,j)/k0) * ((cos(TAU_ANG)*cos(TAU_ANG) + MixRatio*sin(TAU_ANG)*sin(TAU_ANG))/ (c + MixRatio*s));
#else //Isotropic
	tau_d0 = tau_d0 * (k_f_n(n,x,0,j)/k0);
	tau_d1 = tau_d1 * (k_f_n(n,x,1,j)/k0);
#endif
	return c * I_pa_n(n, x, 0, tau_d0, j) + I_pa_n(n, x, 1, tau_d1, j);
}

// Construct the optical depth array, tau[][]
// Given n[], TAU, calculate all tau[][].
void tau_array(double TAU, double tau[][2], const double n[TOTAL_N]) 
{
	int j;
	double k0;
	k0 = k_f_n(n, cos(TAU_ANG), 0, 0);  // OBS_ANG was replaced by TAU_ANG [2010.01.22]
	
	for(j = 0; j < (LEVEL_N - 1); j++) {	
#if TwoD  // use s*s for 2D velocity field, or c*c for 1D velocity field
		tau[j][0] = TAU * (k_f_n(n,cos(OBS_ANG),0,j)/k0) * (sin(TAU_ANG)*sin(TAU_ANG))/(sin(OBS_ANG)*sin(OBS_ANG));
		tau[j][1] = TAU * (k_f_n(n,cos(OBS_ANG),1,j)/k0) * (sin(TAU_ANG)*sin(TAU_ANG))/(sin(OBS_ANG)*sin(OBS_ANG));
#elif OneD
		tau[j][0] = TAU * (k_f_n(n,cos(OBS_ANG),0,j)/k0) * (cos(TAU_ANG)*cos(TAU_ANG))/(cos(OBS_ANG)*cos(OBS_ANG));
		tau[j][1] = TAU * (k_f_n(n,cos(OBS_ANG),1,j)/k0) * (cos(TAU_ANG)*cos(TAU_ANG))/(cos(OBS_ANG)*cos(OBS_ANG));
#elif Mix
		tau[j][0] = TAU * (k_f_n(n,cos(OBS_ANG),0,j)/k0) * ((cos(TAU_ANG)*cos(TAU_ANG) + MixRatio*sin(TAU_ANG)*sin(TAU_ANG))/(cos(OBS_ANG)*cos(OBS_ANG) + MixRatio*sin(OBS_ANG)*sin(OBS_ANG)));
		tau[j][1] = TAU * (k_f_n(n,cos(OBS_ANG),1,j)/k0) * ((cos(TAU_ANG)*cos(TAU_ANG) + MixRatio*sin(TAU_ANG)*sin(TAU_ANG))/(cos(OBS_ANG)*cos(OBS_ANG) + MixRatio*sin(OBS_ANG)*sin(OBS_ANG)));
#else  // Isotropic
		tau[j][0] = TAU * (k_f_n(n,cos(OBS_ANG),0,j)/k0);
		tau[j][1] = TAU * (k_f_n(n,cos(OBS_ANG),1,j)/k0);
#endif
	}
}

// Normalized emerge specific intensity, (intensity divided by F = 2hv^3/c^2)
void I_emerge_n(const double n[], const double tau[][2], double I[][2]) {
	for (int j = 0; j < (LEVEL_N - 1); j++) {
		I[j][0] = (source_f_n(n, cos(OBS_ANG), 0, j + 1) - Br_n[j] / 2)*(1 - exp(-1 * tau[j][0]));
		I[j][1] = (source_f_n(n, cos(OBS_ANG), 1, j + 1) - Br_n[j] / 2)*(1 - exp(-1 * tau[j][1]));
	}
}
