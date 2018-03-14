/*
2018.03.14 Modified by Louis Yang (Kavli IPMU)
2009.04.19 Created by Louis Yang (NTHU)
*/

// Jp denotes J' = j
// dM denotes delta M

#define _CRT_SECURE_NO_WARNINGS    //Disable the SECURE WARNINGS for fopen()

#include <stdio.h>
#include <math.h>
#include "parameters.h"
#include "global_value.h"
#include "integral.h"
#include "physics_function.h"
#include "tools.h"

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
// Normalized source function is defined as S_JJ' / F_JJ' where S_JJ' is the source function 
// for J -> J' = J - 1 and F_JJ' = 2hv_JJ'^3 / c^2.
// Given population n[], angle=cos(\theta), polarization q, higher level J, return the value 
// of normalized source function
double source_f_n(const double n[TOTAL_N], double angle, int q, int J) 
{
	double cos2, sin2;
	int M;
	double sum_1, sum_2; // Summation temps
	int A1_n, A0_n;      // Normalized Einstein A coefficients for q = 1 and 0
	int dA;              // Next changed in A

	switch(q)
	{
	case 1: // Perpendicular, there is no directional dependencies.
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
		cos2 = angle * angle; // angle = cos(angle between magnetic field and line of sight)
		sin2 = fabs(1 - cos2); // sin^2(theta)
				
		M = -J;
		sum_1 = 0;
		sum_2 = 0;
		A1_n = J*(2*J-1);
		dA = 2*J-1;

		while(M < J)
		{
			A0_n = J*J - M*M;
			sum_1 += n[indexN(J, M)] * (A1_n*cos2 + A0_n*sin2);
			sum_2 += (n[indexN(J - 1, M + 1)] * A1_n*cos2 + n[indexN(J - 1, M)] * A0_n*sin2);

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

	switch(q) {
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
// beta_f(tau=0) = 1; beta_f(tau=+inf) = 0; 
double beta_f(double tau) { return (1 - exp(-1*tau)) / tau; }

// Normalized profile-averaged specific intensity (divided by F_JJ')
double I_pa_n(const double n[TOTAL_N], double angle, int q, double tau_d, int j) 
{ //from J -> J', J-1 = J' = j
	double beta; // Escape probability 
	beta = beta_f(tau_d);
	return source_f_n(n, angle, q, j+1) * (1 - beta) + Br_n[j]/2 * beta; 
}

// Calculate stimulated emission coefficient factors R[][]
// Given n[], tau[0], write result to R[][]
// Note: R[][] are dimensionless.
void R_cal(const double n[TOTAL_N], double tau, double R[LEVEL_N-1][2]) {
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
	}
}

// Output integrand of R to file 
// Given n[], tau[0], number of sample points n_sample, write integrand of R to file output_filename
void output_integrand_R(const double n[TOTAL_N], double tau, int n_sample, const char* output_filename) {
	Fn_Param params;
	double c;  // Cos

	FILE *of;  // Output file
	if ((of = fopen(output_filename, "w")) == NULL)	{
		printf("Can not open the file '%s' for debug output\n", output_filename);
		pause();
		return;
	}
	printf("Output integrand of R to file %s\n", output_filename);

	params.n = n;
	params.tau = tau; //the total optical depth tau[0] in this case(computation)
	params.k0 = k_f_n(n, cos(TAU_ANG), 0, 0); //[2010.01.22] OBS_ANG was replaced by TAU_ANG
	
	fprintf(of, "cos, ");
	for (int j = 0; j < (LEVEL_N - 1); j++) { //j = 0 ~ LEVEL_N-2
		fprintf(of, "i_R[%d][0], i_R[%d][1], ", j);
	}
	fprintf(of, "\n");
	for (int i = 0; i < n_sample; i++) {
		c = i * 1.0 / (n_sample - 1);
		fprintf(of, "%.6e, ", c);
		for (int j = 0; j < (LEVEL_N - 1); j++) { //j = 0 ~ LEVEL_N-2
			params.j = j;
			fprintf(of, "%.6e, %.6e, ",
				i_f_r0m(c, &params) * 3.0 / 2.0, i_f_r1m(c, &params) * 3.0 / 4.0);
		}
		fprintf(of, "\n");
	}

	if (of) {
		if (fclose(of)) {
			printf("The file '%s' was not closed\n", output_filename);
		}
	}
}

// Integrand of R_JJ' for dM = 0 (R[J'][0]) without the factor 3/2
double i_f_r0m(double x, void * params) //for R0, integral of I_pa_n[0]*sin^2
{
	static double *n;
	static double cos2, sin2;
	static double tau_d; //tau_d is the tau in this direction(angle) and polariation
	static int j; //J -> j, J-1 = j = J', change at 2009.11.13
	static double k0;
	n = ((Fn_Param *) params)->n;
	j = ((Fn_Param *) params)->j;
	tau_d = ((Fn_Param *) params)->tau;
	k0 = ((Fn_Param *)params)->k0;
	cos2 = x*x; //x = cos
	sin2 = fabs(1-cos2);  //sin^2(\theta)

#if TwoD
	tau_d = tau_d * (k_f_n(n,x,0,j)/k0) * sin(TAU_ANG)*sin(TAU_ANG)/ sin2;
#elif OneD
	tau_d = tau_d * (k_f_n(n,x,0,j)/k0) * (cos(TAU_ANG)*cos(TAU_ANG)/ cos2);
#elif Mix
	tau_d = tau_d * (k_f_n(n,x,0,j)/k0) * ((cos(TAU_ANG)*cos(TAU_ANG) + MixRatio*sin(TAU_ANG)*sin(TAU_ANG))/ (cos2 + MixRatio*sin2));
#else //Isotropic
	tau_d = tau_d * (k_f_n(n, x, 0, j) / k0);
#endif
	return sin2 * I_pa_n(n, x, 0, tau_d, j);
}

// Integrand of R_JJ' for dM = 1 (R[J'][1]) without the factor 3/4
double i_f_r1m(double x, void * params) //for R1, integral of I_pa_n[0]*cos^2 + I_pa_n[1]
{
	static double *n;
	static double cos2, sin2;
	static double tau_d0, tau_d1; //tau_d is the tau in this direction(angle) and polariation
	static int j; //J -> j, J-1 = j = J', change at 2009.11.13
	static double k0;
	n = ((Fn_Param *) params)->n;
	j = ((Fn_Param *) params)->j;
	tau_d0 = ((Fn_Param *) params)->tau;
	tau_d1 = tau_d0;
	k0 =  ((Fn_Param *) params)->k0;
	cos2 = x*x; //x = cos
	sin2 = fabs(1- cos2);

#if TwoD //use s*s for 2D velocity field, or c*c for 1D velocity field
	tau_d0 = tau_d0 * (k_f_n(n,x,0,j)/k0) * sin(TAU_ANG)*sin(TAU_ANG)/ sin2;
	tau_d1 = tau_d1 * (k_f_n(n,x,1,j)/k0) * sin(TAU_ANG)*sin(TAU_ANG)/ sin2;
#elif OneD
	tau_d0 = tau_d0 * (k_f_n(n,x,0,j)/k0) * (cos(TAU_ANG)*cos(TAU_ANG)/ cos2);
	tau_d1 = tau_d1 * (k_f_n(n,x,1,j)/k0) * (cos(TAU_ANG)*cos(TAU_ANG)/ cos2);
#elif Mix
	tau_d0 = tau_d0 * (k_f_n(n,x,0,j)/k0) * ((cos(TAU_ANG)*cos(TAU_ANG) + MixRatio*sin(TAU_ANG)*sin(TAU_ANG))/ (cos2 + MixRatio*sin2));
	tau_d1 = tau_d1 * (k_f_n(n,x,1,j)/k0) * ((cos(TAU_ANG)*cos(TAU_ANG) + MixRatio*sin(TAU_ANG)*sin(TAU_ANG))/ (cos2 + MixRatio*sin2));
#else //Isotropic
	tau_d0 = tau_d0 * (k_f_n(n,x,0,j)/k0);
	tau_d1 = tau_d1 * (k_f_n(n,x,1,j)/k0);
#endif
	return cos2 * I_pa_n(n, x, 0, tau_d0, j) + I_pa_n(n, x, 1, tau_d1, j);
}

// Construct the optical depth array, tau[][]
// Given n[], TAU, calculate all tau[][].
void tau_array(double TAU, double tau[][2], const double n[TOTAL_N]) {
	int j;
	double k0;
	double cos_TAU = cos(TAU_ANG);
	double sin_TAU = sin(TAU_ANG);
	double cos_OBS = cos(OBS_ANG);
	double sin_OBS = sin(OBS_ANG);
	double L_factor;  // LVG characteristic scale length factor, L(OBS_ANG) / L(TAU_ANG)

#if TwoD  // use s*s for 2D velocity field, or c*c for 1D velocity field
	L_factor = (sin_TAU*sin_TAU) / (sin_OBS*sin_OBS);
#elif OneD
	L_factor = (cos_TAU*cos_TAU) / (cos_OBS*cos_OBS);
#elif Mix
	L_factor = (cos_TAU*cos_TAU + MixRatio * sin_TAU*sin_TAU) / (cos_OBS*cos_OBS + MixRatio * sin_OBS*sin_OBS);
#else  // Isotropic
	L_factor = 1.;
#endif

	k0 = k_f_n(n, cos(TAU_ANG), 0, 0);  // OBS_ANG was replaced by TAU_ANG [2010.01.22]
	
	for(j = 0; j < (LEVEL_N - 1); j++) {	
		tau[j][0] = TAU * (k_f_n(n, cos_OBS, 0, j) / k0) * L_factor;
		tau[j][1] = TAU * (k_f_n(n, cos_OBS, 1, j) / k0) * L_factor;
	}
}

// Normalized emerge specific intensity, (intensity divided by F = 2hv^3/c^2)
void I_emerge_n(const double n[], const double tau[][2], double I[][2]) {
	for (int j = 0; j < (LEVEL_N - 1); j++) {
		I[j][0] = (source_f_n(n, cos(OBS_ANG), 0, j + 1) - Br_n[j] / 2)*(1 - exp(-1 * tau[j][0]));
		I[j][1] = (source_f_n(n, cos(OBS_ANG), 1, j + 1) - Br_n[j] / 2)*(1 - exp(-1 * tau[j][1]));
	}
}

// Test source_f_n() for isotropic case
int test_source_f_n_iso() {
	double n[TOTAL_N] = { 0. };
	double S_answ[LEVEL_N - 1] = { 0. };  // Expected answer to the normalized source functions
	double S[6][LEVEL_N - 1];  // Output result
	double err_S[6][LEVEL_N - 1];
	double total_err = 0.;
	//FILE* of;

	printf("Testing source_f_n() for isotropic case:\n");

	// Set test population to be be isotropic with equal spacing.
	// Ex: For LEVEL_N = 3, n_2m = 1, n_1m = 2, n_00 = 3.
	for (int j = 0; j < LEVEL_N; j++) {
		for (int m = 0; m <= j; m++) {
			n[indexN(j, m)] = (LEVEL_N - j); // For isotropic case, n[] should be m independent.
		}
	}

	// Expected answer to the normalized source functions for the test n[]
	// It should be S_JJ' / F_JJ' = n_J / 2 / (n_J' - n_J)
	for (int j = 0; j < LEVEL_N - 1; j++) {
		// Assming n_J' - n_J = 1, S_JJ' / F_JJ' = 0.5 * n_J.
		S_answ[j] = 0.5 * (LEVEL_N - 1 - j);
	}

	for (int j = 0; j < LEVEL_N - 1; j++) {
		S[0][j] = source_f_n(n, 0, 0, j + 1);
		S[1][j] = source_f_n(n, 0, 1, j + 1);
		S[2][j] = source_f_n(n, 1, 0, j + 1);
		S[3][j] = source_f_n(n, 1, 1, j + 1);
		S[4][j] = source_f_n(n, 0.5, 0, j + 1);
		S[5][j] = source_f_n(n, 0.5, 1, j + 1);
	}

	// Print results
	printf("S_answ[]: ");
	for (int j = 0; j < LEVEL_N - 1; j++) {	printf("%.1e ", S_answ[j]);	}
	printf("\n");
	for (int i = 0; i < 6; i++) {
		printf("S[%d]    : ", i);
		for (int j = 0; j < LEVEL_N - 1; j++) {
			printf("%.1e ", S[i][j]);
		}
		printf("\n");
	}

	// Compute Error
	total_err = 0.;
	for (int i = 0; i < 6; i++) {
		total_err += relative_error(S_answ, S[i], LEVEL_N - 1, err_S[i]);
	}

	// Print error
	for (int i = 0; i < 6; i++) {
		printf("error[%d]: ", i);
		for (int j = 0; j < LEVEL_N - 1; j++) {
			printf("%.1e ", err_S[i][j]);
		}
		printf("\n");
	}
		
	if (fabs(total_err) < 0.01) {
		printf("Test on source_f_n() for isotropic pass.\n\n");
		return 1;  // Pass
	}
	else {
		printf("Test on source_f_n() for isotropic fail!\n\n");
		pause();
		return 0;  // Fail
	}
}

#if LEVEL_N == 3
// Answer function for the test_source_f_n_3() for LEVEL_N = 3
// Return answer to S[] = {S00, S01, S10, S11}
void answ_source_f_n_3(const double n[TOTAL_N], double c, double S[4]) {
	double n00 = n[indexN(0, 0)];
	double n10 = n[indexN(1, 0)];
	double n11 = n[indexN(1, 1)];
	double n20 = n[indexN(2, 0)];
	double n21 = n[indexN(2, 1)];
	double n22 = n[indexN(2, 2)];
	double cos2 = c * c;    // cos^2(theta)
	double sin2 = 1 - cos2; // sin^2(theta)
	double sum1, sum2;
	// J -> J' = 1 -> 0
	S[0] = 0.5 * (sin2*n10 + cos2 * n11) / (sin2*(n00 - n10) + cos2 * (n00 - n11));  // q = 0
	S[1] = 0.5 * n11 / (n00 - n11);  // q = 1
	// J -> J' = 2 -> 1
	sum1 = n22 + n21 / 2 + n20 / 6;
	sum2 = n11 - n22 + (n10 - n21) / 2 + (n11 - n20) / 6;
	S[2] = 0.5 * (sin2*(n21 + 2 * n20 / 3) + cos2 * sum1)
		/ (sin2*(n11 - n21 + 2 * (n10 - n20) / 3) + cos2 * sum2);  // q = 0
	S[3] = 0.5 * sum1 / sum2;  // q = 1	
}

// Test source_f_n() for LEVEL_N = 3
#define N_test_source_f_n_3 3
int test_source_f_n_3() {
	double n[TOTAL_N] = { 0. };
	double S_answ[N_test_source_f_n_3][4] = { 0. };  // Expected answer to the normalized source functions
	double S[N_test_source_f_n_3][4];                // Output result
	double err_S[N_test_source_f_n_3][4];            // Relative error
	double angle_list[N_test_source_f_n_3] = { 0, 1, 0.5 };  // selected angle=cos(theta) for testing source_f_n()
	double total_err = 0.;
	//FILE* of;

	printf("Testing source_f_n() for N=3 case:\n");

	// Set test population to be be equal spacing at sublevel level.
	// n00=6, n10=5, n11=4, ...
	printf("test n[]: ");
	for (int i = 0; i < TOTAL_N; i++) {
		n[i] = (TOTAL_N - i);
		printf("%.1e ", n[i]);
	}
	printf("\n");

	// Expected answer to the normalized source functions for the test n[]
	for (int i = 0; i < N_test_source_f_n_3; i++) {
		answ_source_f_n_3(n, angle_list[i], S_answ[i]);
	}

	// Test result
	for (int i = 0; i < N_test_source_f_n_3; i++) {
		S[i][0] = source_f_n(n, angle_list[i], 0, 0+1);
		S[i][1] = source_f_n(n, angle_list[i], 1, 0+1);
		S[i][2] = source_f_n(n, angle_list[i], 0, 1+1);
		S[i][3] = source_f_n(n, angle_list[i], 1, 1+1);
	}

	// Print results
	printf("S_answ[][]: ");
	for (int i = 0; i < N_test_source_f_n_3; i++) {
		for (int j = 0; j < 4; j++) {
			printf("%.1e ", S_answ[i][j]);
		}
		printf("| ");
	}
	printf("\n");
	printf("S[][]     : ");
	for (int i = 0; i < N_test_source_f_n_3; i++) {
		for (int j = 0; j < 4; j++) {
			printf("%.1e ", S[i][j]);
		}
		printf("| ");
	}
	printf("\n");

	// Compute Error
	total_err = 0.;
	for (int i = 0; i < N_test_source_f_n_3; i++) {
		total_err += relative_error(S_answ[i], S[i], 4, err_S[i]);
	}

	// Print error
	printf("error     : ");
	for (int i = 0; i < N_test_source_f_n_3; i++) {
		for (int j = 0; j < 4; j++) {
			printf("%.1e ", err_S[i][j]);
		}
		printf("| ");
	}
	printf("\n");

	if (fabs(total_err) < 0.01) {
		printf("Test on source_f_n() for N=3 pass.\n\n");
		return 1;  // Pass
	}
	else {
		printf("Test on source_f_n() for N=3 fail!\n\n");
		pause();
		return 0;  // Fail
	}
}
#endif