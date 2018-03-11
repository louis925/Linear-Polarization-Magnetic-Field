/*
GSL Liner Algebra Solving Example:
http://www.gnu.org/software/gsl/manual/html_node/Linear-Algebra-Examples.html

Made by Louis
NTHU
2010.01.19
*/

#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <math.h>
#include "physics_function.h"
#include "parameters.h"
#include "global_value.h"
#include "tools.h"
//#define N_loop 50

#include "rate_eq_solv.h"

// Fill a_matrix with A[] and R[], call AR_fill_0 and _1
void rate_eq_fill(double a_matrix[TOTAL_N*TOTAL_N], const double R[LEVEL_N-1][2]); 

// Fill a_matrix with A[] and R[], for M = 0
void AR_fill_0(double a_row[TOTAL_N], const double R[LEVEL_N-1][2], int J);

// Fill a_matrix with A[] and R[], for M > 0
void AR_fill_1(double a_row[TOTAL_N], const double R[LEVEL_N-1][2], int J, int M);

void rate_eq_solve(double n[TOTAL_N], double TAU) {
	double b[TOTAL_N] = {Nt, 0.0};  
	double n_last[TOTAL_N];           // Result n[] in last step for comparing to next step
	double R[LEVEL_N-1][2] = {{0.0}}; // Results of Rate_f_n()
	int converge;                     // Weither n[] converge
	int s, j;
	int l;                            // Number of loops
	int R_sign = 0;                   // Sign of R[][]. 1: R contains negative, 0; all positive

#if 0
	int i,k; //debug use
	FILE *amf; //debug use
#endif

	gsl_permutation *gsl_p;
	gsl_matrix_view gsl_a_m;
	gsl_vector_view gsl_b_v;
	gsl_vector_view gsl_n_v;
	gsl_p = gsl_permutation_alloc (TOTAL_N);
	gsl_a_m = gsl_matrix_view_array (a_matrix, TOTAL_N, TOTAL_N);
	gsl_b_v = gsl_vector_view_array (b, TOTAL_N);
	gsl_n_v = gsl_vector_view_array (n, TOTAL_N);

	l = 0;
	do {
		Rate_f_n_cal(n, TAU, R); // Calculate the R[][] coefficients

#if CUTOFF_R
		for (int i = 0; i < LEVEL_N - 1; i++) {
			R[i][0] = max(R_MIN, R[i][0]);
			R[i][1] = max(R_MIN, R[i][1]);
		}
#endif
			
#if SHOW_R
		//print R[][] for debug*****
		if(l == 0) {
			printf("\nR[][] initial: ");//*****
			print_R(R);
			printf("\n");//*****/
		}
#endif	

		// Copy n[] to n_last[]
		for(j = 0; j < TOTAL_N; j++) { n_last[j] = n[j]; }
		// Copy a_matrix_i[] to a_matrix[]
		for(j = 0; j < TOTAL_N*TOTAL_N; j++) { a_matrix[j] = a_matrix_i[j];	}

		//print_n(n_last);
		
		// Fill rate equations(a_matrix[]) with A[] and R[] ____________________________
		rate_eq_fill(a_matrix, R);
		// ___________________________________________________________________________//

#if 1
		//check a_matrix[]********
		if (TAU == TAU_START && l == 0) {
			printf("[Debug]a_matrix[][] output.\n");
			output_a_matrix(a_matrix, "a_matrix_0[].csv");
		}
#endif	
		// _____________________________________________________________________________
		// Compute n[] by solving a[TOTAL_N][TOTAL_N] x n[TOTAL_N] = b[TOTAL_N] problem
		// using the LU decomposition method
		gsl_linalg_LU_decomp (&gsl_a_m.matrix, gsl_p, &s);
		gsl_linalg_LU_solve (&gsl_a_m.matrix, gsl_p, &gsl_b_v.vector, &gsl_n_v.vector);
		loop_count++;
		// ___________________________________________________________________________//
		
#if 0
		//if(tau[0][0] > 0.018) //debug for LEVEL_N=7 tau=0.01905460718
		//{
			printf("%d:\n", i);//*****
			j=0;//print n[] for debug*****
			while(j < TOTAL_N)
			{
				printf("%.7e, ", n[j]);
				j++;
			}
			printf("\n");//*****/
		//}
#endif	

		// Check whether it converges
		for (j = 0; j < TOTAL_N; j++) {
			if( fabs(n[j] - n_last[j]) < fabs(n[j]) * REL_PREC)	{
				converge = 1; // OK, Converge
			}
			else {
				converge = 0; // Doesn't converge, run again
				j = TOTAL_N;  // Skip checking the rest of j
			}		
		}
		l++;
	}
	while(!converge); // Stop when converge. //l < N_loop
	
	// Check if R[][] are positive
	R_sign = 0;
	for (int i = 0; i < LEVEL_N - 1; i++) {
		if (R[i][0] < 0 || R[i][1] < 0) { R_sign = 1; }
	}
	if (R_sign == 1) { 
		printf("R[][] constains negative values! Solution may not valid.\n"); 
#if 1
		//check a_matrix[]********
		printf("[Debug]a_matrix[][] output.\n");
		output_a_matrix(a_matrix, "a_matrix_n[].csv");
#endif
		getchar();
	}

#if SHOW_R
	printf("R[][] final  : ");
	print_R(R);
	printf("\n");
#endif	

#if 0
	//print n[] for debug*****
	//printf("%d:\n", i);
	printf("\nn[] for debug\n");
	j=0;
	while(j < TOTAL_N)
	{
		printf("%+.3e, ", n[j]);
		j++;
	}
	printf("\n");
#endif

	gsl_permutation_free (gsl_p);	
}

// Initialize n[] using thermal equilibrium at temperature T(K).
void n_initial_cal(double n[TOTAL_N], double T) {
	double e[LEVEL_N];
	double x, et = 0;
	int j, m;
	x = -100.0 * h_CONST * LIGHT_SPEED / k_CONST / T;
	for (j = 0; j < LEVEL_N; j++) {
		e[j] = exp(x*(energy_level[j] - energy_level[0]));
		et += e[j] * (2 * j + 1);
	}
	for (j = 0; j < LEVEL_N; j++) {
		for (m = 0; m <= j; m++) {
			n[indexN(j, m)] = Nt * e[j] / et;
		}
	}
}

// Initialize a_matrix[]
// Fill first row with particle number conservation. Fill the rest rows with collisional excitation rates C[].
void a_matrix_initialize(double a_matrix[TOTAL_N*TOTAL_N]) {
	// a_matrix : a[(J,M)][(J',M')] = a[J(J+1)/2+M][J'(J'+1)/2+M'] := a_matrix[(J(J+1)/2+M)*TOTAL_N + J'(J'+1)/2+M']
	// 2009.11.12 Check OK (for C coeff)
	// Row index: (J, M)
	// Column index: (J', M'), J' also denoted as j (lower case)
	int J, M; // row in the matrix
	int j;    // J'
	int i;    // filling position of a_matrix[i]
	int k;    // upper limit of filling position i
	int i2;   // position of a_matrix for (J', M') = (J, M)
	double cjj;

	// Cleanup the a_matrix
	for (i = 0; i < TOTAL_N*TOTAL_N; i++) {
		a_matrix[i] = 0.;
	}

	// Row 0: particle number conservation --------
	i = 0;
	k = 0;
	for (j = 0; j < LEVEL_N; j++) {
		k += (j + 1); // Add the number of different |M| given j
		a_matrix[i] = 1.0 * PC_FACTOR; //for M'=0
		i++;
		while (i < k) {
			a_matrix[i] = 2.0 * PC_FACTOR; //for M'>0
			i++;
		}
	}	

	// Now the indices i and k should both be at a[1][0]

	// Row 1 to TOTAL_N - 1: Collisional excitation rates C[] -----------	
	for (J = 1; J < LEVEL_N; J++) {
		for (M = 0; M < (J + 1); M++) {
			// For one row, have the same (J, M)------
			// i should be at a_matrix[i] = a[(J,M)][0]

			i2 = i + indexN(J, M); // Position of a_matrix[i2] = a[(J,M)][(J,M)]
			a_matrix[i2] = 0.0; // Coefficient for n[(J, M)]

			//run over all the element in this row, from a[(J,M)][0] to a[(J,M)][TOTAL_N-1]
			j = 0; //j = J', in this loop

			while (j < J)//for each J' < J, from a[(J,M)][J'=0] to a[(J,M)][J'=J-1]
			{
				k += (j + 1);

				cjj = C[indexJJ(J, j)]; //CJJ'
				a_matrix[i2] -= cjj; //JM -> J'M'

				cjj *= (E[indexJJ(J, j)] / (2 * j + 1));
				a_matrix[i] = cjj; //for M'=0
				i++;
				while (i < k) {
					a_matrix[i] = 2.0 * cjj; //for M'>0
					i++;
				}
				j++;
			}
			//Now j = J

			//for J' = J, fill a[][] with zero, unless M = M'
			k += (j + 1);
			while (i < k) {
				if (i != i2) {
					a_matrix[i] = 0.0;
				}
				i++;
			}
			j++;
			//Now j = J+1

			while (j < LEVEL_N)//for each J'=j > J, from a[(J,M)][J'=J+1] to a[(J,M)][J'=LEVEL_N-1]
			{
				k += (j + 1);

				cjj = C[indexJJ(j, J)] / (2 * J + 1);
				a_matrix[i2] -= (cjj * (2 * j + 1) * E[indexJJ(j, J)]);

				a_matrix[i] = cjj; //for M'=0
				i++;
				while (i < k) {
					a_matrix[i] = 2.0 * cjj; //for M'>0
					i++;
				}
				j++;
			}
			//------------------------------------//
		}
	}
	//-------------------------------------------------------//
}

void rate_eq_fill(double a_matrix[TOTAL_N*TOTAL_N], const double R[LEVEL_N - 1][2]) {
	// Fill the rate equations (a_matrix[]) with the contribution from A[] and R[]
	// rates function R[] for each level needed to be calculate first before call this function
	int i;
	int J, M;

	i = TOTAL_N; // Start filling a_matrix from a_matrix[1][0]. Skip the row 0.

	//for J = 1 ~ J = LEVEL_N-1
	for (J = 1; J <= LEVEL_N - 1; J++) {
		M = 0;
		AR_fill_0(&a_matrix[i], R, J);
		i += TOTAL_N;
		M++;

		while (M <= J) //change at 2010/01/19
		{
			AR_fill_1(&a_matrix[i], R, J, M);
			i += TOTAL_N;
			M++;
		}
	}
}

void AR_fill_0(double a_row[TOTAL_N], const double R[LEVEL_N-1][2], int J) //M = 0
{
	double At[6]; //temp. for A coefficients, At[] = {AJM(J-1)M	, AJM(J-1)(M+1), AJM(J-1)(M-1), A(J+1)MJ, A(J+1)(M-1)JM, A(J+1)(M+1)JM}
	int j;
	At[0] = A_coeff(J  ,  0,  0); //AJ0J-10
	At[1] = A_coeff(J  ,  0,  1); //AJ0J-11
	
	j = indexN(J,0);
	a_row[j] -= (  At[0]*(1 + R[J-1][0]) + 2*At[1]*(1 + R[J-1][1]) ); //****R[J]
	if(J < (LEVEL_N - 1))
	{
		At[3] = A_coeff(J+1,  0,  0); //AJ+10J0
		At[4] = A_coeff(J+1, -1,  1); //AJ+1-1J0
		a_row[j] -= ( At[3]*R[J][0] + 2*At[4]*R[J][1] );
	}

	j = indexN(J-1,0);
	a_row[j  ] +=   At[0]*R[J-1][0];
	a_row[j+1] += 2*At[1]*R[J-1][1];
	
	if(J < (LEVEL_N - 1))
	{
		j = indexN(J+1,0);
		a_row[j  ] +=   At[3]*(1 + R[J][0]);
		a_row[j+1] += 2*At[4]*(1 + R[J][1]);
	}
}

void AR_fill_1(double a_row[TOTAL_N], const double R[LEVEL_N-1][2], int J, int M) //M > 0
{
	double At[6]; //temp. for A coefficients, At[] = {AJM(J-1)M	, AJM(J-1)(M+1), AJM(J-1)(M-1), A(J+1)MJ, A(J+1)(M-1)JM, A(J+1)(M+1)JM}
	int j;
	At[0] = A_coeff(J  , M  ,  0); //dM = 0
	At[1] = A_coeff(J  , M  ,  1); //dM = +1
	At[2] = A_coeff(J  , M  , -1); //dM = -1
	
	j = indexN(J,M);
	a_row[j] -= (  At[0] * (1 + R[J-1][0]) + (At[1] + At[2]) * (1 + R[J-1][1]));
	if(J < (LEVEL_N - 1))
	{
		At[3] = A_coeff(J+1, M  ,  0); //dM = 0
		At[4] = A_coeff(J+1, M-1,  1); //dM = +1
		At[5] = A_coeff(J+1, M+1, -1); //dM = -1
		a_row[j] -= ( At[3] *   R[J  ][0]  + (At[4] + At[5]) *      R[J  ][1] );
	}
	
	j = indexN(J-1,M);
	a_row[j  ] += At[0]*R[J-1][0]; //n(J-1,M)
	a_row[j+1] += At[1]*R[J-1][1]; //n(J-1,M+1)
	a_row[j-1] += At[2]*R[J-1][1]; //n(J-1,M-1)
	
	if(J < (LEVEL_N - 1))
	{
		j = indexN(J+1,M);
		a_row[j  ] += At[3]*(1 + R[J][0]); //n(J+1,M)
		a_row[j+1] += At[5]*(1 + R[J][1]); //n(J+1,M+1) //Modified at 2010/01/19
		a_row[j-1] += At[4]*(1 + R[J][1]); //n(J+1,M-1) //Modified at 2010/01/19
	}
}

// Test a_matrix_initialize()
int test_a_matrix_initialize_3() {
	double a_matrix[TOTAL_N*TOTAL_N] = { 11., 11., 11., 11., 11., 11., 11., 11., 11.};
	double C10 = C[indexJJ(1, 0)];
	double C20 = C[indexJJ(2, 0)];
	double C21 = C[indexJJ(2, 1)];
	double E10 = E[indexJJ(1, 0)];
	double E20 = E[indexJJ(2, 0)];
	double E21 = E[indexJJ(2, 1)];
	double a_matrix_ans[TOTAL_N*TOTAL_N] = {
		1 * PC_FACTOR, 1 * PC_FACTOR, 2 * PC_FACTOR, 1 * PC_FACTOR, 2 * PC_FACTOR, 2 * PC_FACTOR,
		C10 * E10, -(C10 + 5 * C21*E21 / 3), 0., C21 / 3, 2 * C21 / 3, 2 * C21 / 3,
		C10 * E10, 0., -(C10 + 5 * C21*E21 / 3), C21 / 3, 2 * C21 / 3, 2 * C21 / 3,
		C20 * E20, C21 * E21 / 3, 2 * C21 * E21 / 3, -(C20 + C21), 0., 0.,
		C20 * E20, C21 * E21 / 3, 2 * C21 * E21 / 3, 0., -(C20 + C21), 0.,
		C20 * E20, C21 * E21 / 3, 2 * C21 * E21 / 3, 0., 0., -(C20 + C21)
	};
	double err_a_matrix[TOTAL_N*TOTAL_N] = { 0. };  // Relative error in a_matrix
	double total_err = 0.;
	a_matrix_initialize(a_matrix);
	output_a_matrix(a_matrix, "a_matrix_test_init[].csv");
	output_a_matrix(a_matrix_ans, "a_matrix_test_init_answ[].csv");
	for (int i = 0; i < TOTAL_N*TOTAL_N; i++) {
		err_a_matrix[i] = a_matrix[i] - a_matrix_ans[i];
		if (err_a_matrix[i] != 0.) {
			err_a_matrix[i] /= (fabs(a_matrix[i]) + fabs(a_matrix_ans[i]));
		}
		total_err += err_a_matrix[i];
	}
	output_a_matrix(err_a_matrix, "a_matrix_test_init_error[].csv");
	if (total_err > 0.01) {
		return 0;
	}
	else {
		return 1;
	}
}
