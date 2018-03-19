#define _CRT_SECURE_NO_WARNINGS    //Disable the SECURE WARNINGS for fopen()

#include <stdio.h>
#include "parameters.h"
#include "physics_function.h"
#include "tools.h"

// Print an array
void print_array(const double A[], int len) {
	for (int i = 0; i < len; i++) {
		printf("%.2e ", A[i]);
	}
}

// Print A[]/v[]^2 horizontally
void print_Av2(const double A[LEVEL_N - 1], const double v[LEVEL_N - 1]) {
	printf("A[]/v[]^2: ");
	for (int i = 0; i < LEVEL_N - 1; i++) {
		printf("%.3e ", A[i] / (v[i] * v[i]));
	}
	printf("\n\n");
}

void output_a_matrix(const double* a_matrix, const char* a_matrix_filename) {
	printf("Output a_matrix[][] to file %s\n", a_matrix_filename);
	FILE *amf;
	if ((amf = fopen(a_matrix_filename, "w")) == NULL)//Open file for debug a_matrix output
	{
		printf("Can not open the file '%s' for debug output\n", a_matrix_filename);
		pause();
	}
	int i = 0;
	int j = 0;
	int k = 0;
	while (j < TOTAL_N)	{
		k += TOTAL_N;
		while (i < k) {
			fprintf(amf, "%.6e", a_matrix[i]);
			fprintf(amf, ",");
			i++;
		}
		fprintf(amf, "\n");
		j++;
	}
	if (amf) {
		if (fclose(amf)) {
			printf("The file '%s' was not closed\n", a_matrix_filename);
		}
	}
}

void print_source_f(const double n[TOTAL_N], const double Br_n[LEVEL_N - 1], const double tau[LEVEL_N - 1][2]) {
	printf("j:     S_f(0)     S_f(1)     Br_n/2     tau[0]     tau[1]\n");
	for (int j = 0; j < (LEVEL_N - 1); j++) {
		printf("%d: % .3e % .3e % .3e % .3e % .3e\n", j,
			source_f_n(n, cos(OBS_ANG), 0, j + 1), source_f_n(n, cos(OBS_ANG), 1, j + 1),
			Br_n[j] / 2, tau[j][0], tau[j][1]);		
	}
}

// Print n[] horizontally
void print_n(const double n[TOTAL_N]) {
	printf("n[]: ");
	for (int i = 0; i < TOTAL_N; i++) {
		printf("%.3e ", n[i]);
	}
	printf("\n");
}

// Print n[] vertically
void print_nv(const double n[TOTAL_N]) {
	for (int i = 0; i < TOTAL_N; i++) {
		printf("%d: %.5e\n", i, n[i]);
	}
}

// Print R[][]
void print_R(const double R[LEVEL_N - 1][2]) {
	for (int i = 0; i < LEVEL_N - 1; i++) {
		printf("%.2e %.2e ", R[i][0], R[i][1]);
	}
}

// Print E[]
void print_E(const double E[TRANS_N]) {
	printf("E[] Ratio of n[]/n[0] in the limit of large optical depth:\n");
	for (int j = 0; j < (LEVEL_N - 1); j++)	{
		printf("%.3e ", E[indexJJ(j + 1, 0)]);		
	}
	printf("\n");
}

void print_I_limit(const double E[TRANS_N], const double Br_n[LEVEL_N - 1], const double v[LEVEL_N - 1]) {
	printf("I_limit[](K) I_t in the limit of large optical depth:\n"); //[2012.11.17]Included the contribution from background radiation
	for (int j = 0; j < (LEVEL_N - 1); j++) {
		printf("%.3e ", (1 / (1 / E[indexJJ(j + 1, j)] - 1) - Br_n[j])*h_CONST*v[j] / k_CONST * 1E9);
	}
	printf("\n");
}

double relative_error(double a, double b) {
	if (a == b) { return 0.; }
	else { return (a - b) / (fabs(a) + fabs(b)); }
}

// Compute relative error between array1 and array2 with length len
// Save relative error in err_array. Return total relative error.
double relative_error_1D(const double array1[], const double array2[], int len, double err_array[]) {
	double total_err = 0.;
	for (int i = 0; i < len; i++) {
		err_array[i] = relative_error(array1[i], array2[i]);
		total_err += fabs(err_array[i]);
	}
	return total_err;
}

// Compute relative error between array1[][2] and array2[][2] with length len
// Save relative error in err_array. Return total relative error.
double relative_error_2D(const double array1[][2], const double array2[][2], int len, double err_array[][2]) {
	double total_err = 0.;
	for (int i = 0; i < len; i++) {
		err_array[i][0] = relative_error(array1[i][0], array2[i][0]);
		err_array[i][1] = relative_error(array1[i][1], array2[i][1]);
		total_err += fabs(err_array[i][0]);
		total_err += fabs(err_array[i][1]);
	}
	return total_err;
}

// Check if test pass or fail
int check_error(char *testname, double total_err, double tolerance) {
	printf("Test on %s ", testname);
	if (fabs(total_err) < tolerance) {
		printf("pass.\n\n");
		return 1;  // Pass
	}
	else {
		printf("fail with error %.2e!\n\n", total_err);
		pause();
		return 0;  // Fail
	}
}

// Set all elements in array A[] to value val
// Size of A[] = len
void erase(double A[], int len, double val) {
	for (int i = 0; i < len; i++) {
		A[i] = val;
	}
}

// Reduce the a_matrix for isotropic case
// Reduce a_matrix by merging the contribution from all the sublevels of each J into one (ex: n1 = n10 + 2*n11)
// a_matrix[]  : ((N + 1)*N) / 2 * ((N + 1)*N) / 2
// a_matrix_red: N * N
void reduce(const double a_matrix[], double a_matrix_red[], int N) {
	int N_row = ((N + 1)*N) / 2;
	for (int J = 0; J < N; J++) {
		int M = 0; // Only deal rows with M = 0. Igonre M != 0 rows.
		for (int j = 0; j < N; j++) {
			double a = 0.;
			for (int m = 0; m <= j; m++) {
				a += a_matrix[indexN(J, M) * N_row + indexN(j, m)];
			}
			a_matrix_red[J * N + j] = a;
		}		
	}
}

// Pause for Windows
void pause() {
#ifdef WINDOWS
	getchar();
#endif
	return;
}


// Test =======================================================
// Test reduce()
int test_reduce() {
	char testname[] = "reduce()";
#define N_test_reduce 3
	int N = N_test_reduce;
	int N_row = ((N_test_reduce + 1)*N_test_reduce) / 2;
	double a_matrix[((N_test_reduce + 1)*N_test_reduce) / 2 * ((N_test_reduce + 1)*N_test_reduce) / 2];
	double a_matrix_red[N_test_reduce * N_test_reduce];
	double a_matrix_red_answ[N_test_reduce * N_test_reduce] = {1, 2, 3, 1, 2, 3, 1, 2, 3 };
	double err_a[N_test_reduce * N_test_reduce];
	double total_err = 0.;

	printf("Testing on %s\n", testname);

	for (int i = 0; i < N_row * N_row; i++) {
		a_matrix[i] = 1.; 
	}

	reduce(a_matrix, a_matrix_red, N);
	total_err = relative_error_1D(a_matrix_red, a_matrix_red_answ, N*N, err_a);

	printf("a_matrix[]: ");
	print_array(a_matrix, N_row * N_row);
	printf("\n");

	printf("a_matrix_red[]     : ");
	print_array(a_matrix_red, N * N);
	printf("\n");

	printf("a_matrix_red_answ[]: ");
	print_array(a_matrix_red_answ, N * N);
	printf("\n");

	return check_error(testname, total_err, 1E-8);
}