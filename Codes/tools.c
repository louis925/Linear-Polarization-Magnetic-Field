#define _CRT_SECURE_NO_WARNINGS    //Disable the SECURE WARNINGS for fopen()

#include <stdio.h>
#include "parameters.h"
#include "physics_function.h"
#include "tools.h"

void output_a_matrix(const double* a_matrix, const char* a_matrix_filename) {
	printf("Output a_matrix[][] to file %s\n", a_matrix_filename);
	FILE *amf;
	if ((amf = fopen(a_matrix_filename, "w")) == NULL)//Open file for debug a_matrix output
	{
		printf("Can not open the file '%s' for debug output\n", a_matrix_filename);
	}
	int i = 0;
	int j = 0;
	int k = 0;
	while (j < TOTAL_N)	{
		k += TOTAL_N;
		while (i < k) {
			fprintf(amf, "%.2e", a_matrix[i]);
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