#pragma once

// Print A[]/v[]^2 horizontally
void print_Av2(const double A[LEVEL_N - 1], const double v[LEVEL_N - 1]);

void output_a_matrix(const double* a_matrix, const char* a_matrix_filename);
void print_source_f(const double n[TOTAL_N], const double Br_n[LEVEL_N - 1], const double tau[LEVEL_N - 1][2]);

// Print n[] horizontally
void print_n(const double n[TOTAL_N]);

// Print n[] vertically
void print_nv(const double n[TOTAL_N]);

// Print R[][]
void print_R(const double R[LEVEL_N - 1][2]);

// Print E[]
void print_E(const double E[TRANS_N]);

// Print I_t in the limit of large optical depth
void print_I_limit(const double E[TRANS_N], const double Br_n[LEVEL_N - 1], const double v[LEVEL_N - 1]);

// Compute relative error between a and b
double relative_error(double a, double b);

// Compute relative error between array1 and array2 with length len
// Save relative error in err_array. Return total relative error.
double relative_error_1D(const double array1[], const double array2[], int len, double err_array[]);

// Compute relative error between array1[][2] and array2[][2] with length len
// Save relative error in err_array. Return total relative error.
double relative_error_2D(const double array1[][2], const double array2[][2], int len, double err_array[][2]);

// Check if the error of the test "testname" is within tolerance
// Return 1: pass, 0: fail
int check_error(char *testname, double total_err, double tolerance);

// Pause for Windows
void pause();