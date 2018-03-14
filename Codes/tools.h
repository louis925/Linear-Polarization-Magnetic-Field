#pragma once
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

// Check relative error between array1 and array2 with length len
// Save relative error in err_array. Return total relative error.
double relative_error(const double array1[], const double array2[], int len, double err_array[]);

// Pause for Windows
void pause();