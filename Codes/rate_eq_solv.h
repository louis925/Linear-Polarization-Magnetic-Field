//In this code, q is defined as polarization
//q = 0 : perpendicular
//q = 1 : parallel

// Solve the transition rate equations at optical depth TAU.
// Write the resulting population to n[].
void rate_eq_solve(double n[TOTAL_N], double TAU);

// Initialize n[] using thermal equilibrium at temperature T(K).
void n_initial_cal(double n[TOTAL_N], double T);

// Initialize a_matrix[]
// Fill first row with particle number conservation. Fill the rest rows with collisional excitation rates C[].
void a_matrix_initialize(double a_matrix[TOTAL_N*TOTAL_N]);//2009.11.12 Check OK (for C coeff)

#if LEVEL_N == 3
// Test a_matrix_initialize() for LEVEL_N = 3
int test_a_matrix_initialize_3();

// Test rate_eq_fill() for LEVEL_N = 3
// Test filling of A[] terms to a_matrix[] without radiation
int test_rate_eq_fill_3();
#endif

// Test rate_eq_solve with only C[] terms and particle conservation (a_matrix_i)
// This is not a test on rate_eq_solve() function but a mathematical check of the solution
int test_solve_a_matrix_i();