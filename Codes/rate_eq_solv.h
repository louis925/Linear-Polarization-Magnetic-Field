//In this code, q is defined as polarization
//q = 0 : perpendicular
//q = 1 : parallel

void rate_eq_solve(double n[TOTAL_N], double TAU);

// Initialize n[] using thermal equilibrium at temperature T(K).
void n_initial_cal(double n[TOTAL_N], double t);

// Initialize a_matrix[]
// Fill first row with particle number conservation. Fill the rest rows with collisional excitation rates C[].
void a_matrix_initialize(double a_matrix[TOTAL_N*TOTAL_N]);//2009.11.12 Check OK (for C coeff)

// Test a_matrix_initialize() for LEVEL_N = 3
int test_a_matrix_initialize_3();

// Test rate_eq_solve with only C[] terms and particle conservation (a_matrix_i)
// This is not a test on rate_eq_solve() function but a mathematical check of the solution
int test_solve_a_matrix_i();