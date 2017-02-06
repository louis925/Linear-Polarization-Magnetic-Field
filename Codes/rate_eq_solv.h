//In this code, q is defined as polarization
//q = 0 : perpendicular
//q = 1 : parallel

void rate_eq_solve(double n[TOTAL_N], double TAU);
//void rate_eq_solve(double n[TOTAL_N], double tau[LEVEL_N-1][2]);
//void rate_eq_solve(double n[TOTAL_N], double tau[LEVEL_N-1][2], double Ri[LEVEL_N-1][2]);
void n_initial_cal(double n[TOTAL_N], double t);
void a_matrix_initialize(double a_matrix[TOTAL_N*TOTAL_N]);//2009.11.12 Check OK (for C coeff)
