/*
2009 NTHU Astrophysics Summer Project-
Linear Polarization from Molecular Cloud in Magnetic Field with Large Velocity Gradient Approximation

First created by Louis Yang (NTHU) on 2009
Last update by Louis Yang (Kavli IPMU) on 2018.03.14

The main purpose of this code is to reproduce the result of Deguchi and Watson (1984) paper and Cortes (2005).
Ref: Cortes (2005) https://arxiv.org/abs/astro-ph/0504258 https://doi.org/10.1086/430815
*/

#define _CRT_SECURE_NO_WARNINGS    //Disable the SECURE WARNINGS for fopen()

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_integration.h>
#include "parameters.h"
#include "physics_function.h"
#include "rate_eq_solv.h"
#include "structure.h"
#include "LAMDA_Data_Reader.h"
#include "coefficients_calculator.h"
#include "tools.h"

// Global Variable --------------------
// for other files(.c .h), please include "global_value.h" to use these variable.
double A[LEVEL_N-1] = {0.0};          // Einstein coefficients for J -> J'=J-1, AJJ' = A[J']
double C[TRANS_N] = {0.0};            // Collisional excitation rates for J -> J', CJJ' = C[(J-1)J/2+J'], C[] = {C10,C20,C21,C30,C31...}
double E[TRANS_N] = {0.0};            // Boltzmann factor, exp(-dEJJ'/kT), for energy difference dE between J and J', EJJ' = E[(J-1)J/2+J']
double F[LEVEL_N-1] = {0.0};          // Flux normalization factor, 2h(vJJ')^3/c^2, FJJ' = F[J']
double v[LEVEL_N-1] = {0.0};          // Frequency (GHz) for J -> J'=J-1, vJJ' = v[J'], GHz = 10^9 Hz
double Br_n[LEVEL_N-1] = {0.0};       // Normalized cosmic blackbody radiation intensity for J -> J'=J-1, Br_JJ' = Br[J']
double energy_level[LEVEL_N] = {0.0}; // Potential energy (cm^-1) at level J (energy_level[J=0] = 0) (in unit of the inverse of wavelength(l), 1/l)
double T;                             // Temperature of the cloud (K)

double a_matrix[TOTAL_N*TOTAL_N];     // A matrix. We use 1D array to represent 2D matrix
double a_matrix_i[TOTAL_N*TOTAL_N];   // Initial a_matrix[]. We use 1D array to represent 2D matrix

gsl_integration_workspace *w;         // GSL integration workspace
gsl_integration_cquad_workspace *ws;  // GSL integration workspace for CQUAD method

unsigned long long int loop_count = 0;        //count of loops
unsigned long long int interval_count = 0;    //count for integral intervals number
// Global Variable ------------------//

void generate_output_filenames(char* file_name, char* file_name_g, int* file_name_index, double T);

int main()
{
	double n[TOTAL_N] = {0.0};
	double n_f[TOTAL_N] = {0.0};
	double tau[LEVEL_N-1][2]; // Optical Depth for each I[j][q]
	double TAU;               // Optical Depth of j=0
	double I[LEVEL_N-1][2];   // Normalized Specific Intensity for each transition radiation,
						      // I[j][q] = I(j+1 -> j)(q = 0/1 for parallel/perpendicular)
	
	int i, j, k;
	double Pt, It, Id, k0, TAUj;
	unsigned long long int loop_count_tmp, interval_count_tmp;  // Total loop and integral interval counts
	
	clock_t calcul_time[N_TAU+1];  //+++++
	time_t current_time;  //for obtaing date
    char* c_time_string;
		
	char file_name[100];
	FILE *fw;  //.csv file and .dem file
//#if GNUPLOT_OUTPUT
	char file_name_g[100];  //GNUPLOT file name
	int file_name_index;
	FILE *fwg;  //.dat file for GNUPlot use
//#endif
	
// Allocate GSL integration workspace
#if GSL_INTEGRAL_CQUAD
	ws = gsl_integration_cquad_workspace_alloc(100);
#else
	w = gsl_integration_workspace_alloc(Gsl_Integ_Space);
#endif


// =========================== Initialization =============================

	printf("Molecular data:      %s\n", MOLE_DATA);
	printf("N of levels:         %d\n", LEVEL_N);

	// Read Coefficients Data _____________________________________________
	if(lamda_data_reader(A, v, C, energy_level, &T)) //read A C v from file 'co.win.dat'
	{
		printf("Cannot read coefficients data.\n");
		pause();
		return 0;
	}
	printf("Finished reading data.\n");

	printf("Gas temperature:     %.0f K\n", T);
	printf("NC:                  %.3e cm^-3\n", NC);
	printf("C[0]:                %.3e s^-1\n", C[0]);
	printf("A[0]:                %.3e s^-1\n", A[0]);
	printf("C/A:                 %lf\n", C[0] / A[0]);
	printf("Observational angle: %.3f PI\n", OBS_ANG / M_PI);
	printf("TAU angle:           %.3f PI\n", TAU_ANG / M_PI);
	printf("Case:                ");
#if TwoD
	printf("2D");
#elif OneD
	printf("1D");
#elif Mix
	printf("Mix with ratio %g", MixRatio);
#elif Isotropic
	printf("Isotropic");
#endif
	printf("\n");
	printf("TAU range:           %lf - %lf [%d steps]\n",
		TAU_START, TAU_START * pow(TAU_INC_RATIO, N_TAU-1), N_TAU);

	// Generate Output file names __________________________________________
	generate_output_filenames(file_name, file_name_g, &file_name_index, T);

	// Open output files ______________________________________________
	if ((fw = fopen(file_name, "w")) == NULL) {
		printf("Cannot open the file %s\n", file_name);
		pause();
		return 0;
	}
#if GNUPLOT_OUTPUT
	//open GNUPLOT output file
	if ((fwg = fopen(file_name_g, "w")) == NULL) {
		printf("Cannot open the file %s\n", file_name_g);
		pause();
		return 0;
	}
#endif	
	
	// Calculate Coefficients ______________________________________________
	coeff_cal(energy_level, v, E, F, Br_n, T);		

	// Initialize a_matrix_i[] _____________________________________________
	a_matrix_initialize(a_matrix_i); //fill a_matrix_i[] with C[]
#if OUTPUT_A_MATRIX_I
	output_a_matrix(a_matrix_i, A_MATRIX_I_FILE);
#endif
	
	printf("\n");

	// Initialize population n[] ___________________________________________
	n_initial_cal(n, TEMP_B);
	//n_initial_cal(n, T);

#if SHOW_NI
	printf("Initial n[]:\n");
	print_nv(n);
	printf("\n");
#endif

#if SHOW_NF
	printf("\nn[] at %fK:\n", T);
	n_initial_cal(n_f, T);
	print_nv(n_f);
#endif

#if SHOW_E
	print_E(E);
	printf("\n");
#endif

#if SHOW_I_LIMIT
	print_I_limit(E, Br_n, v);
	printf("\n");
#endif
	
#if SHOW_S_INIT
	printf("Initial source functions:\n");
	TAU = TAU_START; //set TAU to started tau
	tau_array(TAU, tau, n); // Calculate tau[] for each levels along the line of sight
	print_source_f(n, Br_n, tau);
	printf("\n");
#endif	

#if SHOW_I_INIT
	// Check I_excess from initial n[]
	printf("Initial intensity:\n");
	TAU = TAU_START; //set TAU to started tau
	tau_array(TAU, tau, n); // Calculate tau[] for each levels along the line of sight
	I_emerge_n(n, tau, I);  // Calculate the emitted Intensities
	k0 = k_f_n(n, cos(TAU_ANG), 0, 0);
	printf("k0: %.3e\n", k0);
	printf("j: TAUj, Pt, I[j][0], I[j][1], Id, v[j]\n");
	for (j = 0; j < (LEVEL_N - 1); j++) {
		TAUj = TAU * (k_f_n(n, cos(TAU_ANG), 0, j) / k0);
		Pt = (I[j][1] - I[j][0]) / (I[j][0] + I[j][1]) * 100;  // Fractional polarization
		It = (I[j][0] + I[j][1])*h_CONST*v[j] / k_CONST * 1E9; // Total intensity
		Id = (I[j][1] - I[j][0])*h_CONST*v[j] / k_CONST * 1E9; // Intensity difference

		printf("%d: %.3e, %3f%%, %.3e, %.3e, %.3e, %.3e\n",
			j, TAUj, Pt, I[j][0] * h_CONST*v[j] / k_CONST * 1E9, I[j][1] * h_CONST*v[j] / k_CONST * 1E9, Id, v[j]);
	}
	printf("\n");
#endif
	
#if TEST_A_INIT && LEVEL_N == 3
	test_a_matrix_initialize_3();  // Pass on 2018.03.08
#endif

#if TEST_A_INIT_SOLVE
	test_solve_a_matrix_i();  // Pass on 2018.03.11
#endif

#if TEST_RATE_EQ_FILL && LEVEL_N == 3
	test_rate_eq_fill_3();  // Pass on 2018.03.12
#endif

#if TEST_S_ISO
	test_source_f_n_iso();  // Pass on 2018.03.14
#endif

#if TEST_S_3 && LEVEL_N == 3
	test_source_f_n_3();  // Pass on 2018.03.14
#endif
	printf("========== Initialization done. ==========\n\n");//*****

// ========================= Initialization done =======================//

	
	// Write column name to the output files
	for (i = 0; i < (LEVEL_N - 1); i++) {
		fprintf(fw, "TAU[%d], P%d(%%), I[%d][0](K), I[%d][1](K), Id[%d](K),", i, i, i, i, i);
	}	
	for (i = 0; i < TOTAL_N; i++) {
		fprintf(fw, "n[%d],", i);
	}
	fprintf(fw, "cal_time (%.0es), loop_count, interval_count, 0 = parallel; 1 = perpendicular\n", 1./CLOCKS_PER_SEC);

	
// ==============================================================================
// ========================= Main Calculation ===================================

	loop_count_tmp = loop_count;
	interval_count_tmp = interval_count;
	calcul_time[0] = clock();
	

	// Main Loop ----------------------------------------------------------------
	TAU = TAU_START; //set TAU to started tau
	for (i = 0; i < N_TAU; i++) //N_TAU = number of points in the curve
	{
		printf("%d: %.2e ", i, TAU);

#if SLOW_MODE
		n_initial_cal(n,T); // Use thermal equilibrium to initialize n[] for each TAU
#endif
		
		// Calculate n[] 
		rate_eq_solve(n, TAU);  // Main Calculation

#if SHOW_N
		print_n(n);
#endif
		
		// Calculate I[][] 
		tau_array(TAU, tau, n); // Build tau[] for each levels along the line of sight
		I_emerge_n(n, tau, I);  // Calculate the emitted intensities
		calcul_time[i + 1] = clock(); //+++++

#if SHOW_S
		if (i == 0) print_source_f(n, Br_n, tau);
#endif

		// Write results to the files ----------------------
		// Output TAU, I[], and polarization Pt
		k0 = k_f_n(n, cos(TAU_ANG), 0, 0);
		for(j = 0; j < (LEVEL_N-1); j++) {
			TAUj = TAU * (k_f_n(n, cos(TAU_ANG), 0, j) / k0); // Effective TAU for j
#if 1
			Pt = (I[j][1] - I[j][0])/(I[j][0] + I[j][1])*100; //Fractional polarization
			It = (I[j][0]+I[j][1])*h_CONST*v[j]/k_CONST*1E9; //Total intensity
			Id = (I[j][1]-I[j][0])*h_CONST*v[j]/k_CONST*1E9; //Intensity difference
#else
			Pt = (fabs(I[j][1]) - fabs(I[j][0]))/(fabs(I[j][0]) + fabs(I[j][1]))*100; //Fractional polarization
			It = (fabs(I[j][0])+fabs(I[j][1]))*h_CONST*v[j]/k_CONST*1E9; //Total intensity
			Id = (fabs(I[j][1])-fabs(I[j][0]))*h_CONST*v[j]/k_CONST*1E9; //Intensity difference
#endif

			fprintf(fw, "%.10e, %5f%%, %.10e, %.10e, %.10e,",
				TAUj, Pt, I[j][0]*h_CONST*v[j]/k_CONST*1E9, I[j][1]*h_CONST*v[j]/k_CONST*1E9, Id); //write results to .csv file
#if GNUPLOT_OUTPUT
			fprintf(fwg, "%.10e %5f%% %.10e %.10e %.10e ",
				TAUj, Pt, I[j][0]*h_CONST*v[j]/k_CONST*1E9, It, Id); //write results to .dat file
#endif	
		}

		// Output population n[]
		for (j = 0; j < TOTAL_N; j++) {
			fprintf(fw, "%.10e,", n[j]);
#if GNUPLOT_OUTPUT
			fprintf(fwg, "%.10e ", n[j]);
#endif	
		}
		
		// Output CPU time 
		int d_calcul_time = (int)(calcul_time[i + 1] - calcul_time[i]);
		unsigned long long d_loop_count = loop_count - loop_count_tmp;
		unsigned long long d_interval_count = interval_count - interval_count_tmp;
		fprintf(fw,"%d, %llu, %llu\n", d_calcul_time, d_loop_count, d_interval_count);//+++++
#if GNUPLOT_OUTPUT
		fprintf(fwg,"%d %llu %llu\n", d_calcul_time, d_loop_count, d_interval_count);//+++++
#endif
		
		// Show current progress
		printf("%lfs %lluloops i:%llu %.2lf\n", (double)calcul_time[i+1]/CLOCKS_PER_SEC, d_loop_count, d_interval_count,
			d_interval_count /(double)(d_loop_count));//+++++
		//----------------------------------------------//
		
		loop_count_tmp = loop_count;
		interval_count_tmp = interval_count;
		
		TAU *= TAU_INC_RATIO; //TAU increase	
		printf("\n");
	}
	// Main Loop--------------------------------------------------------------//

	printf("Total loops: %llu    Total integral intervals: %llu\n", loop_count, interval_count);

#if GSL_INTEGRAL_CQUAD
	gsl_integration_cquad_workspace_free(ws);
#else
	gsl_integration_workspace_free (w);  //release GSL integration workspace
#endif

// ========================= Main Calculation done ============================//
// ============================================================================//



	// Write calculation information to the file-----------------------
	fprintf(fw, "N of levels,%d,,", LEVEL_N);
#if TwoD
	fprintf(fw, "2D case,,");
#elif OneD
	fprintf(fw, "1D case,,");
#elif Mix
	fprintf(fw, "Mix case,%g,,", MixRatio);
#elif Isotropic
	fprintf(fw, "Isotropic case,,");
#endif
	fprintf(fw, "MOLE_DATA,%s\n", MOLE_DATA);
	fprintf(fw, "T_Gas,%gK,,T_Background,%gK\n", T, TEMP_B);
	fprintf(fw, "NC,%.10e,,number density of Collisional partner(cm^-3)\n", NC);
	fprintf(fw, "C[0],%.5e,,A[0],%.5e,,C/A,%lf\n", C[0], A[0], C[0]/A[0]);
	fprintf(fw, "OBS_ANG(pi),%lf,,TAU_ANG(pi),%lf,,TAU_Start,%lf,,", OBS_ANG/M_PI, TAU_ANG/M_PI, TAU_START);
	fprintf(fw, "OBS_ANG = observational angle; TAU_ANG = angle for unit of TAU\n");
	fprintf(fw, "\n");

	fprintf(fw, "I_limit[](K)"); //[2012.11.17]Included the contribution from background radiation
	j = 0;
	while(j < (LEVEL_N-1)) {
		fprintf(fw,",%.10e", (1/(1/E[indexJJ(j+1,j)]-1)-Br_n[j])*h_CONST*v[j]/k_CONST*1E9);
		j++;
	}
	fprintf(fw, ",,I_t in the limit of large optical depth\n");
	
	fprintf(fw, "I_last[](K)");
	j = 0;
	while(j < (LEVEL_N-1)) {
		fprintf(fw,",%.10e", (I[j][0]+I[j][1])*h_CONST*v[j]/k_CONST*1E9);
		j++;
	}
	fprintf(fw, ",,Last total excess intensity I_t\n");

	fprintf(fw, "E[]");
	j = 0;
	while(j < (LEVEL_N-1)) {
		fprintf(fw,",%.10e", E[indexJJ(j+1,0)]);
		j++;
	}
	fprintf(fw, ",,Ratio of n[]/n[0] in the limit of large optical depth\n");

	fprintf(fw, "n[]/n[0]");
	j = 1;
	while(j < LEVEL_N) {
		fprintf(fw,",%.10e", n[indexN(j,0)]/n[0]);//ratio of n[]
		j++;
	}
	fprintf(fw, ",,Last n[]/n[0] ratio\n");

	fprintf(fw, "I+B_limit_q_n[]");
	j = 0;
	while(j < (LEVEL_N-1)) {
		fprintf(fw,",%.10e", 1/(1/E[indexJJ(j+1,j)]-1)/2);
		j++;
	}
	fprintf(fw, "\n");

	fprintf(fw, "Br_n[]/2");
	i = 0;
	while(i < (LEVEL_N-1)) {
		fprintf(fw,",%.10e", Br_n[i]/2);
		i++;
	}
	fprintf(fw, ",,Normalized background radiation intensity in one polarization\n");

	fprintf(fw, "S_f_last[0]");
	j = 0;
	while(j < (LEVEL_N-1)) {
		fprintf(fw, ",%.10e", source_f_n(n, cos(OBS_ANG), 0,j+1));
		j++;
	}
	fprintf(fw, ",,Last normailzed source function in 0\n");

	fprintf(fw, "S_f_last[1]");
	j = 0;
	while(j < (LEVEL_N-1)) {
		fprintf(fw, ",%.10e", source_f_n(n, cos(OBS_ANG), 1,j+1));
		j++;
	}
	fprintf(fw, ",,Last normailzed source function in 1\n");
	fprintf(fw, "\n");

	fprintf(fw, "REL_PREC,%.10e,,EpsRel,%.10e,,EpsAbs,%.10e\n", REL_PREC, EpsRel, EpsAbs);
	fprintf(fw, "Gsl_Integ_Space,%d,,Slow Mode,%d\n", Gsl_Integ_Space, SLOW_MODE);
	fprintf(fw, "Time spend(ms),%d\n", (int)calcul_time[N_TAU]);//+++++
	fprintf(fw, "Loops,%llu,,Intervals,%llu\n", loop_count, interval_count);
	current_time = time(NULL);
	c_time_string = ctime(&current_time);
	fprintf(fw, "Version,%s,,Date,%s\n",OUTPUT_VER,c_time_string);
	// write calculation information to the file---------------------//
	
	
	// Close result data files --------------
	if (fw) {
		printf("Result is saved to\n '%s'\n", file_name);
		if (fclose(fw)) {
			printf("Can not close '%s'!\n", file_name);
		}
	}

#if GNUPLOT_OUTPUT
	if (fwg) {
		printf("and\n '%s'\n", file_name_g);
		if (fclose(fwg)) {
			printf("Can not close '%s'!\n", file_name_g);
		}
	}
	

	// Generate GNUPLOT .dem file -------------
	sprintf(file_name + file_name_index, ".dem");
	if( (fwg = fopen( file_name, "w" )) == NULL )//Open file for write in
	{
		printf( "Can not open the file '%s'\n" , file_name);
		pause();
		return 0;
	}
	fprintf(fwg, "set logscale x\n");
	fprintf(fwg, "set title \"Specific Intensity I(K) (T = %gK, C/A = %g)\"\n", T, C[0]/A[0]);
	fprintf(fwg, "plot '%s' using 1:4 title \"1->0\" with lines", file_name_g);
	i = 1;
	while(i < (LEVEL_N-1)) {
		fprintf(fwg, ",\\\n '%s' using %d:%d title \"%d -> %d\" with lines", file_name_g, i*5+1,i*5+4, i+1, i);
		i++;
	}
	fprintf(fwg, "\npause -1\n");
	fprintf(fwg, "set title \"Fractional Linear Polarization P(%%) (T = %gK, C/A = %g)\"\n", T, C[0]/A[0]);
	fprintf(fwg, "plot '%s' using 1:2 title \"P[1 -> 0]\" with lines", file_name_g);
	i = 1;
	while(i < (LEVEL_N-1)) {
		fprintf(fwg, ",\\\n '%s' using %d:%d title \"P[%d -> %d]\" with lines", file_name_g, i*5+1,i*5+2, i+1, i);
		i++;
	}
	fprintf(fwg, "\npause -1\n");
	fprintf(fwg, "reset\n");
	
	if(fwg)	{
		printf("GNUPLOT demo file is saved to\n '%s'\n", file_name);
		if ( fclose(fwg) ) {
			printf( "Can not close '%s'\n" , file_name);
        }
    }
#endif
	// Generate GNUPLOT .dem file -------//
	
	printf("Finished calculation!\n");
	pause();
	return 0;
}

void generate_output_filenames(char* file_name, char* file_name_g, int* file_name_index, double T) {
	int i = 0;
	i += sprintf(file_name + i, OUTPUT_VER); //ex: c3.2
	i += sprintf(file_name + i, "[%d]", LEVEL_N);
#if TwoD
	i += sprintf(file_name + i, "[2D]");
#elif OneD
	i += sprintf(file_name + i, "[1D]");
#elif Mix
	i += sprintf(file_name + i, "[Mix%g]", MixRatio);
#else //Isotropic
	i += sprintf(file_name + i, "[Iso]");
#endif
	i += sprintf(file_name + i, "[%.2e][%gPI][%gK][%gK]", NC, OBS_ANG / M_PI, T, TEMP_B);
	i += sprintf(file_name + i, OUTPUT_FILE_TAG); //ex: test2
												  //#if GNUPLOT_OUTPUT
	(*file_name_index) = i;
	sprintf(file_name_g, "%s.dat", file_name);
	i += sprintf(file_name + i, OUTPUT_TYPE); //ex: .csv
}

