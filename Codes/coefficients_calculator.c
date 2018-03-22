//2012.11.17

#include <stdio.h>
#include <math.h>
#include "parameters.h"
#include "physics_function.h"

#include "coefficients_calculator.h"

// Calculate Boltzmann factor E[], frequency v[], background radiation Br_n[] from energy_level[], temperature T
// energy_level[LEVEL_N]: Potential energy (cm^-1) at level J (energy_level[J=0] = 0) (in unit of the inverse of wavelength(l), 1/l)
// v[LEVEL_N - 1]       : Frequency (GHz) for J -> J'=J-1, vJJ' = v[J'], GHz = 10^9 Hz
// E[TRANS_N]           : Boltzmann factor, exp(-dEJJ'/kT), for energy difference between J and J', dEJJ'. EJJ' = E[(J-1)J/2+J']
// F[LEVEL_N-1]         : Flux normalization factor, 2h(vJJ')^3/c^2, FJJ' = F[J']
// Br_n[LEVEL_N-1]      : Normalized cosmic blackbody radiation intensity for J -> J'=J-1, Br_JJ' = Br[J']
// S_ext_n[LEVEL_N-1]   : Normalized intensity from external source for J -> J'=J-1, S_ext_JJ' / FJJ' = S_ext_n[J']
//                        S_ext_n = (1 - exp(-TAU_ext)) / (exp(hv/kT_ext) - 1)
// T                    : Temperature of the cloud (K)
int coeff_cal(const double *energy_level, double *v, double *E, double *F, double *Br_n, double *S_ext_n, double T) {
	int i;
	int j, jp;
	double x;

#if OUTPUT_E
	FILE *ef; //E[] debug output use
#endif
#if OUTPUT_BR_N
	FILE *brf; //Br_n[] debug output use
#endif

	// E[] Boltzmann factor ====================================
	// E = exp(-hv/kT)
	x = h_CONST / k_CONST * LIGHT_SPEED / T * -100.0;  
	// Note the factor of 100 comes from that energy_level[] are in cm^-1 but LIGHT_SPEED is in m/s
	i = 0; //index of E[]
	j = 1; //upper level
	while(j < LEVEL_N) {
		jp = 0; //lower level
		while(jp < j) {
			E[i] = exp(x*(energy_level[j] - energy_level[jp]));
			i++;
			jp++;
		}
		j++;
	}
	//============================================================//
	
#if OUTPUT_E
	printf("[Debug] E[] output.\n");
	if( (ef = fopen( E_FILE, "w" )) == NULL )//Open file for debug a_matrix output
	{
		printf( "[Debug] Can not open the file '%s' for debug output\n" , E_FILE);
	}
	j = 1;
	while(j < LEVEL_N)
	{
		jp = 0; //lower level
		while(jp < j)
		{
			fprintf(ef, "%.5e,", E[indexJJ(j,jp)]);
			jp++;
		}
		fprintf(ef, "\n");
		j++;
	}
	if(ef)
	{
		if( fclose(ef) )
		{
			printf( "[Debug] The file '%s' was not closed\n" , E_FILE);
		}
	}
#endif


	// v[] frequency ====================================
#if USE_E_LEVEL_FOR_FREQUENCY
	// Use energy level data to calculate frequency of radiation [2012.11.17]
	x = LIGHT_SPEED * 100.0 / 1E9; //the unit of frequency v[] is GHz, GHz = 10^9 Hz
	// Note the factor of 100 comes from that energy_level[] are in cm^-1 but LIGHT_SPEED is in m/s
	for (j = 0; j < (LEVEL_N - 1); j++) {
		v[j] = x*(energy_level[j+1] - energy_level[j]);
	}
#endif

	// Br_n[] normalized intensity of the cosmic background radiation ========================
	// Br_n = 1 / (exp(hv/kTb) - 1) 
	// from Planck's law of black-body radiation: B(T,v) = 2hv^3/c^2 / (exp(hv/kT) - 1)
#if USE_E_LEVEL_FOR_FREQUENCY //[2012.11.17]
	x = LIGHT_SPEED * 100.0 * h_CONST / k_CONST / TEMP_B;
	// Note the factor of 100 comes from that energy_level[] are in cm^-1 but LIGHT_SPEED is in m/s
	for (j = 0; j < (LEVEL_N - 1); j++) {
		Br_n[j] = 1 / (exp(x*(energy_level[j + 1] - energy_level[j])) - 1); //already divided by F = 2hv^3/c^2
	}
#else
	x = h_CONST/k_CONST/TEMP_B*1E9; //the unit of frequency v[] is GHz, GHz = 10^9 Hz
	for (j = 0; j < (LEVEL_N - 1); j++) {
		Br_n[j] = 1 / (exp(x*v[j]) - 1); //already divided by F = 2hv^3/c^2
	}
#endif

#if OUTPUT_BR_N
	printf("[Debug] Br_n[] output.\n");
	if( (brf = fopen( BR_N_FILE, "w" )) == NULL ) {
		printf( "[Debug] Can not open the file '%s' for debug output\n" , BR_N_FILE);
	}
	j = 0;
	while(j < (LEVEL_N-1)) {
		fprintf(brf, "%.5e\n", Br_n[j]);
		j++;
	}
	if(brf) {
		if( fclose(brf) ) {
			printf( "[Debug] The file '%s' was not closed\n" , BR_N_FILE);
		}
	}
#endif

#if SHOW_BR_N
	printf("j: Br_n[]  by level\n");
	x = h_CONST/k_CONST*LIGHT_SPEED/TEMP_B*100.0;
	j = 0;
	while(j < (LEVEL_N - 1)) {
		printf("%d: %.5e %.5e\n", j, Br_n[j], 1/(exp(x*(energy_level[j+1] - energy_level[j]))-1));
		j++;
	}
	printf("\n");
#endif

	// S_ext_n[] normalized radiation from external source ========================
	// S_ext_n = (1 - exp(-TAU_ext)) / (exp(hv/kT_ext) - 1) 
	// from Planck's law of black-body radiation: B(T,v) = 2hv^3/c^2 / (exp(hv/kT) - 1)
#if EXT_SOURCE
#if USE_E_LEVEL_FOR_FREQUENCY
	x = LIGHT_SPEED * 100.0 * h_CONST / k_CONST / TEMP_EXT;
	for (j = 0; j < (LEVEL_N - 1); j++) {
		S_ext_n[j] = (1 - exp(-TAU_EXT)) / (exp(x*(energy_level[j + 1] - energy_level[j])) - 1); //already divided by F = 2hv^3/c^2
	}
#else
	x = h_CONST / k_CONST / TEMP_EXT * 1E9; //the unit of frequency v[] is GHz, GHz = 10^9 Hz
	for (j = 0; j < (LEVEL_N - 1); j++) {
		S_ext_n[j] = (1 - exp(-TAU_EXT)) / (exp(x*v[j]) - 1); //already divided by F = 2hv^3/c^2
	}
#endif
#endif

	printf("Finished coefficient calculation.\n");
	return 0;
}

