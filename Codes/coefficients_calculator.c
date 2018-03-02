//2012.11.17

#include <stdio.h>
#include <math.h>
#include "physics_coeff.h"
#include "physics_function.h"

#include "coefficients_calculator.h"

// Calculate Boltzmann factor E[], frequency v[], background radiation Br_n[] from energy_level[], temperature T
int coeff_cal(const double *energy_level, double *v, double *E, double *F, double *Br_n, double T)
{
	int i;
	int j,jp;
	double x;

#if OUTPUT_E
	FILE *ef; //E[] debug output use
#endif
#if OUTPUT_BR_N
	FILE *brf; //Br_n[] debug output use
#endif

	// E[] Boltzmann factor ====================================
	// E = exp(-hv/kT)
	x = h_CONST/k_CONST*LIGHT_SPEED/T*-100.0;
	i = 0; //index of E[]
	j = 1; //upper level
	while(j < LEVEL_N)
	{
		jp = 0; //lower level
		while(jp < j)
		{
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


#if USE_E_LEVEL_FOR_FREQUENCY
	// v[] frequency ====================================
	// Use energy level data to calculate frequency of radiation [2012.11.17]
	x = LIGHT_SPEED*100.0/1E9; //the unit of frequency v[] is GHz, GHz = 10^9 Hz
	j = 0;
#if SHOW_A_V
	printf(" j: frequence\n");
#endif
	while(j < (LEVEL_N-1))
	{
		v[j] = x*(energy_level[j+1] - energy_level[j]);
#if SHOW_A_V
		printf("%2d: %.3e\n", j, v[j]);
#endif
		j++;
	}
#endif

	// Br_n[] normalized intensity of the cosmic background radiation ========================
	// Br_n = 1/(exp(hv/kTb)-1) 
	// from Planck's law of black-body radiation: B(T,v)= 2hv^3/c^2/(exp(hv/kT)-1)
#if USE_E_LEVEL_FOR_FREQUENCY //[2012.11.17]
	x = LIGHT_SPEED*100.0*h_CONST/k_CONST/TEMP_B;
	j = 0;
	while(j < (LEVEL_N-1))
	{
		Br_n[j] = 1/(exp(x*(energy_level[j+1] - energy_level[j]))-1); //already divided by F = 2hv^3/c^2
		j++;
	}
#else
	x = h_CONST/k_CONST/TEMP_B*1E9; //the unit of frequency v[] is GHz, GHz = 10^9 Hz
	j = 0;
	while(j < (LEVEL_N-1))
	{
		Br_n[j] = 1/(exp(x*v[j])-1); //already divided by F = 2hv^3/c^2
		j++;
	}
#endif

#if OUTPUT_BR_N
	printf("[Debug] Br_n[] output.\n");
	if( (brf = fopen( BR_N_FILE, "w" )) == NULL )//Open file for debug a_matrix output
	{
		printf( "[Debug] Can not open the file '%s' for debug output\n" , BR_N_FILE);
	}
	j = 0;
	while(j < (LEVEL_N-1))
	{
		fprintf(brf, "%.5e\n", Br_n[j]);
		j++;
	}
	if(brf)
	{
		if( fclose(brf) )
		{
			printf( "[Debug] The file '%s' was not closed\n" , BR_N_FILE);
		}
	}
#endif

#if SHOW_BR_N
	printf("j: Br_n[]  by level\n");
	x = h_CONST/k_CONST*LIGHT_SPEED/TEMP_B*100.0;
	j = 0;
	while(j < (LEVEL_N - 1))
	{
		printf("%d: %.5e %.5e\n", j, Br_n[j], 1/(exp(x*(energy_level[j+1] - energy_level[j]))-1));
		j++;
	}
	printf("\n");
#endif

	return 0;
}
