/*
LAMDA Data Reader
Read the Einstein coefficients A, Frequency between each level (v), Collisional coefficients C,
Energy at each level (energy_level), and temperature T
from LAMDA datafiles

References:
Schoier, F.L., van der Tak, F.F.S., van Dishoeck E.F., Black, J.H. 2005, A&A 432, 369-379
http://doi.org/10.1051/0004-6361:20041729
http://www.strw.leidenuniv.nl/~moldata/
http://www.strw.leidenuniv.nl/~moldata/molformat.html

Louis Yang (Kavli IPMU)
2018.03.08
*/
#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include "parameters.h"
#include "LAMDA_Data_Reader.h"

int lamda_data_reader(double *A, double *v, double *C, double *energy_level, double *T)
{
	//'@@@@@...' notes the line of the code for actually storing the coefficients or parameters.

	FILE *file1;
	char file_name_read[] = MOLE_DATA;
	//double A[LEVEL_N-1] = {0}; //Einstein coefficients, AJJ' = A[J'] (J-1 = J')
	//double C[((LEVEL_N-1)*LEVEL_N)/2] = {0}; //CJ->J' = C[(J-1)J/2+J'], C[] = {C10,C20,C21,C30,C31...}
	//double v[LEVEL_N-1] = {0}; //frequency (GHz) from J to J'(vJJ'), vJJ' = v[J'] (J-1 = J')
	//double energy_level[LEVEL_N] = {0}; //energy of each level (cm^-1)(unit of wavelength(l) inverse)(1/l)
	
	int nlev, nlin, ncol, ntemp;
	int i,j;
	int use_nlev, use_nlin, use_ncol, use_temp;
	
	
	//------------------- Open files -------------------------
	if( (file1 = fopen( file_name_read, "r" )) == NULL )//Open file for CO
	{
		printf( "Cannot open file '%s'!\n", file_name_read);
		return 1;
	}

	//========================== Header ============================================
	//C language note: "%*X" means that read the data in X type but don't store it.
	//So "%*[^\n]\n" will jump to next line of the reading file without store this line.
	fscanf(file1, "%*[^\n]\n"); //!MOLECULE
	fscanf(file1, "%*[^\n]\n"); //CO
	fscanf(file1, "%*[^\n]\n"); //!MOLECULAR WEIGHT
	fscanf(file1, "%*[^\n]\n"); //28.0
	
	//======================= Energy Level =========================================
	//get energy_level[]
	fscanf(file1, "%*[^\n]\n"); //!NUMBER OF ENERGY LEVELS
	fscanf(file1, "%d\n", &nlev); //NLEV, (NUMBER OF ENERGY LEVELS)
	use_nlev = LEVEL_N;
	if(use_nlev > nlev)//Check use_nlev
	{
		printf("No Energy level data for the level higher than %d.\n", nlev);
		use_nlin = nlev;
	}

	fscanf(file1, "%*[^\n]\n"); //!LEVEL + ENERGIES(cm^-1) + WEIGHT + J
	
	i = 0;
	while(i < use_nlev)
	{
		fscanf(file1, " %*d %lf %*[^\n]\n", &energy_level[i]); //read, energy of each level (cm^-1)(unit of wavelength(l) inverse)(1/l)
		i++;                                                   //@@@@@@@@@@@@@@@@@@@@@
	}
	while(i < nlev)
	{
		fscanf(file1, "%*[^\n]\n"); //no use
		i++;
	}

	//======================= Radiative Transitions ================================
	//get A[] and v[]
	fscanf(file1, "%*[^\n]\n"); //!NUMBER OF RADIATIVE TRANSITIONS
	fscanf(file1, "%d\n", &nlin); //NLIN, (NUMBER OF RADIATIVE TRANSITIONS) @@@@@@@@@@@@@@@@@@@@@
	fscanf(file1, "%*[^\n]\n"); //!TRANS + UP + LOW + EINSTEINA(s^-1) + FREQ(GHz) + E_u(K)

	use_nlin = use_nlev-1;
	if(use_nlin > nlin)
	{
		printf("No Einstein coefficient data for the level higher than %d.\n", nlin+1);
		use_nlin = nlin;
	}
	i = 0;
	while(i < use_nlin)
	{
#if USE_E_LEVEL_FOR_FREQUENCY
		fscanf(file1, " %*d %*d %*d %le %*[^\n]\n", &A[i]); //[2012.11.17]read A coefficients(s^-1) only
#else
		fscanf(file1, " %*d %*d %*d %le %lf %*[^\n]\n", &A[i], &v[i]); //read A coefficients(s^-1) and Frequency(GHz), GHz = 10^9 Hz
#endif
		i++;                                                           //@@@@@@@@@@@@@@@@@@@@@@@@
	}
	while(i < nlin)
	{
		fscanf(file1, "%*[^\n]\n"); //no use
		i++;
	}

#if SHOW_A_V
#if USE_E_LEVEL_FOR_FREQUENCY
	printf(" j  A[]:\n");
#else
	printf(" j  A[] v[]:\n");
#endif
	i=0;
	while(i < (LEVEL_N-1))
	{
#if USE_E_LEVEL_FOR_FREQUENCY
		printf("%2d: %.3e\n", i, A[i]);
#else
		printf("%2d: %.3e %.3e\n", i, A[i], v[i]);
#endif
		i++;
	}
#endif
	
	//======================= Collisional Coefficients =============================
	//get C[]
	fscanf(file1, "%*[^\n]\n"); //!NUMBER OF COLL PARTNERS
	fscanf(file1, "%*[^\n]\n"); //2
	fscanf(file1, "%*[^\n]\n"); //!COLLISIONS BETWEEN
	fscanf(file1, "%*[^\n]\n"); //2 CO-pH2 from Flower (2001) & Wernli et al. (2006) + extrapolation
	fscanf(file1, "%*[^\n]\n"); //!NUMBER OF COLL TRANS
	fscanf(file1, "%d\n", &ncol);//number of transitions for which collisional data exist (NCOL) 
	if(ncol != (((nlin+1)*nlin)/2)) //Check ncol
	{
		printf("NCOL is not correct. ncol = %d (nlin+1)*nlin/2 = %d\n", ncol,(((nlin+1)*nlin)/2));
	}
	use_ncol = ((LEVEL_N-1)*LEVEL_N)/2;
	if(use_ncol > ncol)//Check use_ncol
	{
		printf("No coll coefficient data for trans higher than %d.\n", ncol);
		use_ncol = ncol;
	}
	
	fscanf(file1, "%*[^\n]\n"); //!NUMBER OF COLL TEMPS
	fscanf(file1, "%d\n", &ntemp);//number of temperatures for which collisional data exist (ntemp)
		
	use_temp = TEMP_SELE;
	if(use_temp > ntemp)//Check use_temp
	{
		printf("No such temperature option (%d).\n", use_temp);
		use_temp = ntemp;
	}
	
	fscanf(file1, "%*[^\n]\n"); //!COLL TEMPS
    j = 1;
	while(j < use_temp)
	{
		fscanf(file1, "%*e"); //neglect the columns before
		j++;
	}
	fscanf(file1, "%le", T); //read temperature T @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	if(use_temp != ntemp)
	{
		fscanf(file1, "%*[^\n]\n"); //neglect the columns after
	}
	else
	{
		fscanf(file1, "\n"); //neglect the columns after
	}
	fscanf(file1, "%*[^\n]\n"); //!TRANS + UP + LOW + COLLRATES(cm^3 s^-1)

	i = 0;
	while(i < use_ncol)
	{
		fscanf(file1, "%*d %*d %*d");
		j = 1;
		while(j < use_temp)
		{
			fscanf(file1, "%*e"); //neglect the columns before
			j++;
		}
		fscanf(file1, "%le", &C[i]); //read Collisional coefficients C[] @@@@@@@@@@@@@@@@@@@@@@
		C[i] *= NC;                  //Multiply by density of the collisional partner
		if(use_temp != ntemp)
		{
			fscanf(file1, "%*[^\n]\n"); //neglect the columns after
		}
		else
		{
			fscanf(file1, "\n");
		}
		i++;
	}

	//------------------- Close files ------------------------
	if(file1)//Close file1
	{
		if ( fclose(file1) )
        {
			printf( "Cannot close the file '%s'.\n", file_name_read);
        }
    }

#if SHOW_C
	printf("C[]:\n");
	i = 0;
	while(i < use_ncol)
	{
		printf("%.2e\n", C[i]);
		i++;
	}
#endif

	return 0;
}
