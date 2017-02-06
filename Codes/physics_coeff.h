//2012.11.20
#include <gsl/gsl_math.h>

#define OUTPUT_VER "c3.6"
#define OUTPUT_FILE_TAG "co" //any additional description of the run
#define OUTPUT_TYPE ".csv"
//output file name will be: OUTPUT_VER[LEVEL_N][Dim][NC NC][OBS_ANG PI][T K][TEMP_B K]OUTPUT_FILE_TAG.OUTPUT_TYPE
//example: c3.4[5][2D][NC 2.1827E1][0.5PI][30K][2.725K]co.csv

#define GNUPLOT_OUTPUT 1            //Creat GNUPLOT demo file or not?

//Molecule data file from LAMDA:
#define MOLE_DATA "co.win.dat"             
//#define MOLE_DATA "sio.win.dat"
//#define MOLE_DATA "co.green.win.dat"
#define USE_E_LEVEL_FOR_FREQUENCY 1   //Use energy level data to calculate frequency of radiation

#define Nt ((LEVEL_N+1)*LEVEL_N)/2                      //total number = n[0] + 2*n[1] + n[2] +..., It will not change the result.
//Number density of Collisional partner(cm^-3):
//#define NC 2.1827273E+03            //co 10,30K E3
//#define NC 2.5725E3                 //co.green 10K E3
//#define NC 3.0012500E+03            //co.green 30K E3
#define NC 225.0                     
//#define NC 5.45370E+02              //co.green 10K Goldreich C/A = 0.212
//#define NC 3.7641975E+04            //sio 10K E4
//#define NC 4.4188406E+04            //sio 30K E4
//#define NC 5.64629630E4             //sio 100K E4
//#define NC 6.62826087E4             //sio 500K E4
                                     
//Velocity Gradient Model Select
#define OneD 0						//cos^2
#define TwoD 0	                    //sin^2
#define Isotropic 0 		        //1, If none of all are chose, Isotropic will be default.
#define Mix 1                       //Mix of OneD and TwoD:
#define MixRatio 0.1	            //MixRatio*sin^2 + cos^2
#if OneD
	#define OBS_ANG 0.49167*M_PI         //Observation angle (2D default: M_PI/2, 1D default: 0.4*M_PI)
	//#define TAU_ANG OBS_ANG              //Use the line of sight as the unit of TAU(optical depth)
	#define TAU_ANG 0.0                 //Use z direction as TAU(Deguchi and Watson's unit)
	#define TAU_START 0.0001               //started TAU, the optical depth of first point in the curve
#else
	#define OBS_ANG 0.5*M_PI             //Observation angle (2D default: M_PI/2, 1D default: 0.4*M_PI)
	#define TAU_ANG OBS_ANG              //Use the line of sight as the unit of TAU(optical depth)
	#define TAU_START 0.001               //started TAU, the optical depth of first point in the curve
#endif

//Calculation range of TAU
#define TAU_N 151                    //number of points in the curve
//#define TAU_N 251                    //number of points in the curve
#define TAU_INC_RATIO 1.09647819614318 //optical depth increase ratio (increase)
//#define TAU_INC_RATIO 0.9120108393559098 //(decrease)

//Number of levels
#define LEVEL_N 8                    //total level number (ex: 2:J=0 ~ 1)
#define TOTAL_N ((LEVEL_N+1)*LEVEL_N)/2 //number of sublevels (size of n[])
#define TRANS_N ((LEVEL_N-1)*LEVEL_N)/2 //total transition number for different J (ex: size of C[])

//Temperature
#define TEMP_SELE 5                  //The collisional temperature column that was selected from the LAMDA data (start from 1)
//#define TEMP_SELE 10                 //sio 100K
//#define TEMP_SELE 31                 //sio 500K
//#define TEMP_B 50.0                 //Cosmic blackbody radiation temperature (data from Wiki 2009.12)
#define TEMP_B 2.725                 //Cosmic blackbody radiation temperature (data from Wiki 2009.12)

//================================================================================================//

//Constant (data from Wiki)
#define h_CONST 6.62606896E-34       //Plank constant
#define LIGHT_SPEED 299792458.0      //Speed of light
#define k_CONST 1.3806504E-23        //Boltzmann constant

#define REL_PREC 1E-6                //n[] precision for rate eq solving, used in rate_eq_solv.c
#define EpsRel 1e-4                  //integral relative precision, used in integral.c
#define EpsAbs 0                     //integral absolute precision, used in integral.c

#define SLOW_MODE 0                  //use the thermal equilibrium population n[] as the initial value for each main loop

#define Gsl_Integ_Space 3000         //Allocated space for GSL Integration

//debug use option-------------------------------
#define SHOW_A_V 0                  //Show the A[] and v[] that read from MOLE_DATA
#define SHOW_C 0                     //Show the C[] that read from MOLE_DATA
#define SHOW_NI 0                    //Show the initial n[] that calculated by n_initial_cal()
#define SHOW_NF 0                    //Show the finial n[] that calculated by n_initial_cal()
#define OUTPUT_A_MATRIX_I 0          //Output the a_matrix_i[] that calculated by a_matrix_initialize() 
                                     //and write it to the file A_MATRIX_I_FILE
#define A_MATRIX_I_FILE "a_matrix_i[].csv" //the output file of a_matrix_i[]
#define OUTPUT_E 0
#define E_FILE "E[].csv"
#define OUTPUT_BR_N 0
#define SHOW_BR_N 0
#define BR_N_FILE "Br_n[].csv"
//debug use option-----------------------------//

//for Visual Stdio 2005 Compiler
#define _CRT_SECURE_NO_WARNINGS 1    //Disable the SECURE WARNINGS for fopen()
