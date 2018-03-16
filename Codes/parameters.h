//2018.03.15
#include <gsl/gsl_math.h>

#define OUTPUT_VER "c3.7"
#define OUTPUT_FILE_TAG "co2009"        // Any additional description of the run
//#define OUTPUT_FILE_TAG "sio"     // Any additional description of the run
#define OUTPUT_TYPE ".csv"
// Output file name will be: OUTPUT_VER[LEVEL_N][Dim][NC NC][OBS_ANG PI][T K][TEMP_B K]OUTPUT_FILE_TAG.OUTPUT_TYPE
// Example: c3.4[5][2D][NC 2.1827E1][0.5PI][30K][2.725K]co.csv

// Molecule data file from LAMDA:
//#define MOLE_DATA "co.win.dat"             
#define MOLE_DATA "co2009.win.dat"
//#define MOLE_DATA "sio.win.dat"
//#define MOLE_DATA "co.green.win.dat"

#define USE_E_LEVEL_FOR_FREQUENCY 1 // Use energy level data to calculate frequency of radiation
#define GNUPLOT_OUTPUT 1            //Creat GNUPLOT demo file or not?
#define WINDOWS                     //Pause at the end of the program

// Number of levels
#define LEVEL_N 2                    // Total number of levels (Ex: 2: J = 0 ~ 1)
#define TOTAL_N ((LEVEL_N+1)*LEVEL_N)/2 // Total number of independent sublevels (size of n[])
#define TRANS_N ((LEVEL_N-1)*LEVEL_N)/2 //total transition number for different J (size of C[])

#define PC_FACTOR 1.              // Particle number conservation factor

//#define Nt PC_FACTOR*((LEVEL_N+1)*LEVEL_N)/2    // = n[0] + 2*n[1] + n[2] +..., It will not change the result.
#define Nt PC_FACTOR*LEVEL_N*LEVEL_N            //total number = n[0] + n[1] + 2*n[2] +..., It will not change the result.
//Number density of collisional partner (cm^-3):
//#define NC 2.1827273E+03            //co 10,30K E3
//#define NC 2.5725E3                 //co.green 10K E3
//#define NC 3.0012500E+03            //co.green 30K E3
//#define NC 225.0                     
//#define NC 5.45370E+02              //co.green 10K Goldreich C/A = 0.212
//#define NC 3.7641975E+04            //sio 10K E4
//#define NC 4.4188406E+04            //sio 30K E4
//#define NC 5.64629630E6             //sio 100K E4
//#define NC 6.62826087E4             //sio 500K E4
#define NC 3E3 // 3E3
                                     
// Velocity Gradient Model Select
#define OneD 0						//cos^2
#define TwoD 1	                    //sin^2
#define Isotropic 0 		        //1, If none of all are chose, Isotropic will be default.
#define Mix 0                       //Mix of OneD and TwoD:
#define MixRatio 0.1	            //MixRatio*sin^2 + cos^2
// OBS_ANG:   Angle between the observational line of sight and the z-axis (direction of B field)
//            (2D default: M_PI/2, 1D default: 0.4*M_PI)
// TAU_ANG:   Angle between the direction for the unit of TAU (optical depth) and the z-axis (direction of B field)
// TAU_START: Initial optical depth. TAU of the first point in the curve.
#if OneD
	#define OBS_ANG 0.49167*M_PI         // Observational angle (2D default: M_PI/2, 1D default: 0.4*M_PI)
	//#define TAU_ANG OBS_ANG            // Use the line of sight as the unit of TAU
	#define TAU_ANG 0.0                  // Use z direction for the unit of TAU (Deguchi and Watson's unit)
	#define TAU_START 0.001
#else
	#define OBS_ANG 0.5*M_PI             // Observational angle (2D default: M_PI/2, 1D default: 0.4*M_PI)
	#define TAU_ANG OBS_ANG              // Use the line of sight for the unit of TAU
	#define TAU_START 0.001
	//#define TAU_START 1000.0
#endif

// Calculation range of TAU
#define N_TAU 151                    //number of points in the curve
//#define N_TAU 251                    //number of points in the curve
#define TAU_INC_RATIO 1.09647819614318 //optical depth increase ratio bewteen each point (increase)
//#define TAU_INC_RATIO 0.9120108393559098 //(decrease)

// Temperature
#define TEMP_SELE 11                 //The collisional temperature column that was selected from the LAMDA data (start from 1) // Before 11
//#define TEMP_SELE 10                 //sio 100K
//#define TEMP_SELE 31                 //sio 500K
//#define TEMP_B 50.0                  //Cosmic blackbody radiation temperature
#define TEMP_B 2.725                 //Cosmic background radiation temperature (K, data from Wiki 2009.12)

// Integral methods
#define GSL_INTEGRAL_QNG 0           // QNG non-adaptive Gauss-Kronrod integration
#define GSL_INTEGRAL_CQUAD 0         // CQUAD doubly-adaptive integration method
// If not using both of them, then in default it uses QAGS adaptive integration with singularities method
#define Gsl_Integ_Space 12000        // Allocated space for GSL Integration if not using CQUAD
#define Gsl_Integ_CQUAD_Space 100    // Allocated space for GSL Integration with CQUAD method

//================================================================================================//

#define REL_PREC 1E-6                //n[] precision for rate eq solving, used in rate_eq_solv.c
#define EpsRel 1e-4                  //integral relative precision, used in integral.c
#define EpsAbs 0.                    //integral absolute precision, used in integral.c

#define SLOW_MODE 0                  //use the thermal equilibrium population n[] as the initial value for each main loop
#define SCALE_A 0                    //Wether to scale a_matrix

//Constant (data from Wiki 2018.03)
#define h_CONST 6.626070040E-34      //Plank constant (J*s)
#define LIGHT_SPEED 299792458.0      //Speed of light (m/s)
#define k_CONST 1.38064852E-23       //Boltzmann constant (J*K^-1)

//debug use option-------------------------------
#define SHOW_A_V 0                   //Show the A[] and v[] that read from MOLE_DATA
#define SHOW_C 1                     //Show the C[] that read from MOLE_DATA
#define SHOW_NI 1                    //Show the initial n[] that calculated by n_initial_cal()
#define SHOW_NF 0                    //Show the finial n[] that calculated by n_initial_cal()
#define OUTPUT_A_MATRIX_I 1          //Output the a_matrix_i[] that calculated by a_matrix_initialize() 
                                     //and write it to the file A_MATRIX_I_FILE
#define A_MATRIX_I_FILE "a_matrix_i[].csv" //the output file of a_matrix_i[]
#define SHOW_E 1
#define SHOW_I_LIMIT 1
#define OUTPUT_E 0
#define E_FILE "E[].csv"
#define OUTPUT_BR_N 0
#define SHOW_BR_N 0
#define BR_N_FILE "Br_n[].csv"
#define SHOW_S 1
#define SHOW_I_INIT 1
#define SHOW_S_INIT 1
#define DEBUG_ZERO_R 0                // Turn off radiation in the rate equation for debugging purpose
#define SHOW_R 1
#define OUTPUT_I_R 1                  // Output integrand of R[][] to file for debugging
#define SHOW_N 1                      // Show n[] at each TAU step
#define CUTOFF_R 0                    // Set the minimum value for R[][] as R_MIN
#define R_MIN 0.
#define ABS_R 0                       // Take absolute value of R[][]
#define TEST_A_INIT 1                 // Test initialization of a_matrix for LEVEL_N == 3
#define TEST_A_INIT_SOLVE 1           // Test rate_eq_solve with only C[] terms and particle conservation (a_matrix_i)
#define TEST_RATE_EQ_FILL 1           // Test rate_eq_fill() for A[] terms without radiation (R[][] = 0)
#define TEST_S_ISO 1                  // Test source_f_n for isotropic cases
#define TEST_S_3 1
//debug use option-----------------------------//

//for Visual Stdio 2005 Compiler
#define _CRT_SECURE_NO_WARNINGS    //Disable the SECURE WARNINGS for fopen()
