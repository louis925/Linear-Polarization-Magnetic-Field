/*
GSL Liner Algebra Solving Example:
http://www.gnu.org/software/gsl/manual/html_node/Linear-Algebra-Examples.html

Made by Louis
NTHU
2010.01.19
*/

#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <math.h>

#include "physics_function.h"
#include "physics_coeff.h"
#include "global_value.h"
//#define N_loop 50

#include "rate_eq_solv.h"

void rate_eq_fill(double a_matrix[TOTAL_N*TOTAL_N], double R[LEVEL_N-1][2]);         //fill a_matrix with A[] and R[], call AR_fill_0 and _1
void AR_fill_0(double a_row[TOTAL_N], double R[LEVEL_N-1][2], int J);        //fill a_matrix with A[] and R[], for M = 0
void AR_fill_1(double a_row[TOTAL_N], double R[LEVEL_N-1][2], int J, int M); //fill a_matrix with A[] and R[], for M > 0

void rate_eq_solve(double n[TOTAL_N], double TAU)
{
	double b[TOTAL_N] = {Nt, 0.0};  
	double n_last[TOTAL_N];           // store the result n[] in last step for comparing to next step
	double R[LEVEL_N-1][2] = {{0.0}}; // results of Rate_f_n()
	int converge;                     // weither n[] converge
	int s, j, l;
	
#if 0
	int i,k; //debug use
	FILE *amf; //debug use
#endif

	gsl_permutation *gsl_p;
	gsl_matrix_view gsl_a_m;
	gsl_vector_view gsl_b_v;
	gsl_vector_view gsl_n_v;
	gsl_p = gsl_permutation_alloc (TOTAL_N);
	gsl_a_m = gsl_matrix_view_array (a_matrix, TOTAL_N, TOTAL_N);
	gsl_b_v = gsl_vector_view_array (b, TOTAL_N);
	gsl_n_v = gsl_vector_view_array (n, TOTAL_N);

	l = 0;
	do {
		Rate_f_n_cal(n, TAU, R); // Calculate the R[] coefficients
			
#if 0
		//print n[] for debug*****
		//if(TAU > 0.018) //debug for LEVEL_N=7 tau=0.01905460718
		//{
			printf("%d:\n", i);//*****
			j=0;
			while(j < TOTAL_N)
			{
				printf("%.5e,", n[j]);
				j++;
			}
			printf("\n");//*****/
		//}
#endif	

		// Copy n[] to n_last[]
		for(j = 0; j < TOTAL_N; j++) { 
			n_last[j] = n[j]; 
		}
		// Copy a_matrix_i[] to a_matrix[]
		for(j = 0; j < TOTAL_N*TOTAL_N; j++) {
			a_matrix[j] = a_matrix_i[j];
			//a_matrix[j] = 0.0;
		}
		
		// Fill rate equations(a_matrix[]) with A[] and R[] ____________________________
		rate_eq_fill(a_matrix, R);
		// ___________________________________________________________________________//

#if 0
		//check a_matrix[]********
		if(TAU == 0.01 && l == 0)
		{
			printf("[Debug]a_matrix[][] output.\n");
			if( (amf = fopen( A_MATRIX_I_FILE, "w" )) == NULL )//Open file for debug a_matrix output
			{
				printf( "Can not open the file '%s' for debug output\n" , A_MATRIX_I_FILE);
			}
			i = 0;
			j = 0;
			k = 0;
			while(j < TOTAL_N)
			{
				k += TOTAL_N;
				while(i < k)
				{
					fprintf(amf, "%.5e",a_matrix[i]);
					fprintf(amf, ",");
					i++;
				}
				fprintf(amf, "\n");
				j++;
			}
			if(amf)
			{
				if( fclose(amf) )
				{
					printf( "The file '%s' was not closed\n" , A_MATRIX_I_FILE);
				}
			}			
		}//*/
#endif	
		// _____________________________________________________________________________
		// Compute n[] by solving a[TOTAL_N][TOTAL_N] x n[TOTAL_N] = b[TOTAL_N] problem
		// using the LU decomposition method
		gsl_linalg_LU_decomp (&gsl_a_m.matrix, gsl_p, &s);
		gsl_linalg_LU_solve (&gsl_a_m.matrix, gsl_p, &gsl_b_v.vector, &gsl_n_v.vector);
		loop_count++;
		// ___________________________________________________________________________//
		
#if 0
		//if(tau[0][0] > 0.018) //debug for LEVEL_N=7 tau=0.01905460718
		//{
			printf("%d:\n", i);//*****
			j=0;//print n[] for debug*****
			while(j < TOTAL_N)
			{
				printf("%.7e, ", n[j]);
				j++;
			}
			printf("\n");//*****/
		//}
#endif	
		
		for (j = 0; j < TOTAL_N; j++) // Checks the accuracy to decide whether it converges
		{
			if( fabs(n[j] - n_last[j]) < fabs(n[j]) * REL_PREC)	{
				converge = 1; // OK, Converge
			}
			else {
				converge = 0; // Doesn't converge, run again
				j = TOTAL_N;  // Skip checking the rest of j
			}		
		}
		l++;
	}
	while(!converge); // Stop when converge. //l < N_loop
	
#if 0
	//print n[] for debug*****
	//printf("%d:\n", i);
	printf("\nn[] for debug\n");
	j=0;
	while(j < TOTAL_N)
	{
		printf("%+.3e, ", n[j]);
		j++;
	}
	printf("\n");
#endif

	gsl_permutation_free (gsl_p);	
}

void n_initial_cal(double n[TOTAL_N], double t)//New version, use thermal equilibrium at t K¡@directly
{
	double e[LEVEL_N];
	double x,et = 0;
	int j,m;
	
	x = h_CONST/k_CONST*LIGHT_SPEED/t*-100.0;
	
	j = 0;
	while(j < LEVEL_N)
	{
		e[j] = exp(x*(energy_level[j] - energy_level[0]));
		et += e[j]*(2*j+1);
		j++;
	}

	j = 0;
	while(j < LEVEL_N)
	{
		m = 0;
		while(m <= j)
		{
			n[indexN(j,m)] = e[j]/et*Nt;
			m++;
		}
		j++;
	}
}//*/

void a_matrix_initialize(double a_matrix[TOTAL_N*TOTAL_N])//2009.11.12 Check OK (for C coeff)
{
	//a_matrix : a[(J,M)][(J',M')] = a[J(J+1)/2+M][J'(J'+1)/2+M'] := a[(J(J+1)/2+M)*TOTAL_N + J'(J'+1)/2+M']
	//int J, M, Mp;
	int i, k, j, i2;   //index for elements in each row, j=J'
	int J, M;          //index for rows in the matrix
	double cjj;
	
	//fill the a[0][] row, particle number conservation--------
	i = 0;
	j = 0;
	k = 0;
	while(j < LEVEL_N)
	{
		k += (j+1);

		a_matrix[i] = 1.0; //for M'=0
		i++;
		while(i < k)
		{
			a_matrix[i] = 2.0; //for M'>0
			i++;
		}
		j++;
	}//------------------------------------------------------//
	//printf("row 1 OK\n");

	//Now the index i,k are at a[1][0]

	//fill collisional excitation rates C[] to a[][]-----------
	J = 1;
	while(J < LEVEL_N)//run over all the a_matrix from a[1][] to a[TOTAL_N-1][]
	{
		
		M = 0;
		while(M < (J+1))
		{
			//for one row, have the same (J,M)------
			i2 = i; //move i2 to the first element in this row, a[(J,M)][0]
			i2 += indexN(J,M); //then move it to point at a[(J,M)][(J,M)]
			a_matrix[i2] = 0.0;

			//run over all the element in this row, from a[(J,M)][0] to a[(J,M)][TOTAL_N-1]
			j = 0; //j = J', in this loop
			
			while(j < J)//for each J' < J, from a[(J,M)][J'=0] to a[(J,M)][J'=J-1]
			{
				k += (j+1);
				
				cjj = C[indexJJ(J,j)]; //CJJ'
				a_matrix[i2] -= cjj;
				cjj *= ( E[indexJJ(J,j)] / (2*j+1) );
				
				a_matrix[i] = cjj; //for M'=0
				i++;
				while(i < k)
				{
					a_matrix[i] = 2.0 * cjj; //for M'>0
					i++;
				}
				j++;
			}
			//Now j = J

			//for J' = J, fill a[][] with zero, unless M = M'
			k += (j+1);
			while(i < k)
			{
				if(i != i2)
				{
					a_matrix[i] = 0.0;
				}
				i++;
			}
			j++;
			//Now j = J+1
			
			
			while(j < LEVEL_N)//for each J'=j > J, from a[(J,M)][J'=J+1] to a[(J,M)][J'=LEVEL_N-1]
			{
				k += (j+1);
				
				cjj = C[indexJJ(j,J)] / (2*J+1);
				a_matrix[i2] -= (cjj * (2*j+1) * E[indexJJ(j,J)]);
				
				a_matrix[i] = cjj; //for M'=0
				i++;
				while(i < k)
				{
					a_matrix[i] = 2.0 * cjj; //for M'>0
					i++;
				}
				j++;
			}
			//------------------------------------//
			
			M++;
		}
		J++;
	}
	//-------------------------------------------------------//
}

void rate_eq_fill(double a_matrix[TOTAL_N*TOTAL_N], double R[LEVEL_N-1][2])
{
	//fill the rate equations (a_matrix[]) with the contribution from A[] and R[]
	//rates function R[] for each level needed to be calculate first before call this function
	int i;
	int J,M;

	i = TOTAL_N;//Now the index i,k are at a[1][0]
	J = 1;
	//for J = 1 ~ J = LEVEL_N-1
	while(J < LEVEL_N)
	{
		M = 0;
		AR_fill_0( &a_matrix[i], R, J);
		i += TOTAL_N;
		M++;

		while(M <= J)//change at 2010/01/19
		{
			AR_fill_1( &a_matrix[i], R, J, M);
			i += TOTAL_N;
			M++;
		}
		J++;
	}
}

void AR_fill_0(double a_row[TOTAL_N], double R[LEVEL_N-1][2], int J) //M = 0
{
	double At[6]; //temp. for A coefficients, At[] = {AJM(J-1)M	, AJM(J-1)(M+1), AJM(J-1)(M-1), A(J+1)MJ, A(J+1)(M-1)JM, A(J+1)(M+1)JM}
	int j;
	At[0] = A_coeff(J  ,  0,  0); //AJ0J-10
	At[1] = A_coeff(J  ,  0,  1); //AJ0J-11
	
	j = indexN(J,0);
	a_row[j] -= (  At[0]*(1 + R[J-1][0]) + 2*At[1]*(1 + R[J-1][1]) ); //****R[J]
	if(J < (LEVEL_N - 1))
	{
		At[3] = A_coeff(J+1,  0,  0); //AJ+10J0
		At[4] = A_coeff(J+1, -1,  1); //AJ+1-1J0
		a_row[j] -= ( At[3]*R[J][0] + 2*At[4]*R[J][1] );
	}

	j = indexN(J-1,0);
	a_row[j  ] +=   At[0]*R[J-1][0];
	a_row[j+1] += 2*At[1]*R[J-1][1];
	
	if(J < (LEVEL_N - 1))
	{
		j = indexN(J+1,0);
		a_row[j  ] +=   At[3]*(1 + R[J][0]);
		a_row[j+1] += 2*At[4]*(1 + R[J][1]);
	}
}

void AR_fill_1(double a_row[TOTAL_N], double R[LEVEL_N-1][2], int J, int M) //M > 0
{
	double At[6]; //temp. for A coefficients, At[] = {AJM(J-1)M	, AJM(J-1)(M+1), AJM(J-1)(M-1), A(J+1)MJ, A(J+1)(M-1)JM, A(J+1)(M+1)JM}
	int j;
	At[0] = A_coeff(J  , M  ,  0); //dM = 0
	At[1] = A_coeff(J  , M  ,  1); //dM = +1
	At[2] = A_coeff(J  , M  , -1); //dM = -1
	
	j = indexN(J,M);
	a_row[j] -= (  At[0] * (1 + R[J-1][0]) + (At[1] + At[2]) * (1 + R[J-1][1]));
	if(J < (LEVEL_N - 1))
	{
		At[3] = A_coeff(J+1, M  ,  0); //dM = 0
		At[4] = A_coeff(J+1, M-1,  1); //dM = +1
		At[5] = A_coeff(J+1, M+1, -1); //dM = -1
		a_row[j] -= ( At[3] *   R[J  ][0]  + (At[4] + At[5]) *      R[J  ][1] );
	}
	
	j = indexN(J-1,M);
	a_row[j  ] += At[0]*R[J-1][0]; //n(J-1,M)
	a_row[j+1] += At[1]*R[J-1][1]; //n(J-1,M+1)
	a_row[j-1] += At[2]*R[J-1][1]; //n(J-1,M-1)
	
	if(J < (LEVEL_N - 1))
	{
		j = indexN(J+1,M);
		a_row[j  ] += At[3]*(1 + R[J][0]); //n(J+1,M)
		a_row[j+1] += At[5]*(1 + R[J][1]); //n(J+1,M+1) //Modified at 2010/01/19
		a_row[j-1] += At[4]*(1 + R[J][1]); //n(J+1,M-1) //Modified at 2010/01/19
	}
}
