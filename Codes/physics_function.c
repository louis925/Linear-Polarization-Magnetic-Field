/*
Made by Louis
NTHU
2009.04.19
*/

//Jp usually means J' = j
//dm usually means delta m

#include <stdio.h>
#include <math.h>
#include "physics_coeff.h"
#include "global_value.h"
#include "integral.h"
#include "physics_function.h"

int indexN(int J, int M)  //index conversion function ,(J,M ) to J(J+1)/2+|M| , for n[]
{
	if(M >= 0)
	{
		return (J*(J+1))/2 + M;
	}
	else
	{
		return (J*(J+1))/2 - M;
	}
}

int indexJJ(int J, int Jp)//index conversion function ,(J,J') to (J-1)J/2+J', for C[], E[]
{
	return (J*(J-1))/2 + Jp;
}

double A_coeff(int J, int M, int dM)//A coefficient for (J,M)->(j=J-1,M+dM)
{
	switch(dM)
	{
	case 0:
		return A[J-1]*(J*J-M*M)/(2*J-1)/J;
		break;
	case -1:
		return A[J-1]*(J+M)*(J-1+M)/2/(2*J-1)/J;
		break;
	case 1:
		return A[J-1]*(J-M)*(J-1-M)/2/(2*J-1)/J;
		break;
	default:
		return 0.0;
	}
}

double source_f_n(double n[TOTAL_N], double angle, int q, int J) //normalized Source function , = source_function/FJJ'
{//normalized source function for J to J' = J-1 //J -> j, J-1 = j = J'  //2009.11.23
	double c,s;
	int M;
	double sum_1, sum_2; //value temp. for sumation
	int A1_n, A0_n;
	int dA;

	switch(q)
	{
	case 1: //perpendicular, there is no direction dependent.
		M = -J;
		sum_1 = 0;
		sum_2 = 0;
		A1_n = J*(2*J-1);
		dA = 2*J-1;

		while(M < (J-1))
		{
			sum_1 += n[indexN(J,M)]*A1_n;
			sum_2 += n[indexN(J-1,M+1)]*A1_n;

			A1_n -= dA;
			dA--;
			M++;
		}
		return sum_1/(sum_2 - sum_1) / 2;
		break;
	case 0: //parallel
		c = angle;        //angle = cos(angle between magnetic field and line of sight)
		c *= c;           //c = cos^2
		s = fabs(1-c);
				
		M = -J;
		sum_1 = 0;
		sum_2 = 0;
		A1_n = J*(2*J-1);
		dA = 2*J-1;

		while(M < J)
		{
			A0_n = J*J - M*M;
			sum_1 += n[indexN(J,M)]*(A1_n*c + A0_n*s);
			sum_2 += (n[indexN(J-1,M+1)]*A1_n*c + n[indexN(J-1, M)]*A0_n*s);

			A1_n -= dA;
			dA--;
			M++;
		}
		return sum_1/(sum_2 - sum_1) / 2;
		break;
	default: //Error
		printf("q is out of range in source_f()!\n");
		return 0;
	}
}


double k_f_n(double n[TOTAL_N], double angle, int q, int j) //normalized Absorption efficients
{//J -> j, J-1 = j = J', change at 2009.11.13
	double c,s;
	int M;
	int J = j+1;
	int A1_n, A0_n;
	int dA;
	double sum;

	switch(q)
	{
	case 1: //perpendicular
		M = -J;
		sum = 0;
		A1_n = J*(2*J-1);
		dA = 2*J-1;

		while(M < (J-1))
		{
			sum += (n[indexN(J-1,M+1)]-n[indexN(J,M)])*A1_n;
			
			A1_n -= dA;
			dA--;
			M++;
		}
		
		return sum / (v[j]*v[j])/(2*J-1)/J*A[j];
		//return 2*sum / (v[j]*v[j])/(2*J-1)/J*A[j];//test, original eq in DW's paper, wrong
		break;
	case 0: //parallel
		c = angle;        //angle = angle between magnetic field and line of sight
		c *= c;           //c = cos^2
		s = fabs(1-c);
		
		M = -J;
		sum = 0;
		A1_n = J*(2*J-1);
		dA = 2*J-1;

		while(M < J)
		{
			A0_n = J*J - M*M;
			sum += ((n[indexN(J-1,M+1)]-n[indexN(J,M)])*A1_n*c + (n[indexN(J-1, M)]-n[indexN(J,M)])*A0_n*s);
			
			A1_n -= dA;
			dA--;
			M++;
		}
		return sum / (v[j]*v[j])/(2*J-1)/J*A[j];//modified at 2009.12.30
		break;
	default: //Error
		printf("q is out of range in k_f()!\n");
		return 0;
	}
}

double beta_f(double tau) //escape probability function
{
	return (1 - exp(-1*tau)) / tau;
}

double I_pa_n(double n[TOTAL_N],double angle, int q, double tau_d, int j) //normalized profile-averaged specific intensity(divided by _F)
{ //from J -> J', J-1 = J' = j
	double beta;
	beta = beta_f(tau_d);
	return source_f_n(n, angle, q, j+1) * (1 - beta) + Br_n[j]/2 * beta; 
}


void Rate_f_n_cal(double n[TOTAL_N], double tau, double R[LEVEL_N-1][2]) //calculation for multi-level normalized rate function
{//This fn doesn't call Rate_f_n()
	int j; //J -> j, J-1 = j = J' //change at 2009.11.13
	double error;
	unsigned int intervals;
	Fn_Param params;

	params.n = n;
	params.tau = tau; //the total optical depth tau[0] in this case(computation)
	params.k0 = k_f_n(n, cos(TAU_ANG),0,0);//[2010.01.22] OBS_ANG was replaced by TAU_ANG
	
	j = 0;
	while(j < (LEVEL_N-1))//j = 0 ~ LEVEL_N-2
	{
		params.j = j;
		R[j][0] = integral(i_f_r0m, &params, &error, &intervals) *3.0 /2.0;
		//R has already integrated over azimuth angle. Thus, it gains a 1/2 factor
		//The unit of R are the same as that in Deguchi and Watson's paper(1984)
		interval_count += intervals;
		
		R[j][1] = integral(i_f_r1m, &params, &error, &intervals) *3.0 /4.0;
		//R has already integrated over azimuth angle. Thus, it gains a 1/2 factor
		interval_count += intervals;
		j++;
	}
}

//[2010.01.22] OBS_ANG was replaced by TAU_ANG
double i_f_r0m(double x, void * params) //for R0, integral of I_pa_n[0]*sin^2
{
	static double *n;
	static double c,s;
	static double tau_d; //tau_d is the tau in this direction(angle) and polariation
	static int j; //J -> j, J-1 = j = J', change at 2009.11.13
	static double k0;
	n = ((Fn_Param *) params)->n;
	j = ((Fn_Param *) params)->j;
	tau_d = ((Fn_Param *) params)->tau;
	k0 =  ((Fn_Param *) params)->k0;
	c = x*x; //x = cos
	s = fabs(1-c);

#if TwoD
	tau_d = tau_d * (k_f_n(n,x,0,j)/k0) * sin(TAU_ANG)*sin(TAU_ANG)/ s;
#elif OneD
	tau_d = tau_d * (k_f_n(n,x,0,j)/k0) * (cos(TAU_ANG)*cos(TAU_ANG)/ c);
#elif Mix
	tau_d = tau_d * (k_f_n(n,x,0,j)/k0) * ((cos(TAU_ANG)*cos(TAU_ANG) + MixRatio*sin(TAU_ANG)*sin(TAU_ANG))/ (c + MixRatio*s));
#else //Isotropic
	tau_d = tau_d * (k_f_n(n,x,0,j)/k0);
#endif
	return s * I_pa_n(n,x,0,tau_d,j);
}

//[2010.01.22] OBS_ANG was replaced by TAU_ANG
double i_f_r1m(double x, void * params) //for R1, integral of I_pa_n[0]*cos^2+I_pa_n[1]
{
	static double *n;
	static double c,s;
	static double tau_d0, tau_d1; //tau_d is the tau in this direction(angle) and polariation
	static int j; //J -> j, J-1 = j = J', change at 2009.11.13
	static double k0;
	n = ((Fn_Param *) params)->n;
	j = ((Fn_Param *) params)->j;
	tau_d0 = ((Fn_Param *) params)->tau;
	tau_d1 = tau_d0;
	k0 =  ((Fn_Param *) params)->k0;
	c = x*x; //x = cos
	s = fabs(1-c);

#if TwoD //use s*s for 2D velocity field, or c*c for 1D velocity field
	tau_d0 = tau_d0 * (k_f_n(n,x,0,j)/k0) * sin(TAU_ANG)*sin(TAU_ANG)/ s;
	tau_d1 = tau_d1 * (k_f_n(n,x,1,j)/k0) * sin(TAU_ANG)*sin(TAU_ANG)/ s;
#elif OneD
	tau_d0 = tau_d0 * (k_f_n(n,x,0,j)/k0) * (cos(TAU_ANG)*cos(TAU_ANG)/ c);
	tau_d1 = tau_d1 * (k_f_n(n,x,1,j)/k0) * (cos(TAU_ANG)*cos(TAU_ANG)/ c);
#elif Mix
	tau_d0 = tau_d0 * (k_f_n(n,x,0,j)/k0) * ((cos(TAU_ANG)*cos(TAU_ANG) + MixRatio*sin(TAU_ANG)*sin(TAU_ANG))/ (c + MixRatio*s));
	tau_d1 = tau_d1 * (k_f_n(n,x,1,j)/k0) * ((cos(TAU_ANG)*cos(TAU_ANG) + MixRatio*sin(TAU_ANG)*sin(TAU_ANG))/ (c + MixRatio*s));
#else //Isotropic
	tau_d0 = tau_d0 * (k_f_n(n,x,0,j)/k0);
	tau_d1 = tau_d1 * (k_f_n(n,x,1,j)/k0);
#endif
	return c * I_pa_n(n, x, 0, tau_d0, j) + I_pa_n(n, x, 1, tau_d1, j);
}


void tau_array(double TAU, double tau[][2], double n[TOTAL_N]) //calculate all tau[][] from tau[0][0]
{
	int j;
	double k0;
	k0 = k_f_n(n, cos(TAU_ANG), 0, 0);//OBS_ANG was replaced by TAU_ANG [2010.01.22]
	j=0;
	while(j < (LEVEL_N - 1))
	{	
#if TwoD //use s*s for 2D velocity field, or c*c for 1D velocity field
		tau[j][0] = TAU * (k_f_n(n,cos(OBS_ANG),0,j)/k0) * (sin(TAU_ANG)*sin(TAU_ANG))/(sin(OBS_ANG)*sin(OBS_ANG));
		tau[j][1] = TAU * (k_f_n(n,cos(OBS_ANG),1,j)/k0) * (sin(TAU_ANG)*sin(TAU_ANG))/(sin(OBS_ANG)*sin(OBS_ANG));
#elif OneD
		tau[j][0] = TAU * (k_f_n(n,cos(OBS_ANG),0,j)/k0) * (cos(TAU_ANG)*cos(TAU_ANG))/(cos(OBS_ANG)*cos(OBS_ANG));
		tau[j][1] = TAU * (k_f_n(n,cos(OBS_ANG),1,j)/k0) * (cos(TAU_ANG)*cos(TAU_ANG))/(cos(OBS_ANG)*cos(OBS_ANG));
#elif Mix
		tau[j][0] = TAU * (k_f_n(n,cos(OBS_ANG),0,j)/k0) * ((cos(TAU_ANG)*cos(TAU_ANG) + MixRatio*sin(TAU_ANG)*sin(TAU_ANG))/(cos(OBS_ANG)*cos(OBS_ANG) + MixRatio*sin(OBS_ANG)*sin(OBS_ANG)));
		tau[j][1] = TAU * (k_f_n(n,cos(OBS_ANG),1,j)/k0) * ((cos(TAU_ANG)*cos(TAU_ANG) + MixRatio*sin(TAU_ANG)*sin(TAU_ANG))/(cos(OBS_ANG)*cos(OBS_ANG) + MixRatio*sin(OBS_ANG)*sin(OBS_ANG)));
#else //Isotropic
		tau[j][0] = TAU * (k_f_n(n,cos(OBS_ANG),0,j)/k0);
		tau[j][1] = TAU * (k_f_n(n,cos(OBS_ANG),1,j)/k0);
#endif

		j++;
	}
}
