#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include"custom_functions.h"

double calc_delta_t(double C,double delta_x,int N,double *a_vec,double *V_vec);

/* Solution data */
const int iter = 50000; // Number of iterations (main loop)
const double C = 0.5; // Courant number
const double GA = 1.4;
const double pN = 0.6784;
const double Cx = 0.2; // Artificial viscosity coefficient
const int print = 1; // Print results and plot graphs?

/* Mesh data */
const int N = 100; // Number of points

int main()
{
	/* Start measuring time */
	clock_t start = clock();
	
	/* Initialize mesh */
	double x_vec[N]; linspace(0,3,N,x_vec);
	double N2 = N, delta_x = 3/(N2-1);

	/* Initialize solution vectors and nozzle shape vector */
	double rho_vec[N];
	double T_vec[N];
	double V_vec[N];
	double a_vec[N];
	double p_vec[N];
	double mach_vec[N];
	double A_vec[N];
	double U1[N];
	double U2[N];
	double U3[N];
	for(int i=0;i<N;i++)
	{
		if(x_vec[i] <= 0.5)
		{
			rho_vec[i] = 1;
			T_vec[i] = 1;
		}
		else if(x_vec[i] > 0.5 && x_vec[i] <= 1.5)
		{
			rho_vec[i] = 1 - 0.366*(x_vec[i]-0.5);
			T_vec[i] = 1 - 0.167*(x_vec[i]-0.5);
		}
		else if(x_vec[i] > 1.5 && x_vec[i] <= 2.1)
		{
			rho_vec[i] = 0.634 - 0.702*(x_vec[i]-1.5);
			T_vec[i] = 0.833 - 0.4908*(x_vec[i]-1.5);
		}
		else
		{
			rho_vec[i] = 0.5892 + 0.10228*(x_vec[i]-2.1);
			T_vec[i] = 0.93968 + 0.0622*(x_vec[i]-2.1);
		}
		
		A_vec[i] = 1 + 2.2*pow((x_vec[i]-1.5),2);
		a_vec[i] = sqrt(T_vec[i]);
		V_vec[i] = 0.59/(rho_vec[i]*A_vec[i]);
		mach_vec[i] = V_vec[i]/a_vec[i];
		p_vec[i] = 1;
		U1[i] = rho_vec[i]*A_vec[i];
		U2[i] = rho_vec[i]*A_vec[i]*V_vec[i];
		U3[i] = rho_vec[i]*(T_vec[i]/(GA-1)+GA/2*pow(V_vec[i],2))*A_vec[i];
	}
	/* Throat position */
	int TH_pos = min_pos_1D(rows(A_vec),A_vec);
	
	/* Vectors for historics (throat) */
	double his_V[iter];
	double his_mach[iter];
	double his_rho[iter];
	double his_T[iter];
	double his_p[iter];
	double his_res_rho[iter];
	double his_res_V[iter];
	double his_res_T[iter];
	
	/* Vectors to store parameters of the equations */
	/* (necessary due to the artificial viscosity mechanism) */
	double diff_U1_F[N];
	double diff_U2_F[N];
	double diff_U3_F[N];
	double U1_bar[N];
	double U2_bar[N];
	double U3_bar[N];
	
	/* Main loop */
	int loop = 0;
	double delta_t;
	double diff_U1_R; // Derivatives (rearward-differences)
	double diff_U2_R;
	double diff_U3_R;
	double diff_U1_AV; // Derivatives (average)
	double diff_U2_AV;
	double diff_U3_AV;
	double rho_bar_old; // Outdated predicted values
	double V_bar_old;
	double T_bar_old;
	double rho_bar[N]; // Predicted values
	double V_bar[N];
	double T_bar[N];
	double p_bar[N];
	double S1; // Artificial viscosity terms
	double S2;
	double S3;
	double S1_bar;
	double S2_bar;
	double S3_bar;
	double F1i; // Terms used to calculate the derivatives
	double F1ip;
	double F2i;
	double F2ip;
	double F3i;
	double F3ip;
	double J2F,J2R;
	double F1_bar;
	double F2_bar;
	double F3_bar;
	double F1_bar_old;
	double F2_bar_old;
	double F3_bar_old;
	for(loop;loop<iter;loop++)
	{
		printf("Iteration %d\n",loop+1);
		
		/* Obtain time step */
		delta_t = calc_delta_t(C,delta_x,N,a_vec,V_vec);
		
		/* Calculate predicted values at point 1 */
		// Calculate relevant terms
		F1i = U2[0];
		F1ip = U2[1];
		F2i = pow(U2[0],2)/U1[0] + (GA-1)/GA*(U3[0]-GA/2*pow(U2[0],2)/U1[0]);
		F2ip = pow(U2[1],2)/U1[1] + (GA-1)/GA*(U3[1]-GA/2*pow(U2[1],2)/U1[1]);
		F3i = GA*U2[0]*U3[0]/U1[0] - GA*(GA-1)/2*pow(U2[0],3)/pow(U1[0],2);
		F3ip = GA*U2[1]*U3[1]/U1[1] - GA*(GA-1)/2*pow(U2[1],3)/pow(U1[1],2);
		J2F = (GA-1)/GA*(U3[0]-GA/2*pow(U2[0],2)/U1[0])*(log(A_vec[1])-log(A_vec[0]))/delta_x;
		
		/* Forward-differences */
		diff_U1_F[0] = -(F1ip - F1i)/delta_x;
		diff_U2_F[0] = -(F2ip - F2i)/delta_x + J2F;
		diff_U3_F[0] = -(F3ip - F3i)/delta_x;
		
		/* Obtain predicted values */
		U1_bar[0] = U1[0] + diff_U1_F[0]*delta_t;
		U2_bar[0] = U2[0] + diff_U2_F[0]*delta_t;
		U3_bar[0] = U3[0] + diff_U3_F[0]*delta_t;

		/* Keep variable values */
		// Calculate relevant terms
		F1_bar_old = U2_bar[0];
		F2_bar_old = pow(U2_bar[0],2)/U1_bar[0] + (GA-1)/GA*(U3_bar[0]-GA/2*pow(U2_bar[0],2)/U1_bar[0]);
		F3_bar_old = GA*U2_bar[0]*U3_bar[0]/U1_bar[0] - GA*(GA-1)/2*pow(U2_bar[0],3)/pow(U1_bar[0],2);

		/* Calculate properties of all internal points */
		/* First part */
		for(int i=1;i<N-1;i++)
		{
			// Calculate relevant terms
			F1i = U2[i];
			F1ip = U2[i+1];
			F2i = pow(U2[i],2)/U1[i] + (GA-1)/GA*(U3[i]-GA/2*pow(U2[i],2)/U1[i]);
			F2ip = pow(U2[i+1],2)/U1[i+1] + (GA-1)/GA*(U3[i+1]-GA/2*pow(U2[i+1],2)/U1[i+1]);
			F3i = GA*U2[i]*U3[i]/U1[i] - GA*(GA-1)/2*pow(U2[i],3)/pow(U1[i],2);
			F3ip = GA*U2[i+1]*U3[i+1]/U1[i+1] - GA*(GA-1)/2*pow(U2[i+1],3)/pow(U1[i+1],2);
			J2F = (GA-1)/GA*(U3[i]-GA/2*pow(U2[i],2)/U1[i])*(log(A_vec[i+1])-log(A_vec[i]))/delta_x;
			
			/* Forward-differences */
			diff_U1_F[i] = -(F1ip - F1i)/delta_x;
			diff_U2_F[i] = -(F2ip - F2i)/delta_x + J2F;
			diff_U3_F[i] = -(F3ip - F3i)/delta_x;
			
			/* Artificial viscosity terms */
			S1 = Cx*fabs(p_vec[i+1]-2*p_vec[i]+p_vec[i-1])/
				(p_vec[i+1]+2*p_vec[i]+p_vec[i-1])*(U1[i+1]-2*U1[i]+U1[i-1]);
			S2 = Cx*fabs(p_vec[i+1]-2*p_vec[i]+p_vec[i-1])/
				(p_vec[i+1]+2*p_vec[i]+p_vec[i-1])*(U2[i+1]-2*U2[i]+U2[i-1]);
			S3 = Cx*fabs(p_vec[i+1]-2*p_vec[i]+p_vec[i-1])/
				(p_vec[i+1]+2*p_vec[i]+p_vec[i-1])*(U3[i+1]-2*U3[i]+U3[i-1]);
			
			/* Obtain predicted values */
			U1_bar[i] = U1[i] + diff_U1_F[i]*delta_t + S1;
			U2_bar[i] = U2[i] + diff_U2_F[i]*delta_t + S2;
			U3_bar[i] = U3[i] + diff_U3_F[i]*delta_t + S3;
		}
		
		/* Estimate pressures */
		U1_bar[N-1] = U1[N-1];
		U2_bar[N-1] = U2[N-1];
		U3_bar[N-1] = U3[N-1];
		for(int i=0;i<N-1;i++)
		{
			rho_bar[i] = U1_bar[i]/A_vec[i];
			V_bar[i] = U2_bar[i]/U1_bar[i];
			T_bar[i] = (GA-1)*(U3_bar[i]/U1_bar[i]-GA/2*pow(V_bar[i],2));
			p_bar[i] = rho_bar[i]*T_bar[i];
		}
		
		/* Second part */
		/* (this separation is necessary due to the U_{i+1} terms in the artificial 
		    viscosity formulas) */
		for(int i=1;i<N-1;i++)
		{
			/* Calculate relevant terms */
			F1_bar = U2_bar[i];
			F2_bar = pow(U2_bar[i],2)/U1_bar[i] + (GA-1)/GA*(U3_bar[i]-GA/2*pow(U2_bar[i],2)/U1_bar[i]);
			F3_bar = GA*U2_bar[i]*U3_bar[i]/U1_bar[i] - GA*(GA-1)/2*pow(U2_bar[i],3)/pow(U1_bar[i],2);
			J2R = (GA-1)/GA*(U3_bar[i]-GA/2*pow(U2_bar[i],2)/U1_bar[i])*(log(A_vec[i])-log(A_vec[i-1]))/delta_x;
			
			/* Rearward-differences */
			diff_U1_R = -(F1_bar - F1_bar_old)/delta_x;
			diff_U2_R = -(F2_bar - F2_bar_old)/delta_x + J2R;
			diff_U3_R = -(F3_bar - F3_bar_old)/delta_x;
			   
			/* Average time derivatives */
			diff_U1_AV = (diff_U1_F[i] + diff_U1_R)/2;
			diff_U2_AV = (diff_U2_F[i] + diff_U2_R)/2;
			diff_U3_AV = (diff_U3_F[i] + diff_U3_R)/2;
			
			/* Artificial viscosity terms */
			S1_bar = Cx*fabs(p_bar[i+1]-2*p_bar[i]+p_bar[i-1])/
				(p_bar[i+1]+2*p_bar[i]+p_bar[i-1])*(U1_bar[i+1]-2*U1_bar[i]+U1_bar[i-1]);
			S2_bar = Cx*fabs(p_bar[i+1]-2*p_bar[i]+p_bar[i-1])/
				(p_bar[i+1]+2*p_bar[i]+p_bar[i-1])*(U2_bar[i+1]-2*U2_bar[i]+U2_bar[i-1]);
			S3_bar = Cx*fabs(p_bar[i+1]-2*p_bar[i]+p_bar[i-1])/
				(p_bar[i+1]+2*p_bar[i]+p_bar[i-1])*(U3_bar[i+1]-2*U3_bar[i]+U3_bar[i-1]);
			
			/* Corrected values */
			U1[i] = U1[i] + diff_U1_AV*delta_t + S1_bar;
			U2[i] = U2[i] + diff_U2_AV*delta_t + S2_bar;
			U3[i] = U3[i] + diff_U3_AV*delta_t + S3_bar;
			
			/* Keep predicted values to use at the next point */
			F1_bar_old = F1_bar;
			F2_bar_old = F2_bar;
			F3_bar_old = F3_bar;
			
			/* Keep residuals at the throat */
			if(i == TH_pos)
			{
				his_res_rho[loop] = fabs(diff_U1_AV);
				his_res_V[loop] = fabs(diff_U2_AV);
				his_res_T[loop] = fabs(diff_U3_AV);
			}
		}
		
		/* Check presence of NaNs */
		/* (throat chosen arbitrarily) */
		if(U1[TH_pos]!=U1[TH_pos] || U2[TH_pos]!=U2[TH_pos] || U3[TH_pos]!=U3[TH_pos])
		{
			puts("Broken solution: presence of NaNs");
			return 0;
		}
		
		/* Recalculate all flow properties (internal points) */
		for(int i=1;i<N-1;i++)
		{
			rho_vec[i] = U1[i]/A_vec[i];
			V_vec[i] = U2[i]/U1[i];
			T_vec[i] = (GA-1)*(U3[i]/U1[i]-GA/2*pow(V_vec[i],2));
			p_vec[i] = rho_vec[i]*T_vec[i];
			a_vec[i] = sqrt(T_vec[i]);
			mach_vec[i] = V_vec[i]/a_vec[i];
		}
		
		/* Apply boundary conditions */
		U2[0] = 2*U2[1] - U2[2];
		U3[0] = U1[0]*(T_vec[0]/(GA-1)+GA/2*pow(V_vec[0],2));
		U1[N-1] = 2*U1[N-2] - U1[N-3];
		U2[N-1] = 2*U2[N-2] - U2[N-3];
		U3[N-1] = pN*A_vec[N-1]/(GA-1) + GA/2*U2[N-1]*V_vec[N-1];
		//
		rho_vec[0] = U1[0]/A_vec[0];
		V_vec[0] = U2[0]/U1[0];
		T_vec[0] = (GA-1)*(U3[0]/U1[0]-GA/2*pow(V_vec[0],2));
		p_vec[0] = rho_vec[0]*T_vec[0];
		a_vec[0] = sqrt(T_vec[0]);
		mach_vec[0] = V_vec[0]/a_vec[0];
		//
		rho_vec[N-1] = U1[N-1]/A_vec[N-1];
		V_vec[N-1] = U2[N-1]/U1[N-1];
		T_vec[N-1] = (GA-1)*(U3[N-1]/U1[N-1]-GA/2*pow(V_vec[N-1],2));
		p_vec[N-1] = rho_vec[N-1]*T_vec[N-1];
		a_vec[N-1] = sqrt(T_vec[N-1]);
		mach_vec[N-1] = V_vec[N-1]/a_vec[N-1];
		
		/* Store values at the throat */
		his_rho[loop] = rho_vec[TH_pos];
		his_T[loop] = T_vec[TH_pos];
		his_p[loop] = p_vec[TH_pos];
		his_V[loop] = V_vec[TH_pos];
		his_mach[loop] = mach_vec[TH_pos];
		
	}
	
	/* Print results to screen */
	clock_t end = clock();
	float seconds = (float)(end - start)/CLOCKS_PER_SEC;
	printf("Time: %.2fs\n",seconds);
	// print_array(rows(mach_vec),mach_vec);
	printf("%.2f\n",mach_vec[TH_pos]);
	
	if(print==1)
	{
		/* Print results to text files */
		/* File 1: air properties along the nozzle */
		/* Columns: x,rho,T,V,a,p,mach,A */
		FILE *f1 = fopen("c7_6_shock_capturing1.res","w");
		if (f1 == NULL)
		{
			printf("Error opening file!\n");
			exit(1);
		}
		for(int i=0;i<N;i++)
		{
			fprintf(f1,"%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n",x_vec[i],
					rho_vec[i],T_vec[i],V_vec[i],a_vec[i],p_vec[i],mach_vec[i],A_vec[i]); // Air properties along the nozzle
		}
		fclose(f1);
		/* File 2: historics of the residuals and air properties at the throat */
		/* Columns: loop,res_rho,res_V,res_T,his_rho,his_T,his_p,his_V,his_mach */
		FILE *f2 = fopen("c7_6_shock_capturing2.res","w");
		if (f2 == NULL)
		{
			printf("Error opening file!\n");
			exit(1);
		}
		for(int i=0;i<iter;i++)
		{
			fprintf(f2,"%d %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f \n",i+1,
					his_res_rho[i],his_res_V[i],his_res_T[i], 						     // Residuals at the throat
					his_rho[i],his_T[i],his_p[i],his_V[i],his_mach[i]);					 // Air properties at the throat
		}
		fclose(f2);
		
		/* Plot air properties */
		FILE *gnuplot_pipe = popen("gnuplot -persistent", "w");
		fprintf(gnuplot_pipe,"%s\n","set terminal wxt 1");
		fprintf(gnuplot_pipe,"%s\n","set grid");
		fprintf(gnuplot_pipe,"%s\n","set style line 1 linetype 1 linewidth 2 pointsize 0 linecolor rgb '#0000ff'");
		fprintf(gnuplot_pipe,"%s\n","set style line 2 linetype 1 linewidth 2 pointsize 0 linecolor rgb '#cc0000'");
		fprintf(gnuplot_pipe,"%s\n","set style line 3 linetype 1 linewidth 2 pointsize 0 linecolor rgb '#00FF00'");
		fprintf(gnuplot_pipe,"%s\n","set style line 4 linetype 1 linewidth 2 pointsize 0 linecolor rgb '#FF00FF'");
		fprintf(gnuplot_pipe,"%s\n","set style line 5 linetype 1 linewidth 2 pointsize 0 linecolor rgb '#FFFF00'");
		fprintf(gnuplot_pipe,"%s\n","set style line 6 linetype 1 linewidth 2 pointsize 0 linecolor rgb '#000000'");
		fprintf(gnuplot_pipe,"%s\n","plot \"c7_6_shock_capturing1.res\" using 1:2 title 'rho' with linespoints linestyle 1,\
										  \"c7_6_shock_capturing1.res\" using 1:3 title 'T' with linespoints linestyle 2,\
										  \"c7_6_shock_capturing1.res\" using 1:4 title 'V' with linespoints linestyle 3,\
										  \"c7_6_shock_capturing1.res\" using 1:5 title 'a' with linespoints linestyle 4,\
										  \"c7_6_shock_capturing1.res\" using 1:6 title 'p' with linespoints linestyle 5,\
										  \"c7_6_shock_capturing1.res\" using 1:7 title 'mach' with linespoints linestyle 6");
		
		/* Plot historics at the throat */
		fprintf(gnuplot_pipe,"%s\n","set terminal wxt 2");
		fprintf(gnuplot_pipe,"%s\n","set grid");
		fprintf(gnuplot_pipe,"%s\n","plot \"c7_6_shock_capturing2.res\" using 1:5 title 'rho' with linespoints linestyle 1,\
										  \"c7_6_shock_capturing2.res\" using 1:6 title 'T' with linespoints linestyle 2,\
										  \"c7_6_shock_capturing2.res\" using 1:8 title 'V' with linespoints linestyle 3,\
										  \"c7_6_shock_capturing2.res\" using 1:7 title 'p' with linespoints linestyle 5,\
										  \"c7_6_shock_capturing2.res\" using 1:9 title 'mach' with linespoints linestyle 6");
		/* Plot residuals */
		fprintf(gnuplot_pipe,"%s\n","set terminal wxt 3");
		fprintf(gnuplot_pipe,"%s\n","set grid");
		fprintf(gnuplot_pipe,"%s\n","plot \"c7_6_shock_capturing2.res\" using 1:2 title 'rho' with linespoints linestyle 1,\
										  \"c7_6_shock_capturing2.res\" using 1:4 title 'T' with linespoints linestyle 2,\
										  \"c7_6_shock_capturing2.res\" using 1:3 title 'V' with linespoints linestyle 3");
	}
	
	return 0;
}

double calc_delta_t(double C,double delta_x,int N,double *a_vec,double *V_vec)
{
	/* Take the minimum value of all possible time steps */
	double vec[N];
	for(int i=0;i<N;i++)
	{
		vec[i] = C*delta_x/(a_vec[i] + V_vec[i]);
	}
	return min_val_1D(rows(vec),vec);
}
