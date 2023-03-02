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

/* Mesh data */
const int N = 1000; // Number of points

int main()
{
	/* Start measuring time */
	clock_t start = clock();
	
	/* Initialize mesh */
	double x_vec[N]; linspace(0,3,N,x_vec);
	double N2 = N, delta_x = 3/(N2-1);

	/* Initialize solution vectors and nozzle shape vector*/
	double rho_vec[N];
	double T_vec[N];
	double V_vec[N];
	double a_vec[N];
	double p_vec[N];
	double mach_vec[N];
	double A_vec[N];
	for(int i=0;i<N;i++)
	{
		rho_vec[i] = 1 - 0.3146*x_vec[i];
		T_vec[i] = 1 - 0.2314*x_vec[i];
		V_vec[i] = (0.1 + 1.09*x_vec[i])*pow(T_vec[i],0.5);
		a_vec[i] = sqrt(T_vec[i]);
		p_vec[i] = 1;
		mach_vec[i] = 0;
		A_vec[i] = 1 + 2.2*pow((x_vec[i]-1.5),2);
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
	
	/* Main loop */
	int loop = 0;
	double delta_t;
	double diff_rho_F; // Derivatives (forward-differences)
	double diff_V_F;
	double diff_T_F;
	double diff_rho_R; // Derivatives (rearward-differences)
	double diff_V_R;
	double diff_T_R;
	double diff_rho_AV; // Derivatives (average)
	double diff_V_AV;
	double diff_T_AV;
	double rho_bar_old; // Outdated predicted values
	double V_bar_old;
	double T_bar_old;
	double rho_bar; // Predicted values
	double V_bar;
	double T_bar;
	for(loop;loop<iter;loop++)
	{
		printf("Iteration %d\n",loop+1);
		
		/* Obtain time step */
		delta_t = calc_delta_t(C,delta_x,N,a_vec,V_vec);
		
		/* Calculate predicted values at point 0 */
		/* Forward differences */
		diff_rho_F = -rho_vec[0]*(V_vec[1] - V_vec[0])/delta_x - 
					 rho_vec[0]*V_vec[0]*(log(A_vec[1])-log(A_vec[0]))/delta_x - 
					 V_vec[0]*(rho_vec[1]-rho_vec[0])/delta_x;
		diff_V_F = -V_vec[0]*(V_vec[1]-V_vec[0])/delta_x - 1/GA*((T_vec[1]-T_vec[0])/delta_x + 
					T_vec[0]/rho_vec[0]*(rho_vec[1]-rho_vec[0])/delta_x);
		diff_T_F = -V_vec[0]*(T_vec[1]-T_vec[0])/delta_x - (GA-1)*T_vec[0]*((V_vec[1]-V_vec[0])/delta_x + 
					V_vec[0]*(log(A_vec[1])-log(A_vec[0]))/delta_x);
		/* Obtain predicted values */
		rho_bar_old = rho_vec[0] + diff_rho_F*delta_t;
		V_bar_old = V_vec[0] + diff_V_F*delta_t;
		T_bar_old = T_vec[0] + diff_T_F*delta_t;
		
		/* Calculate properties of all internal points */
		for(int i=1;i<N-2;i++)
		{
			/* Forward-differences */
			diff_rho_F = -rho_vec[i]*(V_vec[i+1] - V_vec[i])/delta_x -
						 rho_vec[i]*V_vec[i]*(log(A_vec[i+1])-log(A_vec[i]))/delta_x - 
						 V_vec[i]*(rho_vec[i+1]-rho_vec[i])/delta_x;
			diff_V_F = -V_vec[i]*(V_vec[i+1]-V_vec[i])/delta_x - 1/GA*((T_vec[i+1]-T_vec[i])/delta_x +
						T_vec[i]/rho_vec[i]*(rho_vec[i+1]-rho_vec[i])/delta_x);
			diff_T_F = -V_vec[i]*(T_vec[i+1]-T_vec[i])/delta_x - (GA-1)*T_vec[i]*((V_vec[i+1]-V_vec[i])/delta_x +
						V_vec[i]*(log(A_vec[i+1])-log(A_vec[i]))/delta_x);
			/* Obtain predicted values */
			rho_bar = rho_vec[i] + diff_rho_F*delta_t;
			V_bar = V_vec[i] + diff_V_F*delta_t;
			T_bar = T_vec[i] + diff_T_F*delta_t;
			/* Rearward-differences */
			diff_rho_R = -rho_bar*(V_bar-V_bar_old)/delta_x - 
						 rho_bar*V_bar*(log(A_vec[i])-log(A_vec[i-1]))/delta_x - 
						 V_bar*(rho_bar-rho_bar_old)/delta_x;
			diff_V_R = -V_bar*(V_bar-V_bar_old)/delta_x - 1/GA*((T_bar-T_bar_old)/delta_x + 
						T_bar/rho_bar*(rho_bar-rho_bar_old)/delta_x);
			diff_T_R = -V_bar*(T_bar-T_bar_old)/delta_x - (GA-1)*T_vec[i]*((V_bar-V_bar_old)/delta_x + 
						V_bar*(log(A_vec[i])-log(A_vec[i-1]))/delta_x);
			/* Average time derivatives */
			diff_rho_AV = (diff_rho_F + diff_rho_R)/2;
			diff_V_AV = (diff_V_F + diff_V_R)/2;
			diff_T_AV = (diff_T_F + diff_T_R)/2;
			/* Corrected values */
			rho_vec[i] = rho_vec[i] + diff_rho_AV*delta_t;
			V_vec[i] = V_vec[i] + diff_V_AV*delta_t;
			T_vec[i] = T_vec[i] + diff_T_AV*delta_t;
			p_vec[i] = rho_vec[i]*T_vec[i];
			/* Keep predicted values to use at the next point */
			rho_bar_old = rho_bar;
			V_bar_old = V_bar;
			T_bar_old = T_bar;
			/* Keep residuals at the throat */
			if(i == TH_pos)
			{
				his_res_rho[loop] = fabs(diff_rho_AV);
				his_res_V[loop] = fabs(diff_V_AV);
				his_res_T[loop] = fabs(diff_T_AV);
			}
		}
		
		/* Calculate pressure at the end */
		p_vec[N-1] = rho_vec[N-1]*T_vec[N-1];
		
		/* Recalculate a_vec and obtain mach numbers*/
		for(int i=0;i<N;i++)
		{
			a_vec[i] = sqrt(T_vec[i]);
			mach_vec[i] = V_vec[i]/a_vec[i];
		}
		
		/* Apply boundary conditions */
		V_vec[0] = 2*V_vec[1] - V_vec[2];
		V_vec[N-1] = 2*V_vec[N-2] - V_vec[N-3];
		rho_vec[N-1] = 2*rho_vec[N-2] - rho_vec[N-3];
		T_vec[N-1] = 2*T_vec[N-2] - T_vec[N-3];
		
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
	
	/* Print results to text files */
	/* File 1: air proerties along the nozzle */
	/* Columns: x,rho,T,V,a,p,mach,A */
	FILE *f1 = fopen("c7_3_supersonic_isentropic_nozzle1.res","w");
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
	FILE *f2 = fopen("c7_3_supersonic_isentropic_nozzle2.res","w");
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
	fprintf(gnuplot_pipe,"%s\n","plot \"c7_3_supersonic_isentropic_nozzle1.res\" using 1:2 title 'rho' with linespoints linestyle 1,\
									  \"c7_3_supersonic_isentropic_nozzle1.res\" using 1:3 title 'T' with linespoints linestyle 2,\
									  \"c7_3_supersonic_isentropic_nozzle1.res\" using 1:4 title 'V' with linespoints linestyle 3,\
									  \"c7_3_supersonic_isentropic_nozzle1.res\" using 1:5 title 'a' with linespoints linestyle 4,\
									  \"c7_3_supersonic_isentropic_nozzle1.res\" using 1:6 title 'p' with linespoints linestyle 5,\
									  \"c7_3_supersonic_isentropic_nozzle1.res\" using 1:7 title 'mach' with linespoints linestyle 6");
	
	/* Plot historics at the throat */
	fprintf(gnuplot_pipe,"%s\n","set terminal wxt 2");
	fprintf(gnuplot_pipe,"%s\n","set grid");
	fprintf(gnuplot_pipe,"%s\n","plot \"c7_3_supersonic_isentropic_nozzle2.res\" using 1:5 title 'rho' with linespoints linestyle 1,\
									  \"c7_3_supersonic_isentropic_nozzle2.res\" using 1:6 title 'T' with linespoints linestyle 2,\
									  \"c7_3_supersonic_isentropic_nozzle2.res\" using 1:8 title 'V' with linespoints linestyle 3,\
									  \"c7_3_supersonic_isentropic_nozzle2.res\" using 1:7 title 'p' with linespoints linestyle 5,\
									  \"c7_3_supersonic_isentropic_nozzle2.res\" using 1:9 title 'mach' with linespoints linestyle 6");
	/* Plot residuals */
	fprintf(gnuplot_pipe,"%s\n","set terminal wxt 3");
	fprintf(gnuplot_pipe,"%s\n","set grid");
	fprintf(gnuplot_pipe,"%s\n","plot \"c7_3_supersonic_isentropic_nozzle2.res\" using 1:2 title 'rho' with linespoints linestyle 1,\
									  \"c7_3_supersonic_isentropic_nozzle2.res\" using 1:4 title 'T' with linespoints linestyle 2,\
									  \"c7_3_supersonic_isentropic_nozzle2.res\" using 1:3 title 'V' with linespoints linestyle 3");
	
	// set grid
	// set style line 1 linetype 1 linewidth 2 pointsize 0
	// plot "c7_3_supersonic_isentropic_nozzle1.res" using 1:7 w lp with linespoints linestyle 1

	
	// https://stackoverflow.com/questions/14311640/how-to-plot-data-by-c-program
	// https://stackoverflow.com/questions/3521209/making-c-code-plot-a-graph-automatically
	// http://www.gnuplotting.org/plotting-data/
	// https://stackoverflow.com/questions/45802125/plot-only-specific-columns-in-gnuplot
	
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