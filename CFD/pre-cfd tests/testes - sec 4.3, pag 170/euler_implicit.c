#include<stdio.h>
#include<math.h>
#include<string.h>
// #include<stdlib.h>
#include"custom_functions.h"

#define pi 3.141592653589793


/* Declarar funções */
void thomas_algorithm(int nAx,int nAy,double A[nAx][nAy],int N,double r[N],double x[N]);
double heated_bar_analytic(double x,double L,double T_Im,double alpha,double t,int sum_max);

/* Propriedades físicas */
const double L = 1; // [m]
const double T_init = 0; // [°C]
const double T_Im = 100; // [°C]
const double T_0 = 0; // [°C]
const double alpha = 0.0834; // [m^2/s]

/* Dados da malha e da solução */
const int Im = 5; // (seis pontos no total, cinco elementos)
const double delta_x = L/Im;
const double s = 1/20; // (cumpre o requisito s <= 1/2)
const double delta_t = 0.1;
const int iter = 4/delta_t; // Número de iterações

int main()
{
	
	/* Matriz pra guardar as soluções */
	double sol[iter+1][Im+1]; // Os seis pontos na barra, em cada instante de tempo
	memset(sol,0,sizeof(sol[0][0])*cols(sol)*rows(sol));
	
	/* Já se sabe a temperatura da barra inteira em t_0 */
	for(int i=0;i<Im;i++)
	{
		sol[0][i] = T_init;
	}

	/* Também ja se sabe as temperaturas das pontas: */
	for(int i=0;i<iter+1;i++)
	{
		sol[i][0] = T_0;
		sol[i][Im] = T_Im;
	}

	/* Preparar solução analítica */
	double sol_an[iter+1][Im+1];
	memcpy(sol_an, sol, (iter+1)*(Im+1)*sizeof(double));  //print_mat(cols(sol_an),rows(sol_an),*sol_an);putchar('\n');
	double x[Im+1];
	linspace(0,L,Im+1,x); //print_array(rows(x),x);
	double t = 0;
	
	/* Montar a matriz A */
	double S = alpha*delta_t/pow(delta_x,2);
	double A[Im-1][Im-1]; memset(A,0,sizeof(A[0][0])*cols(A)*rows(A));
	/* Primeira linha */
	A[0][0] = 1 + 2*S;
	A[0][1] = -S;
	/* Linhas intermediárias */
	for(int i=1;i<Im-2;i++)
	{
		A[i][i-1] = -S;
		A[i][i] = 1 + 2*S;
		A[i][i+1] = -S;
	}
	/* Última linha */
	A[Im-2][Im-3] = -S;
	A[Im-2][Im-2] = 1 + 2*S;
	
	/* Inicilizar mais algumas coisas */
	int loop;
	double r[Im-1];
	double x_thomas[Im-1];
	double errors[iter+1][Im-1];
	for(loop=1;loop<=iter;loop++)
	{	
		printf("Iteração %d\n",loop);

		/*** Solução numérica ***/
		/* Montar vetor r */
		for(int i=0;i<Im-1;i++)
		{
			r[i] = sol[loop-1][i+1];
		}
		r[0] += S*sol[loop][0];
		r[Im-2] += S*sol[loop][Im];
		thomas_algorithm(cols(A),rows(A),A,Im-1,r,x_thomas);
		/* Inserir resultados na matriz de solução */
		for(int i=0;i<Im-1;i++)
		{
			sol[loop][i+1] = x_thomas[i];
		}
		
		/*** Solução analítica ***/
		for(int node=1;node<Im;node++)
		{
			sol_an[loop][node] = heated_bar_analytic(x[node],L,T_Im,alpha,t,100);
		}
		
		/* Calcular erros entre a solução aproximada e a exata */
		for(int i=1;i<Im-1;i++)
		{
			errors[loop][i-1] = (sol[loop][i]-sol_an[loop][i])/sol[loop][i]*100;
		}
		
		/* Determinar o próximo instante no tempo */
		/* (importante pra solução analítica) */
		t += delta_t;
	}
	
	printf("\nSolução numérica:\n");
	print_mat(cols(sol),rows(sol),*sol);
	printf("\nSolução analítica:\n");
	print_mat(cols(sol_an),rows(sol_an),*sol_an);
	printf("\nErros (%%)\n");
	print_mat(cols(errors),rows(errors),*errors);
}

void thomas_algorithm(int nAx,int nAy,double A[nAx][nAy],int N,double *r,double *x)
{	
	double C[N];
	double D[N];
	
	/* Operações */
	C[0] = A[0][1]/A[0][0];
	D[0] = r[0]/A[0][0];
	
	for(int i=1;i<N;i++)
	{
		if((i-1)<N)
		{
			C[i] = A[i][i+1]/(A[i][i] - A[i][i-1]*C[i-1]);
		}
		
		D[i] = (r[i] - A[i][i-1]*D[i-1])/(A[i][i] - A[i][i-1]*C[i-1]);
	}
	
	x[N-1] = D[N-1];
	for(int i=N-2;i>=0;i--)
	{
		x[i] = D[i] - C[i]*x[i+1];
	}
}

double heated_bar_analytic(double x,double L,double T_Im,double alpha,double t,int sum_max)
{
	double sum=0;
	
	for(int n=1;n<=sum_max;n++)
	{
		sum += pow((-1),n)*2*T_Im/(n*pi)*sin(n*pi*x/L)*exp(-alpha*t*pow((n*pi/L),2));
	}
	
	return x/L*T_Im + sum;
}