#include<stdio.h>
#include"custom_functions.h"

int main()
{
	/* Entradas: matrizes A e r*/
	double A[6][6] = {{1,2,0,0,0,0},
				     {3,4,5,0,0,0},
				     {0,6,3,4,0,0},
  		  		     {0,0,9,1,2,0},
	  			     {0,0,0,3,4,5},
				     {0,0,0,0,1,2}};
	double r[6] = {10,
				  3,
				  5,
				  5,
				  0,
				  20};
	double x[6]; /* Vetor solução */
	
	int N = sizeof(r)/sizeof(r[0]);
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
	
	print_array(sizeof(x)/sizeof(x[0]),x);
	
	
	// printf("%f\n",A[3][4]);
	// printf("Número de elementos: %lu\n",sizeof(A)/sizeof(A[0][0]));
	// printf("Número de linhas:% lu\n",sizeof(A)/sizeof(A[0]));
	// printf("Número de colunas:% lu\n",sizeof(A[0])/sizeof(A[0][0]));
	
	// print_mat(sizeof(A)/sizeof(A[0]),sizeof(A[0])/sizeof(A[0][0]),A);
	// print_mat(sizeof(r)/sizeof(r[0]),sizeof(r[0])/sizeof(r[0][0]),r);
	
	
	
	return 0;
}