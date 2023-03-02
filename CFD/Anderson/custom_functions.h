#include<stdio.h>

#define rows(M) (sizeof(M)/sizeof(M[0])) // For 2D arrays and 1D arrays
#define cols(M) (sizeof(M[0])/sizeof(M[0][0])) // For 2D arrays
#define numel2(M) (sizeof(M)/sizeof(M[0][0])) // For 2D arrays

/*sizeof(A)/sizeof(A[0][0]) encontrar o número de elementos de uma matriz A*/
/*sizeof(A)/sizeof(A[0]) número de linhas */
/*sizeof(A[0])/sizeof(A[0][0]) número de colunas */

void print_array(int n,double *vec)
{
	for(int i=0;i<n;i++)
	{
		printf("%.2f ",vec[i]);
	}
	putchar('\n');
}

// void print_mat(int nx, int ny, const double *mat) ????????
// {
	// /* nx -> number of columns */
	// /* ny -> number of rows */
	// for(int j=0;j<ny;j++)
	// {
		// for(int i=0;i<nx;i++)
		// {
			// printf("%f ",mat[j*ny+i]);
			// printf("%f, %d %d,",mat[j][i],j,i);
		// }
		// putchar('\n');
	// }	
	
	
// }

void print_mat(size_t width, size_t height, const double *A)
{
	for(size_t i = 0; i < height; ++i)
	{
		for(size_t j = 0; j < width; ++j)
		{
			// printf("A[%zu][%zu] = %d\n", i, j, A[i * width + j]);
			printf("%f ",A[i*width + j]);
		}
		putchar('\n');
		
	}
  // https://stackoverflow.com/questions/8715034/print-2d-array-by-calling-a-function-print-array-with-argument-the-2d-array
}

void ones(int nx, int ny, double mat[nx][ny])
{
	for(int i=0;i<nx;i++)
	{
		for(int j=0;j<ny;j++)
		{
			mat[i][j] = 1;
		}
	}
}

void linspace(double p0,double pf,int N,double *vec)
{
	double sep = (pf-p0)/(N-1); // Separation value
	// printf("%f\n\n",sep);
	
	vec[0] = p0;
	for(int i=1;i<N;i++)
	{
		vec[i] = vec[i-1]+sep;
	}
}

double max_val_1D(int n,double *vec)
{
	double tmp = vec[0];
	for(int i=0;i<n;i++)
	{
		if(vec[i]>tmp)
			tmp = vec[i];
		
	}
	return tmp;
}

int max_pos_1D(int n,double *vec)
{
	double tmp = vec[0];
	int pos = 0, i;
	for(i=0;i<n;i++)
	{
		if(vec[i]>tmp)
		{
			tmp = vec[i];
			pos = i;
		}
	}
	return pos;
}

double min_val_1D(int n,double *vec)
{
	double tmp = vec[0];
	for(int i=0;i<n;i++)
	{
		if(vec[i]<tmp)
			tmp = vec[i];
		
	}
	return tmp;
}

int min_pos_1D(int n,double *vec)
{
	double tmp = vec[0];
	int pos = 0, i;
	for(i=0;i<n;i++)
	{
		if(vec[i]<tmp)
		{
			tmp = vec[i];
			pos = i;
		}
	}
	return pos;
}