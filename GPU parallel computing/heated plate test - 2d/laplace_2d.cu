#include<stdio.h>
#include<iostream>

/*
a alterar:
-solucionador não está sendo eficiente-creio que o mesmo processo é repetido em
 múltiplas threads, o que não é necessário
-resolver a questão dos resíduos
*/


void print_2d_array(int Nx,int Ny,int k,double *M){
  for(int j=Ny-1;j>=0;j--){
    for(int i=0;i<Nx;i++)
      printf("%f ",M[k*Nx*Ny+Ny*j+i]);

    putchar('\n');
  }
}

__device__
void dev_print_2d_array(int Nx,int Ny,int k,double *M){
  for(int j=Ny-1;j>=0;j--){
    for(int i=0;i<Nx;i++)
      printf("%f ",M[k*Nx*Ny+Ny*j+i]);

    printf("\n");
  }
}

// void print_2d_array(int Nx,int Ny,int k,double *M){
  // for(int j=0;j<Ny;j++){
    // for(int i=0;i<Nx;i++)
      // printf("%f ",M[k*Nx*Ny+Ny*j+i]);

    // putchar('\n');
  // }
// }

__global__
void initialize(int size,double *M){
	int index = blockIdx.x*blockDim.x + threadIdx.x;
  int stride = blockDim.x*gridDim.x;
  for(int i=index;i<size;i+=stride){
    M[i] = 0;
		// printf("%f ",M[i]);
	}
}

__device__
void solver(int Nx,int Ny,double *sol,int k,double beta){
	for(int j=1;j<Ny-1;j++){
		for(int i=1;i<Nx-1;i++){
			sol[k*Nx*Ny+Ny*j+i] = 1./(2.*(1.+beta*beta))*(sol[k*Nx*Ny+Ny*j+(i+1)] + 
														sol[k*Nx*Ny+Ny*j+(i-1)] + 
														beta*beta*sol[k*Nx*Ny+Ny*(j+1)+i] + 
														beta*beta*sol[k*Nx*Ny+Ny*(j-1)+i]);
		}
	}
}

// Regions are considered as below:
// | [0] [1] |
// | [2] [3] | 
__global__
void physical_boundary_conditions(int Nxs,int Nys,double *sol,double T_left,
																	double T_right,double T_up,double T_down){
	
	int index = blockIdx.x*blockDim.x + threadIdx.x;
  int stride = blockDim.x*gridDim.x;
	for(int j=index;j<Nys-1;j+=stride){
		sol[0*Nxs*Nys+Nys*j+0] = T_left;
		sol[2*Nxs*Nys+Nys*j+0] = T_left;
		sol[1*Nxs*Nys+Nys*j+(Nxs-1)] = T_right;
		sol[3*Nxs*Nys+Nys*j+(Nxs-1)] = T_right;
  }
  for(int i=index;i<Nxs-1;i+=stride){
		sol[2*Nxs*Nys+Nys*0+i] = T_down;
		sol[3*Nxs*Nys+Nys*0+i] = T_down;
		sol[0*Nxs*Nys+Nys*(Nys-1)+i] = T_up;
		sol[1*Nxs*Nys+Nys*(Nys-1)+i] = T_up;
  }
}

__global__
void phantom_boundary_conditions(int Nxs,int Nys,double *sol){
	int index = blockIdx.x*blockDim.x + threadIdx.x;
  int stride = blockDim.x*gridDim.x;
	// Vertical
	for(int j=index;j>0&&j<Nys;j+=stride){
		// Left to right
		sol[3*Nxs*Nys+Nys*j+0] = sol[2*Nxs*Nys+Nys*j+(Nxs-2)];
		sol[1*Nxs*Nys+Nys*(j-1)+0] = sol[0*Nxs*Nys+Nys*(j-1)+(Nxs-2)];
		
		// Right to left
		sol[2*Nxs*Nys+Nys*j+(Nxs-1)] = sol[3*Nxs*Nys+Nys*j+1];
		sol[0*Nxs*Nys+Nys*(j-1)+(Nxs-1)] = sol[1*Nxs*Nys+Nys*(j-1)+1];
	}
	
	// Horizontal
	for(int i=index;i>0&&i<Nxs;i+=stride){
		// Lower to upper
		sol[0*Nxs*Nys+Nys*0+i] = sol[2*Nxs*Nys+Nys*(Nys-2)+i];
		sol[1*Nxs*Nys+Nys*0+(i-1)] = sol[3*Nxs*Nys+Nys*(Nys-2)+(i-1)];
		
		// Upper to lower
		sol[2*Nxs*Nys+Nys*(Nys-1)+i] = sol[0*Nxs*Nys+Nys*1+i];
		sol[3*Nxs*Nys+Nys*(Nys-1)+(i-1)] = sol[1*Nxs*Nys+Nys*1+(i-1)];	
	}
}

__device__
void save_sol_old_loop2(int Nxs,int Nys,double *sol,double *sol_old,int k,
												int j){
	int index = blockIdx.x*blockDim.x + threadIdx.x;
  int stride = blockDim.x*gridDim.x;
	for(int i=index;i>0&&i<Nxs-1;i+=stride)
		sol_old[k*(Nxs-2)*(Nys-2)+(Nys-2)*(j-1)+(i-1)] = sol[k*Nxs*Nys+Nys*j+i];
}

__device__
void save_sol_old(int Nxs,int Nys,double *sol,double *sol_old,int k){
	int index = blockIdx.x*blockDim.x + threadIdx.x;
  int stride = blockDim.x*gridDim.x;
	for(int j=index;j>0&&j<Nys-1;j+=stride)
		save_sol_old_loop2(Nxs,Nys,sol,sol_old,k,j);
}

__device__
void get_residuals_loop2(int Nxs,int Nys,double *sol,double *sol_old,double *res,
									 int k,int j){
	int index = blockIdx.x*blockDim.x + threadIdx.x;
  int stride = blockDim.x*gridDim.x;
	for(int i=index;i>0&&i<Nxs-1;i+=stride)
		res[k*(Nxs-2)*(Nys-2)+(Nys-2)*(j-1)+(i-1)] = fabs(sol[k*Nxs*Nys+Nys*j+i] - 
							 sol_old[k*(Nxs-2)*(Nys-2)+(Nys-2)*(j-1)+(i-1)]);
}

__device__
void get_residuals(int Nxs,int Nys,double *sol,double *sol_old,double *res,
									 int k){
	int index = blockIdx.x*blockDim.x + threadIdx.x;
  int stride = blockDim.x*gridDim.x;
	for(int j=index;j>0&&j<Nys-1;j+=stride)
		get_residuals_loop2(Nxs,Nys,sol,sol_old,res,k,j);
}

double sum_residuals(int Nxs,int Nys,double *res){
	double res_sum = 0.;
	for(int i=0;i<Nxs*Nys*4;i++)
		res_sum += res[i];
	
	return res_sum;
}

__global__ 
void iterate_once(int Nxs,int Nys,double *sol,double *sol_old,double *res,
									double beta){
	int index = blockIdx.x*blockDim.x + threadIdx.x;
  int stride = blockDim.x*gridDim.x;
	// printf("iterate_once: %d %d\n",index,stride);
	for(int k=index;k<4;k+=stride){
		save_sol_old(Nxs,Nys,sol,sol_old,k); // Save current values for convergence checking
		// dev_print_2d_array(Nxs-2,Nys-2,k,sol_old);
		solver(Nxs,Nys,sol,k,beta);
		get_residuals(Nxs,Nys,sol,sol_old,res,k);
	}
}

// __global__
void parallel_solver(int Nxs,int Nys,double *sol,double *sol_old,double *res,
										 int iter,double beta,double eps,int blocksPerGrid,
										 int threadsPerBlock){
	double res_val;
	for(int loop=0;loop<iter;loop++){
		// Calculate one iteration for each region
		iterate_once<<<blocksPerGrid,threadsPerBlock>>>(Nxs,Nys,sol,sol_old,res,beta);
		cudaDeviceSynchronize();
		
		for(int k=0;k<4;k++){
			printf("%d\n",k);
			print_2d_array(Nxs,Nys,k,sol);
			putchar('\n');
			print_2d_array(Nxs-2,Nys-2,k,sol_old);
			putchar('\n');
			print_2d_array(Nxs-2,Nys-2,k,res);
			putchar('\n');
		}
		getchar();
		
		// Check convergence
		res_val = sum_residuals(Nxs-2,Nys-2,res);
		if(res_val <= eps){
			puts("Convergence!");
			break;
		}
		
		printf("Iteration %d | Residuals = %f\n",loop+1,res_val);
		
		// Update phantom boundary conditions
		phantom_boundary_conditions<<<blocksPerGrid,threadsPerBlock>>>(Nxs,Nys,sol);
		cudaDeviceSynchronize();
	}
	
	// cudaFree(res_val);
}

int main(){
	time_t start,end;
	time(&start);
	
	// Problem parameters
  double Lx = 1., Ly = 1.;
  double T_down = 0., T_up = 0.;
  double T_left = 50., T_right = 50.;

  // Mesh parameters
  int Nx = 10, Ny = 10; // Complete domain (even values)
  double deltaX = Lx/((double) Nx);
  double deltaY = Ly/((double) Ny);
  double beta = deltaX/deltaY;
	
	// Numeric solution configuration
  int iter = 50000;
  double eps = 1e-3;
	
	// Initialize arrays
	Nx = (int) (Nx/2)*2; // Ensure values are even
	Ny = (int) (Ny/2)*2;
	int Nxs = Nx/2+1; // Dimensions for the submatrices
	int Nys = Ny/2+1;
	double *sol = new double[4*Nxs*Nys];
	double *sol_old = new double[4*(Nxs-2)*(Nys-2)];
	double *res = new double[4*(Nxs-2)*(Nys-2)];
	cudaMallocManaged(&sol,sizeof(double)*4*Nxs*Nys);
	cudaMallocManaged(&sol_old,sizeof(double)*4*(Nxs-2)*(Nys-2));
	cudaMallocManaged(&res,sizeof(double)*4*(Nxs-2)*(Nys-2));
	
	// Allocate relevant arrays in GPU memory and initialize
	int threadsPerBlock = 256;
	int blocksPerGrid = (4*Nxs*Nys+threadsPerBlock-1)/threadsPerBlock;
	// dim3 threadsPerBlock(threadsPerBlock_,threadsPerBlock_,1);
	// dim3 blocksPerGrid(blocksPerGrid_,blocksPerGrid_,1);
	initialize<<<blocksPerGrid,threadsPerBlock>>>(4*Nxs*Nys,sol);
	initialize<<<blocksPerGrid,threadsPerBlock>>>(4*(Nxs-2)*(Nys-2),sol_old);
	initialize<<<blocksPerGrid,threadsPerBlock>>>(4*(Nxs-2)*(Nys-2),res);
	cudaDeviceSynchronize();
	
  // Insert physical boundary conditions
	physical_boundary_conditions<<<blocksPerGrid,threadsPerBlock>>>(Nxs,Nys,sol,T_left,T_right,T_up,T_down);
	cudaDeviceSynchronize();
	
  // Obtain solution
  parallel_solver(Nxs,Nys,sol,sol_old,res,iter,beta,eps,blocksPerGrid,threadsPerBlock);
	cudaDeviceSynchronize();

  // Print results to file
  FILE *results;
  results = fopen("results.txt","w");
  fprintf(results,"%.8f %.8f ",Lx,Ly);
  for(int i=2;i<Nx;i++)
    fprintf(results,"%.8f ",0.);
  fprintf(results,"\n");
	
	for(int j=Nys-1;j>0;j--){
		for(int i=0;i<Nxs-1;i++)
			fprintf(results,"%.8f ",sol[0*Nxs*Nys+Nys*j+i]);
		for(int i=1;i<Nxs;i++)
			fprintf(results,"%.8f ",sol[1*Nxs*Nys+Nys*j+i]);
		
		fprintf(results,"\n");
	}
	for(int j=Nys-2;j>=0;j--){
		for(int i=0;i<Nxs-1;i++)
			fprintf(results,"%.8f ",sol[2*Nxs*Nys+Nys*j+i]);
		for(int i=1;i<Nxs;i++)
			fprintf(results,"%.8f ",sol[3*Nxs*Nys+Nys*j+i]);
		
		fprintf(results,"\n");
	}
  fclose(results);
	
	// Free memories
	cudaFree(sol);
	cudaFree(sol_old);
	cudaFree(res);
	// free(sol_host);
	// free(sol_old_host);
	// delete[] sol;
	// delete[] sol_old;
	
	time(&end);
	double time_taken = double(end - start);
	printf("\n%f s\n",time_taken);
	
	puts("end");

  return 0;
}