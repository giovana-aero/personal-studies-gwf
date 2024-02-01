#include<stdio.h>
#include<ctime>
// #include<iostream>

/*
v2 difference: variable domain division
*/

void print_2d_array(int Nx,int Ny,int k,double *M);
__global__
void initialize(int size,double *M);
__device__
void solver(int Nx,int Ny,double *sol,int k,double beta);
void physical_boundary_conditions(int Nxs,int Nys,int divX,int divY,double *sol,
																	double T_left,double T_right,double T_up,
																	double T_down);
void vertical_phantom_BC(int Nxs,int Nys,double *sol,int k);
void horizontal_phantom_BC(int Nxs,int Nys,int divX,double *sol,int k);
int element_in_array(int size,int *array,int num);
void phantom_boundary_conditions(int Nxs,int Nys,int divX,int divY,double *sol,
																 int *jump_points);
__device__
void save_sol_old(int Nxs,int Nys,double *sol,double *sol_old,int size);
__device__
void get_residuals(int Nxs,int Nys,int divX,int divY,double *sol,double *sol_old
									 ,double *res);
double sum_residuals(int Nxs,int Nys,double *res);
__global__ 
void iterate_once(int Nxs,int Nys,int divX,int divY,double *sol,double *sol_old,
									double *res,double beta,int size);
void get_jump_points(int divX,int divY,int *jump_points);
void parallel_solver(int Nxs,int Nys,int divX,int divY,double *sol,
										 double *sol_old,double *res,int iter,double beta,
										 double eps,int BpG,int TpB);
void print_results(int Nx,int Ny,int Nxs,int Nys,int divX,int divY,double Lx,
									 double Ly,double *sol);

int main(){
	time_t start,end;
	time(&start);
	
	// Problem parameters
  double Lx = 1., Ly = 1.;
  double T_down = -0., T_up = -100.;
  double T_left = 50., T_right = 50.;

  // Mesh parameters
  int Nx = 402, Ny = 402; // Complete domain (even values)
  double deltaX = Lx/((double) Nx);
  double deltaY = Ly/((double) Ny);
  double beta = deltaX/deltaY;
	
	// Numeric solution configuration
  int iter = 1000000;
  double eps = 1e-3;
	int divX = 10, divY = 10; // Number of domain divisions
	if((Nx-2)%divX != 0){
		puts("Error: (Nx-2) must be divisible by divX");
		return 1;
	}
	if((Ny-2)%divY != 0){
		puts("Error: (Ny-2) must be divisible by divY");
		return 2;
	}
	
	// Allocate relevant arrays in GPU memory
	Nx = (int) (Nx/2)*2; // Ensure values are even
	Ny = (int) (Ny/2)*2;
	int Nxs = (Nx-2)/divX+2; // Dimensions for the submatrices
	int Nys = (Ny-2)/divY+2;
	if(Nxs<3 || Nys<3){
		puts("Error: Nxs and Nys must be equal or larger than 3. Redefine mesh \
				 configurations");
		return 3;
	}
	double *sol = new double[Nxs*Nys*divX*divY];
	double *sol_old = new double[(Nxs-2)*(Nys-2)*divX*divY];
	double *res = new double[(Nxs-2)*(Nys-2)*divX*divY];
	cudaMallocManaged(&sol,sizeof(double)*Nxs*Nys*divX*divY);
	cudaMallocManaged(&sol_old,sizeof(double)*(Nxs-2)*(Nys-2)*divX*divY);
	cudaMallocManaged(&res,sizeof(double)*(Nxs-2)*(Nys-2)*divX*divY);
	
	// Initialize arrays
	int threadsPerBlock = 1024;
	int blocksPerGrid = (Nxs*Nys*divX*divY+threadsPerBlock-1)/threadsPerBlock;
	initialize<<<blocksPerGrid,threadsPerBlock>>>(Nxs*Nys*divX*divY,sol);
	initialize<<<blocksPerGrid,threadsPerBlock>>>((Nxs-2)*(Nys-2)*divX*divY,
						sol_old);
	initialize<<<blocksPerGrid,threadsPerBlock>>>((Nxs-2)*(Nys-2)*divX*divY,res);
	cudaDeviceSynchronize();
	
  // Insert physical boundary conditions
	physical_boundary_conditions(Nxs,Nys,divX,divY,sol,T_left,T_right,T_up,
															 T_down);
	
  // Obtain solution
  parallel_solver(Nxs,Nys,divX,divY,sol,sol_old,res,iter,beta,eps,blocksPerGrid,
									threadsPerBlock);

  // Print results to file
  print_results(Nx,Ny,Nxs,Nys,divX,divY,Lx,Ly,sol);
	
	// Free memories
	cudaFree(sol);
	cudaFree(sol_old);
	cudaFree(res);
	
	time(&end);
	double time_taken = double(end - start);
	printf("\n%f s\n",time_taken);

	puts("end");

  return 0;
}

void print_2d_array(int Nx,int Ny,int k,double *M){
  for(int j=Ny-1;j>=0;j--){
    for(int i=0;i<Nx;i++)
      printf("%f ",M[k*Nx*Ny+Ny*j+i]);

    putchar('\n');
  }
}

// __device__
// void dev_print_2d_array(int Nx,int Ny,int k,double *M){
  // for(int j=Ny-1;j>=0;j--){
    // for(int i=0;i<Nx;i++)
      // printf("%f ",M[k*Nx*Ny+Ny*j+i]);

    // printf("\n");
  // }
// }

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
	}
}

__device__
void solver(int Nx,int Ny,double *sol,int k,double beta){
	if(blockIdx.x*blockDim.x + threadIdx.x == k){
		for(int j=1;j<Ny-1;j++){
			for(int i=1;i<Nx-1;i++)
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
void physical_boundary_conditions(int Nxs,int Nys,int divX,int divY,double *sol,
																	double T_left,double T_right,double T_up,
																	double T_down){
	
	// Left and right
	for(int k1=0,k2=divX-1;k1<=divY*(divX-1);k1+=divX,k2+=divX){
		for(int j=0;j<Nys;j++){
			sol[k1*Nxs*Nys+Nys*j+0] = T_left;
			sol[k2*Nxs*Nys+Nys*j+(Nxs-1)] = T_right;
		}
	}
	
	// Up and down
	for(int k1=0,k2=divX*divY-1;k1<divX;k1++,k2--){
		for(int i=0;i<Nxs;i++){
			sol[k1*Nxs*Nys+Nys*(Nys-1)+i] = T_up;
			sol[k2*Nxs*Nys+Nys*0+i] = T_down;
		}
	}
}

void vertical_phantom_BC(int Nxs,int Nys,double *sol,int k){
	for(int j=1;j<Nys-1;j++){
		// Left to right <-
		sol[(k+1)*Nxs*Nys+Nys*j+0] = sol[k*Nxs*Nys+Nys*j+(Nxs-2)];			
		
		// Right to left
		sol[k*Nxs*Nys+Nys*j+(Nxs-1)] = sol[(k+1)*Nxs*Nys+Nys*j+1];

	}
}

void horizontal_phantom_BC(int Nxs,int Nys,int divX,double *sol,int k){
	for(int i=1;i<Nxs-1;i++){
		// Lower to upper
		sol[k*Nxs*Nys+Nys*0+i] = sol[(k+divX)*Nxs*Nys+Nys*(Nys-2)+i];
		
		// Upper to lower
		sol[(k+divX)*Nxs*Nys+Nys*(Nys-1)+i] = sol[k*Nxs*Nys+Nys*1+i];
	}
}

int element_in_array(int size,int *array,int num){
	for(int i=0;i<size;i++){
		if(array[i] == num)
			return 1;
	}
	return 0;
}

void phantom_boundary_conditions(int Nxs,int Nys,int divX,int divY,double *sol,
																 int *jump_points){
	// Vertical
	for(int k=0;k<divY*divX-1;k++){
		if(element_in_array(divY-1,jump_points,k))
			k++;
		// printf("%d\n",k);
		vertical_phantom_BC(Nxs,Nys,sol,k);
	}
			
	// Horizontal
	for(int k=0;k<divX*(divY-1);k++){
		horizontal_phantom_BC(Nxs,Nys,divX,sol,k);
	}
}

__device__
void save_sol_old(int Nxs,int Nys,double *sol,double *sol_old,int size){
	int index = blockIdx.x*blockDim.x + threadIdx.x;
  int stride = blockDim.x*gridDim.x;
	int k;
	int row;
	int col;
	
	for(int i=index;i<(Nxs-2)*(Nys-2)*size;i+=stride){
		k = i/((Nxs-2)*(Nys-2));
		row = (i - k*(Nxs-2)*(Nys-2))/(Nxs-2);
		col = i - k*(Nxs-2)*(Nys-2) - (Nys-2)*row;
		sol_old[i] = sol[k*Nxs*Nys+Nys*(row+1)+(col+1)];
	}
}

__device__
void get_residuals(int Nxs,int Nys,int divX,int divY,double *sol,double *sol_old
									 ,double *res){
	int index = blockIdx.x*blockDim.x + threadIdx.x;
  int stride = blockDim.x*gridDim.x;
	int k;
	int row;
	int col;
	for(int i=index;i<(Nxs-2)*(Nys-2)*divX*divY;i+=stride){
		k = i/((Nxs-2)*(Nys-2));
		row = (i - k*(Nxs-2)*(Nys-2))/(Nxs-2);
		col = i - k*(Nxs-2)*(Nys-2) - (Nys-2)*row;
		res[i] = fabs(sol[k*Nxs*Nys+Nys*(row+1)+(col+1)] - 
							 sol_old[i]);
	}
}

double sum_residuals(int Nxs,int Nys,int divX,int divY,double *res){
	double res_sum = 0.;
	for(int i=0;i<Nxs*Nys*divX*divY;i++)
		res_sum += res[i];
	
	return res_sum;
}

__global__ 
void iterate_once(int Nxs,int Nys,int divX,int divY,double *sol,double *sol_old,
									double *res,double beta,int size){
	int index = blockIdx.x*blockDim.x + threadIdx.x;
  int stride = blockDim.x*gridDim.x;
	save_sol_old(Nxs,Nys,sol,sol_old,size); // Save values - convergence checking
	for(int k=index;k<size;k+=stride)
		solver(Nxs,Nys,sol,k,beta);
	
	get_residuals(Nxs,Nys,divX,divY,sol,sol_old,res);
}

void get_jump_points(int divX,int divY,int *jump_points){
	for(int i=0;i<divY-1;i++)
		jump_points[i] = (i+1)*divX-1;
}

void parallel_solver(int Nxs,int Nys,int divX,int divY,double *sol,
										 double *sol_old,double *res,int iter,double beta,
										 double eps,int BpG,int TpB){
	double res_val;
	int *jump_points = new int[divY-1];
	cudaMallocManaged(&jump_points,sizeof(int)*(divY-1));
	get_jump_points(divX,divY,jump_points);
	cudaDeviceSynchronize();
	
	for(int loop=0;loop<iter;loop++){
		// Calculate one iteration for each region
		iterate_once<<<BpG,TpB>>>(Nxs,Nys,divX,divY,sol,sol_old,res,beta,divX*divY);
		cudaDeviceSynchronize();
		
		// Check convergence
		res_val = sum_residuals(Nxs-2,Nys-2,divX,divY,res);
		if(loop > 0){
			if(res_val <= eps){
				puts("Convergence!");
				break;
			}
		}
		
		printf("Iteration %d | Residuals = %f\n",loop+1,res_val);
		
		// Update phantom boundary conditions
		phantom_boundary_conditions(Nxs,Nys,divX,divY,sol,jump_points);
	}
	
	cudaFree(jump_points);
}

void print_results(int Nx,int Ny,int Nxs,int Nys,int divX,int divY,double Lx,
									 double Ly,double *sol){
	FILE *results;
  results = fopen("results.txt","w");
  fprintf(results,"%.8f %.8f ",Lx,Ly);
  for(int i=2;i<Nx;i++)
    fprintf(results,"%.8f ",0.);
  fprintf(results,"\n");
	
	int istart,iend;
	int jstart,jend;
	for(int row=0;row<divY;row++){
		
		if(row==0){
			jstart = Nys-1;
			jend = 1;
		}
		else if(row==divY-1){
			jstart = Nys-2;
			jend = 0;
		}
		else{
			jstart = Nys-2;
			jend = 1;
		}
			
		for(int j=jstart;j>=jend;j--){
			for(int col=0;col<divX;col++){
				if(col==0){
					istart = 0;
					iend = Nxs-2;
				}
				else if(col==divX-1){
					istart = 1;
					iend = Nxs-1;
				}
				else{
					istart = 1;
					iend = Nxs-2;
				}
				for(int i=istart;i<=iend;i++)
					fprintf(results,"%.8f ",sol[(row*divY+col)*Nxs*Nys+Nys*j+i]);
				
				// fprintf(results," | ");
			}
			fprintf(results,"\n");
		}
		// fprintf(results,"-------------------------------\n");
	}
	fclose(results);
}