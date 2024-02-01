#include<stdio.h>
#include<math.h>

/*
Pseudo parallelization test: domain is divided into four regions with equal size

v2 difference: usage of "flattened" 3d arrays, convergence checking
*/

void print_2d_array(int Nx,int Ny,int k,double *M);
void initialize(int size,double *M);
void solver(int Nx,int Ny,double *sol,int k,double beta,double deltaX,
						double deltaY);
void physical_boundary_conditions(int Nxs,int Nys,double *sol,double T_left,
																	double T_right,double T_up,double T_down);
void phantom_boundary_conditions(int Nxs,int Nys,double *sol);
void save_sol_old(int Nxs,int Nys,double *sol,double *sol_old,int k);
double get_residuals(int Nxs,int Nys,double *sol,double *sol_old);
void parallel_solver(int Nxs,int Nys,double *sol,double *sol_old,int iter,
										 double beta,double eps);

int main(){
	// Problem parameters
  double Lx = 1., Ly = 1.;
  double T_down = 10., T_up = 20.;
  double T_left = 30., T_right = 40.;

  // Mesh parameters
  int Nx = 10, Ny = 10; // Complete domain (even values)
  double deltaX = Lx/((double) Nx);
  double deltaY = Ly/((double) Ny);
  double beta = deltaX/deltaY;
	
	// Numeric solution configuration
  int iter = 100000;
  double eps = 1e-3;
	
	// Initialize arrays
	Nx = (int) (Nx/2)*2; // Ensure values are even
	Ny = (int) (Ny/2)*2;
	int Nxs = Nx/2+1; // Dimensions for the submatrices
	int Nys = Ny/2+1;
	// double sol[4*Nxs*Nys];
	double *sol = new double[4*Nxs*Nys];
	initialize(4*Nxs*Nys,sol);
	// double sol_old[4*(Nxs-2)*(Nys-2)];
	double *sol_old = new double[4*(Nxs-2)*(Nys-2)];
	initialize(4*(Nxs-2)*(Nys-2),sol_old);

  // Insert physical boundary conditions
	physical_boundary_conditions(Nxs,Nys,sol,T_left,T_right,T_up,T_down);

  // Obtain solution
  parallel_solver(Nxs,Nys,sol,sol_old,iter,beta,eps);
	
	for(int k=0;k<4;k++){
		printf("%d",k);
		print_2d_array(Nxs,Nys,k,sol);
		putchar('\n');
	}

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

	delete[] sol;
	delete[] sol_old;
	puts("end");

  return 0;
}

// void print_2d_array(int Nx,int Ny,int k,double *M){
  // for(int j=Ny-1;j>=0;j--){
    // for(int i=0;i<Nx;i++)
      // printf("%f ",M[k*Nx*Ny+Ny*j+i]);

    // putchar('\n');
  // }
// }

void print_2d_array(int Nx,int Ny,int k,double *M){
  for(int j=0;j<Ny;j++){
    for(int i=0;i<Nx;i++)
      printf("%f ",M[k*Nx*Ny+Ny*j+i]);

    putchar('\n');
  }
}

void initialize(int size,double *M){
  for(int i=0;i<size;i++)
      M[i] = 0;
}

void solver(int Nx,int Ny,double *sol,int k,double beta){
	for(int j=1;j<Ny-1;j++){
		for(int i=1;i<Nx-1;i++){
			// [k*Nx*Ny+Ny*j+i]
			// printf("(%d %d) ",j,i);
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
void physical_boundary_conditions(int Nxs,int Nys,double *sol,double T_left,
																	double T_right,double T_up,double T_down){
	
	for(int j=1;j<Nys-1;j++){
		sol[0*Nxs*Nys+Nys*j+0] = T_left;
		sol[2*Nxs*Nys+Nys*j+0] = T_left;
		sol[1*Nxs*Nys+Nys*j+(Nxs-1)] = T_right;
		sol[3*Nxs*Nys+Nys*j+(Nxs-1)] = T_right;
  }
  for(int i=1;i<Nxs-1;i++){
		sol[2*Nxs*Nys+Nys*0+i] = T_down;
		sol[3*Nxs*Nys+Nys*0+i] = T_down;
		sol[0*Nxs*Nys+Nys*(Nys-1)+i] = T_up;
		sol[1*Nxs*Nys+Nys*(Nys-1)+i] = T_up;
  }
}


void phantom_boundary_conditions(int Nxs,int Nys,double *sol){
	// Vertical
	for(int j=1;j<Nys;j++){
		// Left to right
		sol[3*Nxs*Nys+Nys*j+0] = sol[2*Nxs*Nys+Nys*j+(Nxs-2)];
		sol[1*Nxs*Nys+Nys*(j-1)+0] = sol[0*Nxs*Nys+Nys*(j-1)+(Nxs-2)];
		
		// Right to left
		sol[2*Nxs*Nys+Nys*j+(Nxs-1)] = sol[3*Nxs*Nys+Nys*j+1];
		sol[0*Nxs*Nys+Nys*(j-1)+(Nxs-1)] = sol[1*Nxs*Nys+Nys*(j-1)+1];
	}
	
	// Horizontal
	for(int i=1;i<Nxs;i++){
		// Lower to upper
		sol[0*Nxs*Nys+Nys*0+i] = sol[2*Nxs*Nys+Nys*(Nys-2)+i];
		sol[1*Nxs*Nys+Nys*0+(i-1)] = sol[3*Nxs*Nys+Nys*(Nys-2)+(i-1)];
		
		// Upper to lower
		sol[2*Nxs*Nys+Nys*(Nys-1)+i] = sol[0*Nxs*Nys+Nys*1+i];
		sol[3*Nxs*Nys+Nys*(Nys-1)+(i-1)] = sol[1*Nxs*Nys+Nys*1+(i-1)];	
	}
}

void save_sol_old(int Nxs,int Nys,double *sol,double *sol_old,int k){
	for(int j=1;j<Nys-1;j++){
		for(int i=1;i<Nxs-1;i++)
			sol_old[k*(Nxs-2)*(Nys-2)+(Nys-2)*(j-1)+(i-1)] = sol[k*Nxs*Nys+Nys*j+i];
	}
}

double get_residuals(int Nxs,int Nys,double *sol,double *sol_old){
	double res_sum = 0;
	for(int k=0;k<4;k++){
		for(int j=1;j<Nys-1;j++){
			for(int i=1;i<Nxs-1;i++)
				res_sum += fabs(sol[k*Nxs*Nys+Nys*j+i] - 
									 sol_old[k*(Nxs-2)*(Nys-2)+(Nys-2)*(j-1)+(i-1)]);
		}
	}
	
	return res_sum;
}

void parallel_solver(int Nxs,int Nys,double *sol,double *sol_old,int iter,
										 double beta,double eps){
	double res_val;
	
	for(int loop=0;loop<iter;loop++){
		// Calculate one iteration for each region
		for(int k=0;k<4;k++){
			save_sol_old(Nxs,Nys,sol,sol_old,k); // Save current values for convergence checking
			solver(Nxs,Nys,sol,k,beta);
			// print_2d_array(Nxs,Nys,k,sol);
			// getchar();
		}
		
		// Check convergence
		if((res_val=get_residuals(Nxs,Nys,sol,sol_old)) <= eps){
			puts("Convergence!");
			break;
		}
		
		printf("Iteration %d | Residuals = %f\n",loop+1,res_val);
		
		// Update phantom boundary conditions
		phantom_boundary_conditions(Nxs,Nys,sol);
	}
}