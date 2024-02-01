#include<stdio.h>
#include<math.h>

/*
Pseudo parallelization test: domain is divided into four regions with equal size

to do: implement convergence checking
*/

void print_2d_array(int Nx,int Ny,double **M){
  for(int j=Ny-1;j>=0;j--){
    for(int i=0;i<Nx;i++)
      printf("%f ",M[j][i]);

    putchar('\n');
  }
}

// void print_2d_array(int Nx,int Ny,double **M){
  // for(int j=0;j<Ny;j++){
    // for(int i=0;i<Nx;i++)
      // printf("%f ",M[j][i]);

    // putchar('\n');
  // }
// }

void initialize(int Nx,int Ny,double **M){
	for(int j=0;j<Ny;j++){
    for(int i=0;i<Nx;i++)
      M[j][i] = 0;
  }
}

void solver(int Nx,int Ny,double **sol,double beta){

	for(int j=1;j<Ny-1;j++){
		for(int i=1;i<Nx-1;i++){
			sol[j][i] = 1./(2.*(1.+beta*beta))*(sol[j][i+1] + sol[j][i-1]
							+ beta*beta*sol[j+1][i] + beta*beta*sol[j-1][i]);
		}
	}
}

// Regions are considered as below:
// | [0] [1] |
// | [2] [3] | 
void physical_boundary_conditions(int Nxs,int Nys,double ***sol,double T_left,
																	double T_right,double T_up,double T_down){
	
	for(int j=1;j<Nys-1;j++){
    sol[0][j][0] = T_left;
		sol[2][j][0] = T_left;
    sol[1][j][Nxs-1] = T_right;
		sol[3][j][Nxs-1] = T_right;
  }
  for(int i=1;i<Nxs-1;i++){
    sol[2][0][i] = T_down;
		sol[3][0][i] = T_down;
    sol[0][Nys-1][i] = T_up;
		sol[1][Nys-1][i] = T_up;
  }
}

void phantom_boundary_conditions(int Nxs,int Nys,double ***sol){
	// Vertical
	for(int j=1;j<Nys;j++){
		// Left to right
		sol[3][j][0] = sol[2][j][Nxs-2];
		sol[1][j-1][0] = sol[0][j-1][Nxs-2];
		
		// Right to left
		sol[2][j][Nxs-1] = sol[3][j][1];
		sol[0][j-1][Nxs-1] = sol[1][j-1][1];
	}
	
	// Horizontal
	for(int i=1;i<Nxs;i++){
		// Lower to upper
		sol[0][0][i] = sol[2][Nys-2][i];
		sol[1][0][i-1] = sol[3][Nys-2][i-1];
		
		// Upper to lower
		sol[2][Nys-1][i] = sol[0][1][i];
		sol[3][Nys-1][i-1] = sol[1][1][i-1];
	}
}

void parallel_solver(int Nxs,int Nys,double ***sol,int iter,double beta,
										 double eps){
	
	for(int loop=0;loop<iter;loop++){
		// Calculate one iteration for each region
		for(int k=0;k<4;k++){
			solver(Nxs,Nys,sol[k],beta);
		}
		
		// Update phantom boundary conditions
		phantom_boundary_conditions(Nxs,Nys,sol);
				
		// Check convergence
		
		printf("Iteration %d\n",loop+1);
	}
}

int main(){
	// Problem parameters
  double Lx = 1., Ly = 1.;
  double T_down = 0., T_up = 1000.;
  double T_left = 0., T_right = 500.;

  // Mesh parameters
  int Nx = 500, Ny = 500; // Complete domain (even values)
  double deltaX = Lx/((double) Nx);
  double deltaY = Ly/((double) Ny);
  double beta = deltaX/deltaY;
	
	// Numeric solution configuration
  int iter = 10000;
  double eps = 1e-3;
	
	// Initialize arrays
	Nx = (int) (Nx/2)*2; // Ensure values are even
	Ny = (int) (Ny/2)*2;
	int Nxs = Nx/2+1; // Dimensions for the submatrices
	int Nys = Ny/2+1;
	double ***sol_ptr = new double**[4];
	for(int k=0;k<4;k++){
		sol_ptr[k] = new double*[Nys];
		for(size_t i=0;i<(Nys);i++)
			sol_ptr[k][i] = new double[Nxs];
		initialize(Nxs,Nys,sol_ptr[k]);
	}
	
  // Insert physical boundary conditions
	physical_boundary_conditions(Nxs,Nys,sol_ptr,T_left,T_right,T_up,T_down);

  // Obtain solution
  parallel_solver(Nxs,Nys,sol_ptr,iter,beta,eps);

  // Print results to file
  FILE *results;
  results = fopen("results.txt","w");
  fprintf(results,"%.8f %.8f ",Lx,Ly);
  for(int i=2;i<Nx;i++)
    fprintf(results,"%.8f ",0.);
  fprintf(results,"\n");
	
	for(int j=Nys-1;j>0;j--){
		for(int i=0;i<Nxs-1;i++)
			fprintf(results,"%.8f ",sol_ptr[0][j][i]);
		for(int i=1;i<Nxs;i++)
			fprintf(results,"%.8f ",sol_ptr[1][j][i]);
		
		fprintf(results,"\n");
	}
	for(int j=Nys-2;j>=0;j--){
		for(int i=0;i<Nxs-1;i++)
			fprintf(results,"%.8f ",sol_ptr[2][j][i]);
		for(int i=1;i<Nxs;i++)
			fprintf(results,"%.8f ",sol_ptr[3][j][i]);
		
		fprintf(results,"\n");
	}
  fclose(results);
	
	for(size_t k=0;k<4;k++){
		for(int i=0;i<(Nys);i++)
			delete sol_ptr[k][i];
	}
	delete sol_ptr;

	puts("end");

  return 0;
}