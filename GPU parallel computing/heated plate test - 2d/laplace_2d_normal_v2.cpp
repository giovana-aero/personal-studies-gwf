#include<stdio.h>
#include<math.h>
// #include<iostream>

/* 
v2 differences:
- all 2d arrays are declared dynamically (new)
- sol_old, used to calculate residuals, had its size reduced and now doesn't 
  include boundary conditions (which are unecessary for that purpose)
*/

/*
Técnicas computacionais para dinâmica dos fluidos - Fortuna - Section 7.6

Reminder: the solution matrix has its lines in reverse order in comparison to
the diagrams in the book, and indexes are such that sol[j][i]
*/

// reminder: this already prints reversing the order of the lines
void print_2d_array(int Nx,int Ny,double **M);
void initialize(int Nx,int Ny,double **M);
void solver(int Nx,int Ny,double **sol,int iter,double beta,double eps);

int main(){

  // Problem parameters
  double Lx = 1., Ly = 1.;
  double T_down = 10., T_up = 10.;
  double T_left = 50., T_right = 50.;

  // Mesh parameters
  int Nx = 50, Ny = 50;
  double deltaX = Lx/((double) Nx);
  double deltaY = Ly/((double) Ny);
  double beta = deltaX/deltaY;
	double **sol = new double*[Ny];
	for(size_t i=0;i<Ny;i++)
		sol[i] = new double[Nx];

  // Numeric solution configuration
  int iter = 5000;
  double eps = 1e-3;
	
  // Insert boundary conditions
	initialize(Nx,Ny,sol);
  for(int j=1;j<Ny-1;j++){
    sol[j][0] = T_left;
    sol[j][Nx-1] = T_right;
  }
  for(int i=1;i<Nx-1;i++){
    sol[0][i] = T_down;
    sol[Ny-1][i] = T_up;
  }

  // Obtain solution
  solver(Nx,Ny,sol,iter,beta,eps);

  //
  // print_2d_array(Nx,Ny,sol_ptr);

  // Print results to file
  FILE *results;
  results = fopen("results.txt","w");
  fprintf(results,"%.8f %.8f ",Lx,Ly);
  for(int i=2;i<Nx;i++)
    fprintf(results,"%.8f ",0.);
  fprintf(results,"\n");
  for(int j=Ny-1;j>=0;j--){
    for(int i=0;i<Nx;i++)
      fprintf(results,"%.8f ",sol[j][i]);
    fprintf(results,"\n");
  }
  fclose(results);
	
	for(size_t i=0;i<Ny;i++)
		delete[] sol[i];
	delete[] sol;

  return 0;
}

// reminder: this already prints reversing the order of the lines
void print_2d_array(int Nx,int Ny,double **M){
  for(int j=Ny-1;j>=0;j--){
		// puts("1");
    for(int i=0;i<Nx;i++)
      printf("%f ",M[j][i]);

    putchar('\n');
  }
}

// void print_2d_array(int Nx,int Ny,double *M){
  // for(int j=0;j<Ny-1;j++){
    // for(int i=0;i<Nx;i++)
      // printf("%f ",M[Nx*j+i]);

    // putchar('\n');
  // }
// }

void initialize(int Nx,int Ny,double **M){
	for(int j=0;j<Ny;j++){
    for(int i=0;i<Nx;i++)
      M[j][i] = 0;
  }
}

void solver(int Nx,int Ny,double **sol,int iter,double beta,double eps){
	// putchar('4');
  double res_sum;
	double **sol_old = new double*[Ny-2];
	for(size_t i=0;i<Ny-2;i++)
		sol_old[i] = new double[Nx-2];
	// putchar('5');
	
	initialize(Nx-2,Ny-2,sol_old);
  for(int loop=0;loop<iter;loop++){

    res_sum = 0.;
    for(int j=1;j<Ny-1;j++){
      for(int i=1;i<Nx-1;i++){
        sol_old[j-1][i-1] = sol[j][i];
        sol[j][i] = 1./(2.*(1.+beta*beta))*(sol[j][i+1] + sol[j][i-1]
                + beta*beta*sol[j+1][i] + beta*beta*sol[j-1][i]);
        res_sum += fabs(sol[j][i] - sol_old[j-1][i-1]);
      }
    }
		
    printf("Iteration %d | Residuals = %f\n",loop+1,res_sum);

    // Check convergence
    if(res_sum <= eps){
      puts("Convergence!");
      break;
    }
  }
	
	for(size_t i=0;i<Ny-2;i++){
		// printf("-> %d\n",i);
		delete[] sol_old[i];
	}
	delete[] sol_old;
}