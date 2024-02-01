#include<stdio.h>
#include<math.h>
#include<omp.h>

void initialize(double *sol,int Nx,int Ny){
  int i;

  #pragma omp parallel for shared(sol,Nx,Ny) private(i)
  for(i=0;i<Nx*Ny;i++){
    sol[i] = 0;
  }
}

void boundary_conditions(double *sol,int Nx,int Ny,double T_up,double T_down,
                         double T_left,double T_right){
  int i, j;

  // Up and Down
  #pragma omp parallel for shared(sol,Nx,Ny,T_down,T_up) private(i)
  for(i=0;i<Nx;i++){
    sol[Ny*0+i] = T_down;
    sol[Ny*(Ny-1)+i] = T_up;
    }


  // Left and Right
  #pragma omp parallel for shared(sol,Nx,Ny,T_left,T_right) private(j)
  for(j=1;j<Ny-1;j++){
    sol[Ny*j+0] = T_left;
    sol[Ny*j+(Nx-1)] = T_right;
  }


}

void print_2d_array(double *sol,int Nx,int Ny){
  for(int j=Ny-1;j>=0;j--){
    for(int i=0;i<Nx;i++)
      printf("%f ",sol[Ny*j+i]);
    putchar('\n');
  }
}

void solver(double *sol,double *sol_old,int Nx,int Ny,int iter,double beta,
            double eps){
  double res_sum;
  int i, j;

  for(int loop=0;loop<iter;loop++){
    // puts("1");

    res_sum = 0.;
    #pragma omp parallel for shared(sol,sol_old,Nx,Ny,beta,res_sum) private(i,j)
    for(int j=1;j<Ny-1;j++){
      for(int i=1;i<Nx-1;i++){
        sol_old[(Ny-2)*(j-1)+(i-1)] = sol[Ny*j+i];
        sol[Ny*j+i] = 1./(2.*(1.+beta*beta))*(sol[Ny*j+i+1]+sol[Ny*j+i-1] +
                      beta*beta*sol[Ny*(j+1)+i]+beta*beta*sol[Ny*(j-1)+i]);
        #pragma omp atomic
        res_sum += fabs(sol[Ny*j+i] - sol_old[(Ny-2)*(j-1)+(i-1)]);
      }
    }


    printf("Iteration %d | Residuals = %f\n",loop+1,res_sum);
    // puts("2");
    // Check convergence
    if(res_sum <= eps){
      puts("Convergence!");
      break;
    }
    // puts("3");
  }
  // puts("4");
}

void print_results(double *sol,int Nx,int Ny,double Lx,double Ly){
  FILE *results;
  results = fopen("results.txt","w");
  fprintf(results,"%.8f %.8f ",Lx,Ly);
  for(int i=2;i<Nx;i++)
    fprintf(results,"%.8f ",0.);
  fprintf(results,"\n");

  for(int j=Ny-1;j>=0;j--){
    for(int i=0;i<Nx;i++){
      fprintf(results,"%.8f ",sol[Ny*j+i]);
    }
    fprintf(results,"\n");
  }
}

int main(){

  // Problem parameters
  double Lx = 1., Ly = 1.;
  double T_down = 0., T_up = 100.;
  double T_left = 50., T_right = 50.;

  // Mesh parameters
  int Nx = 10, Ny = 10; // Complete domain (even values)
  double deltaX = Lx/((double) Nx);
  double deltaY = Ly/((double) Ny);
  double beta = deltaX/deltaY;

  // Numeric solution configuration
  int iter = 100000;
  double eps = 1e-3;
  double sol[Nx*Ny];
  double sol_old[(Nx-2)*(Ny-2)];
  double start, end, elapsed;

  // start = omp_get_wtime();

  // Initialize solution array and apply boundary conditions
  initialize(sol,Nx,Ny);
  initialize(sol_old,Nx-2,Ny-2);
  boundary_conditions(sol,Nx,Ny,T_up,T_down,T_left,T_right);

  // Solve and print results
  solver(sol,sol_old,Nx,Ny,iter,beta,eps);
  print_results(sol,Nx,Ny,Lx,Ly);

  // end = omp_get_wtime();
  // elapsed = end - start;
  // printf("\nTime: %.4f s\n",elapsed);

  return 0;

}
