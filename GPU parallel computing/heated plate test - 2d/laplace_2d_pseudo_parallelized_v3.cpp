#include<stdio.h>
#include<math.h>
#include<ctime>

/*
Pseudo parallelization test: domain is divided into four regions with equal size

v3 difference: variable domain division

note: if a majority of the boundary conditions are equal to zero, convergence is 
mistakenly achieved immediately. what could be done to prevent this in an 
elegant manner?
*/

void print_2d_array(int Nx,int Ny,int k,double *M);
void initialize(int size,double *M);
void physical_boundary_conditions(int Nxs,int Nys,int divX,int divY,double *sol,
																	double T_left,double T_right,double T_up,
																	double T_down);
void vertical_phantom_BC(int Nxs,int Nys,double *sol,int k);
void horizontal_phantom_BC(int Nxs,int Nys,int divX,double *sol,int k);
int element_in_array(int size,int *array,int num);
void phantom_boundary_conditions(int Nxs,int Nys,int divX,int divY,double *sol,
																 int *jump_points);
void save_sol_old(int Nxs,int Nys,double *sol,double *sol_old,int k);
double get_residuals(int Nxs,int Nys,int divX,int divY,double *sol,
										 double *sol_old);
void get_jump_points(int divX,int divY,int *jump_points);
void parallel_solver(int Nxs,int Nys,int divX,int divY,double *sol,
										 double *sol_old,int iter,double beta,double eps);
void print_results(int Nx,int Ny,int Nxs,int Nys,int divX,int divY,double Lx,
									 double Ly,double *sol);

int main(){
	time_t start, end;
	start = clock();
	
	// Problem parameters
  double Lx = 1., Ly = 1.;
  double T_down = 0., T_up = 0.;
  double T_left = 0., T_right = 40.;

  // Mesh parameters
  int Nx = 102, Ny = 102; // Complete domain (even values)
  double deltaX = Lx/((double) Nx);
  double deltaY = Ly/((double) Ny);
  double beta = deltaX/deltaY;
	
	// Numeric solution configuration
  int iter = 100000;
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
	
	// Initialize arrays
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
	initialize(4*Nxs*Nys,sol);
	double *sol_old = new double[(Nxs-2)*(Nys-2)*divX*divY];
	initialize(4*(Nxs-2)*(Nys-2),sol_old);
	
  // Insert physical boundary conditions
	physical_boundary_conditions(Nxs,Nys,divX,divY,sol,T_left,T_right,T_up,T_down);
	
  // Obtain solution
  parallel_solver(Nxs,Nys,divX,divY,sol,sol_old,iter,beta,eps);

  // Print results to file
	print_results(Nx,Ny,Nxs,Nys,divX,divY,Lx,Ly,sol);
  
	delete[] sol;
	delete[] sol_old;
	puts("end");
	
	end = clock();
	double time_taken = double(end - start)/double(CLOCKS_PER_SEC);
	printf("\n%f s\n",time_taken);

  return 0;
}

void print_2d_array(int Nx,int Ny,int k,double *M){
  for(int j=Ny-1;j>=0;j--){
    for(int i=0;i<Nx;i++)
      printf("%f ",M[k*Nx*Ny+Ny*j+i]);

    putchar('\n');
  }
}

// void print_2d_array(int Nx,int Ny,int k,double *M){
  // for(int j=0;j<Ny;j++){
    // for(int i=0;i<Nx;i++)
      // printf("%f ",M[k*Nx*Ny+Ny*j+i]);

    // putchar('\n');
  // }
// }

void initialize(int size,double *M){
  for(int i=0;i<size;i++)
      M[i] = 0;
}

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
	for(int j=0;j<Nys;j++){
		// Left to right
		sol[(k+1)*Nxs*Nys+Nys*j+0] = sol[k*Nxs*Nys+Nys*j+(Nxs-2)];			
		
		// Right to left
		sol[k*Nxs*Nys+Nys*j+(Nxs-1)] = sol[(k+1)*Nxs*Nys+Nys*j+1];

	}
}

void horizontal_phantom_BC(int Nxs,int Nys,int divX,double *sol,int k){
	for(int i=1;i<Nxs;i++){
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
		
		vertical_phantom_BC(Nxs,Nys,sol,k);
	}
			
	// Horizontal
	for(int k=0;k<divX*(divY-1);k++)
		horizontal_phantom_BC(Nxs,Nys,divX,sol,k);
}

void save_sol_old(int Nxs,int Nys,double *sol,double *sol_old,int k){
	for(int j=1;j<Nys-1;j++){
		for(int i=1;i<Nxs-1;i++)
			sol_old[k*(Nxs-2)*(Nys-2)+(Nys-2)*(j-1)+(i-1)] = sol[k*Nxs*Nys+Nys*j+i];
	}
}

double get_residuals(int Nxs,int Nys,int divX,int divY,double *sol,double *sol_old){
	double res_sum = 0;
	for(int k=0;k<divX*divY;k++){
		for(int j=1;j<Nys-1;j++){
			for(int i=1;i<Nxs-1;i++)
				res_sum += fabs(sol[k*Nxs*Nys+Nys*j+i] - 
									 sol_old[k*(Nxs-2)*(Nys-2)+(Nys-2)*(j-1)+(i-1)]);
		}
	}
	return res_sum;
}

void get_jump_points(int divX,int divY,int *jump_points){
	for(int i=0;i<divY-1;i++)
		jump_points[i] = (i+1)*divX-1;
}

void parallel_solver(int Nxs,int Nys,int divX,int divY,double *sol,
										 double *sol_old,int iter,double beta,double eps){
	double res_val=1.;
	int *jump_points = new int[divY-1];
	get_jump_points(divX,divY,jump_points);
	
	for(int loop=0;loop<iter;loop++){
		// Calculate one iteration for each region
		for(int k=0;k<divX*divY;k++){
			save_sol_old(Nxs,Nys,sol,sol_old,k); // Save current values for convergence checking
			solver(Nxs,Nys,sol,k,beta);
		}		
		
		// Check convergence
		res_val=get_residuals(Nxs,Nys,divX,divY,sol,sol_old)
		if(loop > 1000){
			if(res_val <= eps){
				puts("Convergence!");
				break;
			}
		}
		printf("Iteration %d | Residuals = %f\n",loop+1,res_val);
		
		// Update phantom boundary conditions
		phantom_boundary_conditions(Nxs,Nys,divX,divY,sol,jump_points);
	}
	
	delete[] jump_points;
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