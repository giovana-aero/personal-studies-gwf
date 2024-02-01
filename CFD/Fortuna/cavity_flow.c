#include<stdio.h>
#include<stdlib.h>
#include<math.h>

/*
Técnicas computacionais para dinâmica dos fluidos - A. de O. Fortuna
Capítulo 6
Escoamento em uma cavidade

Número de Reynolds aqui definido como Re = u0*H/nu, onde H é a altura da
cavidade.
Depois de passar muita raiva procurando erros que não existiam no código,
entendi enfim que esse solucionador quebra com números de Reynolds altos - não
testei muito a fundo realmente, mas se o Reynolds tem mais de três dígitos é
improvável que o solucionador obtenha convergência.

Adicionei também algumas paralelizações com openmp que parecem não ter
eficácia - mas isso também pode ser efeito dos poucos núcleos da CPU que tenho
em mãos no momento. Só não implementei ainda no solucionador PSOR, que é o ponto
mais crítico desse código.

gigi, 01 de fevereiro de 2024, 1941 horas
*/

double **initialize(int Nx,int Ny);
double get_abs_max(double **X,int Nx,int Ny);
double get_deltat(double *deltat,double **umat,double **vmat,int Nx,int Ny,
                  double deltax,double deltay,double nu,double tau);
void save_old_values(double **umat,double **vmat,double **umat_old,
                     double **vmat_old,int Nx,int Ny);
void boundary_conditions(double **umat,double **vmat,int Nx,int Ny,double u0);
double upwind(double ue,double **mat,int i,int j,int op);
double central_differences(double **mat,int i,int j,int op);
double hybrid_scheme(double **mat,double delta,double nu,int i,int j,int op);
void conv_F(double *conv,double **umat,double **vmat,double deltax,
            double deltay,double nu,int i,int j);
void visc_F(double *visc,double **umat,double deltax,double deltay,double nu,
            int i,int j);
void get_F(double **Fmat,double **umat,double **vmat,int Nx,int Ny,
           double deltax,double deltay,double deltat,double nu);
void conv_G(double *conv,double **umat,double **vmat,double deltax,
            double deltay,double nu,int i,int j);
void visc_G(double *visc,double **vmat,double deltax,double deltay,double nu,
            int i,int j);
void get_G(double **Gmat,double **umat,double **vmat,int Nx,int Ny,
           double deltax,double deltay,double deltat,double nu);
void erase_FG_boundaries(double **Fmat,double **Gmat,int Nx,int Ny);
void save_pmat(double **pmat,double **pmat_old,int Nx,int Ny);
void solve_poisson(double **pmat,double **Fmat,double **Gmat,double **pmat_old,
                   int Nx,int Ny,double deltax,double deltay,double deltat,
                   int iterp,double epsp,double rho,double omega,int loop);
void get_uv(double **umat,double **vmat,double **pmat,double **Fmat,
            double **Gmat,int Nx,int Ny,double deltax,double deltay,
            double deltat,double rho);
int check_convergence(double **umat,double **vmat,double **umat_old,
                      double **vmat_old,int Nx,int Ny,double eps,double deltat);
void free_memory(double **X,int Ny);
void unstaggered_umat(double **umat,double **umat_straight,int Nx,int Ny);
void unstaggered_vmat(double **vmat,double **vmat_straight,int Nx,int Ny);
void print_results(double **pmat,double **umat,double **vmat,int Nx,int Ny,
                   double lx,double ly);

int main(){

  // Dados da malha
  double lx = 1.;
  double ly = 1.;
  int Nx = 20; // Número de nós
  int Ny = 20;
  double deltax = lx/((double) Nx);
  double deltay = ly/((double) Ny);

  // Propriedades do fluido
  // Água a 20°C
  // https://www.thermexcel.com/english/tables/eau_atm.htm
  // (Nota posterior: se a intenção for usar água nesse problema em específico,
  // será preciso abaixar drasticamente o número de Reynolds alterando-se os
  // valores da altura da caixa e/ou da velocidade da tampa
  double mu = 10;
  double rho = 1;
  double nu = mu/rho;

  // Configurações da solução geral
  double u0 = 1.; // Velocidade da tampa
  double tau = 1.; // Fator de segurança pro cálculo de deltat
  double deltat;
  int iter = 5000;
  double eps = 1e-6;

  // Configurações do PSOR
  int iterp = 50000;
  double epsp = 1e-8;
  double omega = 1;

  // Matrizes
  double **pmat = initialize(Nx,Ny);
  double **umat = initialize(Nx,Ny);
  double **vmat = initialize(Nx,Ny);
  double **Fmat = initialize(Nx,Ny);
  double **Gmat = initialize(Nx,Ny);
  double **pmat_old = initialize(Nx,Ny);
  double **umat_old = initialize(Nx,Ny);
  double **vmat_old = initialize(Nx,Ny);
  double **umat_straight = initialize(Nx,Ny);
  double **vmat_straight = initialize(Nx,Ny);

  for(int loop=1;loop<=iter;loop++){
    save_old_values(umat,vmat,umat_old,vmat_old,Nx,Ny);

    // Aplicar condições de contorno
    boundary_conditions(umat,vmat,Nx,Ny,u0);
    get_deltat(&deltat,umat,vmat,Nx,Ny,deltax,deltay,nu,tau);

    // Calcular valores de F e G
    get_F(Fmat,umat,vmat,Nx,Ny,deltax,deltay,deltat,nu);
    get_G(Gmat,umat,vmat,Nx,Ny,deltax,deltay,deltat,nu);

    // Zerar os valores de F e G nas fronteiras
    erase_FG_boundaries(Fmat,Gmat,Nx,Ny);

    // Resolver a equação de Poisson
    solve_poisson(pmat,Fmat,Gmat,pmat_old,Nx,Ny,deltax,deltay,deltat,iterp,epsp,
                  rho,omega,loop);

    // Obter novas velocidades
    get_uv(umat,vmat,pmat,Fmat,Gmat,Nx,Ny,deltax,deltay,deltat,rho);

    // printf("Iteration %d\n",loop);

    // Checar convergência
    if(check_convergence(umat,vmat,umat_old,vmat_old,Nx,Ny,eps,deltat)){
      puts("<< Solution convergence >>");
      break;
    }

  }

  // Obter velocidades pra malha não deslocada
  unstaggered_umat(umat,umat_straight,Nx,Ny);
  unstaggered_vmat(vmat,vmat_straight,Nx,Ny);

  // for(int j=0;j<Ny;j++){
  //   for(int i=0;i<Nx;i++){
  //     printf("%.2f ",umat_straight[j][i]);
  //   }
  //   putchar('\n');
  // }
  // putchar('\n');
  // for(int j=0;j<Ny;j++){
  //   for(int i=0;i<Nx;i++){
  //     printf("%.2f ",vmat_straight[j][i]);
  //   }
  //   putchar('\n');
  // }

  // Imprimir resultados
  print_results(pmat,umat_straight,vmat_straight,Nx,Ny,lx,ly);

  // Liberar memória
  free_memory(pmat,Ny);
  free_memory(umat,Ny);
  free_memory(vmat,Ny);
  free_memory(Fmat,Ny);
  free_memory(Gmat,Ny);
  free_memory(umat_old,Ny);
  free_memory(vmat_old,Ny);
  free_memory(umat_straight,Ny);
  free_memory(vmat_straight,Ny);

  return 0;
}

double **initialize(int Nx,int Ny){
  double **X = (double**) malloc(sizeof(double*)*Ny);
  int i, j;
  #pragma omp parallel for shared(X,Nx,Ny) private(i,j)
  for(j=0;j<Ny;j++){
    X[j] = (double*) malloc(sizeof(double)*Nx);
    for(i=0;i<Nx;i++){
      X[j][i] = 0.;
    }
  }

  return X;
}

double get_abs_max(double **X,int Nx,int Ny){
  double answer = 0.;
  for(int j=0;j<Ny;j++){
    for(int i=0;i<Nx;i++){
      if(fabs(X[j][i]) > answer)
        answer = fabs(X[j][i]);
    }
  }

  return answer;
}

double get_deltat(double *deltat,double **umat,double **vmat,int Nx,int Ny,
                  double deltax,double deltay,double nu,double tau){
  double umax = get_abs_max(umat,Nx,Ny);
  double vmax = get_abs_max(vmat,Nx,Ny);

  *deltat = tau/(umax/deltax + vmax/deltay + 2.*nu*(1./(deltax*deltax) +
            1./(deltay*deltay)));
}

void save_old_values(double **umat,double **vmat,double **umat_old,
                     double **vmat_old,int Nx,int Ny){
  int i, j;
  #pragma omp parallel for shared(umat,vmat,umat_old,vmat_old,Nx,Ny)\
          private(i,j)
  for(j=1;j<Ny-1;j++){
    for(i=1;i<Nx-1;i++){
      umat_old[j][i] = umat[j][i];
      vmat_old[j][i] = vmat[j][i];
    }
  }

}

void boundary_conditions(double **umat,double **vmat,int Nx,int Ny,double u0){
  // Paredes inferior e superior
  for(int i=2;i<Nx-1;i++){
    umat[0][i] = -umat[1][i];
    umat[Ny-1][i] = 2*u0 - umat[Ny-2][i];
  }

  // Paredes laterais
  for(int j=2;j<Ny-1;j++){
    vmat[j][0] = -vmat[j][1];
    vmat[j][Nx-1] = -vmat[j][Nx-2];
  }
}

double upwind(double ue,double **mat,int i,int j,int op){
  if(op == 1){ // Direção x
    if(ue >= 0)
      return mat[j][i];
    else
      return mat[j][i+1];
  }
  else{ // Direção y
    if(ue >= 0)
      return mat[j][i];
    else
      return mat[j+1][i];
  }
}

double central_differences(double **mat,int i,int j,int op){
  if(op == 1) // Direção x
    return (mat[j][i] + mat[j][i+1])/2.;
  else // Direção y
    return (mat[j][i] + mat[j+1][i])/2.;
}

double hybrid_scheme(double **mat,double delta,double nu,int i,int j,int op){
  double ue;
  if(op == 1) // Direção x
    ue = (mat[j][i+1] + mat[j][i])/2;
  else
    ue = (mat[j+1][i] + mat[j][i])/2;
  double Pe = ue*delta/nu;
  double FP;

  if(nu<1.9)
    FP = 0.;
  else if(Pe >= 1.9 && Pe < 2)
    FP = (Pe-1.9)/0.1;
  else
    FP = 1.;

  return (1-FP)*central_differences(mat,i,j,op) + FP*upwind(ue,mat,i,j,op);
}

void conv_F(double *conv,double **umat,double **vmat,double deltax,
            double deltay,double nu,int i,int j){
  double ubar1,ubar2,uitp1,uitp2,uitp3,uitp4,vbar1,vbar2;
  // -Velocidade de convecção: definida por média aritmética. Tem a mesma
  //  direção que o delta no denominador da fração
  // -Propriedade transportada: definida por função de interpolação
  ubar1 = (umat[j][i+1] + umat[j][i+2])/2.;       // u_{i+1,j}
  ubar2 = (umat[j][i] + umat[j][i+1])/2.;         // u_{i,j}
  uitp1 = hybrid_scheme(umat,deltax,nu,i+1,j,1);  // u_{i+1,j}
  uitp2 = hybrid_scheme(umat,deltax,nu,i,j,1);    // u_{i,j}
  uitp3 = hybrid_scheme(umat,deltay,nu,i+1,j,2);  // u_{i+1/2,j+1/2}
  uitp4 = hybrid_scheme(umat,deltay,nu,i,j-1,2);  // u_{i+1/2,j-1/2}
  // uitp1 = central_differences(umat,i+1,j,1);  // u_{i+1,j}
  // uitp2 = central_differences(umat,i,j,1);    // u_{i,j}
  // uitp3 = central_differences(umat,i+1,j,2);  // u_{i+1/2,j+1/2}
  // uitp4 = central_differences(umat,i,j-1,2);    // u_{i+1/2,j-1/2}
  vbar1 = (vmat[j+1][i] + vmat[j+1][i+1])/2.;     // v_{i+1/2,j+1/2}
  vbar2 = (vmat[j][i] + vmat[j][i+1])/2.;         // v_{i+1/2,j-1/2}
  *conv = (ubar1*uitp1 - ubar2*uitp2)/deltax +
          (uitp3*vbar1 - uitp4*vbar2)/deltay;
}

// void conv_F(double *conv,double **umat,double **vmat,double deltax,
//             double deltay,double nu,int i,int j){
//   double ubar1,ubar2,uitp1,uitp2,vbar1,vbar2;
//   ubar1 = (umat[j][i+2] + umat[j][i+1])/2;
//   ubar2 = (umat[j][i+1] + umat[j][i])/2;
//   uitp1 = (umat[j][i+1] + umat[j][i])/2;
//   uitp2 = (umat[j][i-1] + umat[j][i])/2;
//   vbar1 = (vmat[j+1][i] + vmat[j+1][i])/2;
//   vbar2 = (vmat[j][i] + vmat[j][i+1])/2;
//   *conv = (ubar1*ubar1 - ubar2*ubar2)/deltax;
//           (uitp1*vbar1 - uitp2*vbar2)/deltay;
//
// }

void visc_F(double *visc,double **umat,double deltax,double deltay,double nu,
            int i,int j){
  *visc = nu*(umat[j][i]-2*umat[j][i+1]+umat[j][i+2])/(deltax*deltax) +
          nu*(umat[j+1][i+1]-2*umat[j][i+1]+umat[j-1][i+1])/(deltay*deltay);
}

void get_F(double **Fmat,double **umat,double **vmat,int Nx,int Ny,
           double deltax,double deltay,double deltat,double nu){
  double conv, visc;
  for(int j=1;j<Ny-1;j++){
    for(int i=1;i<Nx-2;i++){
      visc_F(&visc,umat,deltax,deltay,nu,i,j);
      conv_F(&conv,umat,vmat,deltax,deltay,nu,i,j);
      Fmat[j][i+1] = umat[j][i+1] + deltat*(-conv + visc);
    }
  }
}

void conv_G(double *conv,double **umat,double **vmat,double deltax,
            double deltay,double nu,int i,int j){
  double vbar1,vbar2,vitp1,vitp2,vitp3,vitp4,ubar1,ubar2;
  vbar1 = (vmat[j+1][i] + vmat[j+2][i])/2.;       // v_{i,j+1}
  vbar2 = (vmat[j][i] + vmat[j+1][i])/2.;         // v_{i,j}
  vitp1 = hybrid_scheme(vmat,deltay,nu,i,j+1,2);  // v_{i,j+1}
  vitp2 = hybrid_scheme(vmat,deltay,nu,i,j,2);    // v_{i,j}
  vitp3 = hybrid_scheme(vmat,deltax,nu,i,j+1,1);  // v_{i+1/2,j+1/2}
  vitp4 = hybrid_scheme(vmat,deltax,nu,i-1,j,1);  // v_{i-1/2,j+1/2}
  // vitp1 = central_differences(vmat,i,j+1,2);  // v_{i,j+1}
  // vitp2 = central_differences(vmat,i,j,2);    // v_{i,j}
  // vitp3 = central_differences(vmat,i,j+1,1);  // v_{i+1/2,j+1/2}
  // vitp4 = central_differences(vmat,i-1,j,1);  // v_{i-1/2,j+1/2}
  ubar1 = (umat[j][i+1] + umat[j+1][i+1])/2.;     // v_{i+1/2,j+1/2}
  ubar2 = (umat[j][i] + umat[j+1][i])/2.;         // v_{i-1/2,j+1/2}
  *conv = (vbar1*vitp1 - vbar2*vitp2)/deltay +
          (ubar1*vitp3 - ubar2*vitp4)/deltax;
}

// void conv_G(double *conv,double **umat,double **vmat,double deltax,
//             double deltay,double nu,int i,int j){
//   double vbar1,vbar2,vitp1,vitp2,ubar1,ubar2;
//   vbar1 = (vmat[j+1][i] + vmat[j+2][i])/2.;       // v_{i,j+1}
//   vbar2 = (vmat[j][i] + vmat[j+1][i])/2.;         // v_{i,j}
//   vitp1 = (vmat[j+1][i] + vmat[j+1][i+1])/2.;
//   vitp2 = (vmat[j][i] + vmat[j][i+1])/2.;
//   ubar1 = (umat[j][i+1] + umat[j+1][i+1])/2.;     // v_{i+1/2,j+1/2}
//   ubar2 = (umat[j][i] + umat[j+1][i])/2.;         // v_{i-1/2,j+1/2}
//   *conv = (vbar1*vbar1 - vbar2*vbar2)/deltay +
//           (vitp1*ubar1 - vitp2*ubar2)/deltax;
// }

void visc_G(double *visc,double **vmat,double deltax,double deltay,double nu,
            int i,int j){
  *visc = nu*(vmat[j][i]-2*vmat[j+1][i]+vmat[j+2][i])/(deltay*deltay) +
          nu*(vmat[j+1][i-1]-2*vmat[j+1][i]+vmat[j+1][i+1])/(deltax*deltax);
}

void get_G(double **Gmat,double **umat,double **vmat,int Nx,int Ny,
           double deltax,double deltay,double deltat,double nu){
  double conv, visc;
  for(int j=1;j<Ny-2;j++){
    for(int i=1;i<Nx-1;i++){
      visc_G(&visc,vmat,deltax,deltay,nu,i,j);
      conv_G(&conv,umat,vmat,deltax,deltay,nu,i,j);
      Gmat[j+1][i] = vmat[j+1][i] + deltat*(-conv + visc);
    }
  }
}

void erase_FG_boundaries(double **Fmat,double **Gmat,int Nx,int Ny){
  for(int j=1;j<Ny-1;j++){
    Fmat[j][1] = 0.;
    Fmat[j][Nx-1] = 0.;
  }

  for(int i=1;i<Nx-1;i++){
    Gmat[1][i] = 0.;
    Gmat[Ny-1][i] = 0.;
  }
}

void save_pmat(double **pmat,double **pmat_old,int Nx,int Ny){
  for(int j=0;j<Ny;j++){
    for(int i=0;i<Nx;i++){
      pmat_old[j][i] = pmat[j][i];
    }
  }
}

void solve_poisson(double **pmat,double **Fmat,double **Gmat,double **pmat_old,
                   int Nx,int Ny,double deltax,double deltay,double deltat,
                   int iterp,double epsp,double rho,double omega,int loop){
  // Método PSOR
  double beta = deltax/deltay;
  double fij, p0, Rij, Rsum;
  for(int loopp=1;loopp<=iterp;loopp++){
    printf("%d -- PSOR - Iteration %d\n",loop,loopp);

    // Condições de contorno
    for(int j=1;j<Ny-1;j++){
      pmat[j][0] = pmat[j][1];
      pmat[j][Nx-1] = pmat[j][Nx-2];
    }
    for(int i=1;i<Nx-1;i++){
      pmat[0][i] = pmat[1][i];
      pmat[Ny-1][i] = pmat[Ny-2][i];
    }

    save_pmat(pmat,pmat_old,Nx,Ny);
    for(int j=1;j<Ny-1;j++){
      for(int i=1;i<Nx-1;i++){
        fij = rho/deltat*((Fmat[j][i+1]-Fmat[j][i])/deltax +
              (Gmat[j+1][i]-Gmat[j][i])/deltay);
        pmat[j][i] = (1-omega)*pmat_old[j][i] + omega/(2.*(1.+beta*beta))*
                     (pmat_old[j][i+1] + pmat_old[j][i-1] + beta*beta*pmat_old[j+1][i] +
                     beta*beta*pmat_old[j-1][i] - deltax*deltax*fij);
      }
    }

    // Normalizar
    p0 = pmat[1][1];
    for(int j=1;j<Ny-1;j++){
      for(int i=1;i<Nx-1;i++)
        pmat[j][i] = pmat[j][i] - p0;
    }

    // Checar convergência
    Rsum = 0.;
    for(int j=1;j<Ny-1;j++){
      for(int i=1;i<Nx-1;i++){
        Rij = rho/deltat*((Fmat[j][i+1]-Fmat[j][i])/deltax +
              (Gmat[j+1][i]-Gmat[j][i])/deltay) -
              ((pmat[j][i+1]-2*pmat[j][i]+pmat[j][i-1])/(deltax*deltax) +
              (pmat[j+1][i]-2*pmat[j][i]+pmat[j-1][i])/(deltay*deltay));
        Rsum += Rij*Rij;
      }
    }
    if(sqrt(Rsum) <= epsp){
      puts("-- << PSOR convergence >>");
      break;
    }

    // for(int j=0;j<Ny;j++){
    //   for(int i=0;i<Nx;i++){
    //     printf("%.2f ",pmat[j][i]);
    //   }
    //   putchar('\n');
    // }
    // printf("%f",Rsum);
    // getchar();

  }
}

void get_uv(double **umat,double **vmat,double **pmat,double **Fmat,
            double **Gmat,int Nx,int Ny,double deltax,double deltay,
            double deltat,double rho){
  int i, j;
  #pragma omp parallel for shared(umat,Fmat,deltat,rho,pmat,deltax) private(i,j)
  for(j=1;j<Ny-1;j++){
    for(i=1;i<Nx-2;i++)
      umat[j][i+1] = Fmat[j][i+1] - deltat/rho*(pmat[j][i+1]-pmat[j][i])/deltax;
  }
  #pragma omp parallel for shared(vmat,Gmat,deltat,rho,pmat,deltay) private(i,j)
  for(int j=1;j<Ny-2;j++){
    for(int i=1;i<Nx-1;i++)
      vmat[j+1][i] = Gmat[j+1][i] - deltat/rho*(pmat[j+1][i]-pmat[j][i])/deltay;
  }
}

int check_convergence(double **umat,double **vmat,double **umat_old,
                      double **vmat_old,int Nx,int Ny,double eps,double deltat){
  double Rsum = 0.;
  for(int j=1;j<Ny-1;j++){
    for(int i=1;i<Nx-1;i++){
      Rsum += fabs(umat[j][i]-umat_old[j][i]) + fabs(vmat[j][i]-vmat_old[j][i]);
    }
  }

  if(Rsum/deltat <= eps)
    return 1;
  else
    return 0;
}

void free_memory(double **X,int Ny){
  int j;
  #pragma omp parallel for shared(X,Ny) private(j)
  for(j=0;j<Ny;j++){
    free(X[j]);
  }
}

void unstaggered_umat(double **umat,double **umat_straight,int Nx,int Ny){
  for(int j=0;j<Ny;j++){
    for(int i=0;i<Nx-1;i++)
      umat_straight[j][i] = (umat[j][i] + umat[j][i+1])/2.;
  }
}

void unstaggered_vmat(double **vmat,double **vmat_straight,int Nx,int Ny){
  for(int j=0;j<Ny-1;j++){
    for(int i=0;i<Nx;i++)
      vmat_straight[j][i] = (vmat[j][i] + vmat[j+1][i])/2.;
  }
}

void print_results(double **pmat,double **umat,double **vmat,int Nx,int Ny,
                   double lx,double ly){
  FILE *results1,*results2,*results3;
  results1 = fopen("cavity_flow_p.txt","w");
  results2 = fopen("cavity_flow_u.txt","w");
  results3 = fopen("cavity_flow_v.txt","w");

  fprintf(results1,"%.8f %.8f ",lx,ly);
  for(int i=2;i<Nx;i++)
    fprintf(results1,"%.8f ",0.);
  fprintf(results1,"\n");

  for(int j=0;j<Ny;j++){
    for(int i=0;i<Nx;i++){
      fprintf(results1,"%.8f ",pmat[j][i]);
      fprintf(results2,"%.8f ",umat[j][i]);
      fprintf(results3,"%.8f ",vmat[j][i]);
    }
    fprintf(results1,"\n");
    fprintf(results2,"\n");
    fprintf(results3,"\n");
  }

  fclose(results1);
  fclose(results2);
  fclose(results3);
}
