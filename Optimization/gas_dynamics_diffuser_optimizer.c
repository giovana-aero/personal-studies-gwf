#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define pi 3.14159265358979323846

double cot(double beta);
double deg2rad(double deg);
double rad2deg(double rad);
void calc_delta_s(double *delta_s,double gamma,double cp,double R,double M1);
void calc_p_ratio(double *p_ratio,double delta_s,double R);
void calc_mach(double *M2,double M1,double gamma);
void calc_theta(double *theta,double M1,double beta,double gamma);
void calc_p(double *p2,double p1,double gamma,double M1);
void calc_T(double *T2,double T1,double gamma,double M1);
void calc_v(double *v,double M,double gamma,double R,double T);
void calc_rho(double *rho2,double rho1,double gamma,double M1);
void opt_f(double *p0T,double *beta,double *params);
void gradient(double *grad,double *beta,double *params,double deltaB);

int main(){

  // Thermodynamic constants
  double gamma = 1.4;
  double R = 0.287e3;
  double cp = gamma*R/(gamma-1.);

  // Initial conditions - airflow (gigi - mach 2, altitude 18 km)
  double Minf = 2.;
  double Tinf = 216.66;
  double pinf = 7.5652e3;
  double rhoinf = 1.2164e-1;
  double params[] = {gamma,R,cp,Minf};

  // Initial conditions - airflow (ludimila - mach 3, altitude 22 km)
  // double Minf = 3.;
  // double Tinf = 216.66;
  // double pinf = 4.0420e3;
  // double rhoinf = 6.4995e-2;
  // double params[] = {gamma,R,cp,Minf};

  // Angles [degrees]
  double beta[] = {deg2rad(20.),deg2rad(50.)};
  double theta[2];

  // Flow variables
  double M1, M2, M3, MnA, MnB;
  double p1, p2, p3;
  double T1, T2, T3;
  double vinf, v1, v2, v3; calc_v(&vinf,Minf,gamma,R,Tinf);
  double rho1, rho2, rho3;
  double p0i1, p012, p023, p0T;
  double delta_s;

  // Optimizer configuration
  int iter = 10000;
  double eps = 1e-5;
  double step = 1e-3;
  double grad[2];
  double deltaB = 1e-2;

  // Algorithm (steepest descent)
  for(int i=0;i<iter;i++){
    gradient(grad,beta,params,deltaB);
    beta[0] += step*grad[0];
    beta[1] += step*grad[1];

    opt_f(&p0T,beta,params);
    printf("Iteration %d -> f(%f,%f) = %f | {%f,%f}\n",
              i+1,rad2deg(beta[0]),rad2deg(beta[1]),p0T,grad[0],grad[1]);

    if(fabs(grad[0])<eps && fabs(grad[1])<eps){
      puts("Convergence!");
      break;
    }
  }

  // Get data for the configuration found
  // Fist shockwave
  MnA = Minf*sin(beta[0]);
  calc_delta_s(&delta_s,gamma,cp,R,MnA);
  calc_p_ratio(&p0i1,delta_s,R);
  calc_theta(&theta[0],Minf,beta[0],gamma);
  calc_mach(&MnB,MnA,gamma);
  M1 = MnB/(sin(beta[0]-theta[0]));
  calc_p(&p1,pinf,gamma,MnA);
  calc_T(&T1,Tinf,gamma,MnA);
  calc_v(&v1,M1,gamma,R,T1);
  calc_rho(&rho1,rhoinf,gamma,MnA);

  // Second shockwave
  MnA = M1*sin(beta[1]);
  calc_delta_s(&delta_s,gamma,cp,R,MnA);
  calc_p_ratio(&p012,delta_s,R);
  calc_theta(&theta[1],M1,beta[1],gamma);
  calc_mach(&MnB,MnA,gamma);
  M2 = MnB/(sin(beta[1]-theta[1]));
  calc_p(&p2,p1,gamma,MnA);
  calc_T(&T2,T1,gamma,MnA);
  calc_v(&v2,M2,gamma,R,T2);
  calc_rho(&rho2,rho1,gamma,MnA);

  // Third shockwave
  calc_delta_s(&delta_s,gamma,cp,R,M2);
  calc_p_ratio(&p023,delta_s,R);
  calc_mach(&M3,M2,gamma);
  calc_p(&p3,p2,gamma,M2);
  calc_T(&T3,T2,gamma,M2);
  calc_v(&v3,M3,gamma,R,T3);
  calc_rho(&rho3,rho2,gamma,M2);

  // Results
  printf("\nbeta1 = %.6fº | theta1 = %.6f°\n",
         rad2deg(beta[0]),rad2deg(theta[0]));
  printf("beta2 = %.6f° | theta2 = %.6fº\n",
         rad2deg(beta[1]),rad2deg(theta[1]));
  printf("\nStation          | inf          ");
  printf("| 1            | 2            | 3 \n");
  printf("Mach             | %.6E | %.6E | %.6E | %.6E\n",Minf,M1,M2,M3);
  printf("Pressure [Pa]    | %.6E | %.6E | %.6E | %.6E\n",pinf,p1,p2,p3);
  printf("Temperature [K]  | %.6E | %.6E | %.6E | %.6E\n",Tinf,T1,T2,T3);
  printf("Speed [m/s]      | %.6E | %.6E | %.6E | %.6E\n",vinf,v1,v2,v3);
  printf("Density [kg/m^3] | %.6E | %.6E | %.6E | %.6E\n",
         rhoinf,rho1,rho2,rho3);
  printf("\np0 ratios: p0i1         | p012         | p023          | p0T\n");
  printf("           %.6E | %.6E | %.6E  | %.6E\n",p0i1,p012,p023,p0T);

  return 0;
}

double cot(double beta){
  return cos(beta)/sin(beta);
}

double deg2rad(double deg){
  return deg*pi/180.;
}

double rad2deg(double rad){
 return 180./pi*rad;
}

void calc_delta_s(double *delta_s,double gamma,double cp,double R,double M1){
  *delta_s = cp*log( (1.+2.*gamma/(gamma+1.)*(M1*M1-1.))*(2.+(gamma-1.)*M1*M1)/
             ((gamma+1.)*M1*M1) ) - R*log( 1.+2.*gamma/(gamma+1.)*(M1*M1-1.) );
}

void calc_p_ratio(double *p_ratio,double delta_s,double R){
  *p_ratio = exp(-delta_s/R);
}

void calc_mach(double *M2,double M1,double gamma){
  *M2 = sqrt( (1.+(gamma-1.)/2.*M1*M1)/(gamma*M1*M1-(gamma-1.)/2.) );
}

void calc_theta(double *theta,double M1,double beta,double gamma){
  *theta = atan( 2.*cot(beta)*(M1*M1*pow(sin(beta),2.)-1.)/
           (M1*M1*(gamma+cos(2.*beta))+2.) );
}

void calc_p(double *p2,double p1,double gamma,double M1){
  *p2 = p1*( 1. + 2.*gamma/(gamma+1.)*(M1*M1-1.) );
}

void calc_T(double *T2,double T1,double gamma,double M1){
  *T2 = T1*( (1.+2.*gamma/(gamma+1.)*(M1*M1-1.))*
             (2.+(gamma-1.)*M1*M1)/((gamma+1.)*M1*M1) );
}

void calc_v(double *v,double M,double gamma,double R,double T){
  *v = M*sqrt(gamma*R*T);
}

void calc_rho(double *rho2,double rho1,double gamma,double M1){
  *rho2 = rho1*( (gamma+1.)*M1*M1/(2.+(gamma-1.)*M1*M1) );
}

void opt_f(double *p0T,double *beta,double *params){
  // Angles
  double theta[2];
  double beta1 = *beta; beta++;
  double beta2 = *beta;

  // Airflow data
  // params = [gamma,R,cp,Minf]
  double MnA, MnB;
  double delta_s;
  double M1, M2;
  double p0i1, p012, p023; // Total pressure ratios
  double gamma = *params; params++;
  double R = *params; params++;
  double cp = *params; params++;
  double Minf = *params; params++;

  // Fist shockwave
  MnA = Minf*sin(beta1);
  calc_delta_s(&delta_s,gamma,cp,R,MnA);
  calc_p_ratio(&p0i1,delta_s,R);
  calc_theta(&theta[0],Minf,beta1,gamma);
  calc_mach(&MnB,MnA,gamma);
  M1 = MnB/(sin(beta1-theta[0]));
  if(M1<1){
    puts("Error: M1 < 1. Choose different initial conditions.");
    exit(1);
  }

  // Second shockwave
  MnA = M1*sin(beta2);
  calc_delta_s(&delta_s,gamma,cp,R,MnA);
  calc_p_ratio(&p012,delta_s,R);
  calc_theta(&theta[1],M1,beta2,gamma);
  calc_mach(&MnB,MnA,gamma);
  M2 = MnB/(sin(beta2-theta[1]));
  if(M2<1){
    puts("Error: M2 < 1. Choose different initial conditions.");
    exit(2);
  }

  // Third shockwave
  calc_delta_s(&delta_s,gamma,cp,R,M2);
  calc_p_ratio(&p023,delta_s,R);

  // Total pressure ratio
  *p0T = p0i1*p012*p023;
}

void gradient(double *grad,double *beta,double *params,double deltaB){
  double beta1 = *beta; beta++;
  double beta2 = *beta;
  double tmpB1[2], tmpB2[2];
  double p0T1, p0T2;

  // Calculate grad[0]
  tmpB1[0] = beta1 + deltaB;
  tmpB1[1] = beta2;
  tmpB2[0] = beta1 - deltaB;
  tmpB2[1] = beta2;
  opt_f(&p0T1,tmpB1,params);
  opt_f(&p0T2,tmpB2,params);
  *grad = (p0T1-p0T2)/(2.*deltaB);

  // Calculate grad[1]
  grad++;
  tmpB1[0] = beta1;
  tmpB1[1] = beta2 + deltaB;
  tmpB2[0] = beta1;
  tmpB2[1] = beta2 - deltaB;
  opt_f(&p0T1,tmpB1,params);
  opt_f(&p0T2,tmpB2,params);
  *grad = (p0T1-p0T2)/(2.*deltaB);
}
