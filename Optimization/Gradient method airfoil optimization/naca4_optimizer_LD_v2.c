#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<unistd.h>

#define pi 3.141592653589793

void opt_func(double *LD,double *airfoil_data,double *params,int numPoints);
void gradient(double *grad,double *airfoil_data,double *params,int numPoints,
              double delta_m,double delta_p,double delta_t);
void run_xfoil(double *LD,double Re,double alpha,int iter);
void cosspace(double *vec,double length,int numPoints);
void naca4(double *airfoil_data,int numPoints);

int main(){

  // Airfoil data
 double airfoil_data[] = {2.,4.,12.};
  // double airfoil_data[] = {10.630163,5.714319,4.130613};
  int numPoints = 100; // Total
  naca4(airfoil_data,numPoints); system("mv coordinates.dat original_coordinates.dat");

  // Simulation parameters
  double iter = 50;
  double Re = 1e6;
  double alpha = 5.;
  double params[] = {iter,Re,alpha};

  // Aerodynamic data
  double LD;

  // Algorithm parameters
  int opt_iter = 20;
  double eps = 1e-5;
  double delta_m = 1e-1, delta_p = 1e-1, delta_t = 1e-1;
  double step = 5e-3;
  double grad[3];
  double history[opt_iter];
  double grad_history[opt_iter][3];
  FILE *results;

  // Run algorithm
  for(int i=0;i<opt_iter;i++){
    gradient(grad,airfoil_data,params,numPoints,delta_m,delta_p,delta_t);
    airfoil_data[0] += step*grad[0];
    airfoil_data[1] += step*grad[1];
    airfoil_data[2] += step*grad[2];
    opt_func(&LD,airfoil_data,params,numPoints);
    printf("Iteration %d -> f(%f,%f,%f) = %f | %f %f %f\n",
           i+1,airfoil_data[0],airfoil_data[1],airfoil_data[2],
           LD,grad[0],grad[1],grad[2]);

    history[i] = LD;
    grad_history[i][0] = grad[0];
    grad_history[i][1] = grad[1];
    grad_history[i][2] = grad[2];
    if(fabs(grad[0]) < eps && fabs(grad[1]) < eps && fabs(grad[2]) < eps){
      puts("Convergence!");
      break;
    }
  }

  // putchar('\n');
  // for(int i=0;i<opt_iter;i++){
  //   printf("%f ",history[i]);
  // }
  // putchar('\n');
  // putchar('\n');
  // for(int i=0;i<opt_iter;i++){
  //   printf("%f ",grad_history[i][0]);
  // }
  // putchar('\n');
  // putchar('\n');
  // for(int i=0;i<opt_iter;i++){
  //   printf("%f ",grad_history[i][1]);
  // }
  // putchar('\n');
  // putchar('\n');
  // for(int i=0;i<opt_iter;i++){
  //   printf("%f ",grad_history[i][2]);
  // }
  // putchar('\n');

  results = fopen("naca4_results.txt","w");
  for(int i=0;i<opt_iter;i++){
    fprintf(results,"%.6E %.6E %.6E %.6E\n",history[i],
            grad_history[i][0],grad_history[i][1],grad_history[i][2]);
  }
  fclose(results);

  return 0;
}

void opt_func(double *LD,double *airfoil_data,double *params,int numPoints){

  double iter = (int) params[0];
  double Re = params[1];
  double alpha = params[2];

  numPoints = (int) numPoints/2; // Adapt to the logic of the functions
  naca4(airfoil_data,numPoints);
  run_xfoil(LD,Re,alpha,iter);
}

void gradient(double *grad,double *airfoil_data,double *params,int numPoints,
              double delta_m,double delta_p,double delta_t){
  double airfoil1m[3], airfoil2m[3];
  double airfoil1p[3], airfoil2p[3];
  double airfoil1t[3], airfoil2t[3];
  double LD1m, LD2m, LD1p, LD2p, LD1t, LD2t;
  airfoil1m[0] = airfoil_data[0] + delta_m;
  airfoil1m[1] = airfoil_data[1];
  airfoil1m[2] = airfoil_data[2];
  opt_func(&LD1m,airfoil1m,params,numPoints);
  airfoil2m[0] = airfoil_data[0] - delta_m;
  airfoil2m[1] = airfoil_data[1];
  airfoil2m[2] = airfoil_data[2];
  opt_func(&LD2m,airfoil2m,params,numPoints);
  airfoil1p[0] = airfoil_data[0];
  airfoil1p[1] = airfoil_data[1] + delta_p;
  airfoil1p[2] = airfoil_data[2];
  opt_func(&LD1p,airfoil1p,params,numPoints);
  airfoil2p[0] = airfoil_data[0];
  airfoil2p[1] = airfoil_data[1] - delta_p;
  airfoil2p[2] = airfoil_data[2];
  opt_func(&LD2p,airfoil2p,params,numPoints);
  airfoil1t[0] = airfoil_data[0];
  airfoil1t[1] = airfoil_data[1];
  airfoil1t[2] = airfoil_data[2] + delta_t;
  opt_func(&LD1t,airfoil1t,params,numPoints);
  airfoil2t[0] = airfoil_data[0];
  airfoil2t[1] = airfoil_data[1];
  airfoil2t[2] = airfoil_data[2] - delta_t;
  opt_func(&LD2t,airfoil2t,params,numPoints);

  grad[0] = (LD1m-LD2m)/(2*delta_m);
  grad[1] = (LD1p-LD2p)/(2*delta_p);
  grad[2] = (LD1t-LD2t)/(2*delta_t);
}

void run_xfoil(double *LD,double Re,double alpha,int iter){

  FILE *input, *polar;
  int line_count = 0, ch;
  double CL, CD, v[9];
  input = fopen("xfoil_input.txt","w");

  if(access("polar1.txt", F_OK) == 0){
    system("rm polar1.txt");
  }

  fprintf(input,"LOAD coordinates.dat\n\nOPER\n");
  fprintf(input,"VISC %f\n",Re);
  fprintf(input,"PACC\npolar1.txt\n\n");
  fprintf(input,"ITER %d\nALFA %f\n\nQUIT\n",iter,alpha);
  fclose(input);

  system("xfoil < xfoil_input.txt");
  putchar('\n');

  polar = fopen("polar1.txt","r");
  while((ch = fgetc(polar)) != EOF){
    // printf("%c",ch);
    // putchar(ch);
    if(ch == '\n')
      line_count++;

    if(line_count == 12)
      break;
  }

  fscanf(polar,"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&v[0],&v[1],&v[2],&v[3],&v[4],&v[5],&v[6],&v[7],&v[8]);

  // printf("%f %f\n",v[1],v[2]);

  fclose(polar);

  *LD = v[1]/v[2];

}

void cosspace(double *vec,double length,int numPoints){
    // double theta[numPoints];
    double midPoint = length/2.;
    double angleInc = pi/(((double) numPoints) - 1.);
    double curAngle = angleInc;

    for(int i=1;i<numPoints;i++){
        vec[i] = midPoint*(1.-cos(curAngle));
        curAngle += angleInc;
    }
}

void naca4(double *airfoil_data,int numPoints){
  double m = airfoil_data[0]/100.;
  double p = airfoil_data[1]/10.;
  double t = airfoil_data[2]/100.;
  printf("%f %f %f\n",m,p,t);


  double xVec[numPoints];
  cosspace(xVec,1.,numPoints);

  /* Gerar distribuição de espessura */
  /* (O valor de a4 muda pra -0.1036 no caso de bordos de fuga fechados) */
  double a0 = 0.2969;
  double a1 = -0.1260;
  double a2 = -0.3516;
  double a3 = 0.2843;
  double a4 = -0.1015;
  double yVecThickness[numPoints];
  for(int i=0;i<numPoints;i++){
      yVecThickness[i] = 5.*t*(a0*pow(xVec[i],0.5) + a1*xVec[i] + a2*pow(xVec[i],2.) +
                         a3*pow(xVec[i],3.) + a4*pow(xVec[i],4.));
  }

  double xUpper[numPoints];
  double yUpper[numPoints];
  double xLower[numPoints];
  double yLower[numPoints];
  /* Perfil simétrico */
  if(m==0. && p==0.){
      for(int i=0;i<numPoints;i++){
          xUpper[i] = xVec[i];
          yUpper[i] = yVecThickness[i];
          xLower[i] = xVec[i];
          yLower[i] = -yVecThickness[i];
      }
  }
  /* Perfil assimétrico */
  else{
      double yVecCurvature[numPoints];

      for(int i=0;i<numPoints;i++){
          if(xVec[i] >= 0. && xVec[i] <= p)
              yVecCurvature[i] = m/pow(p,2.)*(2.*p*xVec[i] - pow(xVec[i],2.));

          else if(xVec[i] >= p && xVec[i] <= 1)
              yVecCurvature[i] = m/pow(1.-p,2.)*((1. - 2.*p) + 2.*p*xVec[i] - pow(xVec[i],2.));
      }

      double slope[numPoints];
      double theta[numPoints];
      for(int i=0;i<numPoints;i++){
          if(xVec[i] >= 0. && xVec[i] <= p)
              slope[i] = 2.*m/pow(p,2.)*(p - xVec[i]);

          else if(xVec[i] >= p && xVec[i] <= 1.)
              slope[i] = 2.*m/pow(1.-p,2.)*(p - xVec[i]);

          theta[i] = atan(slope[i]);
      }

      /* Agora, obter as coordenadas por definitivo */
      for(int i=0;i<numPoints;i++){
          xUpper[i] = xVec[i] - yVecThickness[i]*sin(theta[i]);
          yUpper[i] = yVecCurvature[i] + yVecThickness[i]*cos(theta[i]);
          xLower[i] = xVec[i] + yVecThickness[i]*sin(theta[i]);
          yLower[i] = yVecCurvature[i] - yVecThickness[i]*cos(theta[i]);
      }
  }
  
  /* Imprimir arquivo de coordenadas */
  FILE *fp;
  fp = fopen("coordinates.dat","w");

  /* Imprimir extradorso */
  for(int i=numPoints-1;i>=0;i--){
      fprintf(fp,"%f %f\n",xUpper[i],yUpper[i]);
  }
  /* Imprimir intradorso */
  for(int i=1;i<numPoints;i++){
      fprintf(fp,"%f %f\n",xLower[i],yLower[i]);
  }

  fclose(fp);
}
