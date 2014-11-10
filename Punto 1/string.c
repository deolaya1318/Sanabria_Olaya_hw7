#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// This code gives answer to the hyperbolic wave equation

void copy(float *initial, float *final, int n_points);

int main(int argc, char **argv) {

  //First off, initiate every constant, variable, pointer and FILE. The value of rho is given when executing the string.x. The value of rho for this definition of dx and dt cannot be smaller than 1E-5, otherwise the value of c will be bigger than 1, and then this will not give a solution to the given problem.

  float L = 100;
  float T = 40;
  float rho;
  rho = atof(argv[1]);
  float c = sqrt(rho/T);
  FILE *file = fopen("string_rho.dat","w");

  int mult = 5;
  int n_points = 101;
  int n_time = 119 * mult;

  float dx = (float) (L*rho)/(n_points-1);
  float dt = 0.05;
  float r = c*dt/dx;

  float *xt = malloc(sizeof(float) * n_points);
  float *u_initial = malloc(sizeof(float) * n_points);
  float *u_past = malloc(sizeof(float) * n_points);
  float *u_present = malloc(sizeof(float) * n_points);
  float *u_future = malloc(sizeof(float) * n_points);

  int i, j;
  float x, u;
  float future_value;
  
  //Now we initialize the linear space of x from L = 0 to L = 100 and the initial condition of u and its border condition.
  
  for (i=0; i<n_points;i++){
    
    x = (float) i*dx/rho;
    
    if(x <= 0.8*L){
      u = 1.25*x / L;
    }
    else{
      u = 5 - 5*x / L;
    }
    xt[i]=x;
    u_initial[i]=u;
    u_future[i]=0;
    fprintf(file, "%f\t",x);
  }
  fprintf(file, "\n");
  
  u_initial[0] = 0.0;
  u_initial[n_points-1] = 0.0;
  
  //Then we develop the position response of the string at time = 0
  fprintf(file, "%f\t", 0.0);

  for (i=1; i<n_points-1;i++){
    future_value = u_initial[i] + (pow(r,2)/2.0) * (u_initial[i+1] - 2.0 * u_initial[i] + u_initial[i-1]);
    u_future[i] = future_value;
    fprintf(file, "%f\t", future_value);
  }
  fprintf(file, "%f\t", 0.0);
  fprintf(file, "\n");

  copy(u_initial, u_past, n_points);
  copy(u_future, u_present, n_points);

  //Finally we develop the position response of the string from t = 1 to t=120. 
  
  for (j=0; j<n_time;j++){
		
    if(j%mult == 0){
      fprintf(file, "%f\t", 0.0);
    }
    
    for (i=1; i<n_points-1;i++){
      future_value = (2.0 * (1.0-pow(r,2))) * u_present[i] - u_past[i] + pow(r,2)*(u_present[i+1] +  u_present[i-1]);
      u_future[i] = future_value;
      if(j%mult == 0){
	fprintf(file, "%f\t", future_value);
      }
    }
    
    if(j%mult == 0){
      fprintf(file, "%f\t", 0.0);
      fprintf(file, "\n");
    }
    copy(u_present, u_past, n_points);
    copy(u_future, u_present, n_points);
  }

  fclose(file);

  return 0;
}

//It is neccesary to copy the initial pointer repeatedly through the cycle of time, so the function copy is defined

void copy(float *initial, float *final, int n_points)
{
  int i;
  for (i=0;i<n_points;i++){
    final[i]=initial[i];
  }
}
