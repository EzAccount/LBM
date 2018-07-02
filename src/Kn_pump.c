#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
//
// D2Q9 in globals
//
// speed of sound, assuming dx/dt=1;
double e[9][2];            // basic
const int x_size = 202;     // points among x
const int y_size = 52;     // points among y
const double R = 8.31/0.4; //
const double gam = 1.66;   //
double dx = 1;
double dt = 1;
double Pr = 1;
float Kn=0.5;
// Viscosity
//
double mu_powerlow(double temperature){              // \mu = A \cdot T ^ {omega}
  double omega = 0.8;
  double A = 1;       // A = mu_rest \cdot T_rest ^{-omega}

  return A*pow(temperature,omega);
}

double mu_other(double temperature){
    return 0.; //Birdman?
}

double omega(double temperature, double rho)
{
  //double nu = mu_powerlow(temperature)/dx/dx*dt; //*dx*dx/dt;
  return 1./Kn; // report930
  // Check Palabos wiki, models, boundary review.
}
double omega_g(double temperature, double rho)
{
  //double nu = mu_powerlow(temperature)/dx/dx*dt; //*dx*dx/dt;
  return  ((1./omega(temperature, rho)-0.5)/Pr + 0.5);
  // Check Palabos wiki, models, boundary review.
}

//
// Vector type functions :
//

double scalar(double *a, double *b) // dot product in R^2
{
  return a[0]*b[0]+a[1]*b[1];
}

double sq_module(double *a) // dot square
{
  return a[0]*a[0]+a[1]*a[1];
}

void vec_sum(double *a, double *b) // sum in R^2
{
  a[0] += b[0];
  a[1] += b[1];
}

double*  vec_mul(double *a, double b) // vector * scalar
{
  a[0] *= b;
  a[1] *= b;
  return a;
}

//
// LBM functions :
//

double macro_rho(double* f_point){  // mass density as sum of destribution
  double result = 0;
  for (int k = 0; k < 9; ++k)
  {
    result += * (f_point+k);
  }
  return result;
}

double macro_temp(double* g_point, double rho)
{
  double result;
  for(int k=0; k<9; ++k){
    result+= *(g_point+k);
  }
  return result/rho/R;
}
void macro_vel(double * u_point,double * f_point, double rho, double temperature) // speed
{
  double tmp1 = 0;
  double tmp2 = 0;
  for(int k=0; k<9; ++k){
    tmp1 += e[k][0] * *(f_point+k)  / rho;
    tmp2 += e[k][1] * *(f_point+k)  / rho;
  }
  *u_point = tmp1;
  *(u_point + 1) = tmp2;
}

void equalibrum(double * f_eq, double * g_eq,double * u, double rho, double T) // f^eq
{
  double control_sum = 0;
  double c = sqrt(R*T);
  double temp[2];
  temp[0] = *u;
  temp[1] = *(u+1);
  double s[9];
  double w[9]={4/9., 1/9.,1/9.,1/9.,1/9.,1/36.,1/36.,1/36.,1/36.};
  for (int k = 0; k < 9; ++k)
  {
    double tmp;
    float sc = scalar(e[k],temp);
    tmp = w[k] * (3 * sc / c + 4.5 * sc * sc / c / c - 1.5 * sq_module(temp) / c / c);
    *(f_eq+k) = w[k] * rho + rho * tmp;
    *(g_eq+k) = R*T* *(f_eq+k);
  }
}

int main(int argc, char **argv)
{
  //
  // Velocity grid:
  //
  e[0][0] = 0.;
  e[0][1] = 0.;
  e[1][0] = 1.;
  e[1][1] = 0.;
  e[2][0] = 0.;
  e[2][1] = 1.;
  e[3][0] = -1.;
  e[3][1] = 0.;
  e[4][0] = 0.;
  e[4][1] = -1.;
  e[5][0] = 1.;
  e[5][1] = 1.;
  e[6][0] = -1.;
  e[6][1] = 1.;
  e[7][0] = -1.;
  e[7][1] = -1.;
  e[8][0] = 1.;
  e[8][1] = -1.;
  //
  //  General constants
  //
  int random;
  int time = 50;// steps in time
  double a = 0; // a = 0 for Maxwell and a = 1 for mirror
  double T1 = 1.0; // Left temperature
  double T2 = 1.5;
  double P1 = 1.0;
  double P2 = 1.0;
  if (argc > 2){
    sscanf (argv[1],"%d",&time);
    sscanf (argv[2],"%f",&Kn);
  }
  else printf("Standart used\n");
  printf("%f\n", Kn);
  FILE * veldat = fopen("tmp/velx.dat", "w");
  FILE * veldat2 = fopen("tmp/vely.dat", "w");
  FILE * rho = fopen("tmp/rho.dat", "w");
  FILE * Tdat = fopen("tmp/T.dat", "w");
  FILE * debug = fopen("tmp/debug.txt", "w");
  FILE * velvecdat = fopen("tmp/vel.dat", "w");
  FILE * pr = fopen("tmp/P.dat", "w");
  //

  // 1-dim arrays:
  int temp = x_size*y_size*9;
  double f[temp];
  double f_eq[temp];
  double f_temp[temp];
  double g[temp];
  double g_eq[temp];
  double g_temp[temp];
  temp = x_size*y_size;
  double vel[temp*2];
  double T[x_size][y_size];
  double P[x_size][y_size];
  double rho_point[x_size][y_size];
  int x1,x2,y1,y2;
  int h = 5;
  x1 = x_size/4;
  x2 = 3*x_size / 4;
  y1 = (y_size/2) - h;
  y2 = (y_size/2) + h;
  fprintf(debug, "time: %d\n", time);
  fprintf(debug, "Kn: %f\n", Kn);
  fprintf(debug, "x1 :%d\nx2: %d\ny1: %d\ny2: %d \n",x1 ,x2,y1,y2);
  //
  // Data:
  //
  //
  // Initials:
  //
  for (int i = 0; i < x_size*y_size*9; ++i){
    f[i] = 0;
    f_eq[i] = 0;
    f_temp[i] = 0;
    g[i] = 0;
    g_eq[i] = 0;
    g_temp[i] = 0;
  }

  for (int i = 0; i < y_size; ++i){
    double * f_point = f;
    double * feq_point = f_eq;
    double * vel_point = vel;
    double * ftemp_point = f_temp;
    double * g_point  = g;
    double * geq_point = g_eq;
    double * gtemp_point = g_temp;
    for (int j = 0; j <= x1; ++j){
        f_point = f + i*x_size*9 + j*9;
        feq_point = f_eq + i*x_size*9 + j*9;
        ftemp_point = f_temp + i*x_size*9 + j*9;
        g_point = g + i*x_size*9 + j*9;
        geq_point = g_eq + i*x_size*9 + j*9;
        gtemp_point = g_temp + i*x_size*9 + j*9;
        vel_point = vel + i*x_size*2 + j*2;
        *(vel_point) = 0;
        *(vel_point+1) = 0;
        rho_point[j][i] = P2 / T1;//
        P[j][i] = P2;
        T[j][i] = T1;
        equalibrum(feq_point, geq_point, vel_point, rho_point[j][i], T[j][i]);
        *ftemp_point = 0.0;
    }
    for (int j = x2; j < x_size; ++j){
        f_point = f + i*x_size*9 + j*9;
        feq_point = f_eq + i*x_size*9 + j*9;
        ftemp_point = f_temp + i*x_size*9 + j*9;
        g_point = g + i*x_size*9 + j*9;
        geq_point = g_eq + i*x_size*9 + j*9;
        gtemp_point = g_temp + i*x_size*9 + j*9;
        vel_point = vel + i*x_size*2 + j*2;
        *(vel_point) = 0;
        *(vel_point+1) = 0;
        rho_point[j][i] = P2 / T2;//
        P[j][i] = P2;
        T[j][i] = T2;
        equalibrum(feq_point, geq_point, vel_point, rho_point[j][i], T[j][i]);
        *ftemp_point = 0.0;
    }
  } // T1, T2 areas

  for (int i = y1; i<= y2; i++){
    for (int j = x1+1; j < x2; ++j){
        double * f_point = f + i*x_size*9 + j*9;
        double * feq_point = f_eq + i*x_size*9 + j*9;
        double * ftemp_point = f_temp + i*x_size*9 + j*9;
        double * g_point = g + i*x_size*9 + j*9;
        double * geq_point = g_eq + i*x_size*9 + j*9;
        double * gtemp_point = g_temp + i*x_size*9 + j*9;
        double * vel_point = vel + i*x_size*2 + j*2;
        *(vel_point) = 0;
        *(vel_point+1) = 0;
        P[j][i] = P2;
        T[j][i] = T1 + (T2-T1)/(x2-x1)*(j-x1+1);

        rho_point[j][i] = P2 / T[j][i];//
        equalibrum(feq_point, geq_point, vel_point, rho_point[j][i], T[j][i]);
        *ftemp_point = 0.0;
    }
  } // bridge

  // for (int i = 1; i < y_size; ++i){
  //   double * f_point = f;
  //   double * feq_point = f_eq;
  //   double * vel_point = vel;
  //   f_point = f + i*x_size*9;
  //   feq_point = f_eq + i*x_size*9;
  //   vel_point = vel + i*x_size*2;
  //   rho_point[0][i] = P1 / T1;
  //   P[0][i] = P1;
  //   T[0][i] = T1;
  //   equalibrum(feq_point, vel, rho_point[0][i],T[0][i]);
  // }

  for (int i = 0; i < y_size; ++i){
    for (int j = 0; j < x_size; ++j){
      for (int k = 0; k < 9; ++k){
          double * f_point = f + i*x_size*9 + j*9;
          double * feq_point = f_eq + i*x_size*9 + j*9;
          double * g_point = g + i*x_size*9 + j*9;
          double * geq_point = g_eq + i*x_size*9 + j*9;
          //  if (*(feq_point+k) < 0 ) fprintf(debug, "Probably mistake at point: %d %d, f_ %d is %f \n", i, j, k, *(feq_point+k));
          //f[i][j][k] = (1 + temp_tau)*f_eq[0][j][k] - tau * f_eq[i+1][j][k]; // Next in speed grid or ?
          *(f_point  + k) = *(feq_point+k);
          *(g_point  + k) = *(geq_point+k);
        }
    }
  }
  fprintf(debug, "\n");
  fprintf(debug, "Initial parametrs");
  fprintf(debug, "\n");
  fprintf(debug, "X-Speed");
  fprintf(debug, "\n");
  fprintf(debug, "\n");
  for (int j = 0; j < y_size; ++j){
    for (int i =0; i < x_size; ++i) fprintf(debug, "%f ", vel[2*j*x_size+2*i]);
    fprintf(debug, "\n");
  }
  fprintf(debug, "\n");
  fprintf(debug, "Y-Speed");
  fprintf(debug, "\n");
  fprintf(debug, "\n");
  for (int j = 0; j < y_size; ++j){
    for (int i =0; i < x_size; ++i) fprintf(debug, "%f ", vel[2*j*x_size+2*i+1]);
    fprintf(debug, "\n");
  }
  fprintf(debug, "\n");
  fprintf(debug, "Temperature");
  fprintf(debug, "\n");
  fprintf(debug, "\n");
  for (int j = 0; j < y_size; ++j){
    for (int i =0; i < x_size; ++i) fprintf(debug, "%f ", T[i][j]);
    fprintf(debug, "\n");
  }
  fprintf(debug, "\n");
  fprintf(debug, "Rho:");
  fprintf(debug, "\n");
  fprintf(debug, "\n");
  for (int j = 0; j < y_size; ++j){
    for (int i =0; i < x_size; ++i) fprintf(debug, "%f ", rho_point[i][j]);
    fprintf(debug, "\n");
  }
  fprintf(debug, "\n");
  fprintf(debug, "\n");


  //
  // Main loop:
  //

  for (int t = 0; t < time; t++){
    double * ftemp_point = f_temp;
    double * f_point = f;
    double * vel_point;
    double * feq_point;
    double * gtemp_point = g_temp;
    double * g_point = g;
    double * geq_point;
    //
    // inside:
    //
    for (int j = 1; j < y_size-1; ++j){
      for (int i = 1; i < x1; ++i){
        ftemp_point = f_temp + 9*j*x_size + 9*i;
        f_point = f + 9*j*x_size + 9*i;
        *(ftemp_point) = *(f_point);
        *(ftemp_point + 9 + 1) = *(f_point+1);
        *(ftemp_point - 9 * x_size+2) = *(f_point+2);
        *(ftemp_point - 9 + 3) = *(f_point+3);
        *(ftemp_point + 9 * x_size + 4) = *(f_point+4);
        *(ftemp_point - 9 * (x_size - 1) + 5) = *(f_point+5);
        *(ftemp_point - 9 * (x_size + 1) + 6) = *(f_point+6);
        *(ftemp_point + 9 * (x_size - 1) + 7) = *(f_point+7);
        *(ftemp_point + 9 * (x_size + 1) + 8) = *(f_point+8);
        gtemp_point = g_temp + 9*j*x_size + 9*i;
        g_point = g + 9*j*x_size + 9*i;
        *(gtemp_point) = *(g_point);
        *(gtemp_point + 9 + 1) = *(g_point+1);
        *(gtemp_point - 9 * x_size+2) = *(g_point+2);
        *(gtemp_point - 9 + 3) = *(g_point+3);
        *(gtemp_point + 9 * x_size + 4) = *(g_point+4);
        *(gtemp_point - 9 * (x_size - 1) + 5) = *(g_point+5);
        *(gtemp_point - 9 * (x_size + 1) + 6) = *(g_point+6);
        *(gtemp_point + 9 * (x_size - 1) + 7) = *(g_point+7);
        *(gtemp_point + 9 * (x_size + 1) + 8) = *(g_point+8);
      }
      for (int i = x2 + 1; i < x_size-1; ++i){
        ftemp_point = f_temp + 9*j*x_size + 9*i;
        f_point = f + 9*j*x_size + 9*i;
        *(ftemp_point) = *(f_point);
        *(ftemp_point + 9 + 1) = *(f_point+1);
        *(ftemp_point - 9 * x_size+2) = *(f_point+2);
        *(ftemp_point - 9 + 3) = *(f_point+3);
        *(ftemp_point + 9 * x_size + 4) = *(f_point+4);
        *(ftemp_point - 9 * (x_size - 1) + 5) = *(f_point+5);
        *(ftemp_point - 9 * (x_size + 1) + 6) = *(f_point+6);
        *(ftemp_point + 9 * (x_size - 1) + 7) = *(f_point+7);
        *(ftemp_point + 9 * (x_size + 1) + 8) = *(f_point+8);
        gtemp_point = g_temp + 9*j*x_size + 9*i;
        g_point = g + 9*j*x_size + 9*i;
        *(gtemp_point) = *(g_point);
        *(gtemp_point + 9 + 1) = *(g_point+1);
        *(gtemp_point - 9 * x_size+2) = *(g_point+2);
        *(gtemp_point - 9 + 3) = *(g_point+3);
        *(gtemp_point + 9 * x_size + 4) = *(g_point+4);
        *(gtemp_point - 9 * (x_size - 1) + 5) = *(g_point+5);
        *(gtemp_point - 9 * (x_size + 1) + 6) = *(g_point+6);
        *(gtemp_point + 9 * (x_size - 1) + 7) = *(g_point+7);
        *(gtemp_point + 9 * (x_size + 1) + 8) = *(g_point+8);
      }
    }

    for (int i = y1+1; i< y2; i++){
      for (int j = x1; j <= x2; ++j){
        ftemp_point = f_temp + 9*i*x_size + 9*j;
        f_point = f + 9*i*x_size + 9*j;
        *(ftemp_point) = *(f_point);
        *(ftemp_point + 9 + 1) = *(f_point+1);
        *(ftemp_point - 9 * x_size+2) = *(f_point+2);
        *(ftemp_point - 9 + 3) = *(f_point+3);
        *(ftemp_point + 9 * x_size + 4) = *(f_point+4);
        *(ftemp_point - 9 * (x_size - 1) + 5) = *(f_point+5);
        *(ftemp_point - 9 * (x_size + 1) + 6) = *(f_point+6);
        *(ftemp_point + 9 * (x_size - 1) + 7) = *(f_point+7);
        *(ftemp_point + 9 * (x_size + 1) + 8) = *(f_point+8);
        gtemp_point = g_temp + 9*i*x_size + 9*j;
        g_point = g +  9*i*x_size + 9*j;
        *(gtemp_point) = *(g_point);
        *(gtemp_point + 9 + 1) = *(g_point+1);
        *(gtemp_point - 9 * x_size+2) = *(g_point+2);
        *(gtemp_point - 9 + 3) = *(g_point+3);
        *(gtemp_point + 9 * x_size + 4) = *(g_point+4);
        *(gtemp_point - 9 * (x_size - 1) + 5) = *(g_point+5);
        *(gtemp_point - 9 * (x_size + 1) + 6) = *(g_point+6);
        *(gtemp_point + 9 * (x_size - 1) + 7) = *(g_point+7);
        *(gtemp_point + 9 * (x_size + 1) + 8) = *(g_point+8);
      }
    }

    //
    // top boundaries
    //
    for (int i=1; i<x1; ++i){
      ftemp_point = f_temp + i*9;
      f_point = f + i*9;
      vel_point = vel + i * 2;
      // transfer
      *(ftemp_point + 9 * x_size + 4) = *(f_point+4);
      *(ftemp_point + 9 * (x_size - 1) + 7) = *(f_point+7);
      *(ftemp_point + 9 * (x_size + 1) + 8) = *(f_point+8);
      gtemp_point = g_temp + i*9;
      g_point = g + i*9;
      vel_point = vel + i * 2;
      // transfer
      *(gtemp_point + 9 * x_size + 4) = *(g_point+4);
      *(gtemp_point + 9 * (x_size - 1) + 7) = *(g_point+7);
      *(gtemp_point + 9 * (x_size + 1) + 8) = *(g_point+8);
    }
    for (int i=x1+1; i<x2; ++i){
      ftemp_point = f_temp + i*9 + 9*x_size*y1;
      f_point = f + i*9 + 9*x_size*y1;
      vel_point = vel + i * 2 + 2*x_size*y1;
      // transfer
      *(ftemp_point + 9 * x_size + 4) = *(f_point+4);
      *(ftemp_point + 9 * (x_size - 1) + 7) = *(f_point+7);
      *(ftemp_point + 9 * (x_size + 1) + 8) = *(f_point+8);
      gtemp_point = g_temp + i*9 + 9*x_size*y1;
      g_point = g + i*9 + 9*x_size*y1;
      vel_point = vel + i * 2 + 2*x_size*y1;
      // transfer
      *(gtemp_point + 9 * x_size + 4) = *(g_point+4);
      *(gtemp_point + 9 * (x_size - 1) + 7) = *(g_point+7);
      *(gtemp_point + 9 * (x_size + 1) + 8) = *(g_point+8);
    }
    for (int i=x2+1; i<x_size-1; ++i){
      ftemp_point = f_temp + i*9;
      f_point = f + i*9;
      vel_point = vel + i * 2;
      // transfer
      *(ftemp_point + 9 * x_size + 4) = *(f_point+4);
      *(ftemp_point + 9 * (x_size - 1) + 7) = *(f_point+7);
      *(ftemp_point + 9 * (x_size + 1) + 8) = *(f_point+8);
      gtemp_point = g_temp + i*9;
      g_point = g + i*9;
      vel_point = vel + i * 2;
      // transfer
      *(gtemp_point + 9 * x_size + 4) = *(g_point+4);
      *(gtemp_point + 9 * (x_size - 1) + 7) = *(g_point+7);
      *(gtemp_point + 9 * (x_size + 1) + 8) = *(g_point+8);
    }

    // Corners:
    // left-top:
    ftemp_point = f_temp;
    f_point = f;
    *(ftemp_point + 9 * (x_size + 1) + 8) = *(f_point+8);
    *(f_point+8) = *(ftemp_point+6);
    gtemp_point = g_temp;
    g_point = g;
    *(gtemp_point + 9 * (x_size + 1) + 8) = *(g_point+8);
    *(g_point+8) = *(gtemp_point+6);
    // top x1:
    ftemp_point = f_temp + x1*9;
    f_point = f + x1*9;
    *(ftemp_point + 9 * (x_size - 1) + 7) = *(f_point+7);
    *(f_point+7) = *(ftemp_point+5);
    gtemp_point = g_temp + x1*9;
    g_point = g + x1*9;
    *(gtemp_point + 9 * (x_size - 1) + 7) = *(g_point+7);
    *(g_point+7) = *(gtemp_point+5);
    // top x2:
    ftemp_point = f_temp + x2*9;
    f_point = f + x2*9;
    *(ftemp_point + 9 * (x_size + 1) + 8) = *(f_point+8);
    *(f_point+8) = *(ftemp_point+6);
    gtemp_point = g_temp + x2*9;
    g_point = g + x2*9;
    *(gtemp_point + 9 * (x_size + 1) + 8) = *(g_point+8);
    *(g_point+8) = *(gtemp_point+6);

    // right-top
    ftemp_point = f_temp + x_size*9 - 9;
    f_point = f + x_size*9 - 9;
    *(ftemp_point + 9 * (x_size - 1) + 7) = *(f_point+7);
    *(f_point + 7) = *(ftemp_point+5);
    gtemp_point = g_temp + x_size*9 - 9;
    g_point = g + x_size*9 - 9;
    *(gtemp_point + 9 * (x_size - 1) + 7) = *(g_point+7);
    *(g_point + 7) = *(gtemp_point+5);





    //
    // bottom
    //
    for (int i=1; i<x1; ++i){
      ftemp_point = f_temp  +(y_size-1)*x_size*9+ i*9;
      f_point = f +(y_size-1)*x_size*9+ i*9;
      // transfer
      *(ftemp_point - 9 * x_size + 2) = *(f_point+2);
      *(ftemp_point - 9 * (x_size - 1) + 5) = *(f_point+5);
      *(ftemp_point - 9 * (x_size + 1) + 6) = *(f_point+6);
      gtemp_point = g_temp  +(y_size-1)*x_size*9+ i*9;
      g_point = g +(y_size-1)*x_size*9+ i*9;
      // transfer
      *(gtemp_point - 9 * x_size + 2) = *(g_point+2);
      *(gtemp_point - 9 * (x_size - 1) + 5) = *(g_point+5);
      *(gtemp_point - 9 * (x_size + 1) + 6) = *(g_point+6);
    }
    for (int i=x1+1; i<x2; ++i){
      ftemp_point = f_temp  + y2*x_size*9+ i*9;
      f_point = f +y2*x_size*9+ i*9;
      // transfer
      *(ftemp_point - 9 * x_size + 2) = *(f_point+2);
      *(ftemp_point - 9 * (x_size - 1) + 5) = *(f_point+5);
      *(ftemp_point - 9 * (x_size + 1) + 6) = *(f_point+6);
      gtemp_point = g_temp  + y2*x_size*9+ i*9;
      g_point = g +y2*x_size*9+ i*9;
      // transfer
      *(gtemp_point - 9 * x_size + 2) = *(g_point+2);
      *(gtemp_point - 9 * (x_size - 1) + 5) = *(g_point+5);
      *(gtemp_point - 9 * (x_size + 1) + 6) = *(g_point+6);
    }
    for (int i=x2+1; i<x_size-1; ++i){
      ftemp_point = f_temp + i*9 + (y_size-1)*(x_size)*9;
      f_point = f + i*9 + (y_size-1)*(x_size)*9;
      // transfer
      *(ftemp_point - 9 * x_size + 2) = *(f_point+2);
      *(ftemp_point - 9 * (x_size - 1) + 5) = *(f_point+5);
      *(ftemp_point - 9 * (x_size + 1) + 6) = *(f_point+6);
      gtemp_point = g_temp + i*9 + (y_size-1)*(x_size)*9;
      g_point = g + i*9 + (y_size-1)*(x_size)*9;
      // transfer
      *(gtemp_point - 9 * x_size + 2) = *(g_point+2);
      *(gtemp_point - 9 * (x_size - 1) + 5) = *(g_point+5);
      *(gtemp_point - 9 * (x_size + 1) + 6) = *(g_point+6);
    }

    //Corners:
    //left-bot:
    ftemp_point = f_temp + (y_size-1)*x_size*9;
    f_point = f + (y_size-1)*x_size*9;
    *(ftemp_point - 9 * (x_size - 1) + 5) = *(f_point+5);
    gtemp_point = g_temp + (y_size-1)*x_size*9;
    g_point = g + (y_size-1)*x_size*9;
    *(gtemp_point - 9 * (x_size - 1) + 5) = *(g_point+5);
    //right-bot:
    ftemp_point = f_temp + (y_size-1)*x_size*9 + x_size*9 - 9;
    f_point = f + (y_size-1)*x_size*9 + x_size*9 - 9;
    *(ftemp_point - 9 * (x_size + 1) + 6) = *(f_point+6);
    gtemp_point = g_temp + (y_size-1)*x_size*9 + x_size*9 - 9;
    g_point = g + (y_size-1)*x_size*9 + x_size*9 - 9;
    *(gtemp_point - 9 * (x_size + 1) + 6) = *(g_point+6);
    // x2 - bot
    ftemp_point = f_temp + (y_size-1)*x_size*9 + 9*x2;
    f_point = f + (y_size-1)*x_size*9 + 9*x2;
    *(ftemp_point - 9 * (x_size - 1) + 5) = *(f_point+5);
    gtemp_point = g_temp + (y_size-1)*x_size*9 + 9*x2;
    g_point = g + (y_size-1)*x_size*9 + 9*x2;
    *(gtemp_point - 9 * (x_size - 1) + 5) = *(g_point+5);
    // x1 - bot
    ftemp_point = f_temp + (y_size-1)*x_size*9 + 9*x1;
    f_point = f + (y_size-1)*x_size*9 + 9*x1;
    *(ftemp_point - 9 * (x_size + 1) + 6) = *(f_point+6);
    gtemp_point = g_temp + (y_size-1)*x_size*9 + 9*x1;
    g_point = g + (y_size-1)*x_size*9 + 9*x1;
    *(gtemp_point - 9 * (x_size + 1) + 6) = *(g_point+6);

    // inside coreners:
    f_point = f + y1*x_size*9 + x1*9;
    ftemp_point = f_temp + y1*x_size*9 + x1*9;
    *(ftemp_point - 9 + 3) = *(f_point+3);
    *(ftemp_point + 9 * x_size + 4) = *(f_point+4);
    *(ftemp_point - 9 * (x_size + 1) + 6) = *(f_point+6);
    *(ftemp_point + 9 * (x_size - 1) + 7) = *(f_point+7);
    *(ftemp_point + 9 * (x_size + 1) + 8) = *(f_point+8);
    g_point = g + y1*x_size*9 + x1*9;
    gtemp_point = g_temp + y1*x_size*9 + x1*9;
    *(gtemp_point - 9 + 3) = *(g_point+3);
    *(gtemp_point + 9 * x_size + 4) = *(g_point+4);
    *(gtemp_point - 9 * (x_size + 1) + 6) = *(g_point+6);
    *(gtemp_point + 9 * (x_size - 1) + 7) = *(g_point+7);
    *(gtemp_point + 9 * (x_size + 1) + 8) = *(g_point+8);

    f_point = f + y2*x_size*9 + x1*9;
    ftemp_point = f_temp + y2*x_size*9 + x1*9;
    *(ftemp_point - 9 * x_size+2) = *(f_point+2);
    *(ftemp_point - 9 + 3) = *(f_point+3);
    *(ftemp_point - 9 * (x_size - 1) + 5) = *(f_point+5);
    *(ftemp_point - 9 * (x_size + 1) + 6) = *(f_point+6);
    *(ftemp_point + 9 * (x_size - 1) + 7) = *(f_point+7);
    g_point = g + y2*x_size*9 + x1*9;
    gtemp_point = g_temp + y2*x_size*9 + x1*9;
    *(gtemp_point - 9 * x_size+2) = *(g_point+2);
    *(gtemp_point - 9 + 3) = *(g_point+3);
    *(gtemp_point - 9 * (x_size - 1) + 5) = *(g_point+5);
    *(gtemp_point - 9 * (x_size + 1) + 6) = *(g_point+6);
    *(gtemp_point + 9 * (x_size - 1) + 7) = *(g_point+7);

    f_point = f + y1*x_size*9 + x2*9;
    ftemp_point = f_temp + y1*x_size*9 + x2*9;
    *(ftemp_point) = *(f_point);
    *(ftemp_point + 9 + 1) = *(f_point+1);
    *(ftemp_point + 9 * x_size + 4) = *(f_point+4);
    *(ftemp_point - 9 * (x_size - 1) + 5) = *(f_point+5);
    *(ftemp_point + 9 * (x_size - 1) + 7) = *(f_point+7);
    *(ftemp_point + 9 * (x_size + 1) + 8) = *(f_point+8);
    f_point = g + y1*x_size*9 + x2*9;
    ftemp_point = g_temp + y1*x_size*9 + x2*9;
    *(ftemp_point) = *(f_point);
    *(ftemp_point + 9 + 1) = *(f_point+1);
    *(ftemp_point + 9 * x_size + 4) = *(f_point+4);
    *(ftemp_point - 9 * (x_size - 1) + 5) = *(f_point+5);
    *(ftemp_point + 9 * (x_size - 1) + 7) = *(f_point+7);
    *(ftemp_point + 9 * (x_size + 1) + 8) = *(f_point+8);

    f_point = f + y2*x_size*9 + x2*9;
    ftemp_point = f_temp + y2*x_size*9 + x2*9;
    *(ftemp_point) = *(f_point);
    *(ftemp_point + 9 + 1) = *(f_point+1);
    *(ftemp_point - 9 * x_size+2) = *(f_point+2);
    *(ftemp_point - 9 * (x_size - 1) + 5) = *(f_point+5);
    *(ftemp_point - 9 * (x_size + 1) + 6) = *(f_point+6);
    *(ftemp_point + 9 * (x_size + 1) + 8) = *(f_point+8);
    f_point = g + y2*x_size*9 + x2*9;
    ftemp_point = g_temp + y2*x_size*9 + x2*9;
    *(ftemp_point) = *(f_point);
    *(ftemp_point + 9 + 1) = *(f_point+1);
    *(ftemp_point - 9 * x_size+2) = *(f_point+2);
    *(ftemp_point - 9 * (x_size - 1) + 5) = *(f_point+5);
    *(ftemp_point - 9 * (x_size + 1) + 6) = *(f_point+6);
    *(ftemp_point + 9 * (x_size + 1) + 8) = *(f_point+8);

    //
    // left and right
    //
    for (int j = 1; j < y_size-1; ++j){
      ftemp_point = f_temp + j*x_size*9;
      f_point = f + j*x_size*9;
      *ftemp_point = *(f_point);
      *(ftemp_point + 9 + 1) = *(f_point + 1);
      *(ftemp_point - 9 * (x_size - 1) + 5) = *(f_point+5);
      *(ftemp_point + 9 * (x_size + 1) + 8) = *(f_point+8);
      ftemp_point = g_temp + j*x_size*9;
      f_point = g + j*x_size*9;
      *ftemp_point = *(f_point);
      *(ftemp_point + 9 + 1) = *(f_point + 1);
      *(ftemp_point - 9 * (x_size - 1) + 5) = *(f_point+5);
      *(ftemp_point + 9 * (x_size + 1) + 8) = *(f_point+8);
      //
      ftemp_point = f_temp + j*x_size*9+ 9*(x_size-1);
      f_point = f + j*x_size*9+ 9*(x_size-1);
      *(ftemp_point - 9 + 3) = *(f_point+3);
      *(ftemp_point - 9 * (x_size + 1) + 6) = *(f_point + 6);
      *(ftemp_point + 9 * (x_size - 1) + 7) = *(f_point + 7);
      ftemp_point = g_temp + j*x_size*9+ 9*(x_size-1);
      f_point = g + j*x_size*9+ 9*(x_size-1);
      *(ftemp_point - 9 + 3) = *(f_point+3);
      *(ftemp_point - 9 * (x_size + 1) + 6) = *(f_point + 6);
      *(ftemp_point + 9 * (x_size - 1) + 7) = *(f_point + 7);
    }

    // x1 - right:
    for (int j = 1; j < y1; ++j){
      ftemp_point = f_temp + x1*9 + j*x_size*9;
      f_point = f + x1*9 + j*x_size*9;
      *(ftemp_point - 9 + 3) = *(f_point + 3);
      *(ftemp_point - 9 * (x_size + 1) + 6) = *(f_point+6);
      *(ftemp_point + 9 * (x_size - 1) + 7) = *(f_point+7);
      ftemp_point = g_temp + x1*9 + j*x_size*9;
      f_point = g + x1*9 + j*x_size*9;
      *(ftemp_point - 9 + 3) = *(f_point + 3);
      *(ftemp_point - 9 * (x_size + 1) + 6) = *(f_point+6);
      *(ftemp_point + 9 * (x_size - 1) + 7) = *(f_point+7);
    }
    for (int j = y2+1; j < y_size-1; ++j){
      ftemp_point = f_temp + x1*9 + j*x_size*9;
      f_point = f + x1*9 + j*x_size*9;
      *(ftemp_point - 9 + 3) = *(f_point + 3);
      *(ftemp_point - 9 * (x_size + 1) + 6) = *(f_point+6);
      *(ftemp_point + 9 * (x_size - 1) + 7) = *(f_point+7);
      ftemp_point = g_temp + x1*9 + j*x_size*9;
      f_point = g + x1*9 + j*x_size*9;
      *(ftemp_point - 9 + 3) = *(f_point + 3);
      *(ftemp_point - 9 * (x_size + 1) + 6) = *(f_point+6);
      *(ftemp_point + 9 * (x_size - 1) + 7) = *(f_point+7);
    }
    // x2 - left:
    for (int j = 1; j < y1; ++j){
      ftemp_point = f_temp + x2*9 + j*x_size*9;
      f_point = f + x2*9 + j*x_size*9;
      *ftemp_point = *(f_point);
      *(ftemp_point + 9 + 1) = *(f_point + 1);
      *(ftemp_point - 9 * (x_size - 1) + 5) = *(f_point+5);
      *(ftemp_point + 9 * (x_size + 1) + 8) = *(f_point+8);
      ftemp_point = g_temp + x2*9 + j*x_size*9;
      f_point = g + x2*9 + j*x_size*9;
      *ftemp_point = *(f_point);
      *(ftemp_point + 9 + 1) = *(f_point + 1);
      *(ftemp_point - 9 * (x_size - 1) + 5) = *(f_point+5);
      *(ftemp_point + 9 * (x_size + 1) + 8) = *(f_point+8);
    }
    for (int j = y2+1; j < y_size-1; ++j){
      ftemp_point = f_temp + x2*9 + j*x_size*9;
      f_point = f + x2*9 + j*x_size*9;
      *ftemp_point = *(f_point);
      *(ftemp_point + 9 + 1) = *(f_point + 1);
      *(ftemp_point - 9 * (x_size - 1) + 5) = *(f_point+5);
      *(ftemp_point + 9 * (x_size + 1) + 8) = *(f_point+8);
      ftemp_point = g_temp + x2*9 + j*x_size*9;
      f_point = g + x2*9 + j*x_size*9;
      *ftemp_point = *(f_point);
      *(ftemp_point + 9 + 1) = *(f_point + 1);
      *(ftemp_point - 9 * (x_size - 1) + 5) = *(f_point+5);
      *(ftemp_point + 9 * (x_size + 1) + 8) = *(f_point+8);
    }


    //
    // mirroring:
    //
    // left-top:
    ftemp_point = f_temp;
    f_point = f;
    feq_point = f_eq;
    *(f_point+8) = *(ftemp_point+6);
    // x2-top:
    ftemp_point = f_temp + x2*9;
    f_point = f + x2*9;
    feq_point = f_eq + x2*9;
    *(f_point+8) =  *(ftemp_point+6);
    // x1-top
    ftemp_point = f_temp + x1*9;
    f_point = f + x1*9;
    feq_point = f_eq + x1*9;
    *(f_point + 7) = *(ftemp_point+5);
    // right-top
    ftemp_point = f_temp + x_size*9 - 9;
    f_point = f + x_size*9 - 9;
    feq_point = f_eq + x_size*9 - 9;
    *(f_point + 7) =  *(ftemp_point+5);
    // left-bot:
    ftemp_point = f_temp + (y_size-1)*x_size*9;
    f_point = f + (y_size-1)*x_size*9;
    feq_point = f_eq + (y_size-1)*x_size*9;
    *(f_point+5) =  *(ftemp_point+7);
    // x2-bot:
    ftemp_point = f_temp + (y_size-1)*x_size*9 + x2*9;
    f_point = f + (y_size-1)*x_size*9 + x2*9;
    feq_point = f_eq + (y_size-1)*x_size*9 + x2*9;
    *(f_point+5) =  *(ftemp_point+7);
    // x1-bot:
    ftemp_point = f_temp + (y_size-1)*x_size*9 + x1*9;
    f_point = f + (y_size-1)*x_size*9 + x1*9;
    feq_point = f_eq +(y_size-1)*x_size*9 + x1*9;
    *(f_point+6) =  *(ftemp_point+8);
    // right-bot:
    ftemp_point = f_temp + (y_size-2)*x_size*9 + x_size*9 - 9;
    f_point = f + (y_size-1)*x_size*9 + x_size*9 - 9;
    feq_point = f_eq+ (y_size-1)*x_size*9 + x_size*9 - 9;
    *(f_point+6) = *(ftemp_point+8);
    //
    // left and right
    //
    double bal = 1;
    for (int j = 1; j < y_size-1; ++j){
      ftemp_point = f_temp + j*x_size*9;
      f_point = f + j*x_size*9;
      feq_point = f_eq + j*x_size*9;
      bal = (*(ftemp_point+6) + *(ftemp_point+3) + *(ftemp_point+7))/(*(feq_point+5)+*(feq_point+1)+*(feq_point+8));
      *(f_point+5) = a* *(ftemp_point+6);
      *(f_point+1) = a* *(ftemp_point+3);
      *(f_point+8) = a* *(ftemp_point+7);
      *(f_point+5) += (1-a)* *(feq_point+5)*bal;
      *(f_point+1) += (1-a)* *(feq_point+1)*bal;
      *(f_point+8) += (1-a)* *(feq_point+8)*bal;

      ftemp_point = f_temp + j*x_size*9+ 9*(x_size-1);
      f_point = f + j*x_size*9+ 9*(x_size-1);
      feq_point = f_eq + j*x_size*9+ 9*(x_size-1);
      bal =(*(ftemp_point+5) + *(ftemp_point+1) + *(ftemp_point+8))/(*(feq_point+6)+*(feq_point+3)+*(feq_point+7));
      *(f_point+6) = a* *(ftemp_point+5);
      *(f_point+3) = a* *(ftemp_point+1);
      *(f_point+7) = a* *(ftemp_point+8);
      *(f_point+6) += (1-a)* *(feq_point+6)*bal;
      *(f_point+3) += (1-a)* *(feq_point+3)*bal;
      *(f_point+7) += (1-a)* *(feq_point+7)*bal;
    }
    // x1 - right:
    for (int j = 1; j <= y1; ++j){
      ftemp_point = f_temp + x1*9 + j*x_size*9;
      f_point = f + x1*9 + j*x_size*9;
      feq_point = f_eq  + x1*9 + j*x_size*9;
      bal = (*(ftemp_point+5) + *(ftemp_point+1) + *(ftemp_point+8))/(*(feq_point+6)+*(feq_point+3)+*(feq_point+7));
      *(f_point+6) = a* *(ftemp_point+5);
      *(f_point+3) = a* *(ftemp_point+1);
      *(f_point+7) = a* *(ftemp_point+8);
      *(f_point+6) += (1-a)* *(feq_point+6)*bal;
      *(f_point+3) += (1-a)* *(feq_point+3)*bal;
      *(f_point+7) += (1-a)* *(feq_point+7)*bal;
    }
    for (int j = y2; j < y_size-1; ++j){
      ftemp_point = f_temp + x1*9 + j*x_size*9;
      f_point = f + x1*9 + j*x_size*9;
      feq_point = f_eq +x1*9 + j*x_size*9;
      bal = (*(ftemp_point+5) + *(ftemp_point+1) + *(ftemp_point+8))/(*(feq_point+6)+*(feq_point+3)+*(feq_point+7));
      *(f_point+6) = a* *(ftemp_point+5);
      *(f_point+3) = a* *(ftemp_point+1);
      *(f_point+7) = a* *(ftemp_point+8);
      *(f_point+6) +=(1-a) * *(feq_point+6)*bal;
      *(f_point+3) += (1-a) * *(feq_point+3)*bal;
      *(f_point+7) += (1-a) * *(feq_point+7)*bal;
    }
    // x2 - left:
    for (int j = 1; j <= y1; ++j){
      ftemp_point = f_temp + x2*9 + j*x_size*9;
      f_point = f + x2*9 + j*x_size*9;
      feq_point = f_eq + x2*9 + j*x_size*9;
      bal = (*(ftemp_point+6) + *(ftemp_point+3) + *(ftemp_point+7))/(*(feq_point+5)+*(feq_point+1)+*(feq_point+8));
      *(f_point+5) = a**(ftemp_point+6);
      *(f_point+1) = a**(ftemp_point+3);
      *(f_point+8) = a**(ftemp_point+7);
      *(f_point+5) += *(feq_point+5)*(1-a)*bal;
      *(f_point+1) += *(feq_point+1)*(1-a)*bal;
      *(f_point+8) += *(feq_point+8)*(1-a)*bal;
    }
    for (int j = y2; j < y_size-1; ++j){
      ftemp_point = f_temp + x2*9 + j*x_size*9;
      f_point = f + x2*9 + j*x_size*9;
      feq_point = f_eq + x2*9 + j*x_size*9;
      bal =  (*(ftemp_point+6) + *(ftemp_point+3) + *(ftemp_point+7))/(*(feq_point+5)+*(feq_point+1)+*(feq_point+8));
      *(f_point+5) = a**(ftemp_point+6);
      *(f_point+1) = a**(ftemp_point+3);
      *(f_point+8) = a**(ftemp_point+7);
      *(f_point+5) += *(feq_point+5)*(1.-a)*bal;
      *(f_point+1) += *(feq_point+1)*(1-a)*bal;
      *(f_point+8) += *(feq_point+8)*(1-a)*bal;
    }


    //
    // top
    //
    for (int i=1; i<x1; ++i){
      ftemp_point = f_temp + i*9;
      f_point = f + i*9;
      feq_point = f_eq + i*9;
      bal = (1-a)*(*(ftemp_point+2) + *(ftemp_point+5) + *(ftemp_point+6))/(*(feq_point+4)+*(feq_point+8)+*(feq_point+7));
      *(f_point+4)=a**(ftemp_point+2);
      *(f_point+8)=a**(ftemp_point+5);
      *(f_point+7)=a**(ftemp_point+6);
      *(f_point+4)+=bal * *(feq_point+4);
      *(f_point+8)+=bal * *(feq_point+8);
      *(f_point+7)+=bal * *(feq_point+7);
    }
    for (int i=x1; i<=x2; ++i){
      ftemp_point = f_temp + i*9 + 9*x_size*y1; //?
      f_point = f + i*9 + 9*x_size*y1;
      feq_point = f_eq  + i*9 + 9*x_size*y1;
      bal = (1-a)*(*(ftemp_point+2) + *(ftemp_point+5) + *(ftemp_point+6))/(*(feq_point+4)+*(feq_point+8)+*(feq_point+7));
      *(f_point+4)=a**(ftemp_point+2);
      *(f_point+8)=a**(ftemp_point+5);
      *(f_point+7)=a**(ftemp_point+6);
      *(f_point+4)+=bal * *(feq_point+4);
      *(f_point+8)+=bal * *(feq_point+8);
      *(f_point+7)+=bal * *(feq_point+7);
    }
    for (int i=x2+1; i<x_size-1; ++i){
      ftemp_point = f_temp + i*9;
      f_point = f + i*9;
      feq_point = f_eq + i*9;
      bal = (1-a)*(*(ftemp_point+2) + *(ftemp_point+5) + *(ftemp_point+6))/(*(feq_point+4)+*(feq_point+8)+*(feq_point+7));
      *(f_point+4)=a**(ftemp_point+2);
      *(f_point+8)=a**(ftemp_point+5);
      *(f_point+7)=a**(ftemp_point+6);
      *(f_point+4)+=bal * *(feq_point+4);
      *(f_point+8)+=bal * *(feq_point+8);
      *(f_point+7)+=bal * *(feq_point+7);
    }


    //
    // bottom
    //
    for (int i=1; i<x1; ++i){
      ftemp_point = f_temp  +(y_size-1)*x_size*9+ i*9;
      f_point = f +(y_size-1)*x_size*9+ i*9;
      feq_point = f_eq +(y_size-1)*x_size*9+ i*9;
      bal = (1-a)*(*(ftemp_point+4) + *(ftemp_point+8) + *(ftemp_point+7))/(*(feq_point+2)+*(feq_point+5)+*(feq_point+6));
      // transfer
      *(f_point+2) = a**(ftemp_point+4);
      *(f_point+5) = a**(ftemp_point+8);
      *(f_point+6) = a**(ftemp_point+7);
      *(f_point+2)+=bal * *(feq_point+2);
      *(f_point+5)+=bal * *(feq_point+5);
      *(f_point+6)+=bal * *(feq_point+6);
    }
    for (int i=x1; i<=x2; ++i){
      ftemp_point = f_temp  + y2*x_size*9+ i*9;
      f_point = f +y2*x_size*9+ i*9;
      feq_point = f_eq +y2*x_size*9+ i*9;
      bal = (1-a)*(*(ftemp_point+4) + *(ftemp_point+8) + *(ftemp_point+7))/(*(feq_point+2)+*(feq_point+5)+*(feq_point+6));
      // transfer
      *(f_point+2) =a* *(ftemp_point+4);
      *(f_point+5) = a**(ftemp_point+8);
      *(f_point+6) =a**(ftemp_point+7);
      *(f_point+2)+=bal * *(feq_point+2);
      *(f_point+5)+=bal * *(feq_point+5);
      *(f_point+6)+=bal * *(feq_point+6);
    }
    for (int i=x2+1; i<x_size-1; ++i){
      ftemp_point = f_temp + i*9 + (y_size-1)*(x_size)*9;
      f_point = f + i*9 + (y_size-1)*(x_size)*9;
      feq_point = f_eq  + i*9 + (y_size-1)*(x_size)*9;
      bal = (1-a)*(*(ftemp_point+4) + *(ftemp_point+8) + *(ftemp_point+7))/(*(feq_point+2)+*(feq_point+5)+*(feq_point+6));
      *(f_point+2) = a**(ftemp_point+4);
      *(f_point+5) = a**(ftemp_point+8);
      *(f_point+6) =a**(ftemp_point+7);
      *(f_point+2)+=bal * *(feq_point+2);
      *(f_point+5)+=bal * *(feq_point+5);
      *(f_point+6)+=bal * *(feq_point+6);
    }
    // Inside corners:

    // x2-y1:
    ftemp_point = f_temp + x2*9 + y1*x_size*9;
    f_point = f + x2*9+ y1*x_size*9;
    *(f_point+8) = *(ftemp_point+6);
    // x1-y1
    ftemp_point = f_temp + x1*9 +y1*9*x_size;
    f_point = f + x1*9 + y1*9*x_size;
    *(f_point + 7) = *(ftemp_point+5);
    // x2-y2:
    ftemp_point = f_temp + y2*x_size*9 + x2*9;
    f_point = f + y2*x_size*9 + x2*9;
    *(f_point+5) = *(ftemp_point+7);
    // x1-y1:
    ftemp_point = f_temp + y1*x_size*9 + x1*9;
    f_point = f + y1*x_size*9 + x1*9;
    *(f_point+6) = *(ftemp_point+8);

    //
    // Core:
    //
    for (int j = 1; j < y_size-1; ++j){
      for (int i = 1; i < x1; ++i){
        ftemp_point = f_temp + j*x_size*9 + i*9;
        f_point = f + j*x_size*9 + i*9;
        feq_point = f_eq + j*x_size*9 + i*9;
        gtemp_point = g_temp + j*x_size*9 + i*9;
        g_point = g + j*x_size*9 + i*9;
        geq_point = g_eq + j*x_size*9 + i*9;
        vel_point = vel + j*x_size*2 + i*2;
        rho_point[i][j] = macro_rho(ftemp_point);
        macro_vel(vel_point, ftemp_point, rho_point[i][j], T[i][j]);
        T[i][j] = macro_temp(gtemp_point, rho_point[i][j]);
        P[i][j] = rho_point[i][j]*T[i][j];
        equalibrum(feq_point, geq_point, vel_point,rho_point[i][j], T[i][j]);
      }
    }
    for (int j = 1; j < y_size-1; ++j){
      for (int i = x2+1; i < x_size-1; ++i){
        ftemp_point = f_temp + j*x_size*9 + i*9;
        f_point = f + j*x_size*9 + i*9;
        feq_point = f_eq + j*x_size*9 + i*9;
        vel_point = vel + j*x_size*2 + i*2;
        rho_point[i][j] = macro_rho(ftemp_point);
        macro_vel(vel_point, ftemp_point, rho_point[i][j], T[i][j]);
        T[i][j] =macro_temp(gtemp_point, rho_point[i][j]);
        P[i][j] = rho_point[i][j]*T[i][j];
        equalibrum(feq_point, geq_point, vel_point,rho_point[i][j], T[i][j]);
      }
    }
    for (int j = y1+1; j < y2; ++j){
      for (int i = x1; i <= x2; ++i){
        ftemp_point = f_temp + j*x_size*9 + i*9;
        f_point = f + j*x_size*9 + i*9;
        feq_point = f_eq + j*x_size*9 + i*9;
        vel_point = vel + j*x_size*2 + i*2;
        rho_point[i][j] = macro_rho(ftemp_point);
        macro_vel(vel_point, ftemp_point, rho_point[i][j], T[i][j]);
        T[i][j] = macro_temp(gtemp_point, rho_point[i][j]);
        P[i][j] = rho_point[i][j]*T[i][j];
        equalibrum(feq_point, geq_point, vel_point,rho_point[i][j], T[i][j]);
      }
    }

    for (int j = 1; j < y_size-1; ++j){
      for (int i = 1; i < x1; ++i){
        for (int k=0; k<9; ++k)
        {
          f[j*x_size*9 + i*9+k] = f_temp[j*x_size*9 + i*9+k] - omega(T[i][j], rho_point[i][j])*(f_temp[j*x_size*9 + i*9 + k] - f_eq[j*x_size*9 + i*9 + k]);
          g[j*x_size*9 + i*9+k] = g_temp[j*x_size*9 + i*9+k] - omega_g(T[i][j], rho_point[i][j])*(g_temp[j*x_size*9 + i*9 + k] - g_eq[j*x_size*9 + i*9 + k]);

        }
      }
    }



    for (int j = 1; j < y_size-1; ++j){
      for (int i = x2+1; i < x_size-1; ++i){
        for (int k=0; k<9; ++k)
        {
          f[j*x_size*9 + i*9+k] = f_temp[j*x_size*9 + i*9+k] - omega(T[i][j], rho_point[i][j])*(f_temp[j*x_size*9 + i*9 + k] - f_eq[j*x_size*9 + i*9 + k]);
          g[j*x_size*9 + i*9+k] = g_temp[j*x_size*9 + i*9+k] - omega_g(T[i][j], rho_point[i][j])*(g_temp[j*x_size*9 + i*9 + k] - g_eq[j*x_size*9 + i*9 + k]);
        }
      }
    }


    for (int j = y1+1; j < y2; ++j){
      for (int i = x1; i <= x2; ++i){
        for (int k=0; k<9; ++k)
        {
          f[j*x_size*9 + i*9+k] = f_temp[j*x_size*9 + i*9+k] - omega(T[i][j], rho_point[i][j])*(f_temp[j*x_size*9 + i*9 + k] - f_eq[j*x_size*9 + i*9 + k]);
          g[j*x_size*9 + i*9+k] = g_temp[j*x_size*9 + i*9+k] - omega_g(T[i][j], rho_point[i][j])*(g_temp[j*x_size*9 + i*9 + k] - g_eq[j*x_size*9 + i*9 + k]);


        }
      }
    }



    for (int j = 1; j < y_size-1; ++j){
      for (int i = 1; i < x1; ++i){
        f_point = f + j*x_size*9 + i*9;
        vel_point = vel + j*x_size*2 + i*2;
        rho_point[i][j] = macro_rho(f_point);
        macro_vel(vel_point, f_point, rho_point[i][j], T[i][j]);
        T[i][j] = macro_temp(gtemp_point, rho_point[i][j]);
        P[i][j] = rho_point[i][j]*T[i][j];
      }
    }

    for (int j = 1; j < y_size-1; ++j){
      for (int i = x2+1; i < x_size-1; ++i){
        f_point = f + j*x_size*9 + i*9;
        vel_point = vel + j*x_size*2 + i*2;
        rho_point[i][j] = macro_rho(f_point);
        macro_vel(vel_point, f_point, rho_point[i][j], T[i][j]);
        T[i][j] = macro_temp(gtemp_point, rho_point[i][j]);
        P[i][j] = rho_point[i][j]*T[i][j];
      }
    }
    for (int j = y1+1; j < y2; ++j){
      for (int i = x1; i <= x2; ++i){
        f_point = f + j*x_size*9 + i*9;
        vel_point = vel + j*x_size*2 + i*2;
        rho_point[i][j] = macro_rho(f_point);
        macro_vel(vel_point, f_point, rho_point[i][j], T[i][j]);
        T[i][j] =  macro_temp(gtemp_point, rho_point[i][j]);
        P[i][j] = rho_point[i][j]*T[i][j];
      }
    }
  }
  fprintf(debug, "Rho\n");
  fprintf(debug, "\n");
  for (int j=0; j<y_size; ++j){
    for (int i=0; i<x_size; ++i){
      fprintf(debug, "%f ", rho_point[i][j]);
    }
    fprintf(debug, "\n");
  }
  fprintf(debug, "\n\n\n");
  fprintf(debug, "Vel\n");
  for (int j=0; j<y_size; ++j){
    for (int i=0; i<x_size; ++i){
      fprintf(debug, "%f ", vel[2*j*x_size+2*i]);
    }
    fprintf(debug, "\n");
  }
  fprintf(debug, "\n\n\n");
  fprintf(debug, "T\n");
  for (int j=0; j<y_size; ++j){
    for (int i=0; i<x_size; ++i){
      fprintf(debug, "%f ", T[i][j]);
    }
    fprintf(debug, "\n");
  }
  double * f_point;
  double * vel_point;
  double * g_point;
  for (int j = 1; j < y_size-1; ++j){
    for (int i = 1; i < x1; ++i){
      f_point = f + j*x_size*9 + i*9;
      g_point = g + j*x_size*9 + i*9;
      vel_point = vel + j*x_size*2 + i*2;
      rho_point[i][j] = macro_rho(f_point);
      macro_vel(vel_point, f_point, rho_point[i][j], T[i][j]);
      T[i][j] =  macro_temp(g_point, rho_point[i][j]);
      P[i][j] = rho_point[i][j]*T[i][j];
    }
  }
  for (int j = 1; j < y_size-1; ++j){
    for (int i = x2+1; i < x_size-1; ++i){
      f_point = f + j*x_size*9 + i*9;
      g_point = g + j*x_size*9 + i*9;
      vel_point = vel + j*x_size*2 + i*2;
      rho_point[i][j] = macro_rho(f_point);
      macro_vel(vel_point, f_point, rho_point[i][j], T[i][j]);
      T[i][j] = macro_temp(g_point, rho_point[i][j]);
      P[i][j] = rho_point[i][j]*T[i][j];
    }
  }
  for (int j = y1+1; j < y2; ++j){
    for (int i = x1; i <= x2; ++i){
      f_point = f + j*x_size*9 + i*9;
      g_point = g + j*x_size*9 + i*9;
      vel_point = vel + j*x_size*2 + i*2;
      rho_point[i][j] = macro_rho(f_point);
      macro_vel(vel_point, f_point, rho_point[i][j], T[i][j]);
      T[i][j] = macro_temp(g_point, rho_point[i][j]);
      P[i][j] = rho_point[i][j]*T[i][j];
    }
  }
  for (int j=1; j<y_size-1; j++)
    for (int i=1; i<x_size-1; i++)
      fprintf(veldat, "%d %d %f \n", i, j, vel[i*2+j*2*x_size]);
  for (int j=1; j<y_size-1; j++)
    for (int i=1; i<x_size-1; i++)
      fprintf(veldat2, "%d %d %f \n", i, j, vel[i*2+j*2*x_size+1]);
  for (int j=1; j<y_size-1; j++)
    for (int i=1; i<x_size-1; i++)
      fprintf(velvecdat, "%d %d %f %f \n", i, j,vel[i*2+j*2*x_size], vel[i*2+j*2*x_size+1]);
  for (int j=1; j<y_size-1; j++)
    for (int i=1; i<x_size-1; i++)
      fprintf(rho, "%d %d %f \n", i, j, rho_point[i][j]);
  for (int j=1; j<y_size-1; j++)
    for (int i=1; i<x_size-1; i++)
      fprintf(pr, "%d %d %f \n", i, j, rho_point[i][j]*T[i][j]);
  for (int j=1; j<y_size-1; j++)
    for (int i=1; i<x_size-1; i++)
      fprintf(Tdat, "%d %d %f \n", i, j, T[i][j]);
    return 0;
}
