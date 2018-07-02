#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define NUM_COMMANDS 2
//
// D2Q9 in globals
//
 // speed of sound, assuming dx/dt=1;
double e[9][2];         // basic
const int x_size = 12;   // points among x
const int y_size = 102; // points among y
const double R = 8.31/0.4;
const double gam = 1.66;
double dx = 1;
double dt = 1;
//
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
  return  1.6; // report930
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

double macro_temp(double* f_point, double* u_point, double rho)
{
  double result = 0;
  double tmp[2];
  double sq;
  for(int k=0; k<9; ++k){
    tmp[0] = e[k][0] - *(u_point);
    tmp[1] = e[k][1] - *(u_point+1);
    sq = sq_module(tmp);
    result+= sq * *(f_point+k);
  }
  return result/rho*3/2;
}
void macro_vel(double * u_point,double * f_point, double rho, double temperature) // speed
{
  int kappa;
  double tmp1 = 0;
  double tmp2 = 0;
  for(int k=0; k<9; ++k){
    tmp1 += e[k][0] * *(f_point+k)  / rho;
    tmp2 += e[k][1] * *(f_point+k)  / rho;
  }
  *u_point = tmp1;
  *(u_point + 1) = tmp2;
}


void equalibrum(double * f_eq, double * u, double rho, double T) // f^eq
{
  double control_sum = 0;
  double c = sqrt(R);
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
  }

}

int main()
{
  //
  // Velocity grid:
  //
  e[0][0] = 0.;
  e[0][1] = 0.;
  e[1][0] = 1.;;
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
  ////
  int random;
  int time = 8000; // steps in time



  double Tw = 1.0; // Temperature on boundaries *300 K

  double P1 = 1.0;
  double P2 = 1.0;
  FILE * veldat = fopen("data/vel.dat", "w");
  FILE * veldat2 = fopen("data/vel2.dat", "w");
  FILE * rho = fopen("data/rho.dat", "w");
  FILE * Tdat = fopen("data/T.dat", "w");
  FILE * debug = fopen("data/debug.txt", "w");
  //
  // 1-dim arrays:
  int temp = x_size*y_size*9;
  double f[temp];
  double f_eq[temp];
  double f_temp[temp];
  temp = x_size*y_size;
  double vel[temp*2];
  double T[x_size][y_size];
  double P[x_size][y_size];
  double rho_point[x_size][y_size];
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
  }
  for (int i = 0; i < y_size; ++i){
    double * f_point = f;
    double * feq_point = f_eq;
    double * vel_point = vel;
    double * ftemp_point = f_temp;
    for (int j = 0; j < x_size; ++j){
        f_point = f + i*x_size*9 + j*9;
        feq_point = f_eq + i*x_size*9 + j*9;
        ftemp_point = f_temp + i*x_size*9 + j*9;
        vel_point = vel + i*x_size*2 + j*2;
        rho_point[j][i] = P2 / Tw;//
        P[j][i] = P2;
        T[j][i] = Tw;
        equalibrum(feq_point, vel_point, rho_point[j][i], T[j][i]);
        *ftemp_point = 0.0;
    }
  }
  for (int i = 0; i < y_size; ++i){
    double * f_point = f;
    double * feq_point = f_eq;
    double * vel_point = vel;
    f_point = f + i*x_size*9;
    feq_point = f_eq + i*x_size*9;
    vel_point = vel + i*x_size*2;
    rho_point[0][i] = P1 / Tw;
    P[0][i] = P1;
    T[0][i] = Tw;
    equalibrum(feq_point, vel, rho_point[0][i],T[0][i]);
  }

  for (int i = 0; i < y_size; ++i){
    for (int j = 0; j < x_size; ++j){
      for (int k = 0; k < 9; ++k){
          double * f_point = f + i*x_size*9 + j*9;
          double * feq_point = f_eq + i*x_size*9 + j*9;
          //  if (*(feq_point+k) < 0 ) fprintf(debug, "Probably mistake at point: %d %d, f_ %d is %f \n", i, j, k, *(feq_point+k));
          //f[i][j][k] = (1 + temp_tau)*f_eq[0][j][k] - tau * f_eq[i+1][j][k]; // Next in speed grid or ?
          *(f_point  + k) = *(feq_point+k);
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
  /*
  for (int i = 0; i < x_size; ++i){
    for (int j = 0; j < y_size; ++j)
    {
      fprintf(debug, "%d %d\n", i,j);
      for (int k = 0; k < 9; ++k) fprintf(debug, "%f ", f[j*x_size*9+i*9+k]);
      fprintf(debug, "\n");
    }
    fprintf(debug, "\n");
  }*/
  for (int t = 0; t < time; ++t){
    double * ftemp_point = f_temp;
    double * f_point = f;
    double * vel_point;
    double * feq_point;

    if (t%400 == 0) fprintf(debug, "Done: %f \r", t/20.);
    //
    // inside:
    //
    for (int j = 1; j < y_size-1; ++j){
      for (int i = 1; i < x_size-1; ++i){
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
      }
    }

    for (int i=0; i<y_size; ++i){
      f_point =  f + i*x_size*9;
      vel_point = vel + i*x_size*2;
      *vel_point = *(vel_point+2);
      *(vel_point+1) = *(vel_point+3);
      equalibrum(f_point, vel_point, rho_point[0][i], T[0][i]);
    }

    //
    // top boundary
    //
    for (int i=1; i<x_size-1; ++i){

      ftemp_point = f_temp + i*9;
      f_point = f + i*9;
      vel_point = vel + i * 2;
      P[i][0] = P[i][1];
      rho_point[i][0] = P[i][0] / Tw;
      //*vel_point = *(vel_point + x_size*2);
      //*(vel_point+1) = *(vel_point + 1 + x_size*2);
      equalibrum(ftemp_point, vel_point, rho_point[i][0], Tw);

      *ftemp_point = *(f_point);
      *(ftemp_point + 9 + 1) = *(f_point+1);
      *(ftemp_point - 9 + 3) = *(f_point+3);
      *(ftemp_point + 9 * x_size + 4) = *(f_point+4);
      *(ftemp_point + 9 * (x_size - 1) + 7) = *(f_point+7);
      *(ftemp_point + 9 * (x_size + 1) + 8) = *(f_point+8);

    }
    //Corners:
    //left-top:
    // ftemp_point = f;
    // f_point = f;
    // *ftemp_point = *(f_point);
    // *(ftemp_point + 1) = *(f_point+1);
    // *(ftemp_point + 2) = *(f_point+2);
    // *(ftemp_point + 3) = *(f_point+3);
    // *(ftemp_point + 4) = *(f_point+4);
    // *(ftemp_point + 5) = *(f_point+5);
    // *(ftemp_point + 6) = *(f_point+6);
    // *(ftemp_point + 7) = *(f_point+7);
    // *(ftemp_point + 8) = *(f_point+8);
    ftemp_point = f_temp;
    f_point = f;
    *ftemp_point = *(f_point);
    *(ftemp_point + 9 + 1) = *(f_point+1);
    *(ftemp_point + 9 * x_size + 4) = *(f_point+4);
    *(ftemp_point + 9 * (x_size + 1) + 8) = *(f_point+8);
    //right-top:
    // ftemp_point = f + x_size*9 - 9;
    // f_point = f + x_size*9 - 9;
    // *ftemp_point = *(f_point);
    // *(ftemp_point + 1) = *(f_point+1);
    // *(ftemp_point + 2) = *(f_point+2);
    // *(ftemp_point + 3) = *(f_point+3);
    // *(ftemp_point + 4) = *(f_point+4);
    // *(ftemp_point + 5) = *(f_point+5);
    // *(ftemp_point + 6) = *(f_point+6);
    // *(ftemp_point + 7) = *(f_point+7);
    // *(ftemp_point + 8) = *(f_point+8);
    ftemp_point = f_temp + x_size*9 - 9;
    f_point = f + x_size*9 - 9;
    *ftemp_point = *(f_point);
    *(ftemp_point - 9 + 3) = *(f_point+3);
    *(ftemp_point + 9 * x_size + 4) = *(f_point+4);
    *(ftemp_point + 9 * (x_size - 1) + 7) = *(f_point+7);


    //
    // bottom
    //
    for (int i=1; i<x_size-1; ++i){

      ftemp_point = f_temp  +(y_size-1)*x_size*9+ i*9;
      f_point = f + i*9;
      vel_point = vel + (y_size-1)*x_size*2+ i * 2;
      P[i][y_size-1] = P[i][y_size-2];
      rho_point[i][y_size-1] = P[i][y_size-1]  / Tw;
      *vel_point = *(vel_point - x_size*2);
      *(vel_point+1) = *(vel_point + 1 - x_size*2);
      equalibrum(ftemp_point, vel_point, rho_point[i][0], Tw);
      *ftemp_point = *(f_point++);
      *(ftemp_point + 9 +1) = *(f_point+1);
      *(ftemp_point - 9 * x_size+2) = *(f_point+2);
      *(ftemp_point - 9 + 3) = *(f_point+3);
      *(ftemp_point - 9 * (x_size - 1) + 5) = *(f_point+5);
      *(ftemp_point - 9 * (x_size + 1) + 6) = *(f_point+6);

    }

    //Corners:
    //left-bot:


    ftemp_point = f_temp + (y_size-1)*x_size*9;
    f_point = f + (y_size-1)*x_size*9;
    *ftemp_point = *(f_point);
    *(ftemp_point + 9 + 1) = *(f_point+1);
    *(ftemp_point - 9 * x_size + 2) = *(f_point+2);
    *(ftemp_point - 9 * (x_size - 1) + 5) = *(f_point+5);
    //right-bot:


    ftemp_point = f_temp + (y_size-1)*x_size*9 + x_size*9 - 9;
    f_point = f + (y_size-1)*x_size*9 + x_size*9 - 9;
    *ftemp_point = *(f_point);
    *(ftemp_point - 9 + 3) = *(f_point+3);
    *(ftemp_point - 9 * x_size + 2) = *(f_point+2);
    *(ftemp_point - 9 * (x_size + 1) + 6) = *(f_point+6);

    double temp_tau = 0;
    //
    // left and right
    //
    for (int j = 1; j < y_size-1; ++j){
      /*
      ftemp_point = f_temp + j*x_size*9;
      f_point = f + j*x_size*9;
      feq_point = f_eq + j*x_size*9;
      vel_point = vel + j*x_size*2;
      rho_point[0][j] = P1 / gam / T[0][j];
      P[0][j] = P1;
      rho_point[x_size-1][j] = P2 / gam / T[x_size-1][j] ;
      *vel_point = * (vel_point+2);
      *(vel_point+1) = * (vel_point+3);
      *(vel_point+(x_size-1)*2) = *(vel_point+(x_size-2)*2);
      *(vel_point+(x_size-1)*2+1) = *(vel_point+(x_size-2)*2+1);
      for (int k = 0; k < 9; ++k)
      {
        temp_tau = 1/omega(T[0][j]);
        *(f_point+k) = (1 + temp_tau)* *(feq_point+k) - temp_tau * *(feq_point+9+k);
        temp_tau = 1/omega(T[x_size-1][j]);
        *(f_point + 9*(x_size-1) + k) = (1 + temp_tau) * *(feq_point + 9*(x_size-1) + k) - temp_tau * *(feq_point + 9*(x_size-2) + k);
      }
      *ftemp_point = *(f_point);
      *(ftemp_point + 9 +1) = *(f_point+1);
      *(ftemp_point - 9 * x_size + 2) = *(f_point+2);
      *(ftemp_point + 9 * x_size + 4) = *(f_point+4);
      *(ftemp_point - 9 * (x_size - 1) + 5) = *(f_point+5);
      *(ftemp_point + 9 * (x_size + 1) + 8) = *(f_point+8);
      ftemp_point = ftemp_point + 9*(x_size-1);
      f_point = f_point + 9*(x_size-1);
      *(ftemp_point) = *(f_point);
      *(ftemp_point - 9 * x_size+2) = *(f_point+2);
      *(ftemp_point - 9 + 3) = *(f_point+3);
      *(ftemp_point + 9 * x_size + 4) = *(f_point+4);
      *(ftemp_point - 9 * (x_size + 1) + 6) = *(f_point+6);
      *(ftemp_point + 9 * (x_size - 1) + 7) = *(f_point+7);
      */

      /*
      ftemp_point = f_temp + j*x_size*9;
      f_point = f_temp + j*x_size*9 + 9;
      *ftemp_point = *(f_point);
      *(ftemp_point + 1) = *(f_point+1);
      *(ftemp_point + 2) = *(f_point+2);
      *(ftemp_point + 3) = *(f_point+3);
      *(ftemp_point + 4) = *(f_point+4);
      *(ftemp_point + 5) = *(f_point+5);
      *(ftemp_point + 6) = *(f_point+6);
      *(ftemp_point + 7) = *(f_point+7);
      *(ftemp_point + 8) = *(f_point+8);
      ftemp_point = f_temp + j*x_size*9 + 9*(x_size-1);
      f_point = f_temp + j*x_size*9 + 9*(x_size-1)- 9;
      *ftemp_point = *(f_point);
      *(ftemp_point + 1) = *(f_point+1);
      *(ftemp_point + 2) = *(f_point+2);
      *(ftemp_point + 3) = *(f_point+3);
      *(ftemp_point + 4) = *(f_point+4);
      *(ftemp_point + 5) = *(f_point+5);
      *(ftemp_point + 6) = *(f_point+6);
      *(ftemp_point + 7) = *(f_point+7);
      *(ftemp_point + 8) = *(f_point+8);
      */

      // ftemp_point = f + j*x_size*9;
      // f_point = f + j*x_size*9 + 9;
      //
      // *ftemp_point = *(f_point);
      // *(ftemp_point + 1) = *(f_point+1);
      // *(ftemp_point + 2) = *(f_point+2);
      // *(ftemp_point + 3) = *(f_point+3);
      // *(ftemp_point + 4) = *(f_point+4);
      // *(ftemp_point + 5) = *(f_point+5);
      // *(ftemp_point + 6) = *(f_point+6);
      // *(ftemp_point + 7) = *(f_point+7);
      // *(ftemp_point + 8) = *(f_point+8);
      ftemp_point = f_temp + j*x_size*9;
      f_point = f + j*x_size*9;
      *ftemp_point = *(f_point);
      *(ftemp_point + 9 + 1) = *(f_point + 1);
      *(ftemp_point - 9 * x_size + 2) = *(f_point+2);
      *(ftemp_point + 9 * x_size + 4) = *(f_point+4);
      *(ftemp_point - 9 * (x_size - 1) + 5) = *(f_point+5);
      *(ftemp_point + 9 * (x_size + 1) + 8) = *(f_point+8);


      // ftemp_point = f + j*x_size*9 + 9*(x_size-1);
      // f_point = f + j*x_size*9 - 9 + 9*(x_size-1);
      // *ftemp_point = *(f_point);
      // *(ftemp_point + 1) = *(f_point+1);
      // *(ftemp_point + 2) = *(f_point+2);
      // *(ftemp_point + 3) = *(f_point+3);
      // *(ftemp_point + 4) = *(f_point+4);
      // *(ftemp_point + 5) = *(f_point+5);
      // *(ftemp_point + 6) = *(f_point+6);
      // *(ftemp_point + 7) = *(f_point+7);
      // *(ftemp_point + 8) = *(f_point+8);
      ftemp_point = f_temp + j*x_size*9+ 9*(x_size-1);
      f_point = f + j*x_size*9+ 9*(x_size-1);
      *ftemp_point = *(f_point);
      *(ftemp_point - 9 + 3) = *(f_point+3);
      *(ftemp_point - 9 * x_size + 2) = *(f_point+2);
      *(ftemp_point + 9 * x_size + 4) = *(f_point+4);
      *(ftemp_point - 9 * (x_size + 1) + 6) = *(f_point + 6);
      *(ftemp_point + 9 * (x_size - 1) + 7) = *(f_point + 7);
    }



    // Core:
    for (int j = 1; j < y_size-1; ++j){
      for (int i = 1; i < x_size-1; ++i){
        ftemp_point = f_temp + j*x_size*9 + i*9;
        f_point = f + j*x_size*9 + i*9;
        feq_point = f_eq + j*x_size*9 + i*9;
        vel_point = vel + j*x_size*2 + i*2;
        rho_point[i][j] = macro_rho(ftemp_point);
        macro_vel(vel_point, ftemp_point, rho_point[i][j], T[i][j]);
        T[i][j] = macro_temp(ftemp_point,vel_point, rho_point[i][j]);
        P[i][j] = rho_point[i][j]*T[i][j];
        equalibrum(feq_point, vel_point,rho_point[i][j], T[i][j]);
        }
     }
    //
  for (int j = 1; j < y_size-1; ++j){
    for (int i = 1; i < x_size-1; ++i){

      for (int k=0; k<9; ++k)
      {
        f[j*x_size*9 + i*9+k] = f_temp[j*x_size*9 + i*9+k] - omega(T[i][j], rho_point[i][j])*(f_temp[j*x_size*9 + i*9 + k] - f_eq[j*x_size*9 + i*9 + k]);
      }
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







  for (int j=1; j<y_size-1; j++)
    for (int i=1; i<x_size-1; i++)
      fprintf(veldat, "%d %d %f \n", i, j, vel[i*2+j*2*x_size]);
  for (int j=1; j<y_size-1; j++)
    for (int i=1; i<x_size-1; i++)
      fprintf(veldat2, "%d %d %f \n", i, j, vel[i*2+j*2*x_size]);

  for (int j=0; j<y_size; j++)
    for (int i=0; i<x_size; i++)
      fprintf(rho, "%d %d %f \n", i, j, rho_point[i][j]);

  for (int j=1; j<y_size-1; j++)
    for (int i=1; i<x_size-1; i++)
      fprintf(Tdat, "%d %d %f \n", i, j, T[i][j]);

  for (int i=1; i<y_size-1; i++)
      fprintf(debug, "%dx %f \n", i, vel[6*2+i*2*x_size]);
  return 0;
}
