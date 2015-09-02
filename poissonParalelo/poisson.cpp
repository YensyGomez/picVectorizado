#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <complex>
#include <ctime>
#include <cstring>
#include <fstream>

#include <fftw3.h>



#define J_X          129           // Número de puntos de malla X. Recomendado: Del orden 2^n+1
#define J_Y           64           // Número de puntos de malla Y. Recomendado: Del orden 2^n
#define hx  0.9830640


//***********************************************************************
//Cálculo del Potencial Electrostático en cada uno de los puntos de malla
//***********************************************************************
//
using namespace std;


void   poisson2D_dirichletX_periodicY(double *phi, complex<double>  *rho);
//void poisson2D_dirichletX_periodicY_cuda(double *phi, complex<double> *rho); 


int  main(){
  int i;  
  int size2 = J_X * J_Y * sizeof(double);
  int size3 = J_X * J_Y * sizeof(complex<double>);
  double *phi, *x, *y; 
  
  complex<double> *rho; 
  
  rho = (complex<double> *) malloc(size3);
  phi   = (double *) malloc(size2);
  x = (double *) malloc(size2);
  y = (double *) malloc(size2);
  

  ifstream file("rho50000.data");
  for (i  =  0; i < J_X*J_Y; i++) {
    file >> x[i] >> y[i] >> rho[i];  
    
  }
  
  file.close();
  
  //poisson2D_dirichletX_periodicY_cuda(phi, rho); 
  poisson2D_dirichletX_periodicY(phi, rho);
  
  ofstream file2("salida1.data");
  for (int i  =  0; i < J_X; i++) {
    double thisx  =  i * hx;
    for (int j  =  0; j < J_Y; j++) {
      double thisy  =  j * hx;
     file2 << thisx << '\t' << thisy << '\t' << phi[(i * J_Y) + j] << '\n';
    }
    file2 << '\n';
  }
  file2.close();

  free(x);
  free(y);
  free(rho); 
  free(phi);
  
  
  return(0);
}


/*void poisson2D_dirichletX_periodicY_cuda(double *phi, complex<double> *rho) {
  
  int M = J_X - 2, N = J_Y;
  double h = hx;
  double hy = hx;
  double *f;
  //fftw_complex  *f2;
  cufftDoubleComplex *f2_d; 
  
  
  fftw_plan p, p_i; 
  cufftHandle p_y,p_yi ; 
  
  f = (double*) fftw_malloc(sizeof(double)* M);
  cudaMalloc((void**)&f2_d, sizeof(cufftDoubleComplex)*N); 
  
  int size2 = J_X*J_Y*sizeof(cufftDoubleComplex);
  
  p = fftw_plan_r2r_1d(M, f, f, FFTW_RODFT00, FFTW_ESTIMATE);  
  cufftPlan1d(&p_y, N,CUFFT_Z2Z,1);
  p_i = fftw_plan_r2r_1d(M, f, f, FFTW_RODFT00, FFTW_ESTIMATE);
  cufftPlan1d(&p_yi,N,CUFFT_Z2Z,1);
  
  
  // Columnas FFT
  for (int k = 0; k < N; k++) {
    for (int j = 0; j < M; j++)
      f[j] = rho[(j + 1) * N + k].real();
    fftw_execute(p);
    for (int j = 0; j < M; j++)
      rho[(j + 1) * N + k].real() = f[j];
  }
  
  // Filas FFT
  
  for (int j = 0; j < M; j++) {
    for (int k = 0; k < N; k++)
      cudaMemcpy(&f2_d[k], &rho[(j + 1) * N + k],size2, cudaMemcpyHostToDevice);
    cufftExecZ2Z(p_y,f2_d,f2_d,CUFFT_FORWARD);
    for (int k = 0; k < N; k++)
      cudaMemcpy(&rho[(j + 1) * N + k], &f2_d[k],size2, cudaMemcpyDeviceToHost);
  }

  // Resolver en el espacio de Fourier
  complex<double> i(0.0, 1.0);
  double pi = M_PI;
  complex<double> Wy = exp(2.0 * pi * i / double(N));
  complex<double> Wn = 1.0;
  for (int m = 0; m < M; m++) {
    for (int n = 0; n < N; n++) {
      complex<double> denom = h * h * 2.0 + hy * hy * 2.0;
      denom -= hy * hy * (2 * cos((m + 1) * pi / (M + 1))) + h * h * (Wn + 1.0 / Wn);
      if (denom != 0.0)
        rho[(m + 1) * N + n] *= h * h * hy * hy / denom;
      Wn *= Wy;
    }
  }

  // Inversa de las filas
  for (int j = 0; j < M; j++) {
    for (int k = 0; k < N; k++)
      cudaMemcpy(&f2_d[k], &rho[(j + 1) * N + k],size2, cudaMemcpyHostToDevice);
    //fftw_execute(p_yi);
     cufftExecZ2Z(p_yi,f2_d,f2_d,CUFFT_INVERSE);
    for (int k = 0; k < N; k++) {
      cudaMemcpy(&rho[(j + 1) * N + k], &f2_d[k],size2, cudaMemcpyDeviceToHost);
      rho[(j + 1) * N + k] /= double(N); //La transformada debe ser normalizada.
    }
  }

  //Inversa Columnas FFT
  for (int k = 0; k < N; k++) {
    for (int j = 0; j < M; j++)
      f[j]=rho[(j + 1) * N + k].real();
    fftw_execute(p_i);
    for (int j = 0; j < M; j++)
      phi[(j + 1) * N + k] = f[j] / double(2 * (M + 1));
  } // pasar a cuda..
  
  //Aqui es importante pasar el phi de la función para el potencial

  for (int k = 0; k < N; k++) {
    phi[k]=0;
    phi[(J_X - 1) * N + k]=0;
  }

  fftw_destroy_plan(p);
  cufftDestroy(p_y);
  cufftDestroy(p_yi);
  fftw_destroy_plan(p_i);
 
  fftw_free(f);
  cudaFree(f2_d);
 
}*/

void poisson2D_dirichletX_periodicY(double *phi, complex<double> *rho) {
  int M = J_X - 2, N = J_Y;
  double h = hx;
  double hy = hx;
  double *f;
  fftw_complex  *f2;
  fftw_plan p, p_y, p_i, p_yi;
  f = (double*) fftw_malloc(sizeof(double)* M);
  f2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

  p = fftw_plan_r2r_1d(M, f, f, FFTW_RODFT00, FFTW_ESTIMATE);
  p_y = fftw_plan_dft_1d(N, f2, f2, FFTW_FORWARD, FFTW_ESTIMATE);
  p_i = fftw_plan_r2r_1d(M, f, f, FFTW_RODFT00, FFTW_ESTIMATE);
  p_yi = fftw_plan_dft_1d(N, f2, f2, FFTW_BACKWARD, FFTW_ESTIMATE);


  // Columnas FFT
  for (int k = 0; k < N; k++) {
    for (int j = 0; j < M; j++)
      f[j] = rho[(j + 1) * N + k].real();
    fftw_execute(p);
    for (int j = 0; j < M; j++)
      rho[(j + 1) * N + k].real() = f[j];
  }

  // Filas FFT
  for (int j = 0; j < M; j++) {
    for (int k = 0; k < N; k++) 
      memcpy( &f2[k], &rho[(j + 1) * N + k], sizeof( fftw_complex ) );
    fftw_execute(p_y);
    for (int k = 0; k < N; k++)
      memcpy( &rho[(j + 1) * N + k], &f2[k], sizeof( fftw_complex ) );
  }



  // Resolver en el espacio de Fourier
  complex<double> i(0.0, 1.0);
  double pi = M_PI;
  complex<double> Wy = exp(2.0 * pi * i / double(N));
  complex<double> Wn = 1.0;
  for (int m = 0; m < M; m++) {
    for (int n = 0; n < N; n++) {
      complex<double> denom = h * h * 2.0 + hy * hy * 2.0;
      denom -= hy * hy * (2 * cos((m + 1) * pi / (M + 1))) + h * h * (Wn + 1.0 / Wn);
      if (denom != 0.0)
        rho[(m + 1) * N + n] *= h * h * hy * hy / denom;
      Wn *= Wy;
    }
  }

  // Inversa de las filas
  for (int j = 0; j < M; j++) {
    for (int k = 0; k < N; k++)
      memcpy( &f2[k], &rho[(j + 1) * N + k], sizeof( fftw_complex ) );
    fftw_execute(p_yi);
    for (int k = 0; k < N; k++) {
      memcpy( &rho[(j + 1) * N + k], &f2[k], sizeof( fftw_complex ) );
      rho[(j + 1) * N + k] /= double(N); //La transformada debe ser normalizada.
    }
  }

  //Inversa Columnas FFT
  for (int k = 0; k < N; k++) {
    for (int j = 0; j < M; j++)
      f[j]=rho[(j + 1) * N + k].real();
    fftw_execute(p_i);
    for (int j = 0; j < M; j++)
      phi[(j + 1) * N + k] = f[j] / double(2 * (M + 1));
  }

  for (int k = 0; k < N; k++) {
    phi[k]=0;
    phi[(J_X - 1) * N + k]=0;
  }

  fftw_destroy_plan(p);
  fftw_destroy_plan(p_i);
  fftw_destroy_plan(p_y);
  fftw_destroy_plan(p_yi);
  fftw_free(f);
  fftw_free(f2);
}
