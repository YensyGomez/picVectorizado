#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <complex>
#include <ctime>
#include <cstring>
#include <fstream>
#include <fftw3.h>


#define J_X 64
#define J_Y 16
#define hx 0.983064 


using namespace std;

void Poisson2D_DirichletX_PeriodicY_Carlos(double phi[J_X][J_Y], complex<double> rho[J_X][J_Y])
{
    int M=J_X-2,N=J_Y;
    double h = hx;
    double hy = hx;
    double *f;
    fftw_complex  *f2;
    fftw_plan p,p_y,p_i,p_yi;
    f= (double*) fftw_malloc(sizeof(double)* M);
    f2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    p = fftw_plan_r2r_1d(M, f, f, FFTW_RODFT00, FFTW_ESTIMATE);
    p_y = fftw_plan_dft_1d(N, f2, f2, FFTW_FORWARD, FFTW_ESTIMATE);      
    p_i = fftw_plan_r2r_1d(M, f, f, FFTW_RODFT00, FFTW_ESTIMATE);
    p_yi = fftw_plan_dft_1d(N, f2, f2, FFTW_BACKWARD, FFTW_ESTIMATE);


    // Columnas FFT
    for (int k = 0; k < N; k++) {
        for (int j = 0; j < M; j++)
            f[j]=rho[j+1][k].real();
        fftw_execute(p);
        for (int j = 0; j < M; j++)
            rho[j+1][k].real()=f[j];
    }

    // Filas FFT
    for (int j = 0; j < M; j++) {
        for (int k = 0; k < N; k++)
            memcpy( &f2[k], &rho[j+1][k], sizeof( fftw_complex ) ); 
        fftw_execute(p_y);
        for (int k = 0; k < N; k++)
            memcpy( &rho[j+1][k], &f2[k], sizeof( fftw_complex ) );
    }



    // Resolver en el espacio de Fourier
    complex<double> i(0.0, 1.0);
    double pi = M_PI;
    complex<double> Wy = exp(2.0 * pi * i / double(N));
    complex<double> Wn = 1.0;
    for (int m = 0; m < M; m++) {
        for (int n = 0; n < N; n++) {
            complex<double> denom = h*h*2.0+hy*hy*2.0;
            denom -= hy*hy*(2*cos((m+1)*pi/(M+1))) + h*h*(Wn + 1.0 / Wn);
            if (denom != 0.0) 
                rho[m+1][n] *= h*h*hy*hy / denom;
            Wn *= Wy;
        }
    }

   // Inversa de las filas
    for (int j = 0; j < M; j++) {
        for (int k = 0; k < N; k++)
            memcpy( &f2[k], &rho[j+1][k], sizeof( fftw_complex ) );
        fftw_execute(p_yi);
        for (int k = 0; k < N; k++)
        {
            memcpy( &rho[j+1][k], &f2[k], sizeof( fftw_complex ) );
            rho[j+1][k] /= double(N); //La transformada debe ser normalizada.
        }
    }

    //Inversa Columnas FFT
    for (int k = 0; k < N; k++) {
        for (int j = 0; j < M; j++)
            f[j]=rho[j+1][k].real();
        fftw_execute(p_i);
        for (int j = 0; j < M; j++)
            phi[j+1][k]=f[j]/double(2*(M+1));
    }

    for (int k = 0; k < N; k++) 
    {
      phi[0][k]=0;
      phi[J_X-1][k]=0;
    }

    fftw_destroy_plan(p);
    fftw_destroy_plan(p_i);
    fftw_destroy_plan(p_y);
    fftw_destroy_plan(p_yi);
    fftw_free(f); fftw_free(f2);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Poisson2D_DirichletX_PeriodicY_Yen(double *phi, complex<double> *rho)
{
    int M=J_X-2,N=J_Y;
    double h = hx;
    double hy = hx;
    double *f;
    fftw_complex  *f2;
    fftw_plan p,p_y,p_i,p_yi;
    f= (double*) fftw_malloc(sizeof(double)* M);
    f2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    p = fftw_plan_r2r_1d(M, f, f, FFTW_RODFT00, FFTW_ESTIMATE);
    p_y = fftw_plan_dft_1d(N, f2, f2, FFTW_FORWARD, FFTW_ESTIMATE);      
    p_i = fftw_plan_r2r_1d(M, f, f, FFTW_RODFT00, FFTW_ESTIMATE);
    p_yi = fftw_plan_dft_1d(N, f2, f2, FFTW_BACKWARD, FFTW_ESTIMATE);


    // Columnas FFT
    for (int k = 0; k < N; k++) {
        for (int j = 0; j < M; j++)
            f[j]=rho[(j+1)*N+k].real();
        fftw_execute(p);
        for (int j = 0; j < M; j++)
            rho[(j+1)*N+k].real()=f[j];
    }

    // Filas FFT
    for (int j = 0; j < M; j++) {
        for (int k = 0; k < N; k++)
            memcpy( &f2[k], &rho[(j+1)*N+k], sizeof( fftw_complex ) ); 
        fftw_execute(p_y);
        for (int k = 0; k < N; k++)
            memcpy( &rho[(j+1)*N+k], &f2[k], sizeof( fftw_complex ) );
    }



    // Resolver en el espacio de Fourier
    complex<double> i(0.0, 1.0);
    double pi = M_PI;
    complex<double> Wy = exp(2.0 * pi * i / double(N));
    complex<double> Wn = 1.0;
    for (int m = 0; m < M; m++) {
        for (int n = 0; n < N; n++) {
            complex<double> denom = h*h*2.0+hy*hy*2.0;
            denom -= hy*hy*(2*cos((m+1)*pi/(M+1))) + h*h*(Wn + 1.0 / Wn);
            if (denom != 0.0) 
                rho[(m+1)*N+n] *= h*h*hy*hy / denom;
            Wn *= Wy;
        }
    }

   // Inversa de las filas
    for (int j = 0; j < M; j++) {
        for (int k = 0; k < N; k++)
            memcpy( &f2[k], &rho[(j+1)*N+k], sizeof( fftw_complex ) );
        fftw_execute(p_yi);
        for (int k = 0; k < N; k++)
        {
            memcpy( &rho[(j+1)*N+k], &f2[k], sizeof( fftw_complex ) );
            rho[(j+1)*N+k] /= double(N); //La transformada debe ser normalizada.
        }
    }

    //Inversa Columnas FFT
    for (int k = 0; k < N; k++) {
        for (int j = 0; j < M; j++)
            f[j]=rho[(j+1)*N+k].real();
        fftw_execute(p_i);
        for (int j = 0; j < M; j++)
            phi[(j+1)*N+k]=f[j]/double(2*(M+1));
    }

    for (int k = 0; k < N; k++) 
    {
      phi[0*N+k]=0;
      phi[(J_X-1)*N+k]=0;
    }

    fftw_destroy_plan(p);
    fftw_destroy_plan(p_i);
    fftw_destroy_plan(p_y);
    fftw_destroy_plan(p_yi);
    fftw_free(f); fftw_free(f2);
}


int main(){ 

double comodin; 
int size = J_X*J_Y*sizeof(double);
int size1 = J_X*J_Y*sizeof(complex<double>);

double phiCarlos[J_X][J_Y];
double *phiYen; 
complex<double> rhoCarlos[J_X][J_Y];
complex<double>  *rhoYen;

rhoYen = (complex<double> *) malloc(size1);
phiYen = (double *) malloc(size);

FILE *fp, *fp1, *fp2, *fp3;
	fp = fopen ( "./Rho5.data", "r" );        
	if (fp==NULL){
		fputs ("File error leer Rho1",stderr);
		 exit (1);
	}

	for(int i=0; i<J_X*J_Y;i++){
		fscanf(fp,"%lf %lf %lf %lf\n",&comodin, &comodin, &rhoYen[i].real(), &rhoYen[i].imag());
	}

	fclose ( fp );
	
	
	fp3 = fopen ( "./Rho5.data", "r" );        
	if (fp3==NULL){
		fputs ("File error leer Rho2",stderr); exit (1);
	}

	for(int i=0; i<J_X; i++){
		for(int j=0; j<J_Y; j++)
		{
			fscanf(fp3,"%lf %lf %lf %lf\n",&comodin, &comodin, &rhoCarlos[i][j].real(), &rhoCarlos[i][j].imag());
		}
		
	}
	fclose(fp3);
	
	
	
	Poisson2D_DirichletX_PeriodicY_Carlos(phiCarlos, rhoCarlos);
	Poisson2D_DirichletX_PeriodicY_Yen(phiYen,rhoYen);
	
	
  fp1 = fopen ( "potencialYen", "w" );        
	if (fp1==NULL){
		fputs ("File error escribiendo potencial1",stderr);
		 exit (1);
	}

	fp2 = fopen ( "PotencialCarlos", "w" );        
	if (fp2==NULL){
		fputs ("File error escribiendo potencial2",stderr); 
		exit (1);
	}

	for(int i=0; i<J_X*J_Y;i++){
		fprintf(fp1," %lf\n", phiYen[i]);
	}

	for(int i=0; i<J_X; i++){
		for(int j=0; j<J_Y; j++)
		{
			fprintf(fp2," %lf\n", phiCarlos[i][j]);
		}
		
	}
	fclose (fp1);
	fclose (fp2);
	free(phiYen);
	free(rhoYen);

	return 0;
}
