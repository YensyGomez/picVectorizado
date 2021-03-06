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
#define cte_rho 160.722398
#define Factor_carga_e 10
#define n_0 1.000000e-30 




void rhoCarlos(double ne[J_X][J_Y], double ni[J_X][J_Y], double rho[J_X][J_Y]){
 
   for (int j = 0; j < J_X; j++) 
    {
        for (int k = 0; k < J_Y; k++) 
        {
          rho[j][k]= cte_rho* Factor_carga_e*(ni[j][k]- ne[j][k])/n_0;
        }
    }

}


void rhoYen( double *ne, double *ni, double *rho){


for (int j = 0; j < J_X; j++) 
    {
        for (int k = 0; k < J_Y; k++) 
        {
          rho[j*J_Y+k]= cte_rho* Factor_carga_e*(ni[j*J_Y+k]- ne[j*J_Y+k])/n_0;
        }
    }

}


int main(){

	double ne_Carlos[J_X][J_Y];
	double ni_Carlos[J_X][J_Y];
	double rho_Carlos[J_X][J_Y];

	int size = J_X*J_Y*sizeof(double);

	double *rho_Yen, *ne_Yen, *ni_Yen;


	rho_Yen = (double *)malloc(size);
	ne_Yen = (double *)malloc(size);
	ni_Yen = (double *)malloc(size);



	FILE *fp, *fp1, *fp2, *fp3;
	fp = fopen ( "./n5.data", "r" );        
	if (fp==NULL){
		fputs ("File error leer Rho1",stderr);
		 exit (1);
	}

	for(int i=0; i<J_X*J_Y;i++){
		fscanf(fp,"%lf %lf\n", &ni_Yen[i], &ne_Yen[i]);
	}

	fclose ( fp );
	
	
	fp3 = fopen ( "./n5.data", "r" );        
	if (fp3==NULL){
		fputs ("File error leer Rho2",stderr); exit (1);
	}

	for(int i=0; i<J_X; i++){
		for(int j=0; j<J_Y; j++)
		{
			fscanf(fp3,"%lf %lf\n", &ni_Carlos[i][j], &ne_Carlos[i][j]);
		}
		
	}
	fclose(fp3);
	
	
	rhoCarlos(ne_Carlos, ni_Carlos, rho_Carlos);
	rhoYen(ne_Yen, ni_Yen, rho_Yen); 
	
	
  fp1 = fopen ( "rhoYen", "w" );        
	if (fp1==NULL){
		fputs ("File error escribiendo rho1",stderr);
		 exit (1);
	}

	fp2 = fopen ( "rhoCarlos", "w" );        
	if (fp2==NULL){
		fputs ("File error escribiendo rho2",stderr); 
		exit (1);
	}

	for(int i=0; i<J_X*J_Y;i++){
		fprintf(fp1," %lf\n", rho_Yen[i]);
	}

	for(int i=0; i<J_X; i++){
		for(int j=0; j<J_Y; j++)
		{
			fprintf(fp2," %lf\n", rho_Carlos[i][j]);
		}
		
	}
	fclose (fp1);
	fclose (fp2);
	free(rho_Yen);
	free(ne_Yen);
	free(ni_Yen); 

 return 0; 
}
