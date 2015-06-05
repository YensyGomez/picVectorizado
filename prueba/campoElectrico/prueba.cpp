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

void Electric_Field_Carlos (double phi[J_X][J_Y], double E_X[J_X][J_Y], double E_Y[J_X][J_Y])  // Función para calcular el campo eléctrico 
{                                                                                       //en los puntos de malla a partir del potencial. 
  for (int j=1;j<J_X-1;j++)
  {
      for (int k=1;k<J_Y-1;k++)
      {
        E_X[j][k]=(phi[j-1][k]-phi[j+1][k])/(2.*hx);
        E_Y[j][k]=(phi[j][k-1]-phi[j][k+1])/(2.*hx);

        E_X[0][k]=0.0;  //Cero en las fronteras X
        E_Y[0][k]=0.0;
        E_X[J_X-1][k]=0.0; 
        E_Y[J_X-1][k]=0.0;
      }          

      E_X[j][0]=(phi[j-1][0]-phi[j+1][0])/(2.*hx); 
      E_Y[j][0]=(phi[j][J_Y-1]-phi[j][1])/(2.*hx);
      


      E_X[j][J_Y-1]=(phi[j-1][J_Y-1]-phi[j+1][J_Y-1])/(2.*hx);
      E_Y[j][J_Y-1]=(phi[j][J_Y-2]-phi[j][0])/(2.*hx);
  }
}




void Electric_Field_Yen(double *phi, double *E_X, double *E_Y)  // Función para calcular el campo eléctrico 
{                                                                                       //en los puntos de malla a partir del potencial. 
  for (int j=1;j<J_X-1;j++)
  {
      for (int k=1;k<J_Y-1;k++)
      {
        E_X[j*J_Y+k]=(phi[(j-1)*J_Y+k]-phi[(j+1)*J_Y+k])/(2.*hx);
        E_Y[j*J_Y+k]=(phi[j*J_Y+(k-1)]-phi[j*J_Y+(k+1)])/(2.*hx);

        E_X[0*J_Y+k]=0.0;  //Cero en las fronteras X
        E_Y[0*J_Y+k]=0.0;
        E_X[(J_X-1)*J_Y+k]=0.0; 
        E_Y[(J_X-1)*J_Y+k]=0.0;
      }          

      E_X[j*J_Y+0]=(phi[(j-1)*J_Y+0]-phi[((j+1)*J_Y+0)])/(2.*hx); 
      E_Y[j*J_Y+0]=(phi[j*J_Y+(J_Y-1)]-phi[j*J_Y+1])/(2.*hx);
   

      E_X[j*J_Y+(J_Y-1)]=(phi[(j-1)*J_Y+(J_Y-1)]-phi[(j+1)*J_Y+(J_Y-1)])/(2.*hx);
      E_Y[j*J_Y+(J_Y-1)]=(phi[j*J_Y+(J_Y-2)]-phi[j*J_Y+0])/(2.*hx);
  }
}

/*void Electric_Field_Yen (double *phi, double *E_X, double *E_Y)  // Función para calcular el campo eléctrico 
{     

  for(int j= 0; j<J_X; j++){
    for(int k=0; k<J_Y; k++){

    E_X[0*J_Y+k]=0.0; 
    E_X[(J_X-1)*J_Y+k]=0.0; 
    E_Y[0*J_Y+k]=0.0;
    E_Y[(J_X-1)*J_Y+k]=0.0;
   
    }
    
    E_X[j*J_Y+0]=(phi[(j-1)*J_Y+0]-phi[(j+1*J_Y+0)])/(2.*hx); 
    E_Y[j*J_Y+0]=(phi[j*J_Y+(J_Y-1)]-phi[j*J_Y+1])/(2.*hx);
    E_X[j*J_Y+(J_Y-1)]=(phi[(j-1)*J_Y+(J_Y-1)]-phi[(j+1)*J_Y+(J_Y-1)])/(2.*hx);
    E_Y[j*J_Y+(J_Y-1)]=(phi[j*J_Y+(J_Y-2)]-phi[j*J_Y+0])/(2.*hx);
  
  } 
                                                                                  //en los puntos de malla a partir del potencial. 
  for (int j=1;j<J_X-1;j++)
  {
      for (int k=1;k<J_Y-1;k++)
      {
        E_X[j*J_Y+k]=(phi[(j-1)*J_Y+k]-phi[(j+1)*J_Y+k])/(2.*hx);
        E_Y[j*J_Y+k]=(phi[j*J_Y+(k-1)]-phi[j*J_Y+(k+1)])/(2.*hx);

         //Cero en las fronteras X
        
      }                
  }
}*/


using namespace std;

int main(){

	int i, j;
	double comodin;
	double phiC[J_X][J_Y], E_XC[J_X][J_Y], E_YC[J_X][J_Y];
	int size1= J_X*J_Y*sizeof(double);	
	double *E_xYen, *E_yYen, *phiYen;
	E_xYen = (double *) malloc(size1);
 	E_yYen = (double *) malloc(size1);
 	phiYen = (double *) malloc(size1);

	FILE *fp, *fp1, *fp2, *fp3;
	fp = fopen ( "./Poisson5.data", "r" );        
	if (fp==NULL){
		fputs ("File error leer Poisson1",stderr); exit (1);
	}

	for(i=0; i<J_X*J_Y;i++){
		fscanf(fp,"%lf %lf %lf\n",&comodin, &comodin, &phiYen[i]);
	}

	fclose ( fp );

	fp3 = fopen ( "./Poisson5.data", "r" );        
	if (fp3==NULL){
		fputs ("File error Poisson2",stderr); exit (1);
	}

	for(i=0; i<J_X; i++){
		for(j=0; j<J_Y; j++)
		{
			fscanf(fp3,"%lf %lf %lf\n",&comodin, &comodin, &phiC[i][j]);
		}
		
	}
	fclose(fp3);

	Electric_Field_Carlos (phiC, E_XC, E_YC);
	Electric_Field_Yen(phiYen,E_xYen,E_yYen);

	fp1 = fopen ( "EfieldYen", "w" );        
	if (fp1==NULL){
		fputs ("File error escribiendo campo",stderr); exit (1);
	}

	fp2 = fopen ( "EfieldCarlos", "w" );        
	if (fp2==NULL){
		fputs ("File error escribiendo campo2",stderr); exit (1);
	}

	for(i=0; i<J_X*J_Y;i++){
		fprintf(fp1,"%lf %lf %lf\n",E_xYen[i], E_yYen[i], phiYen[i]);
	}

	for(i=0; i<J_X; i++){
		for(j=0; j<J_Y; j++)
		{
			fprintf(fp2,"%lf %lf %lf\n",E_XC[i][j], E_YC[i][j], phiC[i][j]);
		}
		
	}



	fclose (fp1);
	fclose (fp2);

	free(E_xYen);
	free(E_yYen);
	free(phiYen); 
	return(0);
}
