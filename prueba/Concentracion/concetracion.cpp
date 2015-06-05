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
#define max_SPe 10000
#define __d {prueba (__LINE__);}


using namespace std;
void prueba(int d){
        cout<<"LINE: "<< d <<endl;
}

void Concentration_Carlos(double pos[max_SPe][2], double n[J_X][J_Y],int NSP)
{ 
  int j_x,j_y;
  double temp_x,temp_y;
  double jr_x,jr_y;
  for (j_x=0; j_x<J_X; j_x++)
    {
      for (j_y=0; j_y<J_Y; j_y++)
      {
        n[j_x][j_y] = 0.;
      }
    } // Inicializar densidad de carga

  for (int i=0;i<NSP;i++)
    {
       jr_x=pos[i][0]/hx;       // indice (real) de la posición de la superpartícula
       j_x =int(jr_x);       // indice  inferior (entero) de la celda que contiene a la superpartícula
       temp_x = jr_x-j_x;
       jr_y=pos[i][1]/hx;       // indice (real) de la posición de la superpartícula
       j_y =int(jr_y);       // indice  inferior (entero) de la celda que contiene a la superpartícula
       temp_y = jr_y-j_y;

       n[j_x][j_y] += (1.-temp_x)*(1.-temp_y)/(hx*hx*hx);
       n[j_x+1][j_y] += temp_x*(1.-temp_y)/(hx*hx*hx);
       n[j_x][j_y+1] += (1.-temp_x)*temp_y/(hx*hx*hx);
       n[j_x+1][j_y+1] += temp_x*temp_y/(hx*hx*hx);

    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

void Concentration_Yen (double *pos_x, double *pos_y, double *n,int NSP)
{ 
  int j_x,j_y;
  double temp_x,temp_y;
  double jr_x,jr_y;
  for (j_y= 0; j_y<J_Y; j_y++){
      for(j_x=0; j_x<J_X; j_x++){
        n[j_y+j_x*J_Y]=0.0;  
      } 
   }
   
   for (int i=0;i<NSP;i++)
    {
       jr_x=pos_x[i]/hx;       // indice (real) de la posición de la superpartícula
       j_x =int(jr_x);       // indice  inferior (entero) de la celda que contiene a la superpartícula
       temp_x = jr_x-j_x;
       jr_y=pos_y[i]/hx;       // indice (real) de la posición de la superpartícula
       j_y =int(jr_y);       // indice  inferior (entero) de la celda que contiene a la superpartícula
       temp_y = jr_y-j_y;

       n[j_y+j_x*J_Y] += (1.-temp_x)*(1.-temp_y)/(hx*hx*hx);
       n[j_y+(j_x+1)*J_Y] += temp_x*(1.-temp_y)/(hx*hx*hx);
       n[(j_y+1)+j_x*J_Y] += (1.-temp_x)*temp_y/(hx*hx*hx);
       n[(j_y+1)+(j_x+1)*J_Y] += temp_x*temp_y/(hx*hx*hx);

    }
    
  }
    
    
////////////////////////////////////////////////////////////////////////     
int main(){ __d

 //****************************************
  // Inicialización de variables del sistema
  //****************************************
  int i,j;__d
  int size = max_SPe*sizeof(double);__d
  int size1 = J_X*J_Y*sizeof(double);__d

  double  pos_eCarlos[max_SPe][2];__d
  double  neCarlos[J_X][J_Y];__d
  double *pos_e_xYen;__d
  double *pos_e_yYen;__d
  double *neYen; __d
  
  
 pos_e_xYen = (double *) malloc(size);
 pos_e_yYen = (double *) malloc(size);
 neYen = (double *) malloc(size1);
 
 
 
 
 FILE *fp, *fp1, *fp2, *fp3;
	fp = fopen ( "./posicionesXYPrueba", "r" );        
	if (fp==NULL){
		fputs ("File error leer concetracion1",stderr); exit (1);
	}

	for(i=0; i<max_SPe;i++){
		fscanf(fp,"%lf %lf\n",&pos_eCarlos[i][0], &pos_eCarlos[i][1]);
	}

	fclose ( fp );

	fp3 = fopen ( "./posicionesXYPrueba", "r" );        
	if (fp3==NULL){
		fputs ("File error concetracion2",stderr); exit (1);
	}

	for(i=0; i<max_SPe; i++){
		
			fscanf(fp3,"%lf %lf\n",&pos_e_xYen[i], &pos_e_yYen[i]);
	}
	fclose(fp3);
	
	Concentration_Carlos(pos_eCarlos, neCarlos,10);__d
	Concentration_Yen (pos_e_xYen, pos_e_yYen, neYen,10);__d


	fp1 = fopen ( "concetracionYen", "w" );      
	if (fp1==NULL){
		fputs ("File error escribiendo concetracion1",stderr); exit (1);
	}

	fp2 = fopen ( "concetracionCarlos", "w" );      
	if (fp2==NULL){
		fputs ("File error escribiendo concetracion2",stderr); exit (1);
	}

	for(i=0; i<J_X; i++){
		for(j=0; j<J_Y; j++)
		{
			fprintf(fp1," %lf\n",neYen[i*J_Y+j]);
		}
		
	}
	
	for(i=0; i<J_X; i++){
		for(j=0; j<J_Y; j++)
		{
			fprintf(fp2," %lf\n",neCarlos[i][j]);
		}
		
	}



	fclose (fp1);
	fclose (fp2);
  free(pos_e_xYen);
  free(pos_e_yYen);
  free(neYen);
  return 0; 
}
