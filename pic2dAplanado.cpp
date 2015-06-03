
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <complex>
#include <ctime>
#include <cstring>
#include <fstream>
#include <fftw3.h>

#define __d {prueba (__LINE__);}

//Magnitudes físicas

#define m_e         9.10938291e-31  // Masa del Electrón
#define e           1.6021e-19      // Carga del Electrón
#define k_b         1.3806504e-23   // Constante de Boltzmann
#define epsilon_0   8.854187e-12    // Permitividad eléctrica del vacío

#define max_SPe     10000           // Limite (computacional) de Superpartículas electrónicas
#define max_SPi     10000           // Limite (computacional) de Superpartículas iónicas
#define J_X         64         // Número de puntos de malla X. Recomendado: Del orden 2^n+1
#define J_Y         16           // Número de puntos de malla Y. Recomendado: Del orden 2^n

using namespace std;

void prueba(int d){
        cout<<"LINE: "<< d <<endl;
}

double RANDMAX;

void initialize_Particles(double *pos_e_x,double *pos_e_y,double *pos_i_x, double *pos_i_y, double *vel_e_x, double *vel_e_y, double *vel_i_x, double *vel_i_y, int li , int le);
double  create_Velocities_X(double fmax, double vphi);
double  create_Velocities_Y(double fmax, double vphi);
void Concentration (double *pos_x, double *pos_y, double *n,int NSP);
void Poisson2D_DirichletX_PeriodicY(double *phi,complex <double> *rho);
void Electric_Field (double *phi, double *E_X, double *E_Y);
void Motion(double *pos_x, double *pos_y,double *vel_x, double *vel_y, int NSP, int especie, double *E_X, double *E_Y);

int electrones=0;
int Iones = 1;
int X=0, Y=1;

//************************
// Parámetros del sistema
//************************

double  razon_masas = 1.98e5;     // m_i/m_e (Plata)
double  m_i = razon_masas*m_e;    // masa Ión
double  flujo_inicial = 4.5e33;   // Flujo inicial (# part/m^2*s)
double  vflux_i_x = 1e3;  // Componentes de Velocidad de flujo iónico
double  vflux_i_y = 1e3;
double  vflux_i_magnitud = sqrt(vflux_i_x*vflux_i_x+vflux_i_y*vflux_i_y); // Velocidad de flujo iónico (m/s) = sqrt(2*k_b*Te/(M_PI*m_i))
double  vflux_e_x = sqrt(razon_masas)*vflux_i_x;
double  vflux_e_y = sqrt(razon_masas)*vflux_i_y;
double  vflux_e_magnitud = sqrt(vflux_e_x*vflux_e_x+vflux_e_y*vflux_e_y);
double  ni03D = flujo_inicial/vflux_i_magnitud; // Concentración inicial de iones (3D)
double  ne03D=flujo_inicial/vflux_e_x; // Concentración inicial de electrones (3D)
double  Te=M_PI*0.5*m_e*pow(vflux_e_x,2)/k_b;   // Temperatura electrónica inicial (°K)                    
                                                // (vflujo=sqrt(2k_bTe/(pi*me)
double  lambda_D = sqrt(epsilon_0*k_b*Te/(ne03D*pow(e,2)));  //Longitud de Debye
double  om_p = vflux_e_x/lambda_D;                    //Frecuencia del plasma
double  ND = ne03D*pow(lambda_D,3);                          //Parámetro del plasma
int     NTe = 1e5; 
int     NTI = 1e5;                                  //Número de partículas "reales"

//***************************
//Constantes de normalización
//***************************

double  t0=1e-13;                   //Escala de tiempo: Tiempo de vaporización
double  x0=vflux_i_x*t0;   //Escala de longitud: Distancia recorrida en x por un ión en el tiempo t_0
//double  x0=lambda_D;                //Escala de longitud: Longitud de Debye
double  n0=double(NTe)/(x0*x0*x0);


//************************
//Parámetros de simulación
//************************

double  delta_X=lambda_D;   //Paso espacial
double  Lmax_x = (J_X-1)*delta_X;
double  Lmax_y = (J_Y-1)*delta_X;
//double  Lmax[2]={(J_X-1)*delta_X, (J_Y-1)*delta_X}; //Tamaño del espacio de simulación.
int     Factor_carga_e=10, Factor_carga_i=10;       //Número de partículas por superpartícula.
int     k_max_inj;   //Tiempo máximo de inyección
int     K_total;     //Tiempo total
int     Ntv=8;
int     le=0, li=0,kt;
int     NTSPe, NTSPI, max_SPe_dt, max_SPi_dt;
double  L_max_x, L_max_y, vphi_i_x, vphi_i_y, vphi_e_x, vphi_e_y,
        fi_Maxwell_x, fi_Maxwell_y, fe_Maxwell_x, fe_Maxwell_y;

//double  L_max[2], vphi_i[2],vphi_e[2], fi_Maxwell[2],fe_Maxwell[2];
double  T,dt,t_0, ne0_3D,ni0_3D,ni0_1D,ne0_1D;
double  vphi_i_magnitud,vphi_e_magnitud ,vpart,x_0,phi_inic;
double  cte_rho=pow(e*t0,2)/(m_i*epsilon_0*pow(x0,3)); //Normalización de epsilon_0
double  phi0=2.*k_b*Te/(M_PI*e), E0=phi0/x0;
//double  cte_E=t0*e*E0/(vflux_i[X]*m_e),fact_el=-1, fact_i=1./razon_masas;
double  cte_E=razon_masas*x0/(vflux_i_x*t0),fact_el=-1, fact_i=1./razon_masas;


int     total_e_perdidos=0;
int     total_i_perdidos=0;
double  mv2perdidas=0;

FILE    *outPot19,*outEnergia, *outPot0_6, *outPot0_9, *outPot1_5, *outPot3_5, *outPot5_5, *outPot15;
FILE    *outFase_ele[81];
FILE    *outFase_ion[81];
FILE    * outPoisson;

double hx;


int main()
{

  int seed = time (NULL); srand (seed);  // Semilla para generar números aleatorios dependiendo del reloj interno.


//////////////////////////////////////////////////////////////////////////////////////

//******************
  //ARCHIVOS DE SALIDA
  //******************
  outEnergia=fopen("Energia","w");

  
  char buffer[40];
  for(int i = 0; i<=80; i++){
        sprintf(buffer,"./outputs/fase_ele%d", i);
        outFase_ele[i]=fopen(buffer,"w");
  }
  
  for(int i = 0; i<=80; i++){
        sprintf(buffer,"./outputs/fase_ion%d", i);
        outFase_ion[i]=fopen(buffer,"w");
  }


  printf("ni03D=%e \nne03D=%e \nTemperatura electronica=%e eV \nLongitud de Debye=%e  \nFrecuencia del plasma=%e \nTp=%e \nND=%e \nLX=%e \nLY=%e \n",ni03D,ne03D,Te*k_b/(1.602e-19),lambda_D,om_p,2*M_PI/om_p,ND,Lmax_x,Lmax_y);

  printf("cte_E=%e  \ncte_rho=%e  \nTe = %e  \nhx*Ld = %e  \n",cte_E,cte_rho, Te, hx*lambda_D );
  
 int size = max_SPe*sizeof(double);
 int size1 = J_X*J_Y*sizeof(double);
 int size2 = J_X*J_Y*sizeof(complex<double>);
  
 // incializacion variables 
 
 double *pos_e_x, *pos_e_y, *pos_i_x, *pos_i_y;
 double *vel_e_x, *vel_e_y, *vel_i_x, *vel_i_y; 
 double *ne, *ni; 
 
 complex <double> *rho; 
 double *phi, *E_x, *E_y; 
 double E_i, E_e, E_field,E_total,E_perdida;
 
 
 // Asignacion de memoria con malloc
 pos_e_x = (double *) malloc(size);
 pos_e_y = (double *) malloc(size);
 pos_i_x = (double *) malloc(size);
 pos_i_y = (double *) malloc(size);
 vel_e_x = (double *) malloc(size);
 vel_e_y = (double *) malloc(size);
 vel_i_x = (double *) malloc(size);
 vel_i_y = (double *) malloc(size);
 ni = (double *) malloc(size1);
 ne = (double *) malloc(size1);
 E_x = (double *) malloc(size1);
 E_y = (double *) malloc(size1);
 rho = (complex<double> *) malloc(size2);
 phi = (double *) malloc(size1);
 
  //***************************
  // Normalización de variables
  //***************************

  L_max_x=Lmax_x/x0;                      // Longitud región de simulación
  L_max_y=Lmax_y/x0;                      // Longitud región de simulación
  t_0=1;
  x_0=1;
  hx=delta_X/x0;                            // Paso espacial
  double n_0=n0*x0*x0*x0;                   // Densidad de partículas
  dt=1.e-5;                                 // Paso temporal
  ni0_3D=ni03D*pow(x0,3);                   // Concentración de iones inicial 3D 
  ne0_3D=ne03D*pow(x0,3);                   // Concentración de electrones inicial 3D
  vphi_i_x=vflux_i_x/vflux_i_x;    // Velocidad térmica Iónica (X)
  vphi_e_x=vflux_e_x/vflux_i_x;    // Velocidad térmica Electrónica (X)
  vphi_i_y=vflux_i_y/vflux_i_x;    // Velocidad térmica Iónica (Y)
  vphi_e_y=vflux_e_y/vflux_i_x;    // Velocidad térmica Electrónica (Y)
  fi_Maxwell_x=  (2./(M_PI*vphi_i_x));    // Valor Máximo de la función de distribución Semi-Maxwelliana Iónica (X)
  fe_Maxwell_x=  (2./(M_PI*vphi_e_x));    // Valor Máximo de la función de distribución Semi-Maxwelliana elecrónica
  fi_Maxwell_y=  (1./(M_PI*vphi_i_y));    // Valor Máximo de la función de distribución Semi-Maxwelliana Iónica
  fe_Maxwell_y=  (1./(M_PI*vphi_e_y));    // Valor Máximo de la función de distribución Semi-Maxwelliana electrónica
  NTSPe=NTe/Factor_carga_e;
  NTSPI=NTI/Factor_carga_i; // Número total de superpartículas
  
 
   printf("x0^3=%e \nn0i=%e \nlambda/hx=%e \nTemp = %e\n", x0*x0*x0, ni03D, lambda_D/x0,k_b*Te);

  int Kemision=20;  //Pasos para liberar partículas

  double dt_emision=Kemision*dt; //Tiempo para liberar partículas

  max_SPe_dt=NTSPe*dt_emision;   //Número de Superpartículas el. liberadas cada vez.   
  max_SPi_dt=max_SPe_dt;        


  // Ciclo de tiempo

  k_max_inj=t_0/dt;
  K_total=Ntv*k_max_inj;
  int K_total_leapfrog=2*K_total;

  int kk=0;  
 
  initialize_Particles (pos_e_x,pos_e_y,pos_i_x,pos_i_y,vel_e_x,vel_e_y, vel_i_x, vel_i_y, li,le); //Velocidades y posiciones iniciales de las part´iculas (no las libera). 
  

 ofstream init;
 init.open("posicionesXYPrueba");//se escribe un archivo de salida para analizar los datos. la salida corresponde al potencial electrostatico en cada celda conocido como phi.
    for (int i = 0; i < max_SPe; i++){
      init<<pos_e_x[i]<<" "<<pos_i_x[i]<<" "<<pos_e_y[i]<<" "<<pos_i_y[i]<<"\n";
      }
init<<endl;

int test = 0;      
        
  clock_t tiempo0 = clock();
  for (kt=0;kt<=K_total;kt++)
  {
    if (kt%10000==0)
    {
      printf("kt=%d\n",kt);
      printf("le=%d   li=%d \n",le, li );
    }

    if(kt<=k_max_inj && kt==kk)                   // Inyectar superpartículas (i-e)
    {        
      le+=max_SPe_dt;
      li+=max_SPi_dt;
      kk=kk+Kemision;  
    }

    //-----------------------------------------------
    // Calculo de "densidad de carga 2D del plasma"

    Concentration (pos_e_x, pos_e_y, ne, le);           // Calcular concentración de superpartículas electrónicas
    Concentration (pos_i_x, pos_i_y, ni, li);           // Calcular concentración de superpartículas Iónicas
  
    for (int i = 0; i < J_Y; i++) 
    {
        for (int j = 0; j < J_X; j++) 
        {
          rho[i+j*J_Y]= cte_rho* Factor_carga_e*(ni[i+j*J_Y] - ne[i+j*J_Y])/n_0;
        }
    }
   
     // Calcular potencial eléctrico en puntos de malla
    Poisson2D_DirichletX_PeriodicY(phi,rho);

    // Calcular campo eléctrico en puntos de malla
    Electric_Field(phi,E_x, E_y);
    
     // imprimir el potencial electroestatico.
     
     if(kt%50000==0){
       sprintf(buffer,"./outputs/Poisson%d.data", kt);
       ofstream dataFile(buffer);
       for (int j = 0; j < J_X; j++) {
          double thisx = j * hx;
          for (int k = 0; k < J_Y; k++) {
              double thisy = k * hx;
              dataFile << thisx << '\t' << thisy << '\t' << phi[j+k*J_X] << '\n';
          }
          dataFile << '\n';
      }
      dataFile.close(); 
      }
      
      //imprimit la densidad.}
      if(kt%50000==0){
      // Escribir a archivo
      sprintf(buffer,"./outputs/n%d.data", kt);
      ofstream dataFile(buffer);
      for (int j = 0; j < J_X; j++) {
          double thisx = j * hx;
          for (int k = 0; k < J_Y; k++) {
              double thisy = k * hx;
              dataFile << thisx << '\t' << thisy << '\t' << ni[j+k*J_X] << '\t'<< ne[j+k*J_X] << '\t' << E_x[j+k*J_X]<< '\t' << E_y[j+k*J_X] <<'\n';
          }
          dataFile << '\n';
      }
      dataFile.close();
 
      }

    Motion(pos_e_x, pos_e_y, vel_e_x, vel_e_y,le, electrones, E_x,E_y);
    Motion(pos_i_x,pos_i_y,vel_i_x, vel_i_y, li, Iones, E_x, E_y);

  } //Cierre del ciclo principal
  

 free(pos_e_x); 
 free(pos_e_y); 
 free(pos_i_x); 
 free(pos_i_y); 
 free(vel_e_x); 
 free(vel_e_y);   
 free(vel_i_x);
 free(vel_i_y); 
 free(ne);   
 free(ni); 
 free(rho); 
 free(phi);
 free(E_x);
 free(E_y);
 
 
  return 0;   
}
   
double create_Velocities_X(double fmax,double vphi) // función para generar distribución semi-maxwelliana de velocidades de las particulas 
                                             // (Ver pág. 291 Computational Physics Fitzpatrick: Distribution functions--> Rejection Method)
{ 
                                          
  double sigma=vphi;                           // sigma=vflujo=vth    ( "dispersión" de la distribución Maxweliana)
  double vmin= 0. ;                            // Rapidez mínima  
  double vmax= 4.*sigma;                       // Rapidez máxima
  double v,f,f_random;

  
  v=vmin+(vmax-vmin)*double(rand())/double(RAND_MAX); // Calcular valor aleatorio de v uniformemente distribuido en el rango [vmin,vmax]
  f =fmax*exp(-(1.0/M_PI)*pow(v/vphi,2));     //                       
 

  f_random = fmax*double(rand())/double(RAND_MAX);    // Calcular valor aleatorio de f uniformemente distribuido en el rango [0,fmax]

  if (f_random > f) return create_Velocities_X(fmax,vphi);
  else return  v;
}

double create_Velocities_Y(double fmax,double vphi) // función para generar distribución semi-maxwelliana de velocidades de las particulas 
                                             // (Ver pág. 291 Computational Physics Fitzpatrick: Distribution functions--> Rejection Method)
{                                            
  double sigma=vphi;                           // sigma=vflujo=vth    ( "dispersión" de la distribución Maxweliana)
  double vmin= -3.*sigma;                            // Rapidez mínima  
  double vmax= 3.*sigma;                       // Rapidez máxima
  double v,f,f_random;

   
  v=vmin+(vmax-vmin)*double(rand())/double(RAND_MAX); // Calcular valor aleatorio de v uniformemente distribuido en el rango [vmin,vmax]
  f =fmax*exp(-(1.0/M_PI)*pow(v/vphi,2));                      
 
  f_random = fmax*double(rand())/double(RAND_MAX);    // Calcular valor aleatorio de f uniformemente distribuido en el rango [0,fmax]

  if (f_random > f) return create_Velocities_Y(fmax,vphi);
  else return  v;
}



void initialize_Particles (double *pos_e_x, double *pos_e_y, double *pos_i_x, double *pos_i_y, double *vel_e_x , double *vel_e_y,  double *vel_i_x, double *vel_i_y, int le, int li) 
{  
   for (int i=0;i<max_SPe;i++)
   {
     pos_e_x[i+le]=0; 
     vel_e_x[i+le]= create_Velocities_X (fe_Maxwell_x,vphi_e_x);
     pos_e_y[i+le]=L_max_y/2.0;
     vel_e_y[i+le]=create_Velocities_Y(fe_Maxwell_y,vphi_e_y);
     pos_i_x[i+li]=0; 
     vel_i_x[i+li]=create_Velocities_X (fi_Maxwell_x,vphi_i_x);
     pos_i_y[i+li]=L_max_y/2.0;
     vel_i_y[i+li]=create_Velocities_Y (fi_Maxwell_y,vphi_i_y);
   }
}
/***********************************************************************************
Determinación del aporte de carga de cada superpartícula sobre las 4 celdas adyacentes
************************************************************************************/

void Concentration (double *pos_x, double *pos_y, double *n,int NSP)
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
       n[+(j_y+1)+j_x*J_Y] += (1.-temp_x)*temp_y/(hx*hx*hx);
       n[+(j_y+1)+(j_x+1)*J_Y] += temp_x*temp_y/(hx*hx*hx);

    }
   
}



////////////////////////////////////////////////////////////////////////////////////
void Poisson2D_DirichletX_PeriodicY(double *phi, complex<double> *rho)
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Electric_Field(double *phi, double *E_X, double *E_Y)  // Función para calcular el campo eléctrico 
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
//***************************************************************************************************************************************************************


void  Motion(double *pos_x, double *pos_y, double *vel_x, double *vel_y, int NSP, int especie, double *E_X, double *E_Y){
   int j_x,j_y;
   double temp_x,temp_y,Ep_X, Ep_Y,fact;
   double jr_x,jr_y;
   int kk1=0;
   int conteo_perdidas=0;
   
   if(especie==electrones)
    fact=fact_el;                    
    
   if(especie==Iones)
    fact=fact_i;
    
   for (int i=0;i<NSP;i++)
    {
       jr_x=pos_x[i]/hx;     // Índice (real) de la posición de la superpartícula (X)
       j_x =int(jr_x);        // Índice  inferior (entero) de la celda que contiene a la superpartícula (X)
       temp_x = jr_x-double(j_x);
       jr_y=pos_y[i]/hx;     // Índice (real) de la posición de la superpartícula (Y)
       j_y =int(jr_y);        // Índice  inferior (entero) de la celda que contiene a la superpartícula (Y)
       temp_y = jr_y-double(j_y);


       Ep_X=(1-temp_x)*(1-temp_y)*E_X[j_x*J_Y+j_y]+
       temp_x*(1-temp_y)*E_X[(j_x+1)*J_Y+j_y]+
       (1-temp_x)*temp_y*E_X[j_x+(j_y+1)*J_X]+
       temp_x*temp_y*E_X[(j_x+1)+(j_y+1)*J_X];
       
       Ep_Y=(1-temp_x)*(1-temp_y)*E_Y[j_x*J_Y+j_y]+
       temp_x*(1-temp_y)*E_Y[(j_x+1)*J_Y+j_y]+
       (1-temp_x)*temp_y*E_Y[j_x*J_Y+(j_y+1)]+
       temp_x*temp_y*E_Y[(j_x+1)*J_Y+(j_y+1)];
       
       pos_x[i]+=vel_x[i]*dt;
       pos_y[i]+=vel_y[i]*dt;

       vel_x[i]=vel_x[i]+cte_E*fact*Ep_X*dt;
       //vel[i][Y]=vel[i][Y]+cte_E*fact*Factor_carga_e*Ep_Y*dt;


       if(pos_x[i]<0) //Rebote en la pared del material.
       {
          pos_x[i]=-pos_x[i];
          vel_x[i]=-vel_x[i];
       }

       if (pos_x[i]>=L_max_x) //Partícula fuera del espacio de Simulación
       {    
            conteo_perdidas++;
            if(especie==electrones)
            {
              total_e_perdidos++;
              printf("Electron perdido No. %d,  i=%d, kt=%d \n",total_e_perdidos, i ,kt);
              mv2perdidas+=pow( sqrt(vel_x[i]*vel_x[i]+vel_y[i]*vel_y[i]) , 2);
            }
            else
            {
              total_i_perdidos++;
              printf("Ion perdido No. %d,  i=%d, kt=%d \n",total_i_perdidos, i ,kt);
              mv2perdidas+=pow( sqrt(vel_x[i]*vel_x[i]+vel_y[i]*vel_y[i]) , 2)/(razon_masas);
            }
       }



       while(pos_y[i]>L_max_y) //Ciclo en el eje Y.
       {
          pos_y[i]=pos_y[i]-L_max_y;
       }

       while(pos_y[i]<0.0) //Ciclo en el eje Y.
       {

          pos_y[i]=L_max_y+pos_y[i];
       }


       if(pos_x[i]>=0 && pos_x[i]<=L_max_x)
        {
            kk1=kk1+1;
            pos_x[kk1-1]=pos_x[i]; 
            pos_y[kk1-1]=pos_y[i];
            vel_x[kk1-1]=vel_x[i]; 
            vel_y[kk1-1]=vel_y[i];
        }

       //Salida espacio de Fase

       if(kt%10000==0&&especie==electrones)
          {
              fprintf(outFase_ele[kt/10000]," %e   %e  %e  %e  %e \n",kt*dt,pos_x[i],vel_x[i],pos_y[i],vel_y[i]);
              //printf("Fase out\n");
          }
       if(kt % 10000 == 0 && especie == Iones)
          {
              fprintf(outFase_ion[kt/10000]," %e   %e  %e  %e  %e \n", kt*dt ,pos_x[i], vel_x[i], pos_y[i], vel_y[i]);

          }

    }

    if(especie==electrones)
    {
      le=le-conteo_perdidas;
    }
    else
    {
      li=li-conteo_perdidas;
    }
}

//////////////////////////////////////////////////////////////////////////////////











