#include <math.h>
#include <stdio.h>
#include <stdlib.h>

float et(float gamma, float matrix0, float matrix1, float matrix2, float* result);
float Uv(float gamma, float matrix0, float matrix1, float matrix2, float* a,float* b,float* c);
float Fv(float gamma, float arreglo0, float arreglo1, float arreglo2, float *F0,float *F1,float *f2);

float metodo1(float gamma,int tam,float dx, float dt, float re0, float re1, float re2, float re0mas, float re1mas, float re2mas,float re0ceros,float re1ceros,float re2ceros, float *salida0, float *salida1, float *salida2);

float metodo2(float gamma,int tam,float dx, float dt,float re0, float re1, float re2,float retam0, float retam1,float retam2,float retmenos0, float retmenos1,float retmenos2, float *salida0,float *salida1,float *salida2);


float metodo3(float gamma,int tam,float dx, float dt,float re0, float re1, float re2, float re0mas, float re1mas, float re2mas,float retmenos0, float retmenos1,float retmenos2, float *salida0,float *salida1,float *salida2);

int main(int argc, char *argv[])
{ 

float gamma = 1.4;
float t = atof(argv[1]);
int ite =1000;
float dt = t/ite; 
int tam = 1000 ;
float dx = 20/tam;
int i=0.0;
int j=0.0;
int k;
float a,b,c,a1,b1,c1,a2,b2,c2=0.0;


float **matrix0 = malloc(ite*sizeof(float *)); 

float **matrix1 = malloc(ite*sizeof(float *)); 
float **matrix2 = malloc(ite*sizeof(float *));  
for(i=0;i<ite;i++)
    {
    matrix0[i] = (float*)malloc(ite*sizeof(float));
    matrix1[i] = (float*)malloc(ite*sizeof(float));
    matrix2[i] = (float*)malloc(ite*sizeof(float));
    }

//printf("%f\n",matrix0[0][0]) ;

for(i=0;i<tam;i++){
    if (i<=499){                   
        matrix1[i][0]=1.0;        
        matrix2[i][0]=100000.0;}
                                
    if (i>499){                  
        matrix1[i][0]=0.125;     
        matrix2[i][0]=10000.0;}
}

printf("%f\n",matrix0[0][0]) ;

for(i=0;i<tam;i++){
for(j=0;j<ite-1;j++){

if(i==0){

metodo1(gamma,tam,dx,dt,matrix0[i][j], matrix1[i][j],matrix2[i][j], matrix0[i+1][j],matrix1[i+1][j], matrix2[i+1][j],matrix0[0][0],matrix1[0][0],matrix2[0][0], &a,&b,&c);
matrix0[i][j+1] = a;
matrix1[i][j+1] = b;
matrix2[i][j+1] = c;
}
if(i==tam-1){
metodo2(gamma,tam,dx,  dt,matrix0[i][j], matrix1[i][j],matrix2[i][j],matrix0[tam-1][0], matrix1[tam-1][0],matrix2[tam-1][0],matrix0[i-1][j], matrix1[i-1][j],matrix2[i-1][j], &a1,&b1,&c1);
matrix0[i][j+1] = a1;
matrix1[i][j+1] = b1;
matrix2[i][j+1] = c1;

}
else{
metodo3(gamma,tam,dx, dt,matrix0[i][j], matrix1[i][j], matrix2[i][j], matrix0[i+1][j], matrix1[i+1][j], matrix2[i+1][j],matrix0[i-1][j], matrix1[i-1][j],matrix2[i-1][j], &a2,&b2,&c2);

matrix0[i][j+1] = a2;
matrix1[i][j+1] = b2;
matrix2[i][j+1] = c2;

}

}
}

FILE *out;
char filename[100];
sprintf(filename, "estado_%s.dat", t_string);
out = fopen(filename, "w");
for(i=0;i<tam;i++)
{
fprintf(out, "%f %f %f %f\n", dx*i-10.0,matrix0[i][0], matrix1[i][0], matrix2[i][0];
}
fclose(out);
}


return 0;
}

float et(float gamma,float matrix0, float matrix1, float matrix2, float *result)
{

 *result= matrix2/(matrix1*(gamma-1.0)) + (matrix0*matrix0)/2.0;

}

float Uv(float gamma,float matrix0, float matrix1, float matrix2, float *U0,float *U1,float *U2){
    float resultado;
    et(gamma,matrix0,matrix1,matrix2,&resultado);
    *U0 = matrix1;
    *U1 = matrix1*matrix0;
    *U2 = matrix1*resultado;
    
}

float Fv(float gamma,float arreglo0, float arreglo1, float arreglo2, float *F0,float *F1,float *F2){
       
    *F0= arreglo1;
    *F1 = (arreglo1*arreglo1)/arreglo0 + (gamma-1.0)*(arreglo2 - (arreglo1*arreglo1)/(2.0*arreglo0));
    *F2 = (arreglo1/arreglo0)*arreglo2 + (gamma-1.0)*arreglo2*(arreglo1/arreglo0) - (gamma-1.0)*(arreglo1/arreglo0)* (arreglo1*arreglo1)/(2.0*arreglo0);
    }

float metodo3(float gamma,int tam,float dx, float dt,float re0, float re1, float re2, float re0mas, float re1mas, float re2mas,float retmenos0, float retmenos1,float retmenos2, float *salida0, float *salida1,float *salida2){
   
   float U0,U1,U2,Utmenos0,Utmenos1,Utmenos2,U0mas,U1mas,U2mas=0.0;
   float utemppos0,utemppos1,utemppos2,utempneg0,utempneg1,utempneg2=0.0;
   float F0Umas,F1Umas,F2Umas,F0U,F1U,F2U,F0Umenos,F1Umenos,F2Umenos=0.0;
   float Ftemppos0,Ftemppos1,Ftemppos2,Ftempneg0,Ftempneg1,Ftempneg2=0.0;
    
   float unew0,unew1,unew2=0.0;

    

    Uv(gamma,re0, re1,re2, &U0,&U1,&U2);

    Uv(gamma,re0mas, re1mas, re2mas, &U0mas,&U1mas,&U2mas);

    
    Uv(gamma,retmenos0, retmenos1,retmenos2, &Utmenos0,&Utmenos1,&Utmenos2);
    
 
    Fv(gamma,U0mas,U1mas,U2mas, &F0Umas,&F1Umas,&F2Umas);
    Fv(gamma,U0,U1,U2, &F0U,&F1U,&F2U); 
  
    Fv(gamma,Utmenos0,Utmenos1,Utmenos2, &F0Umenos,&F1Umenos,&F2Umenos);
    
    
    
    
            utemppos0 = 0.5*( U0mas+U0  - (dt/dx)*( F0Umas - F0U ) );
            utemppos1 = 0.5*( U1mas+U1  - (dt/dx)*( F1Umas - F1U ) );
            utemppos2 = 0.5*( U2mas+U2  - (dt/dx)*( F2Umas - F2U ) );
 
            utempneg0 = 0.5*( Utmenos0+U0 - (dt/dx)*( F0Umenos - F0U) );
            utempneg1 = 0.5*( Utmenos1+U1 - (dt/dx)*( F0Umenos - F1U) );
            utempneg2 = 0.5*( Utmenos2+U2 - (dt/dx)*( F0Umenos - F2U) );
            Fv(gamma,utemppos0, utemppos1, utemppos2, &Ftemppos0,&Ftemppos1,&Ftemppos2);
            Fv(gamma,utempneg0, utempneg1, utempneg2, &Ftempneg0,&Ftempneg1,&Ftempneg2);

           
    
        unew0 = U0 - (dt/dx)*( Ftemppos0 - Ftempneg0 );
        unew1 = U1 - (dt/dx)*( Ftemppos1 - Ftempneg1 );
        unew2 = U2 - (dt/dx)*( Ftemppos2 - Ftempneg2 );
         *salida1     = unew0;
         *salida0 = unew1/unew0;
         *salida2 = (gamma-1.0)*( unew2 - ( unew0*(unew1/unew0)*(unew1/unew0) )/2.0);
       


  }
    
    

float metodo1(float gamma,int tam,float dx, float dt,float re0, float re1, float re2, float re0mas, float re1mas, float re2mas,float re0ceros,float re1ceros,float re2ceros, float *salida0,float *salida1,float *salida2){
   
   float U0,U1,U2,U0ceros,U1ceros,U2ceros, U0mas,U1mas,U2mas=0.0;
   float utemppos0,utemppos1,utemppos2,utempneg0,utempneg1,utempneg2=0.0;
   float F0Umas,F1Umas,F2Umas,F0U,F1U,F2U,F0Uceros,F1Uceros,F2Uceros=0.0;
   float Ftemppos0,Ftemppos1,Ftemppos2,Ftempneg0,Ftempneg1,Ftempneg2=0.0;
    
   float unew0,unew1,unew2=0.0;

    

    Uv(gamma,re0, re1,re2, &U0,&U1,&U2);
    Uv(gamma,re0ceros, re1ceros, re2ceros, &U0ceros,&U1ceros,&U2ceros);
    Uv(gamma,re0mas, re1mas, re2mas, &U0mas,&U1mas,&U2mas);    
 
    Fv(gamma,U0mas,U1mas,U2mas, &F0Umas,&F1Umas,&F2Umas);
    Fv(gamma,U0,U1,U2, &F0U,&F1U,&F2U);
    Fv(gamma,U0ceros,U1ceros,U2ceros, &F0Uceros,&F1Uceros,&F2Uceros);
    
       
    
            utemppos0= 0.5*( U0mas+U0 - (dt/dx)*( F0Umas- F0U ));
            utemppos1= 0.5*( U0mas+U1 - (dt/dx)*( F1Umas- F1U ));
            utemppos2= 0.5*( U0mas+U2 - (dt/dx)*( F2Umas- F2U ));
   
            utempneg0 = 0.5*( U0ceros+U0 - (dt/dx)*( F0Uceros - F0U  ) );
            utempneg1 = 0.5*( U0ceros+U1 - (dt/dx)*( F1Uceros - F1U  ) );
            utempneg2 = 0.5*( U0ceros+U2 - (dt/dx)*( F2Uceros - F2U  ) );
            Fv(gamma,utemppos0, utemppos1, utemppos2, &Ftemppos0,&Ftemppos1,&Ftemppos2);
            Fv(gamma,utempneg0, utempneg1, utempneg2, &Ftempneg0,&Ftempneg1,&Ftempneg2);
        
   
        unew0 = U0 - (dt/dx)*( Ftemppos0 - Ftempneg0 );
        unew1 = U1 - (dt/dx)*( Ftemppos1 - Ftempneg1 );
        unew2 = U2 - (dt/dx)*( Ftemppos2 - Ftempneg2 );
         *salida1     = unew0;
         *salida0 = unew1/unew0;
         *salida2 = (gamma-1.0)*( unew2 - ( unew0*(unew1/unew0)*(unew1/unew0) )/2.0);
       


  }
    
    

float metodo2(float gamma,int tam,float dx, float dt,float re0, float re1, float re2,float retam0, float retam1,float retam2,float retmenos0, float retmenos1,float retmenos2, float *salida0,float *salida1,float *salida2){
   
   float U0,U1,U2, Utam0,Utam1,Utam2,Utmenos0,Utmenos1,Utmenos2=0.0;
   float utemppos0,utemppos1,utemppos2,utempneg0,utempneg1,utempneg2=0.0;
   float F0U,F1U,F2U,F0Utam,F1Utam,F2Utam,F0Umenos,F1Umenos,F2Umenos=0.0;
   float Ftemppos0,Ftemppos1,Ftemppos2,Ftempneg0,Ftempneg1,Ftempneg2=0.0;
    
   float unew0,unew1,unew2=0.0;

    

    Uv(gamma,re0, re1,re2, &U0,&U1,&U2);
    Uv(gamma,retam0, retam1,retam2, &Utam0,&Utam1,&Utam2);
    Uv(gamma,retmenos0, retmenos1,retmenos2, &Utmenos0,&Utmenos1,&Utmenos2);
    
 
    
    Fv(gamma,U0,U1,U2, &F0U,&F1U,&F2U);   
    Fv(gamma,Utam0,Utam1,Utam2, &F0Utam,&F1Utam,&F2Utam);
    Fv(gamma,Utmenos0,Utmenos1,Utmenos2, &F0Umenos,&F1Umenos,&F2Umenos);
   
        
            utemppos0 = 0.5*( Utam0+U0 - (dt/dx)*( F0Utam - F0U ));
            utemppos1 = 0.5*( Utam1+U1 - (dt/dx)*( F1Utam - F1U ));
            utemppos2 = 0.5*( Utam2+U2 - (dt/dx)*( F2Utam - F2U ));

            utempneg0 = 0.5*( Utmenos0+U0 - (dt/dx)*(F0Umenos - F0U ));
            utempneg1 = 0.5*( Utmenos1+U1 - (dt/dx)*(F1Umenos - F1U ));
            utempneg2 = 0.5*( Utmenos2+U2 - (dt/dx)*(F2Umenos - F2U ));
        
            Fv(gamma,utemppos0, utemppos1, utemppos2, &Ftemppos0,&Ftemppos1,&Ftemppos2);
            Fv(gamma,utempneg0, utempneg1, utempneg2, &Ftempneg0,&Ftempneg1,&Ftempneg2);
         
    
        unew0 = U0 - (dt/dx)*( Ftemppos0 - Ftempneg0 );
        unew1 = U1 - (dt/dx)*( Ftemppos1 - Ftempneg1 );
        unew2 = U2 - (dt/dx)*( Ftemppos2 - Ftempneg2 );
         *salida1     = unew0;
         *salida0 = unew1/unew0;
         *salida2 = (gamma-1.0)*( unew2 - ( unew0*(unew1/unew0)*(unew1/unew0) )/2.0);


}
