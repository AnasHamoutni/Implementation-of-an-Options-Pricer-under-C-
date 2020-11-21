#include<iostream>
#include<stdlib.h>
#include<stdio.h>
#include<conio.h>
#include<math.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<time.h>
#include<gsl/gsl_sf_erf.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_spline.h>
using namespace std;

double implicite(double S, double K, double r, double Vol, double theta)
{

int M=1000,N=100;

//on alloue la memoire pour les matrices et les vecteurs qu'on va utiliser
gsl_matrix * A = gsl_matrix_alloc (N,N);
gsl_matrix * I = gsl_matrix_alloc (N,N);
gsl_matrix * U = gsl_matrix_alloc (N,M+2);
gsl_matrix * F = gsl_matrix_alloc (N,M+2);
gsl_vector * U0 = gsl_vector_alloc (N);
gsl_vector * X = gsl_vector_alloc (N);
gsl_vector * F1 = gsl_vector_alloc (N);
gsl_vector * impl = gsl_vector_alloc (N+2);
//calcule des constantes
double q=2*r/(Vol*Vol);
double a=-1;
double b=1;
double k=0.5*Vol*Vol*theta/(M+1);//pas de discretisation du temps
double h=(b-a)/(N+1);//pas de discretisation de l'espace

    
		
//boucle pour remplir la matrice A
 for(int i=0;i<N;i++)   {	
		for(int j=0;j<N;j++){
			if (i==j){gsl_matrix_set(A,i,j,2);
			}
			else if(i==j+1 || i==j-1){gsl_matrix_set(A,i,j,-1);
			}
			else {gsl_matrix_set(A,i,j,0);
			}
		}}
		

 gsl_matrix_set_identity(I);//declaration de la matrice identité
		
//boucle pour remplir le vecteur U0 
for(int i=0;i<N;i++)   {

gsl_vector_set(U0,i,fmax(exp(0.5*(q+1)*(a+i*h))-exp(0.5*(q-1)*(a+i*h)),0));

 }

//boucle pour remplir la matrice F
for(int j=0;j<M+2;j++)   {
	
		for(int i=0;i<N;i++){
		if (i==N-1){
			
		gsl_matrix_set(F,i,j,exp(0.5*(q+1)*b+0.25*(q+1)*(q+1)*j*k)/(h*h));
		
		}else{
			
		gsl_matrix_set(F,i,j,0);
		
		}
		}
		}

//inserer UO (la premiere colonne de la matrice U) dans U
for(int j=0;j<M+2;j++)   {
		for(int i=0;i<N;i++){
			
			if (j==0){
				
			gsl_matrix_set(U,i,j,gsl_vector_get(U0,i));
			
			}
			else {
				
		    gsl_matrix_set (U, i, j,0);
		    
			}
			}
			}
			
			
	//LES ETAPES DE L ALGORITHME//////////////////////////////		
			
//extraction de F1
for(int i=0;i<N;i++){
	
	gsl_vector_set(F1,i,gsl_matrix_get(F,i,1));
	
	}	
	
gsl_vector_scale(F1,k);	//multiplier F1 par k et sauver le resultat dans le vecteur F1

       
gsl_vector_add(U0,F1);//ajouter F1 a U0

double cste=k/(h*h);

gsl_matrix_scale(A,cste);//multiplier la matrice A par le scalaire cste et sauver le resultat dans la matrice A 

gsl_matrix_add(I,A);//ajouter A a I

//trouver l'inverse de I via la decomposition LU
gsl_permutation *p=gsl_permutation_alloc(N);
int s;
gsl_matrix *inv=gsl_matrix_alloc(N,N);
gsl_linalg_LU_decomp(I,p,&s);
gsl_linalg_LU_invert(I,p,inv);//l'inverse de I est stocké dans inv

for(int j=1;j<=M+1;j++)   {//boucle pour trouver tous les vecteurs de la matrice U

        gsl_blas_dgemv(CblasNoTrans,1,inv,U0,0,X);//trouver le nouveau vecteur de la matrice U et le stocker dans X
        	
		for(int i=0;i<N;i++){
			
			gsl_matrix_set(U,i,j,gsl_vector_get(X,i)); //inserer le nouveau vecteur dans la matrice U
			
			if(j<=M){
				
			gsl_vector_set(F1,i,gsl_matrix_get(F,i,j+1));//extraire le vecteur suivant de la matrice F
			
			}
			
			}
			
			gsl_vector_memcpy(U0,X);//copier X dans U0
			
			gsl_vector_scale(F1,k);//multiplier F1 par k
			
			gsl_vector_add(U0,F1);	//ajouter F1 a U0

			
			}	
		

gsl_vector * UMplus1= gsl_vector_alloc (N+2);	//declaration du vecteur U(M+1), il ne correspond pas au dernier vecteur de U car il n'ont pas la meme dimension. mais le dernier vecteur de U contient les elements de U(M+1) a partir de lindice 1 jusqu'a N, donc on va  ajouter deux elements correspondants aux indices 0 et N+1.

		for(int i=0;i<N+2;i++){
			
		if(i==0){
		
		gsl_vector_set(UMplus1,i,0);
		
		 }else if(i==N+1){
		 	
		 gsl_vector_set(UMplus1,i,exp(0.5*(q+1)*b+0.25*(q+1)*(q+1)*(M+1)*k));
		 
		  }
		  
		else {gsl_vector_set(UMplus1,i,gsl_matrix_get(U,i-1,M+1) );
		
		} 
		
		}
			
		for(int i=0;i<N+2;i++){//boucle pour remplir le vecteur des solutions du schema implicite partir de la transformation du vecteur U(M+1)
			
				gsl_vector_set(impl,i,(gsl_vector_get(UMplus1,i)*K*exp(-0.5*(q-1)*(a+i*h)-0.25*(q+1)*(q+1)*0.5*Vol*Vol*theta)));
				
				}

//interpolation pour trouver la valeur implicite au point desiré log(S/K)
double X0[N+2],Y0[N+2] ,val;    
   
gsl_interp_accel *acc   = gsl_interp_accel_alloc ();

gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline,N+2);  

for(int j=0;j<N+2;j++){//Remplir X0 et Y0
	
      X0[j]=a+j*h;
      
      Y0[j]=gsl_vector_get(impl,j);   
	  
	  }
	  
	    gsl_spline_init (spline,X0,Y0,N+2); 
	    
	    val=gsl_spline_eval(spline,log(S/K),acc);//l'evaluation au point log(S/K)
	    
	     printf("la solution du schema implicite est %f\n\n", val);//afficher la solution implicite au point log(S/K)

                     
/////////////////////on vide la memoire allouée 	
gsl_matrix_free (inv);
gsl_matrix_free (A);
gsl_matrix_free (I);
gsl_vector_free (U0);
gsl_matrix_free (U);
gsl_matrix_free (F);
gsl_vector_free (F1);
gsl_vector_free (X);
gsl_permutation_free (p);
gsl_spline_free (spline);
gsl_interp_accel_free (acc);


return val;//retourner la solution implicite au point log(S/K)

}


double explicite(double S, double K, double r, double Vol, double theta)
{


int M=1000,N=100;

//on alloue la memoire pour les matrices et les vecteurs qu'on va utiliser
gsl_matrix * A = gsl_matrix_alloc (N,N);
gsl_matrix * I = gsl_matrix_alloc (N,N);
gsl_matrix * U = gsl_matrix_alloc (N,M+2);
gsl_matrix * F = gsl_matrix_alloc (N,M+2);
gsl_vector * U0 = gsl_vector_alloc (N);
gsl_vector * X = gsl_vector_alloc (N);
gsl_vector * F0 = gsl_vector_alloc (N);
gsl_vector * expli = gsl_vector_alloc (N+2);
//calcule des constantes
double q=2*r/(Vol*Vol);
double a=-1;
double b=1;
double k=0.5*Vol*Vol*theta/(M+1);//pas de discretisation du temps
double h=(b-a)/(N+1);//pas de discretisation de l'espace

    
		
//boucle pour remplir la matrice A
 for(int i=0;i<N;i++)   {	
 
		for(int j=0;j<N;j++){
			
			if (i==j){gsl_matrix_set(A,i,j,2);
			
			}else if(i==j+1 || i==j-1){
			
			gsl_matrix_set(A,i,j,-1);
			
			}else {gsl_matrix_set(A,i,j,0);
			}
		}}
		

 gsl_matrix_set_identity(I);//declaration de la matrice identité
			
//boucle pour remplir le vecteur U0 
for(int i=0;i<N;i++)   {

gsl_vector_set(U0,i,fmax(exp(0.5*(q+1)*(a+i*h))-exp(0.5*(q-1)*(a+i*h)),0));

 }

//boucle pour remplir la matrice F
for(int j=0;j<M+2;j++)   {
	
		for(int i=0;i<N;i++){
		if (i==N-1){
			
		gsl_matrix_set(F,i,j,exp(0.5*(q+1)*b+0.25*(q+1)*(q+1)*j*k)/(h*h));
		
		}else{
			
		gsl_matrix_set(F,i,j,0);
		
		}
		}
		}

//inserer UO (la premiere colonne de la matrice U) dans U
for(int j=0;j<M+2;j++)   {
		for(int i=0;i<N;i++){
			
			if (j==0){
				
			gsl_matrix_set(U,i,j,gsl_vector_get(U0,i));
			
			}
			else {
				
		    gsl_matrix_set (U, i, j,0);
		    
			}
			}
			}
			
	//LES ETAPES DE L ALGORITHME//////////////////////////////		
			
//extraction de F0
for(int i=0;i<N;i++){
	
	gsl_vector_set(F0,i,gsl_matrix_get(F,i,0));
	
	}	
	
gsl_vector_scale(F0,k);	//multiplier F0 par k et sauver le resultat dans le vecteur F0

double cste=-k/(h*h);

gsl_matrix_scale(A,cste);//multiplier la matrice A par le scalaire cste et sauver le resultat dans la matrice A 

gsl_matrix_add(I,A);//ajouter A a I



for(int j=1;j<M+2;j++)   {//boucle pour trouver tous les vecteurs de la matrice U

        gsl_blas_dgemv(CblasNoTrans,1,I,U0,0,X);///trouver le nouveau vecteur de la matrice U et le stocker dans X
        	
        gsl_vector_add(X,F0);//ajouter F0 a X

		for(int i=0;i<N;i++){
			
			gsl_matrix_set(U,i,j,gsl_vector_get(X,i));//inserer le nouveau vecteur dans la matrice U

			gsl_vector_set(F0,i,gsl_matrix_get(F,i,j));//extraire le vecteur suivant de la matrice F
			
			}
			
			gsl_vector_memcpy(U0,X);//copier X dans U0
			
			gsl_vector_scale(F0,k);//multiplier F0 par k
			
			
			}
			
gsl_vector * UMplus1= gsl_vector_alloc (N+2);	//declaration du vecteur U(M+1), il ne correspond pas au dernier vecteur de U car il n'ont pas la meme dimension. mais le dernier vecteur de U contient les elements de U(M+1) a partir de lindice 1 jusqu'a N, donc on va  ajouter deux elements correspondants aux indices 0 et N+1.

		for(int i=0;i<N+2;i++){
			
		if(i==0){
		
		gsl_vector_set(UMplus1,i,0);
		
		 }else if(i==N+1){
		 	
		 gsl_vector_set(UMplus1,i,exp(0.5*(q+1)*b+0.25*(q+1)*(q+1)*(M+1)*k));
		 
		  }
		  
		else {gsl_vector_set(UMplus1,i,gsl_matrix_get(U,i-1,M+1) );
		
		} 
		
		}
		
		for(int i=0;i<N+2;i++){//boucle pour remplir le vecteur des solutions du schema implicite partir de la transformation du vecteur U(M+1)
			
				gsl_vector_set(expli,i,(gsl_vector_get(UMplus1,i)*K*exp(-0.5*(q-1)*(a+i*h)-0.25*(q+1)*(q+1)*0.5*Vol*Vol*theta)));
				
				}

   //interpolation pour trouver la valeur explicite au point desiré log(S/K)
double X0[N+2],Y0[N+2] ,val;    
   
gsl_interp_accel *acc   = gsl_interp_accel_alloc ();

gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline,N+2);  

for(int j=0;j<N+2;j++){//Remplir X0 et Y0
	
      X0[j]=a+j*h;
      
      Y0[j]=gsl_vector_get(expli,j);   
	  
	  }
	  
	    gsl_spline_init (spline,X0,Y0,N+2); 
	    
	    val=gsl_spline_eval(spline,log(S/K),acc);//l'evaluation au point log(S/K)
	    
	     printf("la solution du schema explicite est %f\n\n", val);//afficher la solution explicite au point log(S/K)

                     
                     
//vider la memoire allouée
gsl_matrix_free (A);
gsl_matrix_free (I);
gsl_vector_free (U0);
gsl_matrix_free (U);
gsl_matrix_free (F);
gsl_vector_free (F0);
gsl_vector_free (X);
gsl_spline_free (spline);
gsl_interp_accel_free (acc);

return val;//retourner la solution explicite au point log(S/K)
}


int main(){
	
	printf("la solution de Crank est %f\n\n ",(implicite(111,100,0.2,0.3,1)+explicite(111,100,0.2,0.3,1))/2);
}

