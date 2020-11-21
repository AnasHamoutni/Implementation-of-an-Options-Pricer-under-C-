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
using namespace std;

void Callexact(double S, double K,double r, double vol, double theta){//la fonction CALL qui va calculer le prix du call selon la formule close de  Black & Scholes.
	if(S>0)
 { double standard_deviation= vol*sqrt(theta);//vol est la volatilité,theta est le temps a la maturité
       double d1=(log(S/K)+r*theta)/standard_deviation+0.5*standard_deviation;//S est le prix du sous jacent, K est le prix d'exercice, r est le taux d'interet instantanné.
       double d2=d1-standard_deviation;
	    double callexact= S*gsl_sf_erf_Q(-d1)-K*exp(-r*theta)*gsl_sf_erf_Q(-d2);
		printf("le prix du call par la formule exacte de BS est: %f\n\n",callexact);
	}
	else
        printf("S est inferieure a 0");	
}


double fct1(double x, double lambda, double sigma, double K)//la fonction qui va calculer le prix du call max(lambda*exp(sigma*x)-K,0)
{
		return fmax(lambda*exp(sigma*x)-K,0);
		
		
}

		double fct2(double x, double lambda, double sigma, double K)//la fonction qui va calculer le prix du put max(K-lambda*exp(sigma*x),0)
{
		return fmax(K-lambda*exp(sigma*x),0);
		
		
}
		


void VA(double S, double K, double r, double Vol, double theta)//la fonction qui va calculer le prix du call en utilisant la simulation MonteCarlo avec des variables antithetiques
{
	int n=1000;// n est le nombre d'iterations
	double G1,G2;
	double sum = 0;
	double sig = 0;
	double lambda;

	gsl_rng *g;//declaration du generateur de nombres aleatoires
	g=gsl_rng_alloc(gsl_rng_mt19937);//allocation memoire du generateur 
	gsl_rng_env_setup();//environnement du generateur 
	gsl_rng_set(g,time(NULL));//initialiser le generateur 
	lambda=S*exp(-0.5*Vol*Vol*theta);//definir lambda

	for(int i=1; i<=n;i++)//boucle pour calculer le prix du call
	{
		G1=fct1(gsl_ran_gaussian(g,1),lambda,Vol*sqrt(theta),K*exp(-r*theta));//calcule de fct1(g) stocké en G1
		G2=fct1(-gsl_ran_gaussian(g,1),lambda,Vol*sqrt(theta),K*exp(-r*theta));//calcule de fct(-g)stocké en G2
		sum=sum+G1+G2;//Stocker la somme dans la variable sum
	}
	
			sum=sum/(2*n);//calculer le prix du call

		printf("le prix du call par la methode des variables antithetiques est: %f\n\n",sum);

	for(int i=1; i<=n;i++)//boucle pour calculer l'ecart type
	{
		G1=fct1(gsl_ran_gaussian(g,1),lambda,Vol*sqrt(theta),K*exp(-r*theta));//calcule de fct1(g) stocké en G1
		G2=fct1(-gsl_ran_gaussian(g,1),lambda,Vol*sqrt(theta),K*exp(-r*theta));//calcule de fct(-g)stocké en G2n
		sig=sig+(G1-sum)*(G1-sum)+(G2-sum)*(G2-sum);//sum est le prix du call deja calculé

	}
			sig=sqrt(sig/(2*n-1));//calcule de l'ecart type stocké dans la variable sig
double a1 = sum - sig*(1.96/sqrt(2*n));//calcule de la borne inferieure de l'intervalle de confiancea 5 %
double a2 = sum + sig*(1.96/sqrt(2*n));//calcule de la borne superieure de l'intervalle de confiancea 5 %

		printf("l'intervalle de confiance a 5 pourcent est IC(VA)=[%f ; %f]\n\n",a2,a1);//afficher l'intervalle de confiancea 5 %

	
}

void VC(double S, double K, double r, double Vol, double theta)//la fonction qui va calculer le prix du call en utilisant la methode des variables de controles.
{
	int n=1000;// n est le nombre d'iterations.
	double G;
	double sum = 0;
	double sig = 0;
	double lambda;
	gsl_rng *g;//declaration du generateur de nombres aleatoires
	g=gsl_rng_alloc(gsl_rng_mt19937);//allocation memoire du generateur 
	gsl_rng_env_setup();//environnement du generateur 
	gsl_rng_set(g,time(NULL));//initialiser le generateur 
	lambda=S*exp(-0.5*Vol*Vol*theta);//definir lambda

	for(int i=1; i<=n;i++)//boucle pour calculer le prix du call 
	{
		G=fct2(gsl_ran_gaussian(g,1),lambda,Vol*sqrt(theta),K*exp(-r*theta));//calcule de fct2(g) stocké en G
		sum=sum+G;
	}
	
	sum= S-K*exp(-r*theta)+sum/(n);//sum/(n) est le prix du put par la simulation de montecarlo
	
		printf("le prix du call par la methode des variables de controles est: %f\n\n",sum);
		
		
for(int i=1; i<=n;i++)//boucle pour calculer l'ecart type
	{
		G=fct2(gsl_ran_gaussian(g,1),lambda,Vol*sqrt(theta),K*exp(-r*theta));
			
		sig=sig+(G-sum)*(G-sum);

	}
			sig=sqrt(sig/n);//calcule de l'ecart type stocké dans la variable sig
double a1 = sum - sig*(1.96/sqrt(n));//calcule de la borne inferieure de l'intervalle de confiancea 5 %
double a2 = sum + sig*(1.96/sqrt(n));//calcule de la borne superieure de l'intervalle de confiancea 5 %

		printf("l'intervalle de confiance a 5 pourcent IC(VC)=[%f ; %f]\n\n",a2,a1);//afficher l'intervalle de confiancea 5 %

}


int main(){
	
Callexact(111,100,0.2,0.3,1);//appeler la fonction Callexact
VA(111,100,0.2,0.3,1);//appeler la fonction VA
VC(111,100,0.2,0.3,1);//appeler la fonction VC

}

