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



int main(){

FILE *stream,*f;//declaration des pointeurs fichiers
int n=10;//initialisation du nombre d'iterations

float sigma=0.1,r=0.1,K=10,B=0,T=10,t=0,Seuler=100,SMilshtein=100,B_0=0,S=100,g=0e0,Sexact;// initialisation des variables,Seuler est qui va contenir le prix de l'actif qu'on va calculer selon le shema deuler , SMilshtein est la variable qui va contenir le prix de l'actif qu'on va calculer selon le shema de Mishtein, Sexact est la variable qui va contenir le prix de l'actif qu'on va calculer selon la formule close de Black & Scholes.les 3 variables sont initialisé a 100
// B et B_0 sont deux mouvements browniens standards qu'on va generer
    gsl_rng *m;//declaration du generateur de nombres aleatoires
	m=gsl_rng_alloc(gsl_rng_mt19937);//allocation memoire du generateur 
	gsl_rng_env_setup();//environnement du generateur 
	gsl_rng_set(m,time(NULL));//initialiser le generateur 
	
if ((stream = fopen( "C:\\Users\\hp\\Documents\\brown.dat","w"))==NULL )//Si on ne reussit pas a creer le fichier brown.dat alors on affiche le message sinon on continue
{
	fprintf(stderr,"cannot open input file.\n");
	return 1;
}
fprintf(stream,"t Sexact Seuler SMilshtein \n");//on ecrit sur la premiere ligne du fichier les titres des colonnes de donnees qu'on va inserer, la premiere colonne correspond au valeurs Sexact, la deuxieme colonne correspond au valeurs Seuler,la troisieme colonne correspond au valeurs SMilshtein
fprintf(stream,"%f %f %f %f \n",t,S,Seuler,SMilshtein);//on remplit les valeurs initiales des trois variables dans la deuxieme ligne
srand(time(NULL));

do //boucle pour calculer les prix selon les trois schema et les remplir dans le fichier
{
	
	t=t+T/n;

    g=g+gsl_ran_gaussian(m,1);//accumulation du g ,la variable gaussienne
    B=sqrt(T/n)*g;
    Sexact=S*exp((r-sigma*sigma/2)*t+sigma*B);//calcule du St de BS
    Seuler=Seuler*(1+r*T/n+sigma*(B-B_0));//calcule du St du schema d'euler
    SMilshtein=SMilshtein*(1+(r-sigma*sigma/2)*T/n+sigma*(B-B_0)+sigma*sigma*(B-B_0)*(B-B_0)/2);//calcule de St selon le schema de SMilshtein 
    B_0=B;//reinitialisation de B_0 qui est le cumule decalé
    fprintf(stream,"%f %f %f %f \n",t,Sexact,Seuler,SMilshtein);//remplissage de la ligne t du fichier

}
while(t<T);//tant que la date est inferieur a T 
fclose(stream);//cloture du fichier
//ouverture du shell et lancement du gnuolot
f=popen("\"D:\\gnuplot2\\bin\\gnuplot.exe\"","w");//gnuplot installé en D
//execution du gnuolot

fprintf(f,"set style data lines\n");
fprintf(f,"plot for [col=2:4] 'C:\\Users\\hp\\Documents\\brown.dat' using 1:col with lines title columnheader \n");//on dessine nos 3 courbes pour les 3 prix calculés, le titre de chaque courbe correspondrait a l'entete de la colonne correspondante a la serie de prix

fflush(f);
printf("press enter\n");
getchar();
pclose(f);



}

