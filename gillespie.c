#include "structs.h"


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include <gsl/gsl_rng.h>					// random number generator functions
#include <gsl/gsl_randist.h>				// random number distributions
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "gillespie.h"
#include "topology.h"

#define SEED	123

double* stochMigration(struct foodweb nicheweb, double* migrationWerte, const double y[] )
{
  int Y = nicheweb.Y;
  int S = nicheweb.S;
  int Tchoice = nicheweb.Tchoice;
  int Rnum = nicheweb.Rnum;
  //gsl_vector *network 	= nicheweb.network;		// Inhalt: A+linksA+Y+linksY+Massen+Trophische_Level = (Rnum+S)²+1+Y²+1+(Rnum+S)+S
  
/*  gsl_vector_view D_view = gsl_vector_subvector(network, (Rnum+S)*(Rnum+S)+1, Y*Y);					// Migrationsmatrix D als Vektor
  gsl_matrix_view ED_mat = gsl_matrix_view_vector(&D_view.vector, Y, Y);								// D als Matrixview
  gsl_matrix *EDmat	 = &ED_mat.matrix;*/	
  
  gsl_vector *c	= gsl_vector_calloc(Y);
  gsl_vector *linkCount	= gsl_vector_calloc(Y);
  double ctemp=0;
  double ctot;
  int l,i;
  
  printf("\n");
  
  // Setze c(l) als (Population in Patch l)/(Population in allen Patches) 
  for(l=0;l<Y;l++)
  {
    for(i=0;i<Rnum+S;i++)
    {
      ctemp += y[l*(Rnum+S)+i];
    }
    gsl_vector_set(c,l,ctemp);
  }
  ctot = gsl_blas_dasum(c);
  gsl_vector_scale(c,1/ctot);
  
//   for(l=0;l<Y;l++)
//   {
//     printf("Eintrag %i von c ist %f\n",l,gsl_vector_get(c,l));
//   }
  
  // Starte Random Number Generator
  srand(SEED);	
  int linkCountTemp, m;

  
	
  gsl_matrix *Dchoice    = SetTopology(Y, Tchoice);	

  //printf("Tchoice ist %i\n",Tchoice);
  
  if( Tchoice == 0 )
  {
    for(i=0; i<3;i++)
    {
      migrationWerte[i]=0;
    }
    printf("Es findet keine Migration statt\n");
  }
  else
  {
    for(l=0;l<Y;l++)
    {
      linkCountTemp=0;
      for(m=0;m<Y;m++)
      {
	if( gsl_matrix_get(Dchoice,l,m)!=0 )
	{
	  linkCountTemp++;
	}
      }
      gsl_vector_set(linkCount,l,linkCountTemp);
    }
    //gsl_vector_scale(linkCount,1/linkCountTot);
  
  
    gsl_vector *a	= gsl_vector_calloc(Y);
    double atot;
    double r,r1,r2;
  
    gsl_vector_memcpy(a,linkCount);
    gsl_vector_mul(a,c);
  
    atot = gsl_blas_dasum(a);
  
    printf("\nBerechne Zeitpunkt, zu dem migriert werden soll\t\t");
    if( atot>0 )
    {
      // rand liefert zufällige Zahl zwischen 0 und INT_MAX
      migrationWerte[0] = -log((double)rand()/INT_MAX) / atot;
    }
    printf("tau: %f\n",migrationWerte[0]);
  
    r = (double)rand()/INT_MAX;
    printf("r ist %f\n",r);
  
    printf("Berechne von welchem Patch aus migriert werden soll\t");
    //printf("r: %f\n",r);
    //printf("atot: %f\n",atot);
    migrationWerte[1] = select_patch(a,atot,r,Y);
    printf("mu: %f\n",migrationWerte[1]);
    int flag=1;
  
    printf("Berechne in welches Patch migriert werden soll\t\t");
    while(flag != 0)
    {
      r1  = (double)rand()/INT_MAX;
      printf("r1 ist %f\n",r1);
      migrationWerte[2] = select_patch(a,atot,r1,Y);
      if(migrationWerte[2]!=migrationWerte[1] && gsl_matrix_get(Dchoice,migrationWerte[2],migrationWerte[1])!=0)
      {
	flag = 0;
      }
    }
    printf("nu: %f\n",migrationWerte[2]);
    
    
    int SpeciesNumber;
    r2 = (double)rand()/INT_MAX;
    printf("r2 ist %f\n",r2);
    int Choice = 0;
    SpeciesNumber = select_species(nicheweb, r2, Choice, y);
    migrationWerte[3] = SpeciesNumber;
    
    printf("SpeciesNumber: %i\n\n", SpeciesNumber);
    printf("Population dieser Spezies ist %f\n",y[SpeciesNumber+Rnum]);
    
    if(SpeciesNumber>S)
    {
      printf("\n\nFehler!!! SpeciesNumber>S \n\n");
      
    }
    
    gsl_vector_free(a);
  }
  

  
  gsl_vector_free(c);
  gsl_vector_free(linkCount);
  gsl_matrix_free(Dchoice);
  
  
  return 0;
}


int select_patch(gsl_vector* a, double atot, double r, int Y)
{
  int i;
  int mu;
  double sum=0;
  
  r = r*atot;
  printf("r*atot in select_patch ist %f\n",r);
  for(i=0;i<Y;i++)
  {
    sum += gsl_vector_get(a,i);
    if( r < sum )
    {
      mu = i;
      break;
    }
  }
  return mu;
}


int select_species(struct foodweb nicheweb, double r, int Choice, const double y[])
{
  int S = nicheweb.S;
  int Y = nicheweb.Y;
  int Rnum = nicheweb.Rnum;
  gsl_vector *network = nicheweb.network;
  gsl_vector_view M_vec = gsl_vector_subvector(network, ((Rnum+S)*(Rnum+S))+1+(Y*Y)+1, (Rnum+S));	// Massenvektor
  gsl_vector *Mvec = &M_vec.vector;
  
  
  
  int i;
  int SpeciesNumber;
  double sum = 0;
  double atot;
  gsl_vector *a = gsl_vector_calloc(S);
  gsl_vector *c = gsl_vector_calloc(S);
  
  
  printf("Berechne, welche Spezies migrieren darf\t\t");
  
  for(i = 0; i< S ;i++)
  {
    gsl_vector_set(a,i,y[Rnum+i]);
  }
  printf("Eintrag 5 von y ist %f\n", y[5]);
  // Nur Abhängigkeit, wie groß Population ist == 0; zusätzlich Massen == 1
  if(Choice == 0)
  {
    atot = gsl_blas_dasum(a);
  }
  else if(Choice == 1)
  {
    double ctot;
    for(i = 0; i < S; i++ )
    {
      gsl_vector_set(c, i , gsl_vector_get(Mvec,i+Rnum));
    }
    ctot = gsl_blas_dasum(c);
    gsl_vector_scale(c,1/ctot);
    gsl_vector_mul(a,c);
    atot = gsl_blas_dasum(a);
  }
    
  // Abhängigkeit
  
  printf("atot ist %f\n",atot);
  
  r = r*atot;

  printf("r*atot ist %f\n",r);
  
  for(i=Rnum;i<S+Rnum;i++)
  {
    sum += gsl_vector_get(a,i);
    if( r < sum )
    {
      SpeciesNumber = i;
      break;
    }
  }
  return SpeciesNumber;
  
}