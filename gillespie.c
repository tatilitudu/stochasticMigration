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

double* stochMigration(struct foodweb nicheweb, struct migration stochastic, const double y[], gsl_rng* rng1, const gsl_rng_type* rng1_T, int migrationEventNumber, gsl_matrix* Dchoice)
{
  int Y = nicheweb.Y;
  int S = nicheweb.S;
  int Tchoice = nicheweb.Tchoice;
  int Rnum = nicheweb.Rnum;
  int Z = nicheweb.Z;
  //gsl_vector *network 	= nicheweb.network;		// Inhalt: A+linksA+Y+linksY+Massen+Trophische_Level = (Rnum+S)²+1+Y²+1+(Rnum+S)+S
  
/*  gsl_vector_view D_view = gsl_vector_subvector(network, (Rnum+S)*(Rnum+S)+1, Y*Y);					// Migrationsmatrix D als Vektor
  gsl_matrix_view ED_mat = gsl_matrix_view_vector(&D_view.vector, Y, Y);								// D als Matrixview
  gsl_matrix *EDmat	 = &ED_mat.matrix;*/	
  
  gsl_vector *c	= gsl_vector_calloc(Y);
  gsl_vector *linkCount	= gsl_vector_calloc(Y);
  double ctemp=0;
  double ctot;
  int l,i;
  
  //printf("\n");
  
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
  //srand(SEED);	
  int linkCountTemp, m;

  
	
  	
  //printf("\n");

//   printf("Tchoice ist %i\n",Tchoice);
//   
//   int k,j;
//   for(k = 0; k< Y;k++)
//   {
//     for(j=0;j<Y;j++)
//     {
//       printf("%f\t",gsl_matrix_get(Dchoice,k,j));
//     }
//     printf("\n");
//     
//   }
  
  
  if( Tchoice == 0 )
  {
    for(i=0; i<3*Z+3;i++)
    {
      gsl_vector_set(nicheweb.migrPara, i, 0);
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
    //double r,r1,r2;
  
    gsl_vector_memcpy(a,linkCount);
    gsl_vector_mul(a,c);
  
    atot = gsl_blas_dasum(a);
  
    //printf("Z ist %i\n", Z);
    
    //printf("Berechne Zeitpunkte, zu den migriert werden soll:\n");
    double tau = choose_time(atot, rng1);
    gsl_vector_set(nicheweb.migrPara, 0, tau);

    //printf("\n");
    //r = (double)rand()/INT_MAX;
    //printf("r ist %f\n",r);
  
    //printf("Berechne von welchem Patch aus migriert werden soll\t");
    //printf("r: %f\n",r);
    //printf("atot: %f\n",atot);
    int mu = select_patch(stochastic, a, atot, rng1, Y, migrationEventNumber, 0);
    //printf("population ist %f\n",gsl_vector_get(stochastic.Biomass_AllMus, migrationEventNumber));
    gsl_vector_set(stochastic.AllMus, migrationEventNumber, mu);
    gsl_vector_set(nicheweb.migrPara, 1, mu);
    //printf("mu: %i\n",mu);
    int flag=1;
  
    //printf("Berechne in welches Patch migriert werden soll\t\t");
    int nu; 
    while(flag != 0)
    {
      //r1  = (double)rand()/INT_MAX;
      //printf("r1 ist %f\n",r1);
      nu = select_patch_random(nicheweb, rng1);
      
      if(nu!= mu  && gsl_matrix_get(Dchoice, nu, mu)!=0)
      {
	gsl_vector_set(nicheweb.migrPara, 2, nu);
	gsl_vector_set(stochastic.AllNus, migrationEventNumber, nu); 
	flag = 0;
      }
    }
    //printf("nu: %i\n", nu);
    
    
    int SpeciesNumber;
    //r2 = (double)rand()/INT_MAX;
    //printf("r2 ist %f\n",r2);
    //int Choice = 0;
    SpeciesNumber = select_species_random(nicheweb, stochastic, rng1);
    gsl_vector_set(stochastic.SpeciesNumbers, migrationEventNumber, SpeciesNumber);
    
    gsl_vector_set(nicheweb.migrPara, 3, SpeciesNumber);
    
    //printf("SpeciesNumber: %i\n\n", SpeciesNumber);
    //printf("Population dieser Spezies ist %f\n",y[SpeciesNumber+Rnum]);
    
    if(SpeciesNumber>S)
    {
      printf("\n\nFehler!!! SpeciesNumber>S \n\n");
      
    }
    
    gsl_vector_free(a);
  }
  

  
  gsl_vector_free(c);
  gsl_vector_free(linkCount);
  
  
  
  return 0;
}

double choose_time(double atot, gsl_rng* rng1)
{
  double tau;
  double r = gsl_rng_uniform_pos(rng1);
  if( atot>0 )
  {
    // rand liefert zufällige Zahl zwischen 0 und INT_MAX
    tau = -log(r)/ atot;
  }
  //printf("tau: %f\t",tau);
    
  return tau;
}


int select_patch(struct migration stochastic, gsl_vector* a, double atot, gsl_rng* rng1, int Y, int migrationEventNumber, int whichPatch)
{
  double r = gsl_rng_uniform_pos(rng1);
  int i;
  int mu;
  double sum=0;
  
  r = r*atot;
  //printf("r*atot in select_patch ist %f\n",r);
  for(i=0;i<Y;i++)
  {
    sum += gsl_vector_get(a,i);
    if( r < sum )
    {
      mu = i;
      break;
    }
  }
  if(whichPatch == 1)
  {
    gsl_vector_set(stochastic.Biomass_AllNus, migrationEventNumber, gsl_vector_get(a,mu));
  }
  else if(whichPatch == 0)
  {
    
    gsl_vector_set(stochastic.Biomass_AllMus, migrationEventNumber, gsl_vector_get(a,mu));
    //printf("population1 ist %f\n",gsl_vector_get(stochastic.Biomass_AllMus, migrationEventNumber));
  }
    
  return mu;
}


int select_patch_random(struct foodweb nicheweb, gsl_rng* rng1)
{
  double r = gsl_rng_uniform_pos(rng1);
  
  int Y = nicheweb.Y;
  int patchNumber;
  float patchNumberFloat;
  
  patchNumberFloat = r*(Y-1);
  //printf("patchNumberFloat ist %f\n", patchNumberFloat);
  
  if(patchNumberFloat>0) patchNumber = (int)(patchNumberFloat + 0.5);
	
  else patchNumber =  (int)(patchNumberFloat - 0.5);
  
  //printf("patchNumber ist %i\n", patchNumber);
  
  return patchNumber;
}


int select_species(struct foodweb nicheweb, struct migration stochastic, gsl_rng* rng1, int Choice, const double y[], int migrationEventNumber)
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
  double r = gsl_rng_uniform_pos(rng1);
  gsl_vector *a = gsl_vector_calloc(S);
  gsl_vector *c = gsl_vector_calloc(S);
  
  
  //printf("Berechne, welche Spezies migrieren darf\t\t\t");
  
  for(i = 0; i< S ;i++)
  {
    gsl_vector_set(a,i,y[Rnum+i]);
  }
  //printf("Eintrag 5 von y ist %f\n", y[5]);
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
  
  //printf("atot ist %f\n",atot);
  
  r = r*atot;

  //printf("r*atot ist %f\n",r);
  
  for(i=0;i<S;i++)
  {
    sum += gsl_vector_get(a,i);
    if( r < sum )
    {
      SpeciesNumber = i;
      break;
    }
  }
  
  gsl_vector_set(stochastic.Biomass_SpeciesNumbers, migrationEventNumber, gsl_vector_get(a,SpeciesNumber));
  
  gsl_vector_free(a);
  gsl_vector_free(c);
  
  return SpeciesNumber;
  
}


int select_species_random(struct foodweb nicheweb, struct migration stochastic, gsl_rng* rng1)
{
  double r = gsl_rng_uniform_pos(rng1);
  
  int S = nicheweb.S;
  int SpeciesNumber;
  float SpeciesNumberFloat;
  
  
  //printf("Berechne, welche Spezies migrieren darf\t\t\t");
  
  
  SpeciesNumberFloat = r*(S-1);
  //printf("SpeciesNumberFloat ist %f\n", SpeciesNumberFloat);
  
  if(SpeciesNumberFloat>0) SpeciesNumber = (int)(SpeciesNumberFloat + 0.5);
	
  else SpeciesNumber =  (int)(SpeciesNumberFloat - 0.5);
  
  //printf("SpeciesNumber ist %i\n", SpeciesNumber);
  
  
  
  return SpeciesNumber;
}