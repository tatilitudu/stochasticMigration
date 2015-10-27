#ifndef STRUCTS_H
#define STRUCTS_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

/*
	Die foodweb Struktur enthält die Konsolenparameter S, B, Rnum, Y, T, M, d, x sowie einen Netzwerk Vektor für die Ergebnisse der Nischennetz-Berechnung.
*/
struct foodweb{

	gsl_vector* network;
	gsl_vector* fixpunkte;				// fix0,1,2, fixp0,1,2, testf0,1,2
	gsl_vector* migrPara;
	
	
	int S;
	int B;
	int Rnum;

	int Y;
	int T;
	int Tchoice;
	
	double d;
	double x;

	int M;
	int Z;
	
};


struct migration{
  
	gsl_vector* SpeciesNumbers;
	gsl_vector* AllMus;
	gsl_vector* AllNus;
	gsl_vector* Biomass_SpeciesNumbers;
	gsl_vector* Biomass_AllMus;
	gsl_vector* Biomass_AllNus;
	
};



struct resource{

	double size;
	double growth;

};

struct data{
      gsl_vector* sini;
      gsl_vector* sfini;
      gsl_vector* bini;
      gsl_vector* bfini;
      gsl_vector* robness;
};
      



#endif
