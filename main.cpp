/* 2015-07-21 14:50:40 

	Projekt: Nischennetz mit Migration mit Robustness Analyse

	Quelldatei: main.c zu Nischennetz_mit_Migration_1.0
	Programmiersprache: C
	Autoren: S. J. Plitzko, M.Hamm
	Version 1.0

###############################################################################################################################
Beschreibung:
Dieses Programm simuliert ein über mehrere Patches ausgedehntes Nahrungsnetz mit dem Nischenmodell. Die Populationsdynamik wird mit einer Holling Typ II Funktion modelliert. Beim Aufrufen können über die Konsole folgende Parameter bestimmt werden:

Anzahl Spezies S
Anzahl basale Spezies B
Allometrie Koeffizient x
Migrationsstärke d
Topologie des Lebensrausm (siehe Set_Topology())
Anzahl der statistischen Wiederholungen L
Anzahl der Patches Y
Migration ja oder nein "M"

Kompilieren mit: gcc -o V1_15_07_13 main.c getargs.c topology.c mignicheweb.c evolveweb.c holling2.c robustness.c -lm -lgsl -lgslcblas -Wall

Eingabe in der Konsole: ./NAME -S X -B X -T X -L X -Y X -d X.X -x X.X -M "X" -R X.X	mit X bzw. X.X numerischer Wert der Variablen (Reihenfolge egal)

Das Ergebnis der Simulation wird in der Datei xxx.out gespeichert. Dabei sind die Daten wie folgt angeordnet (alles eine Zeile): 
S		B		M		x		Y		dpow	T
Rob		Perlok	Perges
Si_ges	Si_TL1	Si_TL2	Si_TL3	Si_TL4	Si_TL>4		Spezies initial
Sf_ges	Sf_TL1	Sf_TL2	Sf_TL3	Sf_TL4	Sf_TL>4		Spezies final
Bi_ges	Bi_TL1	Bi_TL2	Bi_TL3	Bi_TL4	Bi_TL>4		Biomass initial
Bf_ges	Bf_TL1	Bf_TL2	Bf_TL3	Bf_TL4	Bf_TL>4		Biomass final
Sh_ges	Sh_TL1	Sh_TL2	Sh_TL3	Sh_TL4	Sh_TL>4		Spezies Hub
Bh_ges	Bh_TL1	Bh_TL2	Bh_TL3	Bh_TL4	Bh_TL>4		Biomass Hub
Ss_ges	Ss_TL1	Ss_TL2	Ss_TL3	Ss_TL4	Ss_TL>4		
Bs_ges	Bs_TL1	Bs_TL2	Bs_TL3	Bs_TL4	Bs_TL>4
1mit2	2mit3	3mit1
Fixp0	Fixp1	Fixp2	Fixp3	Fixp4	Fixp5	Fixp6	Fixp7

###############################################################################################################################
*/

//--Standardbibliotheken--------------------------------------------------------------

#include <string.h>						// string modification functions
#include <time.h>						// time functions
#include <math.h>						// math functions
#include <stdio.h>						// output functions
#include <stdlib.h>						// standard
#include <limits.h>
#include <gsl/gsl_rng.h>				// random number generator functions
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>			// random number distributions
#include <gsl/gsl_sort_vector.h>		// vector sorting functions
#include <gsl/gsl_odeiv.h>				// differential equation solver
#include <gsl/gsl_errno.h>				// errorhandler


//--Selbstdefinierte Funktionen------------------------------------------------------
#include "structs.h"					// foodweb Struktur

#include "getargs.h"					// getArgs
#include "topology.h"					// SetTopology
#include "mignicheweb.h"				// SetNicheNetwork, LinkElements, CountLinks
#include "evolveweb.h"					// DGL lösen
#include "holling2.h"					// DGL lösen			
#include "robustness.h"					// Robustness Analyse
#include "StandardDeviation.h"				// Standardabweichung berechnen
#include "createOutput.h"				// Output erzeugen

//--Verzeichnis für Ergebnisse-------------------------------------------------------
#define ORT "/home/tatjana/Arbeitsfläche/MichaelasProgramm/stochastischeMigration/ErsteVersuche/"
#define ORT2 "/home/tatjana/Arbeitsfläche/MichaelasProgramm/stochastischeMigration/ErsteVersuche/"
//++START++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int main(int argc, char** argv)
{	

//--Foodweb Struktur mit Standardwerten aufstellen------------------------------------------------------------------------------------------

	gsl_vector* fixpunkte	= gsl_vector_calloc(9);
 	

	struct foodweb nicheweb	= {NULL, fixpunkte, NULL, 18, 3, 1, 5, 0, 0, -7., 0.0, 0, 1};		// Reihenfolge: network, fxpkt, migrPara, AllMus, AllNus, S, B, Rnum, Y, T, Tchoice, d, x, M, Z
	
	struct migration stochastic = {NULL, NULL, NULL, NULL, NULL, NULL};
	
	struct resource res = {500.0, 0.0};											// Resource: Größe, Wachstum
	
	
	
	
//--Konsoleneingabe-------------------------------------------------------------------------------------------------------------------------
	
	int L = 5;	// Statistik
	int i = 0;	// Counter
	
	int checksum = getArgs(argc, argv, &(nicheweb.S), &(nicheweb.B), &(nicheweb.T), &(nicheweb.d), &L, &(nicheweb.Y), &(nicheweb.x), &(nicheweb.M), &(res.size), &(nicheweb.Z));	

		if (checksum != 10 && checksum!=(int)(argc-1)/2) 	// Alles gesetzt?									
		 {	
			printf("Bitte gültige Eingabe für Parameter machen!\nProgramm wird beendet.\n");		
			return(0);		
		 }

	int length			= ((nicheweb.Rnum+nicheweb.S)*(nicheweb.S+nicheweb.Rnum)+1+nicheweb.Y*nicheweb.Y+1+(nicheweb.Rnum+nicheweb.S)+nicheweb.S);	// Länge des Rückabewerts
	nicheweb.network = gsl_vector_calloc(length);
	
	printf("Z = %i\n",nicheweb.Z);
	nicheweb.migrPara = gsl_vector_calloc(7); // Reihenfolge: tau, mu, nu, SpeciesNumber, momentanes t, ymigr, migrationEventNumber	 
	stochastic.SpeciesNumbers = gsl_vector_calloc(L*nicheweb.Z);
	stochastic.AllMus = gsl_vector_calloc(nicheweb.Z);
	stochastic.AllNus = gsl_vector_calloc(nicheweb.Z);
	stochastic.Biomass_SpeciesNumbers = gsl_vector_calloc(nicheweb.Z);
	stochastic.Biomass_AllMus = gsl_vector_calloc(nicheweb.Z);
	stochastic.Biomass_AllNus = gsl_vector_calloc(nicheweb.Z);
	
//--Zufallszahlengenerator initialisieren--------------------------------------------------------------------------------

		const gsl_rng_type *rng1_T;											// ****
		gsl_rng *rng1;   													// initialize random number generator
		gsl_rng_env_setup();   												// ermöglicht Konsolenparameter
		rng1_T = gsl_rng_default;   										// default random number generator (so called mt19937)
		//gsl_rng_default_seed = 0;											// default seed for rng
		gsl_rng_default_seed = ((unsigned)time(NULL));						// random starting seed for rng
		rng1 = gsl_rng_alloc(rng1_T);
		
		
		
		
//--Struct initialisieren für patchweise Ausgabe----------------------------------------------------------------------------------------		 

	struct data patchwise[nicheweb.Y];
	for(i=0; i<nicheweb.Y; i++)
	{
	  gsl_vector* sini  = gsl_vector_calloc(6);
	  gsl_vector* sfini  = gsl_vector_calloc(6);
	  gsl_vector* bini  = gsl_vector_calloc(6);
	  gsl_vector* bfini  = gsl_vector_calloc(6);
	  gsl_vector* robness  = gsl_vector_calloc(2);
	  
	  //struct data tempo = {sini,sfini,bini,bfini,robness};
	  struct data temp = {sini,sfini,bini,bfini,robness};
	  patchwise[i] = temp;
	}
	
	//printf("test");
//--Initialisierungen---------------------------------------------------------------------------------------------------
	nicheweb.Tchoice = nicheweb.T;
	nicheweb.T = 0;
	nicheweb.d = nicheweb.d/10;
	res.size = res.size/10;
	//int len	= ((nicheweb.Rnum+nicheweb.S)*(nicheweb.S+nicheweb.Rnum)+1+nicheweb.Y*nicheweb.Y+1+(nicheweb.Rnum+nicheweb.S)+nicheweb.S);	// Länge des Rückabewerts

	gsl_vector *populationFIN 	= gsl_vector_calloc((nicheweb.Rnum + nicheweb.S)*(nicheweb.Y)*5 + (nicheweb.S) + 3);				// Gleiche Länge wie Rückgabe von evolveNetwork
	gsl_vector *robustness		= gsl_vector_calloc(63);
	gsl_vector *resultEvolveWeb	= gsl_vector_calloc((nicheweb.Rnum+nicheweb.S)*nicheweb.Y*5 + 3 + nicheweb.S); 				// y[Simulation], y0, ymax, ymin, yavg, fixp, TL
	gsl_vector *resultRobustness 	= gsl_vector_calloc(63);
	gsl_matrix *D			= gsl_matrix_calloc(nicheweb.Y,nicheweb.Y);
	
	gsl_vector *robustnesstemp	= gsl_vector_calloc(63);
	gsl_vector *meanSquOfDataAll 	= gsl_vector_calloc(63);
	gsl_vector *meanSquOfDataAlltemp = gsl_vector_calloc(63);
	gsl_vector *standardDeviationAll = gsl_vector_calloc(63);
	
	gsl_vector *meanOfData	= gsl_vector_calloc((6*4+2)*nicheweb.Y);
	gsl_vector *meanOfDatatemp = gsl_vector_calloc((6*4+2)*nicheweb.Y);
	gsl_vector *meanSquOfData = gsl_vector_calloc((6*4+2)*nicheweb.Y);
	
	gsl_vector *standardDeviation = gsl_vector_calloc((6*4+2)*nicheweb.Y);
	
	gsl_vector_set_zero(robustness);
	gsl_vector_set_zero(meanOfData);
	gsl_vector_set_zero(meanSquOfData); 
	gsl_vector_set_zero(nicheweb.migrPara);
	gsl_vector_set_zero(meanSquOfDataAll);
	
	double SpeciesNumber[L*nicheweb.Z][2]; 
	double AllMu[L*nicheweb.Z][2];
	double AllNu[L*nicheweb.Z][2];
	
	double ymigr = 0;
	double mu = 0;
	double nu = 0;
	double ymigrtemp;
	double ymigrSqu = 0;
	double ymigrDeviation;

//--Simulation---------------------------------------------------------------------------------------------------	
	D    = SetTopology(nicheweb.Y, nicheweb.T, D);						// migration matrix
	for(i = 0; i < L; i++)																							
	 { 	
// 		const gsl_rng_type *rng1_T;											// ****
// 		gsl_rng *rng1;   													// initialize random number generator
// 		gsl_rng_env_setup();   												// ermöglicht Konsolenparameter
// 		rng1_T = gsl_rng_default;   										// default random number generator (so called mt19937)
// 		gsl_rng_default_seed = 0;											// default seed for rng
// 		//gsl_rng_default_seed = ((unsigned)time(NULL));						// random starting seed for rng
// 		rng1 = gsl_rng_alloc(rng1_T);
		
		
		printf("\nStarte Durchlauf L = %i\n", i);
//--Starte Simulation-----------------------------------------------------------------------------------------------			
		SetNicheNetwork(nicheweb, res, D, rng1, rng1_T);
		//printf("erstes Element ist %f\n",gsl_vector_get(nicheweb.network,0));
		gsl_vector_set_zero(resultEvolveWeb);
		populationFIN	 = EvolveNetwork(nicheweb, stochastic, rng1, rng1_T, resultEvolveWeb);
			
		gsl_vector_set_zero(resultRobustness);										
		gsl_vector_memcpy(robustnesstemp, EvaluateRobustness(populationFIN, nicheweb, patchwise, resultRobustness));	// Robustness Analyse
		
		
//--Standardabweichung für Mittelung vorbereiten-----------------------------------------------------------------------------------------		
		determineMean(robustnesstemp, 63, robustness);
		determineMeanSqu(robustnesstemp, 63, meanSquOfDataAll);
		
//--Ausgabewerte----------------------------------------------------------------------------------------------------------		
		ymigrtemp = gsl_vector_get(nicheweb.migrPara, 5);
		for(int j= 0; j<nicheweb.Z; j++)
		{
		  AllMu[i*nicheweb.Z+j][0] = gsl_vector_get(stochastic.AllMus, j);
		  AllNu[i*nicheweb.Z+j][0] = gsl_vector_get(stochastic.AllNus, j);
		  SpeciesNumber[i*nicheweb.Z+j][0] = gsl_vector_get(stochastic.SpeciesNumbers,j);
		  
		  AllMu[i*nicheweb.Z+j][1] = gsl_vector_get(stochastic.Biomass_AllMus, j);
		  AllNu[i*nicheweb.Z+j][1] = gsl_vector_get(stochastic.Biomass_AllNus, j);
		  SpeciesNumber[i*nicheweb.Z+j][1] = gsl_vector_get(stochastic.Biomass_SpeciesNumbers,j);
		}
		//printf("SpeciesNumber ist %f\n",SpeciesNumber[i]);
		ymigr += ymigrtemp;
		
		ymigrSqu += (ymigrtemp*ymigrtemp); 
		
		
//--Mittelwert und Vorbereitungen für Standardabweichung für die patchweise Ausgabe berechnen--------------------------------		
		linkElements(patchwise, nicheweb.Y, meanOfDatatemp);
		
		determineMean(meanOfDatatemp, (6*4+2)*nicheweb.Y, meanOfData);
		determineMeanSqu(meanOfDatatemp, (6*4+2)*nicheweb.Y, meanSquOfData);
		
		printf("\nBeende Durchlauf L = %i\n", i);
	 }


//-- Standardabweichung berechnen--------------------------------------------------------------------------------------

	ymigrSqu = ymigrSqu/L;
	ymigr = ymigr/L;
	ymigrDeviation = sqrt(ymigrSqu - ymigr*ymigr);

//-- Für patchweise Ausgabe-------------------------------------------------------------------------------------------
	standardDeviation = determineStandardDeviation((6*4+2)*nicheweb.Y, meanOfData, meanSquOfData, L, standardDeviation);
	standardDeviationAll = determineStandardDeviation(63, robustness, meanSquOfDataAll, L, standardDeviationAll);
	 //printf("der 3. Eintrag in standardDeviationAll ist %f\n", gsl_vector_get(standardDeviationAll,3));
// 	 printf("S ist %f\n", gsl_vector_get(robustness,3));
// 	 printf("Standardabweichung von S ist %f\n", gsl_vector_get(standardDeviationAll,3));
// 	 printf("meanOfDataSqu ist %f\n", gsl_vector_get(meanOfDataSquAll,3));
// 	 printf("meanSquOfData ist %f\n", gsl_vector_get(meanSquOfDataAll,3));
	 
	 
	 
	printf("L=%i\tspeciesini=%f\tspeciesfinal=%f\n", L, gsl_vector_get(robustness, 3)/L, gsl_vector_get(robustness, 9)/L);
	

	
	
//--Abspeichern in File-------------------------------------------------------------------------------------	
	char aims[255] = ORT;

	
	createOutputGeneral(nicheweb, res, aims, robustness, standardDeviationAll, L, mu, nu, ymigr, ymigrDeviation);													// Datei schließen

	
    
	
//--Daten patchweise abspeichern----------------------------------------------------------------------	
	printf("population ist %f\n",gsl_vector_get(stochastic.Biomass_AllMus,0));
     
     for(int l = 0 ; l< nicheweb.Y; l++)
     {
       //char name[100];
       char aims2[255] = ORT2;
       
       createOutputPatchwise(nicheweb, res, aims2, meanOfData, standardDeviation, L, l);
       
      }
      
      if(nicheweb.Tchoice != 0)
      {
//--Ausgewählte Spezies rausschreiben, die migrieren darf---------------------------------------------------------------------------
	
	char aims3[255] = ORT;
	
	createOutputSpeciesNumber(nicheweb, res, aims3, SpeciesNumber, L);

      
//--Ausgewählte Verbindung rausschreiben, über die migriert werden darf---------------------------------------------------------------------------
	
	char aims4[255] = ORT;
	
	createOutputPatchlink(nicheweb, res, aims4, AllMu, AllNu, L);
      }
      
      printf("\nSimulation abgespeichert\n\n");
	
//--free----------------------------------------------------------------------------------------------------------------  
	free(nicheweb.network);
	
	gsl_vector_free(fixpunkte);
	
	for(i=0; i<nicheweb.Y; i++)
	{
	  gsl_vector_free(patchwise[i].sini);
	  gsl_vector_free(patchwise[i].sfini);
	  gsl_vector_free(patchwise[i].bini);
	  gsl_vector_free(patchwise[i].bfini);
	  gsl_vector_free(patchwise[i].robness);

	}
	

	gsl_vector_free(nicheweb.migrPara);
	gsl_vector_free(stochastic.AllMus);
	gsl_vector_free(stochastic.AllNus);
	gsl_vector_free(stochastic.SpeciesNumbers);
	gsl_vector_free(stochastic.Biomass_AllMus);
	gsl_vector_free(stochastic.Biomass_AllNus);
	gsl_vector_free(stochastic.Biomass_SpeciesNumbers);
	gsl_vector_free(populationFIN);
	gsl_vector_free(robustness);	
	gsl_vector_free(robustnesstemp);
	gsl_vector_free(meanOfData);
	gsl_vector_free(meanOfDatatemp);
	gsl_vector_free(meanSquOfData);
	gsl_vector_free(standardDeviation);
	gsl_vector_free(standardDeviationAll);
	gsl_vector_free(meanSquOfDataAll);
	gsl_vector_free(meanSquOfDataAlltemp);
	gsl_rng_free(rng1);
	
	return(0);

}




















