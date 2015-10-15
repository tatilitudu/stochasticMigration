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
#include "StandardDeviation.h"					// Standardabweichung berechnen

//--Verzeichnis für Ergebnisse-------------------------------------------------------
#define ORT "/home/tatjana/Arbeitsfläche/MichaelasProgramm/stochastischeMigration/ErsteVersuche/"
#define ORT2 "/home/tatjana/Arbeitsfläche/MichaelasProgramm/stochastischeMigration/ErsteVersuche/"
//++START++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int main(int argc, char** argv)
{	

//--Foodweb Struktur mit Standardwerten aufstellen------------------------------------------------------------------------------------------

	gsl_vector* fixpunkte	= gsl_vector_calloc(9);
 	

	struct foodweb nicheweb	= {NULL, fixpunkte, NULL, 18, 3, 1, 5, 0, 0, -7., 0.0, 0, 1};		// Reihenfolge: network, fxpkt, S, B, Rnum, Y, T, Tchoice, d, x, M, Z
	
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

	printf("Z = %i\n",nicheweb.Z);
	nicheweb.migrPara = gsl_vector_calloc(nicheweb.Z + 6);	 
		 
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
//--Simulation---------------------------------------------------------------------------------------------------
	nicheweb.Tchoice = nicheweb.T;
	nicheweb.T = 0;
	//int len	= ((nicheweb.Rnum+nicheweb.S)*(nicheweb.S+nicheweb.Rnum)+1+nicheweb.Y*nicheweb.Y+1+(nicheweb.Rnum+nicheweb.S)+nicheweb.S);	// Länge des Rückabewerts

	gsl_vector *populationFIN 	= gsl_vector_calloc((nicheweb.Rnum + nicheweb.S)*(nicheweb.Y)*5 + (nicheweb.S) + 3);				// Gleiche Länge wie Rückgabe von evolveNetwork
	gsl_vector *robustness		= gsl_vector_calloc(63);

	gsl_vector *robustnesstemp	= gsl_vector_calloc(63);
	gsl_vector *meanOfDataSquAll 	= gsl_vector_calloc(63);
	gsl_vector *meanSquOfDataAll 	= gsl_vector_calloc(63);
	gsl_vector *meanSquOfDataAlltemp = gsl_vector_calloc(63);
	
	gsl_vector *standardDeviationAll = gsl_vector_calloc(63);
	
	gsl_vector *meanOfData	= gsl_vector_calloc((6*4+2)*nicheweb.Y);
	gsl_vector *meanOfDatatemp = gsl_vector_calloc((6*4+2)*nicheweb.Y);
	gsl_vector *meanSquOfData = gsl_vector_calloc((6*4+2)*nicheweb.Y);
	gsl_vector *meanSquOfDatatemp = gsl_vector_calloc((6*4+2)*nicheweb.Y);
	gsl_vector *meanOfDataSqu = gsl_vector_calloc((6*4+2)*nicheweb.Y);
	gsl_vector *standardDeviation = gsl_vector_calloc((6*4+2)*nicheweb.Y);
	
	gsl_vector_set_zero(robustness);
	gsl_vector_set_zero(meanOfData);
	gsl_vector_set_zero(meanSquOfData); 
	gsl_vector_set_zero(meanSquOfDatatemp);
	gsl_vector_set_zero(nicheweb.migrPara);
	gsl_vector_set_zero(meanOfDataSquAll);
	gsl_vector_set_zero(meanSquOfDataAll);
	
	double SpeciesNumber[L]; 
	
	double ymigr = 0;
	double mu = 0;
	double nu = 0;
	double ymigrtemp;
	double ymigrSqu = 0;
	double ymigrDeviation;
	//printf("test branch bla 2");
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
		nicheweb.network = SetNicheNetwork(nicheweb, res, rng1, rng1_T);
		populationFIN	 = EvolveNetwork(nicheweb, rng1, rng1_T);
			
												
		gsl_vector_memcpy(robustnesstemp, EvaluateRobustness(populationFIN, nicheweb, patchwise));	// Robustness Analyse
		
// 		int length = 0;
// 		length = NumberOfItemsOfObject(robustnesstemp, length);
// 		printf("length ist %i\n",length);
		
//--Standardabweichung für Mittelung vorbereiten-----------------------------------------------------------------------------------------		
		determineMean(robustnesstemp, 63, robustness);
		determineMeanSqu(robustnesstemp, 63, meanSquOfDataAll);
		
//--Ausgabewerte----------------------------------------------------------------------------------------------------------		
		ymigrtemp = gsl_vector_get(nicheweb.migrPara, nicheweb.Z+4);
		mu += gsl_vector_get(nicheweb.migrPara, nicheweb.Z);
		nu += gsl_vector_get(nicheweb.migrPara, nicheweb.Z+1);
		SpeciesNumber[i] = gsl_vector_get(nicheweb.migrPara,nicheweb.Z+2);
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
	standardDeviation = determineStandardDeviation((6*4+2)*nicheweb.Y, meanOfData, meanSquOfData, L);
	standardDeviationAll = determineStandardDeviation(63, robustness, meanSquOfDataAll, L);
	 //printf("der 3. Eintrag in standardDeviationAll ist %f\n", gsl_vector_get(standardDeviationAll,3));
// 	 printf("S ist %f\n", gsl_vector_get(robustness,3));
// 	 printf("Standardabweichung von S ist %f\n", gsl_vector_get(standardDeviationAll,3));
// 	 printf("meanOfDataSqu ist %f\n", gsl_vector_get(meanOfDataSquAll,3));
// 	 printf("meanSquOfData ist %f\n", gsl_vector_get(meanSquOfDataAll,3));
	 
	 
	 
	printf("L=%i\tspeciesini=%f\tspeciesfinal=%f\n", L, gsl_vector_get(robustness, 3)/L, gsl_vector_get(robustness, 9)/L);
	

	
	
//--Abspeichern in File-------------------------------------------------------------------------------------	
	
	FILE *statistics;				// neuer FILE pointer
    char aims[255] = ORT;			// Ausgewähltes Verzeichnis
    char buffers[100];				// Speicher für Dateiname

    sprintf(buffers,"S%dB%d_M%d_x%1.1fY%dd%2.1fT%dL%dRSize%3.1f.out",nicheweb.S,nicheweb.B,nicheweb.M,nicheweb.x,nicheweb.Y,nicheweb.d,nicheweb.Tchoice,L,res.size);		
	// sprintf: schreibt eine Zeichenkette in den Speicherbereich von buffers

    statistics = fopen(strcat(aims, buffers),"w");											// strcat: klebt zwei Strings aneinander (buffers an aims) -> Pfad+Name
    // fopen(*filename, "w") erzeugt eine neue Datei in die geschrieben werden kann. Existiert schon eine Datei dieses Namens wird diese überschrieben.


      fprintf(statistics,"RSize\tS\tB\tM\tx\tY\tdpow\tT\tRob\tPerlok\tPerges\tSi_ges\tSi_TL1\tSi_TL2\tSi_TL3\tSi_TL4\tSi_TL>4\tSf_ges\tSf_TL1\tSf_TL2\tSf_TL3\tSf_TL4\tSf_TL>4\tBi_ges\tBi_TL1\tBi_TL2\tBi_TL3\tBi_TL4\tBi_TL>4\tBf_ges\tBf_TL1\tBf_TL2\tBf_TL3\tBf_TL4\tBf_TL>4\tSh_ges\tSh_TL1\tSh_TL2\tSh_TL3\tSh_TL4\tSh_TL>4\tBh_ges\tBh_TL1\tBh_TL2\tBh_TL3\tBh_TL4\tBh_TL>4\tSs_ges\tSs_TL1\tSs_TL2\tSs_TL3\tSs_TL4\tSs_TL>4\tBs_ges\tBs_TL1\tBs_TL2\tBs_TL3\tBs_TL4\tBs_TL>4\t1mit2\t2mit3\t3mit1\tFixp0\tFixp1\tFixp2\tFixp3\tFixp4\tFixp5\tFixp6\tFixp7\tRob2\tydotMigration\tmu\tnu\tTchoice\n");

    fclose(statistics);
    statistics = fopen(aims,"a");												// fopen(*filename, a): schreibt am Ende der Datei weiter		

//--Daten in Datei schreiben---------------------------------------------------------------------------------------------------------------------------------				
    fprintf(statistics,"%5.1f\t%d\t%d\t%d\t%2.1f\t%d\t%2.1f\t%d\t", res.size, nicheweb.S, nicheweb.B, nicheweb.M, nicheweb.x, nicheweb.Y, nicheweb.d, nicheweb.T);		// Konsolenparameter im Namen

    for(i=0; i<51; i++)															// Rob, Perlok, Perges; +48 Elemente Spezies Informationen
      {
	gsl_vector_set(robustness, i, gsl_vector_get(robustness, i)/L);			// Scale (1/L)
        fprintf(statistics,"%5.3f\t", gsl_vector_get(robustness, i));			// %: Charakter; number.number: min Anzahl an Stellen die gedruckt werden; f: dezimal float
      }

	// muss hier skaliert werden oder nicht?
    for(i=51; i<62; i++)	//1mit2... Fixp1...7
      {
	fprintf(statistics,"%5.0f\t", gsl_vector_get(robustness, i));
      }

    for(i=62; i<63; i++)	//Rob2 -> regio 
      {
	gsl_vector_set(robustness, i, gsl_vector_get(robustness, i)/L);
        fprintf(statistics,"%5.3f\t", gsl_vector_get(robustness, i));
      }

      /* ymigr wurde oben schon durch L geteilt */
      fprintf(statistics,"%7.6f\t", ymigr);
      fprintf(statistics, "%5.2f\t",mu/L);
      fprintf(statistics, "%5.2f\t",nu/L);
      fprintf(statistics, "%d\t",nicheweb.Tchoice);
    
    
    fprintf(statistics,"\n");
    
//--Standardabweichung hinzufügen-------------------------------------------------------------------------------------------------    
    for(i=0; i<8; i++)
    {
	  fprintf(statistics,"%d\t",0);
    }
    
    
    for(i = 0 ; i< 63; i++)
    {
      fprintf(statistics, "%5.3f\t",gsl_vector_get(standardDeviationAll,i));
    }
    
    fprintf(statistics,"%7.6f\t", ymigrDeviation);
    fprintf(statistics, "%d\t",0);
    fprintf(statistics, "%d\t",0);
    
    
    fclose(statistics);															// Datei schließen

	
    
	
//--Daten patchweise abspeichern----------------------------------------------------------------------	
	
     //char name[100];
     //char aims2[255] = ORT2;
     //strcat(aims2, name);
     
     //printf("aims2 %s\n", aims2);
     
     for(int l = 0 ; l< nicheweb.Y; l++)
     {
       char name[100];
       char aims2[255] = ORT2;
       //printf("aims2 %s\n", aims2);
       FILE *statForPatchl;
       sprintf(name,"Patch_l%dS%dB%d_M%d_x%1.1fY%dd%2.1fT%dL%dRSize%3.1f.out",l,nicheweb.S,nicheweb.B,nicheweb.M,nicheweb.x,nicheweb.Y,nicheweb.d,nicheweb.Tchoice,L,res.size);		
	// sprintf: schreibt eine Zeichenkette in den Speicherbereich von buffers
	
	statForPatchl = fopen(strcat(aims2, name),"w");											// strcat: klebt zwei Strings aneinander (buffers an aims) -> Pfad+Name
	// fopen(*filename, "w") erzeugt eine neue Datei in die geschrieben werden kann. Existiert schon eine Datei dieses Namens wird diese überschrieben.

	//printf("aims2 unten %s\n", aims2);
	
	fprintf(statForPatchl,"l\tRSize\tS\tB\tM\tx\tY\tdpow\tT\tSi_ges\tSi_TL1\tSi_TL2\tSi_TL3\tSi_TL4\tSi_TL>4\tSf_ges\tSf_TL1\tSf_TL2\tSf_TL3\tSf_TL4\tSf_TL>4\tBi_ges\tBi_TL1\tBi_TL2\tBi_TL3\tBi_TL4\tBi_TL>4\tBf_ges\tBf_TL1\tBf_TL2\tBf_TL3\tBf_TL4\tBf_TL>4\tRob\tRob2\n");
    
	fprintf(statForPatchl,"%d\t%5.1f\t%d\t%d\t%d\t%2.1f\t%d\t%2.1f\t%d\t",l, res.size, nicheweb.S, nicheweb.B, nicheweb.M, nicheweb.x, nicheweb.Y, nicheweb.d, nicheweb.T);		// Konsolenparameter im Namen

	//tempo = patchwise[l];
 
	for(i = (0 + (6*4+2)*l); i< ((6*2)+(6*4+2)*l); i++)
	 {
	   fprintf(statForPatchl,"%5.3f\t",gsl_vector_get(meanOfData,i));      
	 }
	 
	for(i = ((6*2)+(6*4+2)*l); i< ((6*4+2)+(6*4+2)*l); i++)
	 {
	   fprintf(statForPatchl,"%6.5f\t",gsl_vector_get(meanOfData,i));      
	 }
	 
	 fprintf(statForPatchl,"\n");
	
	/* Standardabweichung für alle Werte dazuspeichern */
	for(i=0; i<9; i++)
	{
	  fprintf(statForPatchl,"%d\t",0);
	}
	  
	for(i = (0 + (6*4+2)*l); i< ((6*2)+(6*4+2)*l); i++)
	 {
	   fprintf(statForPatchl,"%5.3f\t",gsl_vector_get(standardDeviation,i));      
	 }
	 
	for(i = ((6*2)+(6*4+2)*l); i< ((6*4+2)+(6*4+2)*l); i++)
	 {
	   fprintf(statForPatchl,"%6.5f\t",gsl_vector_get(standardDeviation,i));      
	 }
	 
	 printf("\n Es wird Datei für Patch %i erstellt\n",l);
      }
      
//--Ausgewählte Spezies rausschreiben, die migrieren darf---------------------------------------------------------------------------
	FILE *SpeciesNumbers;
	char aims3[255] = ORT;
	char buffers3[100];
	
	printf("\n Es wird Ausgabe erstellt, welche Spezies zum Migrieren ausgewählt wurden\n");
	
	sprintf(buffers3,"SpeciesNumbersS%dB%d_M%d_x%1.1fY%dd%2.1fT%dL%dRSize%5.1f.out",nicheweb.S,nicheweb.B,nicheweb.M,nicheweb.x,nicheweb.Y,nicheweb.d,nicheweb.T,L,res.size);
      
	// sprintf: schreibt eine Zeichenkette in den Speicherbereich von buffers

	SpeciesNumbers = fopen(strcat(aims3, buffers3),"w");											// strcat: klebt zwei Strings aneinander (buffers an aims) -> Pfad+Name
	// fopen(*filename, "w") erzeugt eine neue Datei in die geschrieben werden kann. Existiert schon eine Datei dieses Namens wird diese überschrieben.

	for(i = 0; i<L; i++)
	{
	  fprintf(SpeciesNumbers,"%5.1f\t",SpeciesNumber[i]);
	}
	
	fprintf(SpeciesNumbers,"\n");
	fclose(SpeciesNumbers);
	
	
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
	gsl_vector_free(populationFIN);
	gsl_vector_free(robustness);	
	gsl_vector_free(meanOfData);
	gsl_vector_free(meanOfDatatemp);
	gsl_vector_free(meanSquOfData);
	gsl_vector_free(meanOfDataSqu);
	gsl_vector_free(meanSquOfDatatemp);
	gsl_vector_free(standardDeviation);
	gsl_vector_free(standardDeviationAll);
	gsl_vector_free(meanOfDataSquAll);
	gsl_vector_free(meanSquOfDataAll);
	gsl_vector_free(meanSquOfDataAlltemp);
	gsl_rng_free(rng1);
	
	return(0);

}




















