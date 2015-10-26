#include "structs.h"



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>						// string modification functions

#include <gsl/gsl_rng.h>					// random number generator functions
#include <gsl/gsl_randist.h>				// random number distributions
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>



int createOutputGeneral(struct foodweb nicheweb, struct resource res, char* aims, gsl_vector* robustness, gsl_vector* standardDeviationAll, int L, double mu, double nu, double ymigr, double ymigrDeviation)
{
  int i;
  FILE *statistics;
  
  

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
    
    
    fclose(statistics);		
    
    return 0;
    
}


int createOutputPatchwise(struct foodweb nicheweb, struct resource res, char* aims, gsl_vector* meanOfData, gsl_vector* standardDeviation, int L, int l)
{
  int i;
  char name[100];
  
  FILE *statForPatchl;
  
  sprintf(name,"Patch_l%dS%dB%d_M%d_x%1.1fY%dd%2.1fT%dL%dRSize%3.1f.out",l,nicheweb.S,nicheweb.B,nicheweb.M,nicheweb.x,nicheweb.Y,nicheweb.d,nicheweb.Tchoice,L,res.size);		
  // sprintf: schreibt eine Zeichenkette in den Speicherbereich von buffers
  
  statForPatchl = fopen(strcat(aims, name),"w");											// strcat: klebt zwei Strings aneinander (buffers an aims) -> Pfad+Name
  // fopen(*filename, "w") erzeugt eine neue Datei in die geschrieben werden kann. Existiert schon eine Datei dieses Namens wird diese überschrieben.

  //printf("aims2 unten %s\n", aims2);
	
  fprintf(statForPatchl,"l\tRSize\tS\tB\tM\tx\tY\tdpow\tT\tSi_ges\tSi_TL1\tSi_TL2\tSi_TL3\tSi_TL4\tSi_TL>4\tSf_ges\tSf_TL1\tSf_TL2\tSf_TL3\tSf_TL4\tSf_TL>4\tBi_ges\tBi_TL1\tBi_TL2\tBi_TL3\tBi_TL4\tBi_TL>4\tBf_ges\tBf_TL1\tBf_TL2\tBf_TL3\tBf_TL4\tBf_TL>4\tRob\tRob2\n");
    
  fprintf(statForPatchl,"%d\t%5.1f\t%d\t%d\t%d\t%2.1f\t%d\t%2.1f\t%d\t",l, res.size, nicheweb.S, nicheweb.B, nicheweb.M, nicheweb.x, nicheweb.Y, nicheweb.d, nicheweb.T);		// Konsolenparameter im Namen

  //tempo = patchwise[l];
 
  for(i = (0 + (6*4+2)*l); i< ((6*2)+(6*4+2)*l); i++)
  {
    fprintf(statForPatchl,"%5.3f\t",gsl_vector_get(meanOfData,i)/L);      
    
  }
	 
  for(i = ((6*2)+(6*4+2)*l); i< ((6*4+2)+(6*4+2)*l); i++)
  {
    fprintf(statForPatchl,"%6.5f\t",gsl_vector_get(meanOfData,i)/L);      
    
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
  fclose(statForPatchl);	
  
  return 0;
  
}


int createOutputSpeciesNumber(struct foodweb nicheweb, struct resource res, char* aims, double SpeciesNumber[], int L)
{
  int i;
  
  FILE *SpeciesNumbers;
  char buffers3[100];
	
  printf("\n Es wird Ausgabe erstellt, welche Spezies zum Migrieren ausgewählt wurden\n");
	
  sprintf(buffers3,"SpeciesNumbersS%dB%d_M%d_x%1.1fY%dd%2.1fT%dL%dRSize%5.1f.out",nicheweb.S,nicheweb.B,nicheweb.M,nicheweb.x,nicheweb.Y,nicheweb.d,nicheweb.T,L,res.size);
      
  // sprintf: schreibt eine Zeichenkette in den Speicherbereich von buffers

  SpeciesNumbers = fopen(strcat(aims, buffers3),"w");											// strcat: klebt zwei Strings aneinander (buffers an aims) -> Pfad+Name
  // fopen(*filename, "w") erzeugt eine neue Datei in die geschrieben werden kann. Existiert schon eine Datei dieses Namens wird diese überschrieben.

  for(i = 0; i<nicheweb.Z*L; i++)
  {
    fprintf(SpeciesNumbers,"%5.1f\t",SpeciesNumber[i]);
  }
	
  fprintf(SpeciesNumbers,"\n");
  fclose(SpeciesNumbers);
  
  return 0;
  
}



int createOutputPatchlink(struct foodweb nicheweb, struct resource res, char* aims, double AllMu[], double AllNu[], int L)
{
  int i;
  
  FILE *Patchlink;
  char buffers4[100];
	
  printf("\n Es wird Ausgabe erstellt, welcher Link zum Migrieren ausgewählt wurde\n");
	
  sprintf(buffers4,"PatchlinkS%dB%d_M%d_x%1.1fY%dd%2.1fT%dL%dRSize%5.1f.out",nicheweb.S,nicheweb.B,nicheweb.M,nicheweb.x,nicheweb.Y,nicheweb.d,nicheweb.T,L,res.size);
      
  // sprintf: schreibt eine Zeichenkette in den Speicherbereich von buffers

  Patchlink = fopen(strcat(aims, buffers4),"w");											// strcat: klebt zwei Strings aneinander (buffers an aims) -> Pfad+Name
  // fopen(*filename, "w") erzeugt eine neue Datei in die geschrieben werden kann. Existiert schon eine Datei dieses Namens wird diese überschrieben.

  for(i = 0; i<nicheweb.Z*L; i++)
  {
    fprintf(Patchlink,"%5.1f\t",AllMu[i]);
  }
	
  fprintf(Patchlink,"\n");
  
  for(i = 0; i<nicheweb.Z*L; i++)
  {
    fprintf(Patchlink,"%5.1f\t",AllNu[i]);
  }
	
  fprintf(Patchlink,"\n");
  
  fclose(Patchlink);
  
  return 0;
  
}