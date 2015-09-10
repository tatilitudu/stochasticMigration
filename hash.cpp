#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>

#include "hash.h"
#include "structs.h"					// foodweb Struktur

//--Verzeichnis für Ergebnisse-------------------------------------------------------
#define ORT "/home/tatjana/Arbeitsfläche/MichaelasProgramm/stochastischeMigration/ErsteVersuche/"

using namespace std;

hash::hash()
{
  int i;
  for(i = 0 ; i< tableSize; i++)
  {
    HashTable[i] = new item;
    HashTable[i]->l = 0;
    HashTable[i]->speciesini = NULL;
    HashTable[i]->speciesfini = NULL;
    HashTable[i]->biomassini = NULL;
    HashTable[i]->biomassfini = NULL;
    HashTable[i]->next = NULL;
  }

}

 void hash::AddItem(int l, double* speciesini, double* speciesfini, double* biomassini, double* biomassfini)
{
  int index = Hash(l);
  
  if(HashTable[index]->l == 0)
  {
    HashTable[index]->l = l;
    HashTable[index]->speciesini = speciesini;
    HashTable[index]->speciesfini = speciesfini;
    HashTable[index]->biomassini = biomassini;
    HashTable[index]->biomassfini = biomassfini;
  }
  else
  {
    item* Ptr = HashTable[index];
    item* n = new item;
    n->l = l;
    n->speciesini = speciesini;
    n->speciesfini = speciesfini;
    n->biomassini = biomassini;
    n->biomassfini = biomassfini;
    n->next = NULL;
    while(Ptr->next != NULL)
    {
      Ptr = Ptr->next;
    }
    Ptr->next = n;
  }
  
}


int hash::NumberOfItemsInIndex(int index)
{
  int count = 0;
  
  if(HashTable[index]->l == 0)
  {
    return count;
  }
  else
  {
    count++;
    item* Ptr = HashTable[index];
    while(Ptr->next != NULL)
    {
      count++;
      Ptr = Ptr->next;
    }
  }
  return count;
}


void hash::PrintTable()
{
  int number;
  int i;
  
  for(i = 0; i< tableSize; i++)
  {
    number = NumberOfItemsInIndex(i);
    printf("--------------------------\n");
    printf("index = %i\n",i);
    cout << HashTable[i]->l<< endl;
    cout << HashTable[i]->speciesini[0] << endl;
    cout << HashTable[i]->speciesfini[0] << endl;
    cout << HashTable[i]->biomassini[0] << endl;
    cout << HashTable[i]->biomassfini[0] << endl;
    printf("# of items = %i\n",number);
    printf("--------------------------\n");
  }
}

 void hash::PrintItemsInIndex(int index)
 {
   item* Ptr = HashTable[index];
   
   if(Ptr->l == 0)
   {
     cout << "index = " << index << " is empty" << endl;
   }
   else
   {
     cout << "index " << index << " contains the following items \n";
     
     while(Ptr != NULL)
     {
       cout << "--------------------\n";
       cout << Ptr->l << endl;
       cout << Ptr->speciesini << endl;
       cout << Ptr->speciesfini << endl;
       cout << Ptr->biomassini << endl;
       cout << Ptr->biomassfini << endl;
       cout << "--------------------\n";
       Ptr = Ptr->next;
     }
     
   }
   
   
 }

 
 void hash::PrintData(struct foodweb nicheweb, struct resource res, int l, int L)
 {

   item* Ptr = HashTable[l];
   
   if(Ptr->l == 0)
   {
     printf("\nBucket ist leer - es wird keine Datei erstellt\n\n");
   }
   else
   {
     FILE *statForPatchl;
     char aims[255] = ORT;
     char buffers[100];
     
     sprintf(buffers,"Patch_l%iS%dB%d_M%d_x%1.1fY%dd%2.1fT%dL%dRSize%3.1f.out",l,nicheweb.S,nicheweb.B,nicheweb.M,nicheweb.x,nicheweb.Y,nicheweb.d,nicheweb.T,L,res.size);		
	// sprintf: schreibt eine Zeichenkette in den Speicherbereich von buffers

    statForPatchl = fopen(strcat(aims, buffers),"w");											// strcat: klebt zwei Strings aneinander (buffers an aims) -> Pfad+Name
    // fopen(*filename, "w") erzeugt eine neue Datei in die geschrieben werden kann. Existiert schon eine Datei dieses Namens wird diese überschrieben.

    fprintf(statForPatchl,"RSize\tS\tB\tM\tx\tY\tdpow\tT\tRob\tRob2\tSi_ges\tSi_TL1\tSi_TL2\tSi_TL3\tSi_TL4\tSi_TL>4\tSf_ges\tSf_TL1\tSf_TL2\tSf_TL3\tSf_TL4\tSf_TL>4\tBi_ges\tBi_TL1\tBi_TL2\tBi_TL3\tBi_TL4\tBi_TL>4\tBf_ges\tBf_TL1\tBf_TL2\tBf_TL3\tBf_TL4\tBf_TL>4\n");
    
    fprintf(statForPatchl,"%5.1f\t%d\t%d\t%d\t%2.1f\t%d\t%2.1f\t%d\t", res.size, nicheweb.S, nicheweb.B, nicheweb.M, nicheweb.x, nicheweb.Y, nicheweb.d, nicheweb.T);		// Konsolenparameter im Namen

    int i;
    for(i = 0; i<6; i++)
    {
      fprintf(statForPatchl,"%5.3f\t",HashTable[l]->speciesini[i]);
      
    }
    printf("\n Es wird Datei für Patch %i erstellt\n\n",l);
  }
 }
   

 int hash::Hash(int l)
 {
   //int hash = 0;
   int index;
   
   //index = key.length();
   
//    cout << "key[0] = " << (int)key[0] <<endl;
//    cout << "key[1] = " << (int)key[1] <<endl;
//    cout << "key[2] = " << (int)key[2] <<endl;
//    cout << "key[3] = " << (int)key[3] <<endl;
  
  
//   unsigned int i;
//   for(i=0; i< key.length(); i++)
//   {
//    hash = hash + (int)key[i]; 
//    //cout << "hash = " << hash << endl;
//   }
  
  index = l;
  
   return index;
 }