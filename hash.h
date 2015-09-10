#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>

#include "structs.h"					// foodweb Struktur

using namespace std;


#ifndef HASH_H
#define HASH_H

class hash{
 private:
   static const int tableSize = 4;	// erlaubt die Größe direkt hier zu definieren und zu bestimmen; 
					//es kann auch nur int genommen und später definiert werden
  
  struct item{
    int l;
    double* speciesini;
    double* speciesfini;
    double* biomassini;
    double* biomassfini;
    item* next;
  };
  
  item* HashTable[tableSize];
  
 public:
  hash();
  int Hash(int l);
  void AddItem(int l, double* speciesini, double* speciesfini, double* biomassini, double* biomassfini);
  int NumberOfItemsInIndex(int index);
  void PrintTable();
  void PrintData(struct foodweb nicheweb, struct resource res, int l, int L);
  void PrintItemsInIndex(int index);
  
  
};

#endif