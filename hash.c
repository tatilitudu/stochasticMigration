#include <stdlib.h>
#include <iostream>
#include <string>

#include "hash.h"

using namespace std;

hash::hash()
{
  int i;
  for(i = 0 ; i< tableSize; i++)
  {
    HashTable[i] = new item;
    HashTable[i]->name = "empty";
    HashTable[i]->drink = 0;
    HashTable[i]->next = NULL;
  }

}

void hash::AddItem(string name, double drink)
{
  int index = Hash(name);
  
  if(HashTable[index]->name == "empty")
  {
    HashTable[index]->name = name;
    HashTable[index]->drink = drink;
  }
  else
  {
    item* Ptr = HashTable[index];
    item* n = new item;
    n->name = name;
    n->drink = drink;
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
  
  if(HashTable[index]->name == "empty")
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
    cout << "--------------------------\n";
    cout << "index = " << i << endl;
    cout << HashTable[i]->name << endl;
    cout << HashTable[i]->drink << endl;
    cout << "# of items = " << number << endl;
    cout << "--------------------------\n";
  }
}

 void hash::PrintItemsInIndex(int index)
 {
   item* Ptr = HashTable[index];
   
   if(Ptr->name == "empty")
   {
     cout << "index = " << index << " is empty" << endl;
   }
   else
   {
     cout << "index " << index << " contains the following items \n";
     
     while(Ptr != NULL)
     {
       cout << "--------------------\n";
       cout << Ptr->name << endl;
       cout << Ptr->drink << endl;
       cout << "--------------------\n";
       Ptr = Ptr->next;
     }
     
   }
   
   
 }


 int hash::Hash(string key)
 {
   int hash = 0;
   int index;
   
   //index = key.length();
   
//    cout << "key[0] = " << (int)key[0] <<endl;
//    cout << "key[1] = " << (int)key[1] <<endl;
//    cout << "key[2] = " << (int)key[2] <<endl;
//    cout << "key[3] = " << (int)key[3] <<endl;
  
  
  int i;
  for(i=0; i< key.length(); i++)
  {
   hash = hash + (int)key[i]; 
   //cout << "hash = " << hash << endl;
  }
  
  index = hash % tableSize;
  
   return index;
 }