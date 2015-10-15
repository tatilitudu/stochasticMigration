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




int NumberOfItemsOfObject(gsl_vector* object, int count)
{
  printf("Wert ist %f\n",gsl_vector_get(object, count+1));
  while( gsl_vector_get(object,count+1) != false)
  {
    printf("in while-schleife");
    count++;
  }
  return count;
}


gsl_vector* determineMean(gsl_vector* object, int length, gsl_vector* mean)
{
  
  gsl_vector_add(mean, object);
  
  return 0;
  
}


gsl_vector* determineMeanSqu(gsl_vector* object, int length, gsl_vector* meanSqu)
{
  gsl_vector* meanSqutemp = gsl_vector_calloc(length);
  
  gsl_vector_memcpy(meanSqutemp, object);
  gsl_vector_mul(meanSqutemp, object);
  
  gsl_vector_add(meanSqu, meanSqutemp);
  
  return 0;
  
}
  


gsl_vector* determineStandardDeviation(int length, gsl_vector* meanOfData, gsl_vector* meanSquOfData, int L )
{
  gsl_vector* standardDeviation = gsl_vector_calloc(length);
  gsl_vector* meanOfDataSqu = gsl_vector_calloc(length);
  
  printf("meanSquOfData ist oben %f\n", gsl_vector_get(meanSquOfData,3));
  gsl_vector_memcpy(meanOfDataSqu, meanOfData);
  printf("meanOfData ist %f\n", gsl_vector_get(meanOfDataSqu,3));
  
  gsl_vector_mul(meanOfDataSqu, meanOfData);

  printf("meanOfDataSqu ist oben %f\n", gsl_vector_get(meanOfDataSqu,3));
  printf("L ist %i\n",L);
  for(int i = 0 ; i < length; i++)
  {
    gsl_vector_set(meanSquOfData,i,gsl_vector_get(meanSquOfData,i)/L);
    gsl_vector_set(meanOfDataSqu,i,gsl_vector_get(meanOfDataSqu,i)/(L*L));
    if(i==3)
    {
      printf("meanOfDataSqu ist %f\n", gsl_vector_get(meanOfDataSqu,i));
      printf("meanSquOfData ist %f\n", gsl_vector_get(meanSquOfData,i));
    }
    gsl_vector_set(standardDeviation, i, sqrt(gsl_vector_get(meanSquOfData,i) - gsl_vector_get(meanOfDataSqu,i)));
  }
  
  return standardDeviation;
  
}
  
  
  
int linkElements(struct data patchwise[], int Y, gsl_vector* meanOfDatatemp)
{
  for(int l = 0; l< Y ; l++)
  {
    for( int j = 0; j < 6; j++)
    {
      gsl_vector_set(meanOfDatatemp, j+(4*6+2)*l,  gsl_vector_get(patchwise[l].sini,j));
      gsl_vector_set(meanOfDatatemp, 6+j+(4*6+2)*l, gsl_vector_get(patchwise[l].sfini,j));
      gsl_vector_set(meanOfDatatemp, 2*6+j+(4*6+2)*l, gsl_vector_get(patchwise[l].bini,j));
      gsl_vector_set(meanOfDatatemp, 3*6+j+(4*6+2)*l, gsl_vector_get(patchwise[l].bfini,j));
    }
    gsl_vector_set(meanOfDatatemp, 4*6+0+(4*6+2)*l,  gsl_vector_get(patchwise[l].robness,0));
    gsl_vector_set(meanOfDatatemp, 4*6+1+(4*6+2)*l,  gsl_vector_get(patchwise[l].robness,1));
  }
  
  return 0;
}
  
  
  
  