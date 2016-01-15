#ifndef STANDARDDEVIATION_H
#define STANDARDDEVIATION_H

#include <gsl/gsl_vector.h>

int NumberOfItemsOfObject(gsl_vector* object, int count);	
gsl_vector* determineMean(gsl_vector* object, int length, gsl_vector* mean);
gsl_vector* determineMeanSqu(gsl_vector* object, int length, gsl_vector* meanSqu);
gsl_vector* determineStandardDeviation(int length, gsl_vector* meanOfData, gsl_vector* meanSquOfData, int L, gsl_vector* );
int linkElements(struct data patchwise[], int Y, gsl_vector* meanOfDatatemp);

#endif
