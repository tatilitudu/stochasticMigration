#ifndef ROBUSTNESS_H
#define ROBUSTNESS_H

#include <gsl/gsl_vector.h>

gsl_vector *EvaluateRobustness(gsl_vector *, struct foodweb nicheweb, struct data patchwise[], gsl_vector*);	

#endif
