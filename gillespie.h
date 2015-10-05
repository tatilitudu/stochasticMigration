#ifndef GILLESPIE_H
#define GILLESPIE_H

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


double* stochMigration(struct foodweb nicheweb, double*, const double[]);
int select_patch(gsl_vector*, double, double, int);
int select_species(struct foodweb nicheweb, double r, int Choice, const double y[]);

#endif