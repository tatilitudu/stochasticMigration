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


double* stochMigration(struct foodweb nicheweb, const double[], gsl_rng* rng1, const gsl_rng_type* rng1_T);
int select_patch(gsl_vector*, double, gsl_rng* rng1, int);
int select_species(struct foodweb nicheweb, gsl_rng* rng1, int, const double[]);
double choose_time(double, gsl_rng* rng1);

#endif