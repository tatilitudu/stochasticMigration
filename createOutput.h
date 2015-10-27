#ifndef CREATEOUTPUT_H
#define CREATEOUTPUT_H

#include <gsl/gsl_vector.h>

int createOutputGeneral(struct foodweb nicheweb, struct resource res, char* aims, gsl_vector* robustness, gsl_vector* standardDeviationAll, int L, double mu, double nu, double ymigr, double ymigrDeviation);
int createOutputPatchwise(struct foodweb nicheweb, struct resource res, char* aims, gsl_vector* robustness, gsl_vector* standardDeviationAll, int L, int l);
int createOutputSpeciesNumber(struct foodweb nicheweb, struct resource res, char* aims, double SpeciesNumber[][2], int L);
int createOutputPatchlink(struct foodweb nicheweb, struct resource res, char* aims, double AllMu[][2], double AllNu[][2], int L);

#endif