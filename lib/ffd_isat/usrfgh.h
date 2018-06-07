#include"ffd_ISAT.h"
#include "io.h"

#ifndef USRFGH_H_INCLUDED
#define USRFGH_H_INCLUDED

void USRFGH(int need[], int *nx, double x[], int *nf, int *nh, int iusr[], double rusr[], double f[], double** g, double h[]);
void numericalDifferentiation (double g[][nf_SIZE]);

#endif