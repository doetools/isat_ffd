
#include "main_gpu.h"
#include "common.h"
#include "io.h"

#ifndef FFD_ISAT_H_INCLUDED
#define FFD_ISAT_H_INCLUDED

void ffd_ISAT (int need[], double x[], double f[], double g[][nf_SIZE], void *p);
int read_existing();

#endif
