#include "common.h"
#include "io.h"
#include "ffd_ISAT.h"

#ifndef MPC_H_INCLUDED
#define MPC_H_INCLUDED

void update_digit();
void nDemArrEva(int dimension, double xStep[]);
void evaluate ();
void randomCall (int nCall, int useNormalDistribution, int useUnboundedTest);
void accuracyTest ();
void writeRecord ();
double getRandom (int dimension, int useNormalDistribution, int useUnboundedTest);
double randNormal (double mu, double sigma);
void binaryTrain ();
double my_round(double x, int digits);
void findHash();
void addHash();


#endif

//#include "io.h"



