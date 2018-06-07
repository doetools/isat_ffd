#ifndef COMMON_H_INCLUDED

#include<stdio.h>
#include <stdlib.h>
#include <windows.h>
#include <time.h>
#include <math.h>

#include <string.h>  /* strcpy */
#include "uthash.h"
 
/******************************************************************************
| Set input and output array size
|
| Note: usually only change nx_SIZE and nf_SIZE
******************************************************************************/
#define READ_FFD_RESULT 0
#define nx_SIZE 2                   // Dimension of x
#define nf_SIZE 1                   // Dimension of f
#define nh_SIZE 1                   // Dimension of h, Set to 1 if h(x) is not required to initialize pointer ha

//#define totalTimeSteps 8000               // Total simulation steps of simulation
//#define skipSteps 0               // Time steps to skiped befere record simulation result at each step

/****************************************************************************
| fixme
****************************************************************************/
typedef enum{MPC_WARNING, MPC_ERROR, MPC_NORMAL, MPC_NEW} MPC_MSG_TYPE;

/****************************************************************************
| Hash Table struct
****************************************************************************/
typedef struct {
  double x[nx_SIZE];
} hashKey; //hash key

typedef struct {
  hashKey key;           /* we'll use this field as the key */
  /* ... other data ... */
  UT_hash_handle hh; /* makes this structure hashable */
} hashStruct; //hash unit



typedef struct{
  const int *nx;
  const int *nf;
  const int *ng;
  double *x;
  double *f;
  double *g;
} ffdIO;


char logMsg[1000];

#define COMMON_H_INCLUDED
#endif



