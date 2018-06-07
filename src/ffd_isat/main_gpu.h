///////////////////////////////////////////////////////////////////////////////
///
/// \file   main.c
///
/// \brief  main entrance of parallel FFD program
///
/// \author Wei Tian
///         University of Miami, Schneider Electric
///         w.tian@umiami.edu, Wei.Tian@Schneider-Electric.com
///         Thomas Sevilla
///         University of Miami
///         t.sevilla@umiami.edu
///
/// \date   7/07/2017
///
///\ This codes have been extensively validated and for the results please refer
///   to the journal paper:
///   Tian, Wei, Thomas Alonso Sevilla, and Wangda Zuo. "A systematic 
///   evaluation of accelerating indoor airflow simulations using cross-platform parallel computing." 
///   Journal of Building Performance Simulation 10.3 (2017): 243-255.
///
///\  All RIGHTS RESERVED.
///////////////////////////////////////////////////////////////////////////////
#ifndef _MAIN_GPU_H
#define _MAIN_GPU_H

#include "timing.h"
#include "data_writer.h"
#include "initialization.h"

///////////////////////////////////////////////////////////////////////////////
/// Main routine of FFD
///
///\para coupled simulation Integer to identify the simulation type
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int ffd_prep(int cosimulation);

///////////////////////////////////////////////////////////////////////////////
/// Allocate memory for variables
///
///\param para Pointer to FFD parameters
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
int allocate_memory(PARA_DATA *para);

#endif
