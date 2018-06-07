///////////////////////////////////////////////////////////////////////////////
///
/// \file   utility.c
///
/// \brief  Some frequently used functions for FFD
///
/// \author Wangda Zuo, Ana Cohen
///         University of Miami
///         W.Zuo@miami.edu
///         Purdue University
///         Mingang Jin, Qingyan Chen
///         Jin55@purdue.edu, YanChen@purdue.edu
///         Wei Tian
///         University of Miami, Schneider Electric
///         w.tian@umiami.edu, Wei.Tian@Schneider-Electric.com
///
/// \date   7/05/2017
/// \add: add a function min_distance to calculate the distance of
///					a fluid cell to the nearest solid boundary condition, which
///					is to be used by Chen's zero equation turbulence model
///
///////////////////////////////////////////////////////////////////////////////

#ifndef _UTILITY_H
#define _UTILITY_H


#include "data_structure.h"


///////////////////////////////////////////////////////////////////////////////
/// Write the log file
///
///\param message Pointer the message
///\param msg_type Type of message
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
void ffd_log(char *message, FFD_MSG_TYPE msg_type);

///////////////////////////////////////////////////////////////////////////////
/// Free memory for BINDEX
///
///\param BINDEX Pointer to the boundary index
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
void free_index(int **BINDEX);

///////////////////////////////////////////////////////////////////////////////
/// Free memory for FFD simulation variables
///
///\param var Pointer to FFD simulation variables
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
void free_data(REAL **var);

REAL average_volume(PARA_DATA *para, REAL **var, REAL *psi);
REAL fluid_volume(PARA_DATA *para, REAL **var);
REAL vol(PARA_DATA *para, REAL **var, int i, int j, int k);

///////////////////////////////////////////////////////////////////////////////
/// Calculate time averaged value
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int average_time(PARA_DATA *para, REAL **var);

///////////////////////////////////////////////////////////////////////////////
/// Check the minimum to the solid boudaries
/// cells used in calculation of zero equation tuebulence model
///
///\param para Pointer to FFD parameters
///\param var Pointer to the variable
///\param
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int min_distance(PARA_DATA *para, REAL **var, int **BINDEX);

///////////////////////////////////////////////////////////////////////////////
/// Check flow rates through all the tiles
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to the boundary index
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int check_tile_flowrate(PARA_DATA *para, REAL **var, int **BINDEX);

REAL vol_inflow(PARA_DATA *para, REAL **var, int **BINDEX);

///////////////////////////////////////////////////////////////////////////////
/// Check flow rates at inlets when t=0
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to the boundary index
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
REAL initial_inflows(PARA_DATA *para, REAL **var, int **BINDEX);

int assign_tile_velocity(PARA_DATA *para, REAL **var, int **BINDEX);

///////////////////////////////////////////////////////////////////////////////////////
/// Calculate the flow rates through the tiles using corrected pressure
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///\p_corr corrected pressure
///\return 0 if no error occurred
///\ Wei Tian
///\ 09/05/2017
//////////////////////////////////////////////////////////////////////////////////////
REAL pressure_correction(PARA_DATA *para, REAL **var, int **BINDEX, REAL p_corr);

REAL area_xy(PARA_DATA *para, REAL **var, int i, int j, int k);

REAL area_yz(PARA_DATA *para, REAL **var, int i, int j, int k);

REAL area_zx(PARA_DATA *para, REAL **var, int i, int j, int k);

REAL length_x(PARA_DATA *para, REAL **var, int i, int j, int k);

REAL length_y(PARA_DATA *para, REAL **var, int i, int j, int k);

REAL length_z(PARA_DATA *para, REAL **var, int i, int j, int k);

///////////////////////////////////////////////////////////////////////////////
/// Calculate the area of inlet or outlet of rack
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///\param A Pointer to the array of area
///
///\return 0 if no error occurred
///\ Version 1.0
///\ Wei Tian, 1-20-2018, Wei.Tian@Schneider-Electric.com
///////////////////////////////////////////////////////////////////////////////
int rack_fluid_area(PARA_DATA *para, REAL **var, int **BINDEX);

///////////////////////////////////////////////////////////////////////////////////////
/// The black box model of rack, which treat the rack as a box with inlet outlet and heat dissipation
/// The temperature stratification of inlet temperature is kept in the outlet temperature
/// The velocity at inlet and outlet is the same
/// The inlet of rack is treated as outlet for the DC room while the outlet of rack is treated as inlet for DC room
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///\return 0 if no error occurred
///\ Version 1.0
///\ Wei Tian, 1-20-2018, Wei.Tian@Schneider-Electric.com
//////////////////////////////////////////////////////////////////////////////////////
int rack_model_black_box(PARA_DATA *para, REAL **var, int **BINDEX);

#endif
