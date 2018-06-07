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

#include "utility.h"
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
//char msg[1000];
FILE *file_log;

///////////////////////////////////////////////////////////////////////////////
/// Check the residual of equation
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param psi Pointer to the variable
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
REAL check_residual(PARA_DATA *para, REAL **var, REAL *x) {
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int i, j, k;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL *aw = var[AW], *ae = var[AE], *as = var[AS], *an = var[AN];
  REAL *ap = var[AP], *ab = var[AB], *af = var[AF], *b = var[B];
  REAL tmp, residual = 0.0;

  FOR_EACH_CELL
    tmp = ap[IX(i,j,k)]*x[IX(i,j,k)]
        - ae[IX(i,j,k)]*x[IX(i+1,j,k)] - aw[IX(i,j,k)]*x[IX(i-1,j,k)]
        - an[IX(i,j,k)]*x[IX(i,j+1,k)] - as[IX(i,j,k)]*x[IX(i,j-1,k)]
        - af[IX(i,j,k)]*x[IX(i,j,k+1)] - ab[IX(i,j,k)]*x[IX(i,j,k-1)]
        - b[IX(i,j,k)];
    residual += tmp * tmp;
  END_FOR

  return residual / (imax*jmax*kmax);

}// End of check_residual( )

///////////////////////////////////////////////////////////////////////////////
/// Write the log file
///
///\param message Pointer the message
///\param msg_type Type of message
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
void ffd_log(char *message, FFD_MSG_TYPE msg_type) {
  if(msg_type==FFD_NEW) {
    if((file_log=fopen("log_gpu.ffd","w+"))==NULL) {
        fprintf(stderr, "Error:can not open error file!\n");
        exit(1);
    }
  }
  else if((file_log=fopen("log_gpu.ffd","a+"))==NULL) {
    fprintf(stderr,"Error:can not open error file!\n");
    exit(1);
  }

  switch(msg_type) {
    case FFD_WARNING:
      fprintf(file_log, "WARNING in %s\n", message);
      break;
    case FFD_ERROR:
      fprintf(file_log, "ERROR in %s\n", message);
      break;
    // Normal log
    default:
      fprintf(file_log, "%s\n", message);
  }
  fclose(file_log);
} // End of ffd_log()

///////////////////////////////////////////////////////////////////////////////
/// Calculate time averaged value
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int average_time(PARA_DATA *para, REAL **var) {
  int i, j, k;
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int kmax = para->geom->kmax;
  int IMAX = imax + 2, IJMAX = (imax + 2)*(jmax + 2);
  int step = para->mytime->step_mean;

  FOR_ALL_CELL
  var[VXM][IX(i, j, k)] = var[VXM][IX(i, j, k)] / step;
  var[VYM][IX(i, j, k)] = var[VYM][IX(i, j, k)] / step;
  var[VZM][IX(i, j, k)] = var[VZM][IX(i, j, k)] / step;
  var[TEMPM][IX(i, j, k)] = var[TEMPM][IX(i, j, k)] / step;
  var[MVM][IX(i, j, k)] = var[MVM][IX(i, j, k)] / step;
  END_FOR

  // Wall surfaces
    for (i = 0; i<para->bc->nb_wall; i++)
      para->bc->temHeaMean[i] = para->bc->temHeaMean[i] / step;

  // Fluid ports
  for (i = 0; i<para->bc->nb_port; i++) {
    para->bc->TPortMean[i] = para->bc->TPortMean[i] / step;
    para->bc->velPortMean[i] = para->bc->velPortMean[i] / step;

    for (j = 0; j<para->bc->nb_Xi; j++)
      para->bc->XiPortMean[i][j] = para->bc->XiPortMean[i][j] / step;
    for (j = 0; j<para->bc->nb_C; j++)
      para->bc->CPortMean[i][j] = para->bc->CPortMean[i][j] / step;
  }

  return 0;
} // End of average_time()

///////////////////////////////////////////////////////////////////////////////
/// Free memory for BINDEX
///
///\param BINDEX Pointer to the boundary index
///
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
void free_index(int **BINDEX) {
  if(BINDEX[0]) free(BINDEX[0]);
  if(BINDEX[1]) free(BINDEX[1]);
  if(BINDEX[2]) free(BINDEX[2]);
  if (BINDEX[3]) free(BINDEX[3]);
  if (BINDEX[4]) free(BINDEX[4]);  
} // End of free_index ()

///////////////////////////////////////////////////////////////////////////////
/// Free memory for FFD simulation variables
///
///\param var Pointer to FFD simulation variables
///
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
void free_data(REAL **var) {
	int nb_var = C2BC + 1;
	int i;
	for (i = 0; i < nb_var; i++) {
			if (var[i]) free(var[i]);
	}
  
	/*
  if(var[X]) free(var[X]);
  if(var[Y]) free(var[Y]);
  if(var[Z]) free(var[Z]);
  if(var[VX]) free(var[VX]);
  if(var[VY]) free(var[VY]);
  if(var[VZ]) free(var[VZ]);
  if(var[VXS]) free(var[VXS]);
  if(var[VYS]) free(var[VYS]);
  if(var[VZS]) free(var[VZS]);
  if(var[VXM]) free(var[VXM]);
  if(var[VYM]) free(var[VYM]);
  if(var[VZM]) free(var[VZM]);
  if(var[TEMP]) free(var[TEMP]);
  if(var[TEMPM]) free(var[TEMPM]);
  if(var[TEMPS]) free(var[TEMPS]);
  if(var[IP]) free(var[IP]);
  if(var[TMP1]) free(var[TMP1]);
  if(var[TMP2]) free(var[TMP2]);
  if(var[TMP3]) free(var[TMP3]);
  if(var[AP]) free(var[AP]);
  if(var[AN]) free(var[AN]);
  if(var[AS]) free(var[AS]);
  if(var[AE]) free(var[AE]);
  if(var[AW]) free(var[AW]);
  if(var[AF]) free(var[AF]);
  if(var[AB]) free(var[AB]);
  if(var[B])  free(var[B]);
  if(var[GX])  free(var[GX]);
  if(var[GY])  free(var[GY]);
  if(var[GZ])  free(var[GZ]);
  if(var[AP0])  free(var[AP0]);
  if(var[PP])  free(var[PP]);
  if(var[FLAGP])  free(var[FLAGP]);
  if(var[FLAGU])  free(var[FLAGU]);
  if(var[FLAGV])  free(var[FLAGV]);
  if(var[FLAGW])  free(var[FLAGW]);
  if(var[LOCMIN])  free(var[LOCMIN]);
  if(var[LOCMAX])  free(var[LOCMAX]);
  if(var[VXBC])  free(var[VXBC]);
  if(var[VYBC])  free(var[VYBC]);
  if(var[VZBC])  free(var[VZBC]);
  if(var[TEMPBC])  free(var[TEMPBC]);
  if(var[Xi1])  free(var[Xi1]);
  if(var[Xi2])  free(var[Xi2]);
  if(var[Xi1BC])  free(var[Xi1BC]);
  if(var[Xi2BC])  free(var[Xi2BC]);
  if(var[C1])  free(var[C1]);
  if(var[C2])  free(var[C2]);
  if(var[C1BC])  free(var[C1BC]);
  if(var[C2BC])  free(var[C2BC]);
  if(var[QFLUXBC])  free(var[QFLUXBC]);
  if(var[QFLUX])  free(var[QFLUX]);
	*/
} // End of free_data()

  ///////////////////////////////////////////////////////////////////////////////
  /// Calculate volume weighted averaged value of psi in a space
  ///
  /// The average is weighted by volume of each cell
  ///
  ///\param para Pointer to FFD parameters
  ///\param var Pointer to FFD simulation variables
  ///\param psi Pointer to the variable
  ///
  ///\return Volume weighted average
  ///////////////////////////////////////////////////////////////////////////////
REAL average_volume(PARA_DATA *para, REAL **var, REAL *psi) {
	int imax = para->geom->imax, jmax = para->geom->jmax;
	int kmax = para->geom->kmax;
	int i, j, k;
	int IMAX = imax + 2, IJMAX = (imax + 2)*(jmax + 2);
	REAL tmp1 = 0, tmp2 = 0, tmp3 = 0;

	if (para->geom->volFlu == 0)
		return 0;
	else {
		FOR_EACH_CELL
			if (var[FLAGP][IX(i, j, k)] == FLUID) {
				//tianwei: temporary add to calculate the half room average temperature
				if (k>kmax / 2) continue;
				tmp1 = vol(para, var, i, j, k);
				tmp2 += psi[IX(i, j, k)] * tmp1;
				tmp3 += tmp1;

			}
			else
				continue;
		END_FOR

			//return tmp2 / para->geom->volFlu;
			//printf("the volume is %f\n", tmp3);
			return tmp2 / tmp3;
	}

}// End of average_volume( )

 ///////////////////////////////////////////////////////////////////////////////
 /// Calculate total fluid volume in the space
 ///
 ///\param para Pointer to FFD parameters
 ///\param var Pointer to FFD simulation variables
 ///
 ///\return Volume weighted average
 ///////////////////////////////////////////////////////////////////////////////
REAL fluid_volume(PARA_DATA *para, REAL **var) {
	int imax = para->geom->imax, jmax = para->geom->jmax;
	int kmax = para->geom->kmax;
	int i, j, k;
	int IMAX = imax + 2, IJMAX = (imax + 2)*(jmax + 2);
	REAL V = 0;

	FOR_EACH_CELL
		if (var[FLAGP][IX(i, j, k)] == FLUID) {
			V += vol(para, var, i, j, k);
		}
		else
			continue;
	END_FOR

		return V;
}// End of fluid_volume( )

 ///////////////////////////////////////////////////////////////////////////////
 /// Calculate the volume of control volume (i,j,k)
 ///
 ///\param para Pointer to FFD parameters
 ///\param var Pointer to FFD simulation variables
 ///\param i I-index of the control volume
 ///\param j J-index of the control volume
 ///\param K K-index of the control volume
 ///
 ///\return Volume
 ///////////////////////////////////////////////////////////////////////////////
REAL vol(PARA_DATA *para, REAL **var, int i, int j, int k) {

	return area_xy(para, var, i, j, k)
		* length_z(para, var, i, j, k);
} // End of vol()


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
int min_distance(PARA_DATA *para, REAL **var, int **BINDEX) {
		int imax = para->geom->imax, jmax = para->geom->jmax;
		int kmax = para->geom->kmax;
		int i, j, k;
		REAL *x = var[X], *y = var[Y], *z = var[Z];
		int it, i_bc, j_bc, k_bc;
		int index = para->geom->index;
		int IMAX = imax + 2, IJMAX = (imax + 2)*(jmax + 2);
		REAL tmp = 1e12; // define a large number
		REAL lx, ly, lz, l;

		// Loop all the fluid cells
		FOR_EACH_CELL
				// initiate the tmp variable after every iteragion
				tmp = 1e12;
		// pass if the cell is not fluid
		if (var[FLAGP][IX(i, j, k)] >= 0) continue;
		// go through all the boudnary conditions and find minimal
		for (it = 0; it < index; it++) {
				i_bc = BINDEX[0][it];
				j_bc = BINDEX[1][it];
				k_bc = BINDEX[2][it];
				if (var[FLAGP][IX(i_bc, j_bc, k_bc)] == INLET || var[FLAGP][IX(i_bc, j_bc, k_bc)] == OUTLET) continue;
				// caculate the distance in each dimension and find the Euler distance l
				lx = fabs(x[IX(i, j, k)] - x[IX(i_bc, j_bc, k_bc)]);
				ly = fabs(y[IX(i, j, k)] - y[IX(i_bc, j_bc, k_bc)]);
				lz = fabs(z[IX(i, j, k)] - z[IX(i_bc, j_bc, k_bc)]);
				l = sqrt(lx*lx + ly*ly + lz*lz);
				// store the minimal value during the looping to tmp
				if (l < tmp) {
						tmp = l;
				}
		}
		// store the minimal value associated with (i,j,k) to global var
		var[MIN_DISTANCE][IX(i, j, k)] = tmp;

		//printf("Distance [%d, %d, %d] is %f\n", i, j, k, tmp);


		END_FOR
				return 0;

}// End of min_distance( )

///////////////////////////////////////////////////////////////////////////////
	/// Check flow rates through all the tiles
	///
	///\param para Pointer to FFD parameters
	///\param var Pointer to FFD simulation variables
	///\param BINDEX Pointer to the boundary index
	///
	///\return 0 if no error occurred
	///////////////////////////////////////////////////////////////////////////////
FILE *FILE_TILE_FLOW;
int check_tile_flowrate(PARA_DATA *para, REAL **var, int **BINDEX) {
		int i, j, k;
		int imax = para->geom->imax, jmax = para->geom->jmax;
		int kmax = para->geom->kmax;
		int IMAX = imax + 2, IJMAX = (imax + 2)*(jmax + 2);
		int it, id;
		int index = para->geom->index;
		REAL *flagp = var[FLAGP];
		int nb_ports = para->bc->nb_port;
		REAL A=0.0, V_tmp=0.0;
		REAL Axy=0.0, Ayz=0.0, Azx=0.0;
		REAL *QPort = para->bc->QPort;
		REAL *u = var[VX], *v = var[VY], *w = var[VZ];
		int put_X = para->geom->tile_putX, put_Y =para->geom->tile_putY, put_Z = para->geom->tile_putZ;

		// Set the value of para->bc->QPort =0
		for (i = 0; i < para->bc->nb_port; i++) {
				QPort[i] = 0.0;
		}

		// Loop all the boundary cells and calculate the flow rates at tiles
		for (it = 0; it < index; it++) {
				i = BINDEX[0][it];
				j = BINDEX[1][it];
				k = BINDEX[2][it];
				id = BINDEX[4][it];

				if (flagp[IX(i, j, k)] == TILE || flagp[IX(i, j, k)] == OUTLET) {
					// West or East Boundary
					if (put_X) {
							A = area_yz(para, var, i, j, k);

							if (i > 0)
									V_tmp = u[IX(i - 1, j, k)];
							else
									V_tmp = u[IX(i, j, k)];

							QPort[id] += V_tmp*A;
					}
					// South and North Boundary
					if (put_Y) {
							A = area_zx(para, var, i, j, k);

							if (j > 0)
									V_tmp = v[IX(i, j - 1, k)];
							else
									V_tmp = v[IX(i, j, k)];

							QPort[id] += V_tmp*A;
					}
					// Ceiling and Floor Boundary
					if (put_Z) {
							A = area_xy(para, var, i, j, k);

							if (k >0 )
									V_tmp = w[IX(i, j, k - 1)];
							else
									V_tmp = w[IX(i, j, k)];

							QPort[id] += V_tmp*A;
					}
				}
			else if (flagp[IX(i, j, k)] == INLET) {
				Ayz = area_yz(para, var, i, j, k);
				Azx = area_zx(para, var, i, j, k);
				Axy = area_xy(para, var, i, j, k);

				if (i > 0)
						V_tmp = u[IX(i - 1, j, k)];
				else
						V_tmp = u[IX(i, j, k)];

				QPort[id] += V_tmp*Ayz;

				if (j > 0)
						V_tmp = v[IX(i, j - 1, k)];
				else
						V_tmp = v[IX(i, j, k)];

				QPort[id] += V_tmp*Azx;

				if (k >0 )
						V_tmp = w[IX(i, j, k - 1)];
				else
						V_tmp = w[IX(i, j, k)];

				QPort[id] += V_tmp*Axy;
			}
		} //end of for

		// Write the results into files
		if (para->mytime->step_current == 0) {
				// create a new .dat file
				if ((FILE_TILE_FLOW = fopen("tile_flowrates.dat", "w+")) == NULL) {
						fprintf(stderr, "Error:can not open error file!\n");
						exit(1);
				}
				fprintf(FILE_TILE_FLOW, "Time\t");
				for (i = 0; i < para->bc->nb_port; i++) {
						fprintf(FILE_TILE_FLOW, "%s\t", para->bc->portName[i]);
				}
				fprintf(FILE_TILE_FLOW, "\n");
				fprintf(FILE_TILE_FLOW, "%.4f\t", para->mytime->t);
				for (i = 0; i < para->bc->nb_port; i++) {
						fprintf(FILE_TILE_FLOW, "%f\t", QPort[i] * 2118.88);
				}
				fprintf(FILE_TILE_FLOW, "\n");
				fclose(FILE_TILE_FLOW);
		}
		else {
				if ((FILE_TILE_FLOW = fopen("tile_flowrates.dat", "a+")) == NULL) {
						fprintf(stderr, "Error:can not open error file!\n");
						exit(1);
				}
				fprintf(FILE_TILE_FLOW, "%.4f\t", para->mytime->t);
				for (i = 0; i < para->bc->nb_port; i++) {
						fprintf(FILE_TILE_FLOW, "%f\t", QPort[i] * 2118.88);
				}
				fprintf(FILE_TILE_FLOW, "\n");

				fclose(FILE_TILE_FLOW);
		}
		// Output to screen or txt files
		for (i = 0; i < para->bc->nb_port; i++) {
				printf("%s--->>>>>>>>>>%f CFM\t\n", para->bc->portName[i],QPort[i]* 2118.88);
		}
		return 0;
}

///////////////////////////////////////////////////////////////////////////////
/// Check flow rates at inlets when t=0
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to the boundary index
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
REAL initial_inflows(PARA_DATA *para, REAL **var, int **BINDEX) {
		int i, j, k;
		int imax = para->geom->imax, jmax = para->geom->jmax;
		int kmax = para->geom->kmax;
		int IMAX = imax + 2, IJMAX = (imax + 2)*(jmax + 2);
		int it, id;
		int index = para->geom->index;
		REAL *flagp = var[FLAGP];
		int nb_ports = para->bc->nb_port;
		REAL A = 0.0, V_tmp = 0.0;
		REAL *QPort = para->bc->QPort;
		REAL *u = var[VXBC], *v = var[VYBC], *w = var[VZBC];
		REAL inflow = 0.0;

		// Loop all the boundary cells and calculate the flow rates at tiles
		for (it = 0; it < index; it++) {
				i = BINDEX[0][it];
				j = BINDEX[1][it];
				k = BINDEX[2][it];
				id = BINDEX[4][it];
				if (flagp[IX(i, j, k)] == INLET) {
						// West or East Boundary
						if (i == 0 || i == imax + 1) {
								A = area_yz(para, var, i, j, k);

								if (i == 0)
										V_tmp = u[IX(i, j, k)];
								else
										V_tmp = -1*u[IX(i, j, k)];

								QPort[id] += V_tmp*A;
						}
						// South and North Boundary
						if (j == 0 || j == jmax + 1) {
								A = area_zx(para, var, i, j, k);

								if (j == 0)
										V_tmp = v[IX(i, j, k)];
								else
										V_tmp = -1*v[IX(i, j, k)];

								inflow += V_tmp*A;
						}
						// Ceiling and Floor Boundary
						if (k == 0 || k == kmax + 1) {
								A = area_xy(para, var, i, j, k);

								if (k == 0)
										V_tmp = w[IX(i, j, k)];
								else
										V_tmp = -1*w[IX(i, j, k)];

								inflow += V_tmp*A;
						}
				}
		} //end of for

		return inflow;
}

///////////////////////////////////////////////////////////////////////////////
 /// Check the volumetric inflow rate
 ///
 ///\param para Pointer to FFD parameters
 ///\param var Pointer to FFD simulation variables
 ///\param psi Pointer to the variable
 ///\param BINDEX Pointer to the boundary index
 ///
 ///\return 0 if no error occurred
 ///////////////////////////////////////////////////////////////////////////////
REAL vol_inflow(PARA_DATA *para, REAL **var, int **BINDEX) {
		int i, j, k;
		int it;
		int imax = para->geom->imax, jmax = para->geom->jmax;
		int kmax = para->geom->kmax;
		int index = para->geom->index;
		int IMAX = imax + 2, IJMAX = (imax + 2)*(jmax + 2);
		REAL *gx = var[GX], *gy = var[GY], *gz = var[GZ];
		REAL *u = var[VX], *v = var[VY], *w = var[VZ];
		REAL mass_in = 0;
		REAL *flagp = var[FLAGP];

		/*---------------------------------------------------------------------------
		| Compute the total inflow
		---------------------------------------------------------------------------*/
		for (it = 0; it<index; it++) {
				i = BINDEX[0][it];
				j = BINDEX[1][it];
				k = BINDEX[2][it];

				if (flagp[IX(i, j, k)] == 0) {
						if (i == 0) mass_in += u[IX(i, j, k)] * (gy[IX(i, j, k)]
								- gy[IX(i, j - 1, k)])* (gz[IX(i, j, k)] - gz[IX(i, j, k - 1)]);

						if (i == imax + 1) mass_in += (-u[IX(i, j, k)])*(gy[IX(i, j, k)]
								- gy[IX(i, j - 1, k)])* (gz[IX(i, j, k)] - gz[IX(i, j, k - 1)]);

						if (j == 0) mass_in += v[IX(i, j, k)] * (gx[IX(i, j, k)]
								- gx[IX(i - 1, j, k)])* (gz[IX(i, j, k)] - gz[IX(i, j, k - 1)]);

						if (j == jmax + 1) mass_in += (-v[IX(i, j, k)])*(gx[IX(i, j, k)]
								- gx[IX(i - 1, j, k)])* (gz[IX(i, j, k)] - gz[IX(i, j, k - 1)]);

						if (k == 0) mass_in += w[IX(i, j, k)] * (gx[IX(i, j, k)]
								- gx[IX(i - 1, j, k)])* (gy[IX(i, j, k)] - gy[IX(i, j - 1, k)]);

						if (k == kmax + 1) mass_in += (-w[IX(i, j, k)])*(gx[IX(i, j, k)]
								- gx[IX(i - 1, j, k)])* (gy[IX(i, j, k)] - gy[IX(i, j - 1, k)]);
				}
		}
		return mass_in;
} // End of vol_inflow()

///////////////////////////////////////////////////////////////////////////////
/// Calculate the XY area of control volume (i,j,k)
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param i I-index of the control volume
///\param j J-index of the control volume
///\param K K-index of the control volume
///
///\return Area of XY surface
///////////////////////////////////////////////////////////////////////////////
REAL area_xy(PARA_DATA *para, REAL **var, int i, int j, int k) {
  return length_x(para, var, i, j, k)
       * length_y(para, var, i, j, k);
} // End of area_xy()

///////////////////////////////////////////////////////////////////////////////
/// Calculate the YZ area of control volume (i,j,k)
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param i I-index of the control volume
///\param j J-index of the control volume
///\param K K-index of the control volume
///
///\return Area of YZ surface
///////////////////////////////////////////////////////////////////////////////
REAL area_yz(PARA_DATA *para, REAL **var, int i, int j, int k) {
  return length_y(para, var, i, j, k)
       * length_z(para, var, i, j, k);
} // End of area_yz();

///////////////////////////////////////////////////////////////////////////////
/// Calculate the ZX area of control volume (i,j,k)
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param i I-index of the control volume
///\param j J-index of the control volume
///\param K K-index of the control volume
///
///\return Area of ZX surface
///////////////////////////////////////////////////////////////////////////////
REAL area_zx(PARA_DATA *para, REAL **var, int i, int j, int k) {
  return length_z(para, var, i, j, k)
       * length_x(para, var, i, j, k);
} // End of area_zx()

///////////////////////////////////////////////////////////////////////////////
/// Calculate the X-length of control volume (i,j,k)
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param i I-index of the control volume
///\param j J-index of the control volume
///\param K K-index of the control volume
///
///\return Length in X-direction
///////////////////////////////////////////////////////////////////////////////
REAL length_x(PARA_DATA *para, REAL **var, int i, int j, int k) {
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);

  if(i==0)
    return 0;
  else
    return (REAL) fabs(var[GX][IX(i,j,k)]-var[GX][IX(i-1,j,k)]);
} // End of length_x()

///////////////////////////////////////////////////////////////////////////////
/// Calculate the Y-length of control volume (i,j,k)
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param i I-index of the control volume
///\param j J-index of the control volume
///\param K K-index of the control volume
///
///\return Length in Y-direction
///////////////////////////////////////////////////////////////////////////////
REAL length_y(PARA_DATA *para, REAL **var, int i, int j, int k) {
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);

  if(j==0)
    return 0;
  else
    return (REAL) fabs(var[GY][IX(i,j,k)]-var[GY][IX(i,j-1,k)]);
} // End of length_y()

///////////////////////////////////////////////////////////////////////////////
/// Calculate the Z-length of control volume (i,j,k)
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param i I-index of the control volume
///\param j J-index of the control volume
///\param K K-index of the control volume
///
///\return Length in Z-direction
///////////////////////////////////////////////////////////////////////////////
REAL length_z(PARA_DATA *para, REAL **var, int i, int j, int k) {
  int imax = para->geom->imax, jmax = para->geom->jmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);

  if(k==0)
    return 0;
  else
    return (REAL) fabs(var[GZ][IX(i,j,k)]-var[GZ][IX(i,j,k-1)]);
} // End of length_z()

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
int rack_fluid_area(PARA_DATA *para, REAL **var, int **BINDEX) {
	int i, j, k, it, id, obj_type;
	REAL *flagp = var[FLAGP];
	int index = para->geom->index, imax = para->geom->imax,
		jmax = para->geom->jmax, kmax = para->geom->kmax;
	int IMAX = imax + 2, IJMAX = (imax + 2)*(jmax + 2);
	REAL axy, ayz, azx;

	// Initialize the area
	for (i = 0; i<para->bc->nb_rack; i++) {
		para->bc->RackArea[i] = 0.0;
	}

	// Loop all the boundary cells
	for (it = 0; it<index; it++) {
		i = BINDEX[0][it];
		j = BINDEX[1][it];
		k = BINDEX[2][it];
		id = BINDEX[4][it];
		obj_type = BINDEX[5][it];

		// calculate the area
		axy = area_xy(para, var, i, j, k);
		ayz = area_yz(para, var, i, j, k);
		azx = area_zx(para, var, i, j, k);

		// If it is rack cell and it is a rack inlet boundary
		if (obj_type == RACK && flagp[IX(i, j, k)] == RACK_INLET) {
			if (para->bc->RackDir[id] == 1 || para->bc->RackDir[id] == -1) {
				para->bc->RackArea[id] += ayz;
			}
			else if (para->bc->RackDir[id] == 2 || para->bc->RackDir[id] == -2) {
				para->bc->RackArea[id] += azx;
			}
			else if (para->bc->RackDir[id] == 3 || para->bc->RackDir[id] == -3) {
				para->bc->RackArea[id] += axy;
			}
			else {
				ffd_log("rack_fluid_area(): fail to detect the flow direction of the rack", FFD_ERROR);
				return 1;
			}
		} //end of if (obj_type == RACK && flagp[IX(i,j,k)]==RACK_INLET)
	} //end of for(it=0; it<index; it++)
			// return 0
	return 0;
}

///////////////////////////////////////////////////////////////////////////////////////
/// Assign the the velocity for the tiles after determining the pressure correction
/// A bisec method is used to solve the non-linear equations.
/// For more insights of the solver, refer to http://cims.nyu.edu/~donev/Teaching/NMI-Fall2010/Lecture6.handout.pdf
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///\param BINDEX Pointer to boundary index
///\return 0 if no error occurred
///\ Wei Tian
///\ 09/05/2017
//////////////////////////////////////////////////////////////////////////////////////
int assign_tile_velocity(PARA_DATA *para, REAL **var, int **BINDEX) {
	REAL *u = var[VX], *v = var[VY], *w = var[VZ];
	int i, j, k;
	int imax = para->geom->imax, jmax = para->geom->jmax;
	int kmax = para->geom->kmax;
	int IMAX = imax + 2, IJMAX = (imax + 2)*(jmax + 2);
	int it;
	int index = para->geom->index;
	REAL *flagp = var[FLAGP];
	REAL total_resistance_initil = 0.000000001;
	REAL in_flowrate = 0.0;
	REAL P_ave = 0.0;
	REAL *p = var[IP]; // pressure
	REAL p_corr1 = 0, p_corr2 = 0, p_corr3 = 0; // the correction for the pressure that is applied to the whole fluid field
	REAL epsilon1 = 1.0, epsilon2 = 1.0, epsilon3 = 1.0;
	int NEXT = 1;
	REAL corrected_flow = 0.0;
	REAL axy=0.0, ayz=0.0, azx=0.0;
	REAL rho = para->prob->rho;
	// if it is the begin of the simulation, then do an estimation
	if (para->mytime->step_current == 0) {
			for (it = 0; it < index; it++) {
					i = BINDEX[0][it];
					j = BINDEX[1][it];
					k = BINDEX[2][it];
					if (flagp[IX(i, j, k)] == TILE) {
							total_resistance_initil += 1 / pow(var[TILE_RESI_BC][IX(i, j, k)], 0.5);
					}
			} //end of for
			// find the volumetric inflow
			// since it is before the start of calculation, calling inflow function is not working.
			in_flowrate = initial_inflows(para, var, BINDEX);
			//printf("the inflow rate is %f\n", in_flowrate);
			//getchar();
			//in_flowrate = 0.472;
			// caclulate the P_ave
			P_ave = pow(in_flowrate, 2) / pow(total_resistance_initil, 2);
			// Calculate the velocity for each tile
			for (it = 0; it < index; it++) {
					i = BINDEX[0][it];
					j = BINDEX[1][it];
					k = BINDEX[2][it];
					axy = area_xy(para, var, i, j, k);
					ayz = area_yz(para, var, i, j, k);
					azx = area_zx(para, var, i, j, k);
					if (flagp[IX(i, j, k)] == TILE) {
							if (i == imax + 1) {
									var[TILE_FLOW_BC][IX(i, j, k)] = pow(P_ave / var[TILE_RESI_BC][IX(i, j, k)], 0.5)/ayz;
							}
							if (i == 0) {
									var[TILE_FLOW_BC][IX(i, j, k)] = -1*pow(P_ave / var[TILE_RESI_BC][IX(i, j, k)], 0.5) / ayz;
							}
							if (j == jmax + 1) {
									var[TILE_FLOW_BC][IX(i, j, k)] = pow(P_ave / var[TILE_RESI_BC][IX(i, j, k)], 0.5) / azx;
							}
							if (j == 0) {
									var[TILE_FLOW_BC][IX(i, j, k)] = -1 * pow(P_ave / var[TILE_RESI_BC][IX(i, j, k)], 0.5) / azx;
							}
							if (k == kmax + 1) {
									var[TILE_FLOW_BC][IX(i, j, k)] = pow(P_ave / var[TILE_RESI_BC][IX(i, j, k)], 0.5) / axy;
							}
							if (k == 0) {
									var[TILE_FLOW_BC][IX(i, j, k)] = -1 * pow(P_ave / var[TILE_RESI_BC][IX(i, j, k)], 0.5) / axy;
							}
							//printf("revised speed is--->>> %f\n", var[TILE_FLOW_BC][IX(i, j, k)]);
					}
			} //end of for

	}
	// otherwise, correct the pressure to meet the mass flow rate
	else {
			// initial values for p_corr with two rough assumptions
			p_corr1 = 0;
			epsilon1 = pressure_correction(para, var, BINDEX, p_corr1);
			//printf("***********************************\n");
			//printf("epsilon1 is %f\n", epsilon1);
			//getchar();
			/*current pressure is too large to provide surplus outflows*/
			if (epsilon1 >  1e-6 )
					p_corr2 = -1000;
			/*current pressure is too small to provide enough outflows*/
			else if (epsilon1 < -1* 1e-6 )
					p_corr2 = +1000;
			/*current pressure is good and correction is not needed*/
			else {
					p_corr2 = 0;
					p_corr3 = 0;
					NEXT = 0;
			}


			if (NEXT){
					// evaluate the guessed corrected pressure good enough or not
					epsilon2 = pressure_correction(para, var, BINDEX, p_corr2);
					//printf("epsilon1 and epsilon2 is %f %f\n", epsilon1, epsilon2);
					//getchar();

					if (epsilon2*epsilon1 < 0.0)
							p_corr3 = 0.5*(p_corr1 + p_corr2);
					else {
							sprintf(msg, "assign_tile_velocity(): epsilon 2, epsilon1, p_corr2 are %f %f %f", epsilon2, epsilon1, p_corr2);
							ffd_log(msg, FFD_NORMAL);
							ffd_log("assign_tile_velocity(): the inital value is not right", FFD_ERROR);
							return 1;
					}
			}

			// find p_corr
			while (NEXT) {
					epsilon3 = pressure_correction(para, var, BINDEX, p_corr3);
					if (fabs(epsilon3) < 1e-5) {
							NEXT = 0;
							//printf("epsilon3 and p_corr is %f %f\n", epsilon3, p_corr3);
					}

					if (epsilon3*epsilon1 < 0.0) {
							p_corr2 = p_corr3;
							epsilon2 = epsilon3;
							p_corr3 = 0.5*(p_corr1 + p_corr3);
					}
					else {
							p_corr1 = p_corr3;
							epsilon1 = epsilon3;
							p_corr3 = 0.5*(p_corr2 + p_corr3);
					}

			}
			//printf("***********************************\n");
			//printf("epsilon3 is %f\n",epsilon3);
			//getchar();
			// correct all the pressures
			FOR_ALL_CELL
					p[IX(i, j, k)] += p_corr3;
			END_FOR

			// update the flow rates at the tiles
					for (it = 0; it < index; it++) {
							i = BINDEX[0][it];
							j = BINDEX[1][it];
							k = BINDEX[2][it];
							if (flagp[IX(i, j, k)] == TILE) {
									if (i==imax+1 || j==imax+1 || k==kmax+1) {
											if (p[IX(i, j, k)] >0)
													var[TILE_FLOW_BC][IX(i, j, k)] = pow((p[IX(i, j, k)]) / (var[TILE_RESI_BC][IX(i, j, k)] * rho), 0.5);
											else
													var[TILE_FLOW_BC][IX(i, j, k)] = -1*pow(fabs((p[IX(i, j, k)]) / (var[TILE_RESI_BC][IX(i, j, k)] * rho)), 0.5);
									}
									else {
											if (p[IX(i, j, k)] >0)
													var[TILE_FLOW_BC][IX(i, j, k)] = -1*pow((p[IX(i, j, k)]) / (var[TILE_RESI_BC][IX(i, j, k)] * rho), 0.5);
											else
													var[TILE_FLOW_BC][IX(i, j, k)] = pow(fabs((p[IX(i, j, k)]) / (var[TILE_RESI_BC][IX(i, j, k)] * rho)), 0.5);
									}
									//printf("***********************************\n");
									//printf("revised speed is %f, P is %f\n", var[TILE_RESI_BC][IX(i, j, k)], p[IX(i, j, k)]);
							}

					} //end of for
			//getchar();
	}// end of if (para->mytime->step_current == 0)

	// Output the tile flow rates information
	if (check_tile_flowrate(para, var, BINDEX) != 0) {
			ffd_log("assign_tile_velocity: can not output the flow rates at tiles", FFD_ERROR);
	}
	return 0;
}

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
REAL pressure_correction(PARA_DATA *para, REAL **var, int **BINDEX, REAL p_corr) {
	int i, j, k;
	int imax = para->geom->imax, jmax = para->geom->jmax;
	int kmax = para->geom->kmax;
	int IMAX = imax + 2, IJMAX = (imax + 2)*(jmax + 2);
	int it;
	int index = para->geom->index;
	REAL *flagp = var[FLAGP];
	REAL *p = var[IP]; // pressure
	REAL epsilon=1e6, corrected_flow=0.0, in_flowrate=0.0;
	REAL axy, ayz, azx;
	REAL rho = para->prob->rho;

	in_flowrate = vol_inflow(para, var, BINDEX);
	//printf("the inflow is %f\n", in_flowrate);
	for (it = 0; it < index; it++) {
			i = BINDEX[0][it];
			j = BINDEX[1][it];
			k = BINDEX[2][it];
			axy = area_xy(para, var, i, j, k);
			ayz = area_yz(para, var, i, j, k);
			azx = area_zx(para, var, i, j, k);
			if (flagp[IX(i, j, k)] == TILE) {
					//printf("the coefficient is %f\n", var[TILE_RESI_BC][IX(i, j, k)]);
					if (i==imax+1 || i==0){
							if ((p[IX(i, j, k)] + p_corr)>0)
									corrected_flow += pow((p[IX(i, j, k)] + p_corr) / (var[TILE_RESI_BC][IX(i, j, k)]*rho), 0.5)*ayz;
							else
									corrected_flow -= pow(fabs((p[IX(i, j, k)] + p_corr) / (var[TILE_RESI_BC][IX(i, j, k)] * rho)), 0.5)*ayz;
					}
					if (j == jmax + 1 || j == 0) {
							if ((p[IX(i, j, k)] + p_corr)>0)
									corrected_flow += pow((p[IX(i, j, k)] + p_corr) / (var[TILE_RESI_BC][IX(i, j, k)] * rho), 0.5)*azx;
							else
									corrected_flow -= pow(fabs((p[IX(i, j, k)] + p_corr) / (var[TILE_RESI_BC][IX(i, j, k)] * rho)), 0.5)*azx;
					}
					if (k == kmax + 1 || k == 0) {
							if ((p[IX(i, j, k)] + p_corr)>0)
									corrected_flow += pow((p[IX(i, j, k)] + p_corr) / (var[TILE_RESI_BC][IX(i, j, k)] * rho), 0.5)*axy;
							else
									corrected_flow -= pow(fabs((p[IX(i, j, k)] + p_corr) / (var[TILE_RESI_BC][IX(i, j, k)] * rho)), 0.5)*axy;
					}

			}
	} //end of for

	epsilon = corrected_flow - in_flowrate;
	//printf("difference is %f, p corr is %f\n", epsilon, p_corr);
	//getchar();
	return epsilon;
}

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
int rack_model_black_box(PARA_DATA *para, REAL **var, int **BINDEX) {
	int i, j, k, it, id, obj_type;
	int iin, jin, kin;
	REAL *flagp = var[FLAGP];
  int index= para->geom->index, imax = para->geom->imax,
      jmax = para->geom->jmax, kmax = para->geom->kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  REAL axy, ayz, azx;
  REAL mDot_Cp, Q_dot;

  // Loop all the boundary cells
	for(it=0; it<index; it++) {
    i = BINDEX[0][it];
    j = BINDEX[1][it];
    k = BINDEX[2][it];
    id = BINDEX[4][it];
    obj_type = BINDEX[5][it];

    // calculate the area
    axy = area_xy(para, var, i, j, k);
    ayz = area_yz(para, var, i, j, k);
    azx = area_zx(para, var, i, j, k);

    // If it is rack cell and it is a rack inlet boundary
    if (obj_type == RACK) {
    	// Assign velocity to the inlet of rack
    	if (flagp[IX(i,j,k)]==RACK_INLET){
    		if (para->bc->RackDir[id] == 1 || para->bc->RackDir[id] == -1) {
    			var[VXBC][IX(i,j,k)] = para->bc->RackFlowRate[id]/para->bc->RackArea[id]*sign(para->bc->RackDir[id])/*direction*/;
    			var[VYBC][IX(i,j,k)] = 0.0;
    			var[VZBC][IX(i,j,k)] = 0.0;
    			// Assign the adjacent fluid cell temperature to the inlet of rack
    			var[TEMPBC][IX(i,j,k)] = var[TEMP][IX(i-sign(para->bc->RackDir[id]),j,k)];
    		}
    		else if (para->bc->RackDir[id] == 2 || para->bc->RackDir[id] == -2) {
    			var[VYBC][IX(i,j,k)] = para->bc->RackFlowRate[id]/para->bc->RackArea[id]*sign(para->bc->RackDir[id])/*direction*/;
    			var[VXBC][IX(i,j,k)] = 0.0;
    			var[VZBC][IX(i,j,k)] = 0.0;
    			// Assign the adjacent fluid cell temperature to the inlet of rack
    			var[TEMPBC][IX(i,j,k)] = var[TEMP][IX(i,j-sign(para->bc->RackDir[id]),k)];
    		}
    		else if (para->bc->RackDir[id] == 3 || para->bc->RackDir[id] == -3) {
    			var[VZBC][IX(i,j,k)] = para->bc->RackFlowRate[id]/para->bc->RackArea[id]*sign(para->bc->RackDir[id])/*direction*/;
    			var[VXBC][IX(i,j,k)] = 0.0;
    			var[VYBC][IX(i,j,k)] = 0.0;
    			// Assign the adjacent fluid cell temperature to the inlet of rack
    			var[TEMPBC][IX(i,j,k)] = var[TEMP][IX(i,j,k-sign(para->bc->RackDir[id]))];
    		}
    		else {
      		ffd_log("rack_model_black_box(): fail to detect the flow direction of the rack", FFD_ERROR);
      		return 1;
    		}

    	}
    	// Assign velocity and temperature to outlet of rack
    	else if (flagp[IX(i,j,k)]==RACK_OUTLET) {
    		if (para->bc->RackDir[id] == 1 || para->bc->RackDir[id] == -1) {
    			var[VXBC][IX(i,j,k)] = para->bc->RackFlowRate[id]/para->bc->RackArea[id]*sign(para->bc->RackDir[id])/*direction*/;
    			var[VYBC][IX(i,j,k)] = 0.0;
    			var[VZBC][IX(i,j,k)] = 0.0;
    			// Calculate the temperature at the outlet of rack
    			if (k==0){// This is to eliminate the divide by zero scenario
    				ayz = area_yz(para, var, i, j, k+1);
    			}
					iin = i - sign(para->bc->RackDir[id])*para->bc->RackMap[id][0];
					jin = j - sign(para->bc->RackDir[id])*para->bc->RackMap[id][1];
					kin = k - sign(para->bc->RackDir[id])*para->bc->RackMap[id][2];
					Q_dot = para->bc->HeatDiss[id]*ayz/para->bc->RackArea[id]; // heat dissipation by area
					mDot_Cp = para->prob->rho*var[VXBC][IX(i,j,k)]*ayz*para->prob->Cp; // mass flow rate multiply Cp
					var[TEMPBC][IX(i,j,k)] = var[TEMPBC][IX(iin,jin,kin)] + sign(para->bc->RackDir[id])*Q_dot/mDot_Cp;
					printf("rack_model_black_box(): temperature at outlet of rack [%d, %d, %d]: %f\n",i,j,k,var[TEMPBC][IX(i,j,k)]);
					printf("rack_model_black_box(): velocity at outlet of rack [%d, %d, %d]: %f\n",i,j,k,var[VXBC][IX(i,j,k)]);
    		}
    		else if (para->bc->RackDir[id] == 2 || para->bc->RackDir[id] == -2) {
    			var[VYBC][IX(i,j,k)] = para->bc->RackFlowRate[id]/para->bc->RackArea[id]*sign(para->bc->RackDir[id])/*direction*/;
    			var[VXBC][IX(i,j,k)] = 0.0;
    			var[VZBC][IX(i,j,k)] = 0.0;
    			// Calculate the temperature at the outlet of rack
    			iin = i - sign(para->bc->RackDir[id])*para->bc->RackMap[id][0];
					jin = j - sign(para->bc->RackDir[id])*para->bc->RackMap[id][1];
					kin = k - sign(para->bc->RackDir[id])*para->bc->RackMap[id][2];
					Q_dot = para->bc->HeatDiss[id]*azx/para->bc->RackArea[id]; // heat dissipation by area
					mDot_Cp = para->prob->rho*var[VYBC][IX(i,j,k)]*azx*para->prob->Cp; // mass flow rate multiply Cp
    			var[TEMPBC][IX(i,j,k)] = var[TEMPBC][IX(iin,jin,kin)] + Q_dot/mDot_Cp;
    		}
    		else if (para->bc->RackDir[id] == 3 || para->bc->RackDir[id] == -3) {
    			var[VZBC][IX(i,j,k)] = para->bc->RackFlowRate[id]/para->bc->RackArea[id]*sign(para->bc->RackDir[id])/*direction*/;
    			var[VXBC][IX(i,j,k)] = 0.0;
    			var[VYBC][IX(i,j,k)] = 0.0;
    			// Calculate the temperature at the outlet of rack
    			iin = i - sign(para->bc->RackDir[id])*para->bc->RackMap[id][0];
					jin = j - sign(para->bc->RackDir[id])*para->bc->RackMap[id][1];
					kin = k - sign(para->bc->RackDir[id])*para->bc->RackMap[id][2];
					Q_dot = para->bc->HeatDiss[id]*axy/para->bc->RackArea[id]; // heat dissipation by area
					mDot_Cp = para->prob->rho*var[VZBC][IX(i,j,k)]*axy*para->prob->Cp; // mass flow rate multiply Cp
    			var[TEMPBC][IX(i,j,k)] = var[TEMPBC][IX(iin,jin,kin)] + Q_dot/mDot_Cp;
    		}
    		else {
      		ffd_log("rack_model_black_box(): fail to detect the flow direction of the rack", FFD_ERROR);
      		return 1;
    		}
    	} //end of else if (flagp[IX(i,j,k)]==RACK_OUTLET)
    	// Pass internal rack cells
    	else {
    		continue;
    	} //end of else
    } //end of if (obj_type == RACK)
	} //end of for(it=0; it<index; it++)

	// return 0
	return 0;
}
