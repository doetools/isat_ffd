///////////////////////////////////////////////////////////////////////////////
/// \Overview of the functions used in the kernels (shown in the next section)
///
/// \Disclaimer:the codes are of absolute NO warranties in any form. Use at own risks.
///
/// \author: Wei Tian, w.tian@umaimi.edu
/// 
/// \date:   08/08/2017
///
/// \The author acknowledges Tom Sevilla for providing assistance in debugging the codes.
///
/// \Recent updates:
///
/// \add chen's zero equation model and convective heat transfer coefficient
///
///
/// \ALL RIGHTS RESERVED @2017
///////////////////////////////////////////////////////////////////////////////
/// func. 0: REAL length_x(__global PARA_DATA_SIMP *para, __global REAL *var, int *ip, int *jp, int *kp)
/// func. 1: REAL length_y(__global PARA_DATA_SIMP *para, __global REAL *var, int *ip, int *jp, int *kp)
/// func. 2: REAL length_z(__global PARA_DATA_SIMP *para, __global REAL *var, int *ip, int *jp, int *kp)
/// func. 3: REAL area_xy(__global PARA_DATA_SIMP *para, __global REAL *var, int *ip, int *jp, int *kp)
/// func. 4: REAL area_zx(__global PARA_DATA_SIMP *para, __global REAL *var, int *ip, int *jp, int *kp)
/// func. 5: REAL area_yz(__global PARA_DATA_SIMP *para, __global REAL *var, int *ip, int *jp, int *kp)
/// func. 6: REAL interpolation(__global PARA_DATA_SIMP *para, __global REAL *d0, REAL *x_1p, REAL *y_1p, REAL *z_1p,
///                   int *pp, int *qp, int *rp)
/// func. 7: void set_x_location(__global PARA_DATA_SIMP *para, __global REAL *var, __global REAL *flag, __global REAL *x, REAL *u0p, 
///                    int *ip, int *jp, int *kp, 
///                    REAL *OL, int *OC, int *LOC, int *COOD)
/// func. 8: void set_y_location(__global PARA_DATA_SIMP *para, __global REAL *var, __global REAL *flag, __global REAL *y, REAL *v0p, 
///                    int *ip, int *jp, int *kp,
///                    REAL *OL, int *OC, int *LOC, int *COOD)
/// func. 9: void set_z_location(__global PARA_DATA_SIMP *para, __global REAL *var, __global REAL *flag, __global REAL *z, REAL *w0p,
///                    int *ip, int *jp, int *kp,
///                    REAL *OL, int *OC, int *LOC, int *COOD)
/// func. 10: REAL nu_t_chen_zero_equ(__global PARA_DATA_SIMP *para, __global REAL *var, int *ip, int *jp, int *kp)
/// func. 11: REAL h_coef (__global PARA_DATA_SIMP *para, __global REAL *var, int *ip, int *jp, int *kp, int *dp)
/// func. 12: 
///////////////////////////////////////////////////////////////////////////////

//#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#include"data_structure_gpu.h"

///////////////////////////////////////////////////////////////////////////////
/// check lenth X
/// 6/7/2015
/// Wei Tian
///////////////////////////////////////////////////////////////////////////////
REAL length_x(__global PARA_DATA_SIMP *para, __global REAL *var, int *ip, int *jp, int *kp) {
  int i = ip[0],j = jp[0],k = kp[0];
  int imax = para->geom.imax, jmax = para->geom.jmax, kmax = para->geom.kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int size = (imax + 2)*(jmax + 2)*(kmax + 2);
  if(i==0)
    return 0;
  else
    return (REAL) fabs(var[GX*size+IX(i,j,k)]-var[GX*size+IX(i-1,j,k)]);
}

///////////////////////////////////////////////////////////////////////////////
/// check lenth Y
/// 6/7/2015
/// Wei Tian
///////////////////////////////////////////////////////////////////////////////
REAL length_y(__global PARA_DATA_SIMP *para, __global REAL *var, int *ip, int *jp, int *kp) {
  int i = ip[0],j = jp[0],k = kp[0];
  int imax = para->geom.imax, jmax = para->geom.jmax, kmax = para->geom.kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int size = (imax + 2)*(jmax + 2)*(kmax + 2);

  if(j==0)
    return 0;
  else
    return (REAL) fabs(var[GY*size+IX(i,j,k)]-var[GY*size+IX(i,j-1,k)]);
} // End of length_y()

///////////////////////////////////////////////////////////////////////////////
/// check lenth Z
/// 6/7/2015
/// Wei Tian
///////////////////////////////////////////////////////////////////////////////
REAL length_z(__global PARA_DATA_SIMP *para, __global REAL *var, int *ip, int *jp, int *kp) {
  int i = ip[0],j = jp[0],k = kp[0];
  int imax = para->geom.imax, jmax = para->geom.jmax, kmax = para->geom.kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int size = (imax + 2)*(jmax + 2)*(kmax + 2);

  if(k==0)
    return 0;
  else
    return (REAL) fabs(var[GZ*size+IX(i,j,k)]-var[GZ*size+IX(i,j,k-1)]);
} // End of length_z()

///////////////////////////////////////////////////////////////////////////////
/// check area_xy area
/// 6/7/2015
/// Wei Tian
///////////////////////////////////////////////////////////////////////////////
REAL area_xy(__global PARA_DATA_SIMP *para, __global REAL *var, int *ip, int *jp, int *kp) {
  return length_x(para, var, ip, jp, kp)
       * length_y(para, var, ip, jp, kp);
}

///////////////////////////////////////////////////////////////////////////////
/// check area_zx area
/// 6/7/2015
/// Wei Tian
///////////////////////////////////////////////////////////////////////////////

REAL area_zx(__global PARA_DATA_SIMP *para, __global REAL *var, int *ip, int *jp, int *kp) {
  return length_z(para, var, ip, jp, kp)
       * length_x(para, var, ip, jp, kp);
} // End of area_zx()

///////////////////////////////////////////////////////////////////////////////
/// check area_yz area
/// 6/7/2015
/// Wei Tian
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
REAL area_yz(__global PARA_DATA_SIMP *para, __global REAL *var, int *ip, int *jp, int *kp) {
  return length_y(para, var, ip, jp, kp)
       * length_z(para, var, ip, jp, kp);
} // End of area_yz();

///////////////////////////////////////////////////////////////////////////////
/// Bilinear interpolation for advection
/// 6/17/2015
/// Wei Tian
///////////////////////////////////////////////////////////////////////////////
REAL interpolation(__global PARA_DATA_SIMP *para, __global REAL *d0, REAL *x_1p, REAL *y_1p, REAL *z_1p,
                   int *pp, int *qp, int *rp) {
  int imax = para->geom.imax, jmax = para->geom.jmax, kmax = para->geom.kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int size = (imax+2)*(jmax+2)*(kmax+2);
  int p = pp[0], q = qp[0],r = rp[0];
  REAL x_1 = x_1p[0],y_1 = y_1p[0],z_1 = z_1p[0];
  REAL x_0, y_0, z_0;
  REAL tmp0, tmp1;
  REAL d000, d010, d100, d110;
  REAL d001, d011, d101, d111;
  
  // assign the coefficients
  d000 = d0[IX(p,q,r)];
  d010 = d0[IX(p,q+1,r)];
  d100 = d0[IX(p+1,q,r)];
  d110 = d0[IX(p+1,q+1,r)];
  d001 = d0[IX(p,q,r+1)];
  d011 = d0[IX(p,q+1,r+1)];
  d101 = d0[IX(p+1,q,r+1)];
  d111 = d0[IX(p+1,q+1,r+1)];
  
  // interpolation
  x_0 = (REAL) 1.0 - x_1;
  y_0 = (REAL) 1.0 - y_1;
  z_0 = (REAL) 1.0 - z_1;

  tmp0 = x_0*(y_0*d000+y_1*d010) + x_1*(y_0*d100+y_1*d110);
  tmp1 = x_0*(y_0*d001+y_1*d011) + x_1*(y_0*d101+y_1*d111);

  return z_0*tmp0+z_1*tmp1;
 
} // End of interpolation()


///////////////////////////////////////////////////////////////////////////////
/// sub_function for advection of vectors: Find x location
/// 6/16/2015
// Wei Tian
///////////////////////////////////////////////////////////////////////////////
void set_x_location(__global PARA_DATA_SIMP *para, __global REAL *var, __global REAL *flag, __global REAL *x, REAL *u0p, 
                    int *ip, int *jp, int *kp, 
                    REAL *OL, int *OC, int *LOC, int *COOD) {
    int i = ip[0],j = jp[0],k = kp[0];
  REAL u0 = u0p[0];
  int imax = para->geom.imax, jmax = para->geom.jmax;
  int kmax = para->geom.kmax;
  int size = (imax + 2)*(jmax + 2)*(kmax + 2);
  int IMAX = imax + 2, IJMAX = (imax + 2)*(jmax + 2);

/****************************************************************************
  | If the previous location is equal to current position
  | stop the process (COOD[X] = 0) 
  ****************************************************************************/
  if(OL[X]==x[IX(OC[X],OC[Y],OC[Z])]) 
    COOD[X]=0;  
 /****************************************************************************
  | Otherwise, if previous location is on the west of the current position
  ****************************************************************************/
  else if(OL[X]<x[IX(OC[X],OC[Y],OC[Z])]) {
    // If donot reach the boundary yet, move to west 
    if(OC[X]>0) 
      OC[X] -=1;

    // If the previous position is on the east of new location, stop the process
    if(OL[X]>=x[IX(OC[X],OC[Y],OC[Z])]) 

      COOD[X]=0; 

    // If the new position is solid 
    if(flag[IX(OC[X],OC[Y],OC[Z])]==1) {
      // Use the east cell for new location
      OL[X] = x[IX(OC[X]+1,OC[Y],OC[Z])]; 
      OC[X] +=1;
      // Hit the boundary
      LOC[X] = 0; 
      // Stop the trace process
      COOD[X] = 0; 
    } // End of if() for solid

    // If the new position is inlet or outlet
    if(flag[IX(OC[X],OC[Y],OC[Z])]==0||flag[IX(OC[X],OC[Y],OC[Z])]==2) {
      // Use new position
      OL[X] = x[IX(OC[X],OC[Y],OC[Z])];
      // use east cell for coordinate
      OC[X] += 1;
      // Hit the boundary
      LOC[X] = 0; 
      // Stop the trace process
      COOD[X]=0; 
    } // End of if() for inlet or outlet
  } // End of if() for previous position is on the west of new position 
  /****************************************************************************
  | Otherwise, if previous location is on the east of the current position
  ****************************************************************************/
  else {
    // If not at the east boundary
    if(OC[X]<=imax) 
      // Move to east
      OC[X] +=1;

    // If the previous position is  on the west of new position
    if(OL[X]<=x[IX(OC[X],OC[Y],OC[Z])]) 
      // Stop the trace process
      COOD[X]=0;

    // If the cell is solid
    if(flag[IX(OC[X],OC[Y],OC[Z])]==1) {
      // Use west cell
      OL[X] = x[IX(OC[X]-1,OC[Y],OC[Z])]; 
      OC[X] -= 1;
      // Hit the boundary
      LOC[X] = 0;
      // Stop the trace process
      COOD[X] = 0;
    } // End of if() for solid

    // If the new position is inlet or outlet 
    if(flag[IX(OC[X],OC[Y],OC[Z])]==0||flag[IX(OC[X],OC[Y],OC[Z])]==2) {
      // Use the current cell for previous location
      OL[X] = x[IX(OC[X],OC[Y],OC[Z])];
      // Use the west cell for coordinate
      OC[X] -=1; 
      // Hit the boundary
      LOC[X] = 0;
      // Stop the trace process
      COOD[X]=0;
    } // End of if() for inlet or outlet
  } // End of if() for previous position is on the east of new position
} // End of set_x_location()

///////////////////////////////////////////////////////////////////////////////
/// sub_function for advection of vectors: Find y location
/// 6/16/2015
// Wei Tian
///////////////////////////////////////////////////////////////////////////////
void set_y_location(__global PARA_DATA_SIMP *para, __global REAL *var, __global REAL *flag, __global REAL *y, REAL *v0p, 
                    int *ip, int *jp, int *kp,
                    REAL *OL, int *OC, int *LOC, int *COOD) {
    int i = ip[0],j = jp[0],k = kp[0];
  REAL v0 = v0p[0];
  int imax = para->geom.imax, jmax = para->geom.jmax;
  int kmax = para->geom.kmax;
  int size = (imax + 2)*(jmax + 2)*(kmax + 2);
  int IMAX = imax + 2, IJMAX = (imax + 2)*(jmax + 2);

  /****************************************************************************
  | If the previous location is equal to current position,
  | stop the process (COOD[X] = 0) 
  ****************************************************************************/
  if(OL[Y]==y[IX(OC[X],OC[Y],OC[Z])]) 
    COOD[Y] = 0;
  /****************************************************************************
  | Otherwise, if previous location is on the south of the current position
  ****************************************************************************/
  else if(OL[Y]<y[IX(OC[X],OC[Y],OC[Z])]) {
    // If donot reach the boundary yet
    if(OC[Y]>0) 
      OC[Y] -= 1;

    // If the previous position is on the north of new location 
    if(OL[Y]>=y[IX(OC[X],OC[Y],OC[Z])]) 
      // Stop the process 
      COOD[Y] = 0;

    // If the new position is solid 
    if(flag[IX(OC[X],OC[Y],OC[Z])]==1) {
      // Use the north cell for new location
      OL[Y] = y[IX(OC[X],OC[Y]+1,OC[Z])]; 
      OC[Y] += 1; 
      // Hit the boundary
      LOC[Y] = 0;
      // Stop the trace process
      COOD[Y] = 0;
    } // End of if() for solid

    // If the new position is inlet or outlet
    if(flag[IX(OC[X],OC[Y],OC[Z])]==0||flag[IX(OC[X],OC[Y],OC[Z])]==2) {
      // Use new position
      OL[Y] = y[IX(OC[X],OC[Y],OC[Z])];
      // Use north cell for coordinate
      OC[Y] += 1; 
      // Hit the boundary
      LOC[Y] = 0;
      // Stop the trace process
      COOD[Y] = 0;
    } // End of if() for inlet or outlet
  } // End of if() for previous position is on the south of new position 
  /****************************************************************************
  | Otherwise, if previous location is on the north of the current position
  ****************************************************************************/
  else {
    // If not at the north boundary
    if(OC[Y]<=jmax) 
      // Move to north
      OC[Y] +=1;

    // If the previous position is on the south of new position
    if(OL[Y]<=y[IX(OC[X],OC[Y],OC[Z])]) 
      // Stop the trace process
      COOD[Y] = 0;

    // If the cell is solid
    if(flag[IX(OC[X],OC[Y],OC[Z])]==1) {
      // Use south cell
      OL[Y] = y[IX(OC[X],OC[Y]-1,OC[Z])]; 
      OC[Y] -= 1;
      // Hit the boundary
      LOC[Y] = 0;
      // Stop the trace process
      COOD[Y] = 0;
    } // End of if() for solid

    // If the new position is inlet or outlet 
    if(flag[IX(OC[X],OC[Y],OC[Z])]==0||flag[IX(OC[X],OC[Y],OC[Z])]==2) {
      // Use the current cell for previous location
      OL[Y] = y[IX(OC[X],OC[Y],OC[Z])];
      // Use the south cell for coordinate
      OC[Y] -= 1;
      // Hit the boundary
      LOC[Y]=0;
      // Stop the trace process
      COOD[Y]=0;
    } // End of if() for inlet or outlet
  } // End of if() for previous position is on the east of new position
} // End of set_y_location()

///////////////////////////////////////////////////////////////////////////////
/// sub_function for advection of vectors: Find z location
/// 6/16/2015
// Wei Tian
///////////////////////////////////////////////////////////////////////////////
void set_z_location(__global PARA_DATA_SIMP *para, __global REAL *var, __global REAL *flag, __global REAL *z, REAL *w0p,
                    int *ip, int *jp, int *kp,
                    REAL *OL, int *OC, int *LOC, int *COOD) {
    int i = ip[0],j = jp[0],k = kp[0];
  REAL w0 = w0p[0];
  int imax = para->geom.imax, jmax = para->geom.jmax;
  int kmax = para->geom.kmax;
  int size = (imax + 2)*(jmax + 2)*(kmax + 2);
  int IMAX = imax + 2, IJMAX = (imax + 2)*(jmax + 2);

  /****************************************************************************
  | If the previous location is equal to current position,
  | stop the process (COOD[Z] = 0) 
  ****************************************************************************/
  if(OL[Z]==z[IX(OC[X],OC[Y],OC[Z])]) 
    COOD[Z] = 0;
  /****************************************************************************
  | Otherwise, if previous location is on the floor of the current position
  ****************************************************************************/
  else if(OL[Z]<z[IX(OC[X],OC[Y],OC[Z])]) {
    // If donot reach the boundary yet
    if(OC[Z]>0) 
      OC[Z] -= 1;
  
    // If the previous position is on the ceiling of new location 
    if(OL[Z]>=z[IX(OC[X],OC[Y],OC[Z])]) 
      // Stop the process 
      COOD[Z] = 0;
    
    // If the new position is solid 
    if(flag[IX(OC[X],OC[Y],OC[Z])]==1) {
      // Use the ceiling cell for new location
      OL[Z] = z[IX(OC[X],OC[Y],OC[Z]+1)];
      OC[Z] += 1;
      // Hit the boundary
      LOC[Z] = 0;
      // Stop the trace process
      COOD[Z] = 0;
    } // End of if() for solid

    // If the new position is inlet or outlet
    if(flag[IX(OC[X],OC[Y],OC[Z])]==0||flag[IX(OC[X],OC[Y],OC[Z])]==2) {
      // Use new position
      OL[Z] = z[IX(OC[X],OC[Y],OC[Z])];
      // Use ceiling cell for coordinate
      OC[Z] += 1;
      // Hit the boundary
      LOC[Z] = 0;
      // Stop the trace process
      COOD[Z] = 0;
    } // End of if() for inlet or outlet
  } // End of if() for previous position is on the floor of new position 
  /****************************************************************************
  | Otherwise, if previous location is on the ceiling of the current position
  -***************************************************************************/
  else {
    // If not at the ceiling boundary
    if(OC[Z]<=kmax) 
      // Move to ceiling
      OC[Z] += 1;

    // If the previous position is on the floor of new position
    if(OL[Z] <=z[IX(OC[X],OC[Y],OC[Z])]) 
      // Stop the trace process
      COOD[Z] = 0;

    // If the cell is solid
    if(flag[IX(OC[X],OC[Y],OC[Z])]==1) {
      // Use floor cell
      OL[Z] = z[IX(OC[X],OC[Y],OC[Z]-1)]; 
      OC[Z] -= 1;
      // Hit the boundary
      LOC[Z] = 0;
      // Stop the trace process
      COOD[Z] = 0;
    } // End of if() for solid

    // If the new position is inlet or outlet 
    if(flag[IX(OC[X],OC[Y],OC[Z])]==0||flag[IX(OC[X],OC[Y],OC[Z])]==2) {
      // Use the current cell for previous location
      OL[Z]=z[IX(OC[X],OC[Y],OC[Z])];
      // Use the floor cell for coordinate
      OC[Z] -= 1;
      // Hit the boundary
      LOC[Z]=0;
      // Stop the trace process
      COOD[Z]=0;
    } // End of if() for inlet or outlet
  } // End of if() for previous position is on the east of new position
} // End of set_z_location()

///////////////////////////////////////////////////////////////////////////////
/// functions to calculate the turbulent viscosity/thermal diffusivity using Chen's zero equation 
/// 8/08/2017
/// Wei Tian
///////////////////////////////////////////////////////////////////////////////
REAL nu_t_chen_zero_equ(__global PARA_DATA_SIMP *para, __global REAL *var, int *ip, int *jp, int *kp){

  int i = ip[0],j = jp[0],k = kp[0];
  int imax = para->geom.imax, jmax = para->geom.jmax;
  int kmax = para->geom.kmax;
  int size = (imax + 2)*(jmax + 2)*(kmax + 2);
  int IMAX = imax + 2, IJMAX = (imax + 2)*(jmax + 2);
	REAL l= 0.0;
	REAL jim_a = 0.0185;
	REAL coeff = para->prob.chen_a;
	REAL nu_t;

	__global float *u = &var[VX*size];
	__global float *v = &var[VY*size];
	__global float *w = &var[VZ*size];
	__global float *flagp = &var[FLAGP*size];

	l = var[MIN_DISTANCE*size + IX(i, j, k)];

  if (flagp[IX(i - 1, j, k)] >= 0 || flagp[IX(i + 1, j, k)] >= 0 ||
				flagp[IX(i, j - 1, k)] >= 0 || flagp[IX(i, j + 1, k)] >= 0 ||
				flagp[IX(i, j, k - 1)] >= 0 || flagp[IX(i, j, k + 1)] >= 0) {
				//if the cell is adjacent to solid boundaries, assign a otherwise coffecient
				coeff = jim_a;
	}
  else {
				// if the cell is not adjacent to solid boundaries, assign a standard coffecient
				coeff = para->prob.chen_a;
		}

  nu_t = coeff * l
				* (REAL)sqrt(u[IX(i,j,k)]*u[IX(i,j,k)]
										+v[IX(i,j,k)]*v[IX(i,j,k)]
										+w[IX(i,j,k)]*w[IX(i,j,k)] );

  return nu_t;
}


///////////////////////////////////////////////////////////////////////////////
/// functions to calculate convective heat transfer coefficient
/// 8/14/2017
/// Wei Tian
///////////////////////////////////////////////////////////////////////////////
REAL h_coef (__global PARA_DATA_SIMP *para, __global REAL *var, int *ip, int *jp, int *kp, int *dp) {
	REAL D = dp[0];
	REAL nu = para->prob.nu;
  REAL rhoCp = para->prob.rho * para->prob.Cp;
  REAL coef_h=para->prob.coeff_h;
	REAL h, kapa;

  if (para->prob.tur_model == LAM) {
		return coef_h*rhoCp;
	}
	else if (para->prob.tur_model == CONSTANT) {
		return coef_h*rhoCp*101;
	}
	else if (para->prob.tur_model == CHEN) {
		kapa = nu + nu_t_chen_zero_equ(para, var, ip, jp, kp);
		h = (rhoCp* para->prob.alpha / D) *(kapa / nu);
		return h;
	}
	else
		return coef_h*rhoCp;

}// end of h_coef























///////////////////////////////////////////////////////////////////////////////
/// \brief the functions of kernels. 
///
/// \Disclaimer:the codes are of absolute NO warranties in any form. Use at own risks.
///
/// \define all the kernels that run on the multi-core devices, like CPUs, GPUs
///
/// \author: Wei Tian, w.tian@umiami.edu
/// \date:   08/08/2017
///
/// \ The author acknowledges Tom Sevilla for providing assistance in debugging the codes.
///
/// \Recent updates:
///
/// \add implicit scheme to the advection terms
///
/// \advection may not be vectorized perfectly. Need improvement.
///
/// \use existing parallel libraries for solving the quation. The libray should be
///  independent of the platforms, or, at least should be compatible to NVIDIA platform.
///
/// \ ALL RIGHTS RESERVED @2017
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/// 0: void adve_VX(__global PARA_DATA_SIMP *para, __global REAL *var, __global  int *BINDEX)
/// 1: void adve_VY(__global PARA_DATA_SIMP *para, __global REAL *var, __global  int *BINDEX)
/// 2: void adve_VZ(__global PARA_DATA_SIMP *para, __global REAL *var, __global  int *BINDEX)
/// 3: void adve_T(__global PARA_DATA_SIMP *para, __global REAL *var, __global  int *BINDEX)
/// 4: void diff_VX(__global PARA_DATA_SIMP *para, __global REAL *var,  __global  int *BINDEX)
/// 5: void diff_VY(__global PARA_DATA_SIMP *para, __global REAL *var, __global  int *BINDEX)
/// 6: void diff_VZ(__global PARA_DATA_SIMP *para, __global REAL *var,  __global  int *BINDEX)
/// 7: void diff_T(__global PARA_DATA_SIMP *para, __global REAL *var,  __global  int *BINDEX) 
/// 8: void project(__global PARA_DATA_SIMP *para, __global REAL *var, __global int *BINDEX)
/// 9: void project_velo_corr(__global PARA_DATA_SIMP *para, __global REAL *var, __global int *BINDEX)
/// 10: void Ax_BSolver(__global PARA_DATA_SIMP *para, __global REAL *var,__global int *index, __global int *WF)
/// 11: void Ax_BSolver_P(__global PARA_DATA_SIMP *para, __global REAL *var)
/// 12: void ap_coeff(__global PARA_DATA_SIMP *para, __global REAL *var, __global int *var_type)
/// 13: set_bnd_T(__global PARA_DATA_SIMP *para, __global REAL *var, __global REAL *psi, __global int *BINDEX)
/// 14: set_bnd_pressure(__global PARA_DATA_SIMP *para, __global REAL *var,__global int *BINDEX)
/// 15: set_bnd_VX(__global PARA_DATA_SIMP *para, __global REAL *var, __global int *WF, __global int *BINDEX)
/// 16: set_bnd_VY(__global PARA_DATA_SIMP *para, __global REAL *var, __global int *WF, __global int *BINDEX)
/// 17: set_bnd_VZ(__global PARA_DATA_SIMP *para, __global REAL *var, __global int *WF, __global int *BINDEX)
/// 18: adjust_velocity(PARA_DATA *para, REAL **var, int **BINDEX) 
/// 19: mass_conservation(__global PARA_DATA_SIMP *para, __global REAL *var, __global int *BINDEX) 
/// 20: void reset_time_averaged_data(__global PARA_DATA_SIMP *para, __global REAL *var)
/// 21: void add_time_averaged_data(__global PARA_DATA_SIMP *para, __global REAL *var)
/// 22: void time_averaged(__global PARA_DATA_SIMP *para, __global REAL *var)
/// 23: void Ax_BSolver_upd(__global PARA_DATA_SIMP *para, __global REAL *var,__global int *index, __global int *WF)
/// 24: void Ax_BSolver_P_upd(__global PARA_DATA_SIMP *para, __global REAL *var)
/// 25: void adve_VX_im(__global PARA_DATA_SIMP *para, __global REAL *var, __global  int *BINDEX)
/// 26: void adve_VY_im(__global PARA_DATA_SIMP *para, __global REAL *var, __global  int *BINDEX)
/// 27: void adve_VZ_im(__global PARA_DATA_SIMP *para, __global REAL *var, __global  int *BINDEX)
/// 28: void adve_T_im(__global PARA_DATA_SIMP *para, __global REAL *var, __global  int *BINDEX)
/// 29: void ap_im_coeff(__global PARA_DATA_SIMP *para, __global REAL *var, __global int *var_type)
/// 30: set_bnd_T_im(__global PARA_DATA_SIMP *para, __global REAL *var, __global int *BINDEX)
/// 31: set_bnd_VX_im(__global PARA_DATA_SIMP *para, __global REAL *var, __global int *BINDEX)
/// 32: set_bnd_VY_im(__global PARA_DATA_SIMP *para, __global REAL *var, __global int *BINDEX)
/// 33: set_bnd_VZ_im(__global PARA_DATA_SIMP *para, __global REAL *var, __global int *BINDEX)
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
/// kenerl setting diffusion coefficient for VX
/// 6/7/2015
///////////////////////////////////////////////////////////////////////////////
__kernel void diff_VX(__global PARA_DATA_SIMP *para, __global REAL *var, __global  int *BINDEX) {

  int i = get_global_id(0);
  int j = get_global_id(1);
  int k = get_global_id(2);
  int bar =0;
  int imax = para->geom.imax, jmax = para->geom.jmax;
  int kmax = para->geom.kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int size=(imax+2)*(jmax+2)*(kmax+2);

  __global float *psi = &var[VX*size];
  __global float *psi0 = &var[TMP1*size];
  __global float *aw = &var[AW*size];
  __global float *ae = &var[AE*size];
  __global float *as = &var[AS*size];
  __global float *an = &var[AN*size];
  __global float *af = &var[AF*size];
  __global float *ab = &var[AB*size]; 
  __global float *ap = &var[AP*size];
  __global float *ap0 = &var[AP0*size]; 
  __global float *b = &var[B*size];
  __global float *gx = &var[GX*size];
  __global float *x = &var[X*size];
  __global float *y = &var[Y*size];
  __global float*z = &var[Z*size];
  __global float *gy = &var[GY*size];
  __global float *gz = &var[GZ*size];
  __global float *pp = &var[PP*size];
  __global float *Temp = &var[TEMP*size];
  __global float *flagp = &var[FLAGP*size];
   
  REAL dxe, dxw, dyn, dys, dzf, dzb, Dx, Dy, Dz;
  REAL dt = para->mytime.dt, beta = para->prob.beta;
  REAL Temp_Buoyancy = para->prob.Temp_Buoyancy;
  REAL gravx = para->prob.gravx, gravy = para->prob.gravy,
       gravz = para->prob.gravz;
  REAL kapa=para->prob.nu;
  REAL coef_CONSTANT = 120;

  // FOR_U_CELL
  if (i>0 && i<imax && j>0 && j<jmax+1 && k>0 && k<kmax+1) {

				if (para->prob.tur_model==CHEN) {
						kapa = nu_t_chen_zero_equ(para, var, &i, &j, &k)+ para->prob.nu;
				}


        dxe = gx[IX(i+1,j  ,k)] - gx[IX(i  ,j,k)];
        dxw = gx[IX(i  ,j  ,k)] - gx[IX(i-1,j,k)];
        dyn =  y[IX(i  ,j+1,k)] -  y[IX(i  ,j,k)];
        dys =  y[IX(i,j,k)] -      y[IX(i,j-1,k)];
        dzf =  z[IX(i,j,k+1)] -    z[IX(i,j,k)];
        dzb =  z[IX(i,j,k)] -      z[IX(i,j,k-1)];
        Dx =   x[IX(i+1,j,k)] -    x[IX(i,j,k)];
        Dy =  gy[IX(i,j,k)] -     gy[IX(i,j-1,k)];
        Dz =  gz[IX(i,j,k)] -     gz[IX(i,j,k-1)];

        aw[IX(i,j,k)] = kapa*Dy*Dz/dxw;
        ae[IX(i,j,k)] = kapa*Dy*Dz/dxe;
        an[IX(i,j,k)] = kapa*Dx*Dz/dyn;
        as[IX(i,j,k)] = kapa*Dx*Dz/dys;
        af[IX(i,j,k)] = kapa*Dx*Dy/dzf;
        ab[IX(i,j,k)] = kapa*Dx*Dy/dzb;
        ap0[IX(i,j,k)] = Dx*Dy*Dz/dt;
        b[IX(i,j,k)] = psi0[IX(i,j,k)]*ap0[IX(i,j,k)]
                     - beta*gravx*(Temp[IX(i,j,k)]-Temp_Buoyancy)*Dx*Dy*Dz
                     + (pp[IX(i,j,k)]-pp[IX(i+1,j,k)])*Dy*Dz;

  } //end of FOR_U_CELL
}// end of __kernel void diff_VX()

///////////////////////////////////////////////////////////////////////////////
/// kenerl setting advection coefficient for VX
/// 8/7/2015
///////////////////////////////////////////////////////////////////////////////
__kernel void adve_VX_im(__global PARA_DATA_SIMP *para, __global REAL *var, __global  int *BINDEX) {
  int i = get_global_id(0);
  int j = get_global_id(1);
  int k = get_global_id(2);
  int imax = para->geom.imax, jmax = para->geom.jmax;
  int kmax = para->geom.kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int size=(imax+2)*(jmax+2)*(kmax+2);

  __global float *psi = &var[TMP1*size];
  __global float *psi0 = &var[VX*size];
  __global float *aw = &var[AW*size];
  __global float *ae = &var[AE*size];
  __global float *as = &var[AS*size];
  __global float *an = &var[AN*size];
  __global float *af = &var[AF*size];
  __global float *ab = &var[AB*size]; 
  __global float *ap0 = &var[AP0*size]; 
  __global float *b = &var[B*size];
  __global float *x = &var[X*size];
  __global float *y = &var[Y*size];
  __global float *z = &var[Z*size];
	__global float *gx = &var[GX*size];
  __global float *gy = &var[GY*size];
  __global float *gz = &var[GZ*size];
	__global float *u = &var[VX*size];
	__global float *v = &var[VY*size];
	__global float *w = &var[VZ*size];
   
  REAL dxe, dxw, dyn, dys, dzf, dzb, Dx, Dy, Dz;
	REAL uw, ue, vs, vn, wb, wf;
	REAL Fw, Fe, Fs, Fn, Fb, Ff;
  REAL dt = para->mytime.dt;

	if (i>0 && i<imax && j>0 && j<jmax+1 && k>0 && k<kmax+1) {
				// define the dimensions
				dxe = gx[IX(i + 1, j, k)] - gx[IX(i, j, k)];
				dxw = gx[IX(i, j, k)] - gx[IX(i - 1, j, k)];
				Dx = x[IX(i + 1, j, k)] - x[IX(i, j, k)];
				Dy = gy[IX(i, j, k)] - gy[IX(i, j - 1, k)];
				Dz = gz[IX(i, j, k)] - gz[IX(i, j, k - 1)];

				// define the velocity at the surface
				uw = 0.5*(u[IX(i - 1, j, k)] + u[IX(i, j, k)]);
				ue = 0.5*(u[IX(i, j, k)] + u[IX(i + 1, j, k)]);
				vs = 0.5*(v[IX(i, j-1, k)] + v[IX(i + 1, j-1, k)]);
				vn = 0.5*(v[IX(i, j, k)] + v[IX(i + 1, j, k)]);
				wb = 0.5*(w[IX(i, j, k-1)] + w[IX(i+1, j, k-1)]);
				wf = 0.5*(w[IX(i, j, k)] + w[IX(i + 1, j, k)]);

				// define the flow rate at the surface
				Fw = uw * Dy*Dz;
				Fe = ue * Dy*Dz;
				Fs = vs * Dx*Dz;
				Fn = vn * Dx*Dz;
				Fb = wb * Dx*Dy;
				Ff = wf * Dx*Dy;

				// define the coefficient for calculation
				aw[IX(i, j, k)] = max (Fw, 0);
				ae[IX(i, j, k)] = max (-Fe, 0);
				as[IX(i, j, k)] = max(Fs, 0);
				an[IX(i, j, k)] = max (-Fn, 0);
				ab[IX(i, j, k)] = max(Fb, 0);
				af[IX(i, j, k)] = max (-Ff, 0);
				ap0[IX(i, j, k)] = Dx * Dy * Dz / dt;
				b[IX(i, j, k)] = psi0[IX(i, j, k)] * ap0[IX(i, j, k)];

	} //end of if
}

///////////////////////////////////////////////////////////////////////////////
/// kenerl setting diffusion coefficient for VY
/// 6/7/2015
///////////////////////////////////////////////////////////////////////////////

__kernel void diff_VY(__global PARA_DATA_SIMP *para, __global REAL *var,  __global  int *BINDEX) {
  int i = get_global_id(0);
  int j = get_global_id(1);
  int k = get_global_id(2);
  int imax = para->geom.imax, jmax = para->geom.jmax;
  int kmax = para->geom.kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int size=(imax+2)*(jmax+2)*(kmax+2);

  __global float *psi = &var[VY*size];
  __global float *psi0 = &var[TMP2*size];

  __global float *aw = &var[AW*size];
  __global float *ae = &var[AE*size];
  __global float *as = &var[AS*size];
  __global float *an = &var[AN*size];
  __global float *af = &var[AF*size];
  __global float *ab = &var[AB*size]; 
  __global float *ap = &var[AP*size];
  __global float *ap0 = &var[AP0*size]; 
  __global float *b = &var[B*size];
  __global float *gx = &var[GX*size];
  __global float *x = &var[X*size];
  __global float *y = &var[Y*size];
  __global float *z = &var[Z*size];
  __global float *gy = &var[GY*size];
  __global float *gz = &var[GZ*size];
  __global float *pp = &var[PP*size];
  __global float *Temp = &var[TEMP*size];
  __global float *flagp = &var[FLAGP*size];
   
  REAL dxe, dxw, dyn, dys, dzf, dzb, Dx, Dy, Dz;
  REAL dt = para->mytime.dt, beta = para->prob.beta;
  REAL Temp_Buoyancy = para->prob.Temp_Buoyancy;
  REAL gravx = para->prob.gravx, gravy = para->prob.gravy,
       gravz = para->prob.gravz;
  REAL kapa=para->prob.nu;
  REAL coef_CONSTANT = 120;

  // FOR_V_CELL
  if (i>0 && i<imax+1 && j>0 && j<jmax && k>0 && k<kmax+1) {

				if (para->prob.tur_model==CHEN) {
						kapa = nu_t_chen_zero_equ(para, var, &i, &j, &k)+ para->prob.nu;
				}
        dxe = x[IX(i+1,j,k)] - x[IX(i,j,k)];
        dxw = x[IX(i,j,k)] - x[IX(i-1,j,k)];
        dyn = gy[IX(i,j+1,k)] - gy[IX(i,j,k)];
        dys = gy[IX(i,j,k)] - gy[IX(i,j-1,k)];
        dzf = z[IX(i,j,k+1)] - z[IX(i,j,k)];
        dzb = z[IX(i,j,k)] - z[IX(i,j,k-1)];
        Dx = gx[IX(i,j,k)] - gx[IX(i-1,j,k)];
        Dy = y[IX(i,j+1,k)] - y[IX(i,j,k)];
        Dz = gz[IX(i,j,k)] - gz[IX(i,j,k-1)];

        aw[IX(i,j,k)] = kapa*Dy*Dz/dxw;
        ae[IX(i,j,k)] = kapa*Dy*Dz/dxe;
        an[IX(i,j,k)] = kapa*Dx*Dz/dyn;
        as[IX(i,j,k)] = kapa*Dx*Dz/dys;
        af[IX(i,j,k)] = kapa*Dx*Dy/dzf;
        ab[IX(i,j,k)] = kapa*Dx*Dy/dzb;
        ap0[IX(i,j,k)] = Dx*Dy*Dz/dt;
        b[IX(i,j,k)] = psi0[IX(i,j,k)]*ap0[IX(i,j,k)]
                     - beta*gravy*(Temp[IX(i,j,k)]-Temp_Buoyancy)*Dx*Dy*Dz
                     + (pp[IX(i,j,k)]-pp[IX(i ,j+1,k)])*Dx*Dz;

  } //end of FOR_V_CELL
} //end of __kernel void diff_VY()

///////////////////////////////////////////////////////////////////////////////
/// kenerl setting advection coefficient for VY
/// 8/7/2015
///////////////////////////////////////////////////////////////////////////////
__kernel void adve_VY_im(__global PARA_DATA_SIMP *para, __global REAL *var, __global  int *BINDEX) {
  int i = get_global_id(0);
  int j = get_global_id(1);
  int k = get_global_id(2);
  int imax = para->geom.imax, jmax = para->geom.jmax;
  int kmax = para->geom.kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int size=(imax+2)*(jmax+2)*(kmax+2);

  __global float *psi = &var[TMP2*size];
  __global float *psi0 = &var[VY*size];
  __global float *aw = &var[AW*size];
  __global float *ae = &var[AE*size];
  __global float *as = &var[AS*size];
  __global float *an = &var[AN*size];
  __global float *af = &var[AF*size];
  __global float *ab = &var[AB*size]; 
  __global float *ap0 = &var[AP0*size]; 
  __global float *b = &var[B*size];
  __global float *x = &var[X*size];
  __global float *y = &var[Y*size];
  __global float *z = &var[Z*size];
	__global float *gx = &var[GX*size];
  __global float *gy = &var[GY*size];
  __global float *gz = &var[GZ*size];
	__global float *u = &var[VX*size];
	__global float *v = &var[VY*size];
	__global float *w = &var[VZ*size];
   
  REAL dxe, dxw, dyn, dys, dzf, dzb, Dx, Dy, Dz;
	REAL uw, ue, vs, vn, wb, wf;
	REAL Fw, Fe, Fs, Fn, Fb, Ff;
  REAL dt = para->mytime.dt;

	if (i>0 && i<imax+1 && j>0 && j<jmax && k>0 && k<kmax+1) {
				// define the dimensions
				dyn = gy[IX(i, j + 1, k)] - gy[IX(i, j, k)];
				dys = gy[IX(i, j, k)] - gy[IX(i, j - 1, k)];
				Dx = gx[IX(i, j, k)] - gx[IX(i - 1, j, k)];
				Dy = y[IX(i, j + 1, k)] - y[IX(i, j, k)];
				Dz = gz[IX(i, j, k)] - gz[IX(i, j, k - 1)];

				// define the velocity at the surface
				uw = 0.5*(u[IX(i - 1, j, k)] + u[IX(i - 1, j + 1, k)]);
				ue = 0.5*(u[IX(i, j, k)] + u[IX(i, j + 1, k)]);
				vs = 0.5*(v[IX(i, j - 1, k)] + v[IX(i, j, k)]);
				vn = 0.5*(v[IX(i, j, k)] + v[IX(i, j + 1, k)]);
				wb = 0.5*(w[IX(i, j, k - 1)] + w[IX(i, j + 1, k - 1)]);
				wf = 0.5*(w[IX(i, j, k)] + w[IX(i, j +1, k)]);

				// define the flow rate at the surface
				Fw = uw * Dy*Dz;
				Fe = ue * Dy*Dz;
				Fs = vs * Dx*Dz;
				Fn = vn * Dx*Dz;
				Fb = wb * Dx*Dy;
				Ff = wf * Dx*Dy;

				// define the coefficient for calculation
				aw[IX(i, j, k)] = max(Fw, 0);
				ae[IX(i, j, k)] = max(-Fe, 0);
				as[IX(i, j, k)] = max(Fs, 0);
				an[IX(i, j, k)] = max(-Fn, 0);
				ab[IX(i, j, k)] = max(Fb, 0);
				af[IX(i, j, k)] = max(-Ff, 0);
				ap0[IX(i, j, k)] = Dx * Dy * Dz / dt;
				b[IX(i, j, k)] = psi0[IX(i, j, k)] * ap0[IX(i, j, k)];

	} //end of FOR_V_CELL
}

///////////////////////////////////////////////////////////////////////////////
/// kenerl setting diffusion coefficient for VZ
/// 6/7/2015
///////////////////////////////////////////////////////////////////////////////

__kernel void diff_VZ(__global PARA_DATA_SIMP *para, __global REAL *var, __global  int *BINDEX) {

  int i = get_global_id(0);
  int j = get_global_id(1);
  int k = get_global_id(2);
  int imax = para->geom.imax, jmax = para->geom.jmax;
  int kmax = para->geom.kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int size=(imax+2)*(jmax+2)*(kmax+2);

  __global float *psi = &var[VZ*size];
  __global float *psi0 = &var[TMP3*size];

  __global float *aw = &var[AW*size];
  __global float *ae = &var[AE*size];
  __global float *as = &var[AS*size];
  __global float *an = &var[AN*size];
  __global float *af = &var[AF*size];
  __global float *ab = &var[AB*size]; 
  __global float *ap = &var[AP*size];
  __global float *ap0 = &var[AP0*size]; 
  __global float *b = &var[B*size];
  __global float *gx = &var[GX*size];
  __global float *x = &var[X*size];
  __global float *y = &var[Y*size];
  __global float*z = &var[Z*size];
  __global float *gy = &var[GY*size];
  __global float *gz = &var[GZ*size];
  __global float *pp = &var[PP*size];
  __global float *Temp = &var[TEMP*size];
  __global float *flagp = &var[FLAGP*size];
   
  REAL dxe, dxw, dyn, dys, dzf, dzb, Dx, Dy, Dz;
  REAL dt = para->mytime.dt, beta = para->prob.beta;
  REAL Temp_Buoyancy = para->prob.Temp_Buoyancy;
  REAL gravx = para->prob.gravx, gravy = para->prob.gravy,
       gravz = para->prob.gravz;
  REAL kapa=para->prob.nu;
  REAL coef_CONSTANT = 120;

  // FOR_W_CELL
  if (i>0 && i<imax+1 && j>0 && j<jmax+1 && k>0 && k<kmax) {

				if (para->prob.tur_model==CHEN) {
						kapa = nu_t_chen_zero_equ(para, var, &i, &j, &k)+ para->prob.nu;
				}

        dxe = x[IX(i+1,j,k)] - x[IX(i,j,k)];
        dxw = x[IX(i,j,k)] - x[IX(i-1,j,k)];
        dyn = y[IX(i,j+1,k)] - y[IX(i,j,k)];
        dys = y[IX(i,j,k)] - y[IX(i,j-1,k)];
        dzf = gz[IX(i,j,k+1)] - gz[IX(i,j,k)];
        dzb = gz[IX(i,j,k)] - gz[IX(i,j,k-1)];
        Dx = gx[IX(i,j,k)] - gx[IX(i-1,j,k)];
        Dy = gy[IX(i,j,k)] - gy[IX(i,j-1,k)];
        Dz = z[IX(i,j,k+1)] - z[IX(i,j,k)];

        aw[IX(i,j,k)] = kapa*Dy*Dz/dxw;
        ae[IX(i,j,k)] = kapa*Dy*Dz/dxe;
        an[IX(i,j,k)] = kapa*Dx*Dz/dyn;
        as[IX(i,j,k)] = kapa*Dx*Dz/dys;
        af[IX(i,j,k)] = kapa*Dx*Dy/dzf;
        ab[IX(i,j,k)] = kapa*Dx*Dy/dzb;
        ap0[IX(i,j,k)] = Dx*Dy*Dz/dt;
        b[IX(i,j,k)] = psi0[IX(i,j,k)]*ap0[IX(i,j,k)]
                     - beta*gravz*(Temp[IX(i,j,k)]-Temp_Buoyancy)*Dx*Dy*Dz
                     + (pp[IX(i,j,k)]-pp[IX(i ,j,k+1)])*Dy*Dx;
        
  } //end of FOR_W_CELL
}// end of __kernel void diff_VZ()

///////////////////////////////////////////////////////////////////////////////
/// kenerl setting advection coefficient for VZ
/// 8/7/2015
///////////////////////////////////////////////////////////////////////////////
__kernel void adve_VZ_im(__global PARA_DATA_SIMP *para, __global REAL *var, __global  int *BINDEX) {
  int i = get_global_id(0);
  int j = get_global_id(1);
  int k = get_global_id(2);
  int imax = para->geom.imax, jmax = para->geom.jmax;
  int kmax = para->geom.kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int size=(imax+2)*(jmax+2)*(kmax+2);

  __global float *psi = &var[TMP3*size];
  __global float *psi0 = &var[VZ*size];
  __global float *aw = &var[AW*size];
  __global float *ae = &var[AE*size];
  __global float *as = &var[AS*size];
  __global float *an = &var[AN*size];
  __global float *af = &var[AF*size];
  __global float *ab = &var[AB*size]; 
  __global float *ap0 = &var[AP0*size]; 
  __global float *b = &var[B*size];
  __global float *x = &var[X*size];
  __global float *y = &var[Y*size];
  __global float *z = &var[Z*size];
	__global float *gx = &var[GX*size];
  __global float *gy = &var[GY*size];
  __global float *gz = &var[GZ*size];
	__global float *u = &var[VX*size];
	__global float *v = &var[VY*size];
	__global float *w = &var[VZ*size];
   
  REAL dxe, dxw, dyn, dys, dzf, dzb, Dx, Dy, Dz;
	REAL uw, ue, vs, vn, wb, wf;
	REAL Fw, Fe, Fs, Fn, Fb, Ff;
  REAL dt = para->mytime.dt;

  // FOR_W_CELL
  if (i>0 && i<imax+1 && j>0 && j<jmax+1 && k>0 && k<kmax) {
	  // define the dimensions
				dzf = gz[IX(i, j, k + 1)] - gz[IX(i, j, k)];
				dzb = gz[IX(i, j, k)] - gz[IX(i, j, k - 1)];
				Dx = gx[IX(i, j, k)] - gx[IX(i - 1, j, k)];
				Dy = gy[IX(i, j, k)] - gy[IX(i, j - 1, k)];
				Dz = z[IX(i, j, k + 1)] - z[IX(i, j, k)];

				// define the velocity at the surface
				uw = 0.5*(u[IX(i - 1, j, k)] + u[IX(i - 1, j, k + 1)]);
				ue = 0.5*(u[IX(i, j, k)] + u[IX(i, j, k + 1)]);
				vs = 0.5*(v[IX(i, j - 1, k)] + v[IX(i, j - 1, k + 1)]);
				vn = 0.5*(v[IX(i, j, k)] + v[IX(i, j, k + 1)]);
				wb = 0.5*(w[IX(i, j, k - 1)] + w[IX(i, j, k)]);
				wf = 0.5*(w[IX(i, j, k)] + w[IX(i, j, k + 1)]);


				// define the flow rate at the surface
				Fw = uw * Dy*Dz;
				Fe = ue * Dy*Dz;
				Fs = vs * Dx*Dz;
				Fn = vn * Dx*Dz;
				Fb = wb * Dx*Dy;
				Ff = wf * Dx*Dy;

				// define the coefficient for calculation
				aw[IX(i, j, k)] = max(Fw, 0);
				ae[IX(i, j, k)] = max(-Fe, 0);
				as[IX(i, j, k)] = max(Fs, 0);
				an[IX(i, j, k)] = max(-Fn, 0);
				ab[IX(i, j, k)] = max(Fb, 0);
				af[IX(i, j, k)] = max(-Ff, 0);
				ap0[IX(i, j, k)] = Dx * Dy * Dz / dt;
				b[IX(i, j, k)] = psi0[IX(i, j, k)] * ap0[IX(i, j, k)];

	} //end of FOR_W_CELL
}

///////////////////////////////////////////////////////////////////////////////
/// kenerl setting diffusion coefficient for T
/// 6/7/2015
///////////////////////////////////////////////////////////////////////////////

__kernel void diff_T(__global PARA_DATA_SIMP *para, __global REAL *var, __global  int *BINDEX) {
  //printf("running diff\n");
  //printf("running diff\n");
  int i = get_global_id(0);
  int j = get_global_id(1);
  int k = get_global_id(2);
  int imax = para->geom.imax, jmax = para->geom.jmax;
  int kmax = para->geom.kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int size=(imax+2)*(jmax+2)*(kmax+2);

  __global float *psi = &var[TEMP*size];
  __global float *psi0 = &var[TMP1*size];

  __global float *aw = &var[AW*size];
  __global float *ae = &var[AE*size];
  __global float *as = &var[AS*size];
  __global float *an = &var[AN*size];
  __global float *af = &var[AF*size];
  __global float *ab = &var[AB*size]; 
  __global float *ap = &var[AP*size];
  __global float *ap0 = &var[AP0*size]; 
  __global float *b = &var[B*size];
  __global float *gx = &var[GX*size];
  __global float *x = &var[X*size];
  __global float *y = &var[Y*size];
  __global float*z = &var[Z*size];
  __global float *gy = &var[GY*size];
  __global float *gz = &var[GZ*size];
  __global float *pp = &var[PP*size];
  __global float *Temp = &var[TEMP*size];
  __global float *flagp = &var[FLAGP*size];
   
  REAL dxe, dxw, dyn, dys, dzf, dzb, Dx, Dy, Dz;
  REAL dt = para->mytime.dt, beta = para->prob.beta;
  REAL Temp_Buoyancy = para->prob.Temp_Buoyancy;
  REAL gravx = para->prob.gravx, gravy = para->prob.gravy,
       gravz = para->prob.gravz;
  REAL kapa;
  REAL coef_CONSTANT = 120;
	REAL kapaE, kapaW, kapaN, kapaS, kapaF, kapaB;
	int i_plus = 0, i_minus = 0, j_plus = 0, j_minus = 0;
	int k_plus = 0, k_minus = 0;
  
 // FOR_ALL_CELL
 if (i>0 && i<imax+1 && j>0 && j<jmax+1 && k>0 && k<kmax+1) {

        dxe = x[IX(i+1,j,k)] - x[IX(i,j,k)];
        dxw = x[IX(i,j,k)] - x[IX(i-1,j,k)];
        dyn = y[IX(i,j+1,k)] - y[IX(i,j,k)];
        dys = y[IX(i,j,k)] - y[IX(i,j-1,k)];
        dzf = z[IX(i,j,k+1)] - z[IX(i,j,k)];
        dzb = z[IX(i,j,k)] - z[IX(i,j,k-1)];
        Dx = gx[IX(i,j,k)] - gx[IX(i-1,j,k)];
        Dy = gy[IX(i,j,k)] - gy[IX(i,j-1,k)];
        Dz = gz[IX(i,j,k)] - gz[IX(i,j,k-1)];

				if (para->prob.tur_model==CHEN) {
						i_plus = i+1;
						i_minus = i-1;
						j_plus = j+1;
						j_minus = j-1;
						k_plus = k+1;
						k_minus = k-1;
						kapa = nu_t_chen_zero_equ(para, var, &i, &j, &k)+ para->prob.nu;
						kapaE = nu_t_chen_zero_equ(para, var, &i_plus, &j, &k) + para->prob.nu;
						kapaW = nu_t_chen_zero_equ(para, var, &i_minus, &j, &k) + para->prob.nu;
						kapaN = nu_t_chen_zero_equ(para, var, &i, &j_plus, &k) + para->prob.nu;
						kapaS = nu_t_chen_zero_equ(para, var, &i, &j_minus, &k) + para->prob.nu;
						kapaF = nu_t_chen_zero_equ(para, var, &i, &j, &k_plus) + para->prob.nu;
						kapaB = nu_t_chen_zero_equ(para, var, &i, &j, &k_minus) + para->prob.nu;
						aw[IX(i, j, k)] = 0.5*(kapa + kapaW)*Dy*Dz / dxw;
						ae[IX(i, j, k)] = 0.5*(kapa + kapaE)*Dy*Dz / dxe;
						an[IX(i, j, k)] = 0.5*(kapa + kapaN)*Dx*Dz / dyn;
						as[IX(i, j, k)] = 0.5*(kapa + kapaS)*Dx*Dz / dys;
						af[IX(i, j, k)] = 0.5*(kapa + kapaF)*Dx*Dy / dzf;
						ab[IX(i, j, k)] = 0.5*(kapa + kapaB)*Dx*Dy / dzb;
				}

				else {
						aw[IX(i,j,k)] = kapa*Dy*Dz/dxw;
						ae[IX(i,j,k)] = kapa*Dy*Dz/dxe;
						an[IX(i,j,k)] = kapa*Dx*Dz/dyn;
						as[IX(i,j,k)] = kapa*Dx*Dz/dys;
						af[IX(i,j,k)] = kapa*Dx*Dy/dzf;
						ab[IX(i,j,k)] = kapa*Dx*Dy/dzb;
				}

        ap0[IX(i,j,k)] = Dx*Dy*Dz/dt;
        b[IX(i,j,k)] = psi0[IX(i,j,k)]*ap0[IX(i,j,k)];

  } // end of FOR_ALL_CELL
} // __kernel void diff_T

///////////////////////////////////////////////////////////////////////////////
/// kenerl setting advection coefficient for T
/// 8/7/2015
///////////////////////////////////////////////////////////////////////////////
__kernel void adve_T_im(__global PARA_DATA_SIMP *para, __global REAL *var, __global  int *BINDEX) {
  int i = get_global_id(0);
  int j = get_global_id(1);
  int k = get_global_id(2);
  int imax = para->geom.imax, jmax = para->geom.jmax;
  int kmax = para->geom.kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int size=(imax+2)*(jmax+2)*(kmax+2);

  __global float *psi = &var[TMP1*size];
  __global float *psi0 = &var[TEMP*size];
  __global float *aw = &var[AW*size];
  __global float *ae = &var[AE*size];
  __global float *as = &var[AS*size];
  __global float *an = &var[AN*size];
  __global float *af = &var[AF*size];
  __global float *ab = &var[AB*size]; 
  __global float *ap0 = &var[AP0*size]; 
  __global float *b = &var[B*size];
  __global float *x = &var[X*size];
  __global float *y = &var[Y*size];
  __global float *z = &var[Z*size];
	__global float *gx = &var[GX*size];
  __global float *gy = &var[GY*size];
  __global float *gz = &var[GZ*size];
	__global float *u = &var[VX*size];
	__global float *v = &var[VY*size];
	__global float *w = &var[VZ*size];
   
  REAL dxe, dxw, dyn, dys, dzf, dzb, Dx, Dy, Dz;
	REAL uw, ue, vs, vn, wb, wf;
	REAL Fw, Fe, Fs, Fn, Fb, Ff;
  REAL dt = para->mytime.dt;

	 if (i>0 && i<imax+1 && j>0 && j<jmax+1 && k>0 && k<kmax+1) {
	      // define the dimensions
				Dx = gx[IX(i, j, k)] - gx[IX(i - 1, j, k)];
				Dy = gy[IX(i, j, k)] - gy[IX(i, j - 1, k)];
				Dz = gz[IX(i, j, k)] - gz[IX(i, j, k - 1)];

				// define the velocity at the surface
				uw = u[IX(i - 1, j, k)];
				ue = u[IX(i, j, k )];
				vs = v[IX(i, j - 1, k)];
				vn = v[IX(i, j, k)];
				wb = w[IX(i, j, k - 1)];
				wf = w[IX(i, j, k)];


				// define the flow rate at the surface
				Fw = uw * Dy*Dz;
				Fe = ue * Dy*Dz;
				Fs = vs * Dx*Dz;
				Fn = vn * Dx*Dz;
				Fb = wb * Dx*Dy;
				Ff = wf * Dx*Dy;

				// define the coefficient for calculation
				aw[IX(i, j, k)] = max(Fw, 0);
				ae[IX(i, j, k)] = max(-Fe, 0);
				as[IX(i, j, k)] = max(Fs, 0);
				an[IX(i, j, k)] = max(-Fn, 0);
				ab[IX(i, j, k)] = max(Fb, 0);
				af[IX(i, j, k)] = max(-Ff, 0);
				ap0[IX(i, j, k)] = Dx * Dy * Dz / dt;
				b[IX(i, j, k)] = psi0[IX(i, j, k)] * ap0[IX(i, j, k)];

	 }// end of FOR_ALL_CELL
}

///////////////////////////////////////////////////////////////////////////////
/// kenerl of advection for VX
/// 6/16/2015
///////////////////////////////////////////////////////////////////////////////
__kernel void adve_VX(__global PARA_DATA_SIMP *para, __global REAL *var, __global  int *BINDEX) {
  int i = get_global_id(0);
  int j = get_global_id(1);
  int k = get_global_id(2);
  int itmax = 20000;
  int it;
  int imax = para->geom.imax, jmax = para->geom.jmax;
  int kmax = para->geom.kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int size=(imax+2)*(jmax+2)*(kmax+2);
  REAL x_1, y_1, z_1;
  REAL dt = para->mytime.dt; 
  REAL u0, v0, w0;
  int  COOD[3], LOC[3];
  REAL OL[3];
  int  OC[3];

  // global variables
  __global REAL *x = &var[X*size], *y = &var[Y*size],  *z = &var[Z*size]; 
  __global REAL *gx = &var[GX*size]; 
  __global REAL *u = &var[VX*size], *v = &var[VY*size], *w = &var[VZ*size];
  __global REAL *flagu = &var[FLAGU*size], *flagp = &var[FLAGP*size];

  // Acustomized variables
  // In advection, VX is used to calcalute new value to update temporal TMP1
  // Then TMP1 is used in diffusion to update the VX
  // Therefore in diffusion the order of d and d0 is reversed
  // note by wei, 6/17/2015
  // checkpoint: adve_VX()

  __global REAL *d = &var[TMP1*size];
  __global REAL *d0 = &var[VX*size];

  // FOR_U_CELL and the cell is of fluid type
  if (i>0 && i<imax && j>0 && j<jmax+1 && k>0 && k<kmax+1 && flagu[IX(i,j,k)] < 0) {
      u0 = u[IX(i,j,k)];
    v0 = (REAL) 0.5 
      * ((v[IX(i,  j,k)]+v[IX(i,  j-1,k)])*( x[IX(i+1,j,k)]-gx[IX(i,j,k)])
        +(v[IX(i+1,j,k)]+v[IX(i+1,j-1,k)])*(gx[IX(i,  j,k)]- x[IX(i,j,k)])) 
      / (x[IX(i+1,j,k)]-x[IX(i,j,k)]);
    w0 = (REAL) 0.5 
      * ((w[IX(i,  j,k)]+w[IX(i  ,j, k-1)])*( x[IX(i+1,j,k)]-gx[IX(i,j,k)])
        +(w[IX(i+1,j,k)]+w[IX(i+1,j, k-1)])*(gx[IX(i,  j,k)]- x[IX(i,j,k)]))
      / (x[IX(i+1,j,k)]-x[IX(i,j,k)]); 
    // Find the location at previous time step
    OL[X] =gx[IX(i,j,k)] - u0*dt;
    OL[Y] = y[IX(i,j,k)] - v0*dt;
    OL[Z] = z[IX(i,j,k)] - w0*dt;
    // Initialize the coordinates of previous step 
    OC[X] = i; 
    OC[Y] = j; 
    OC[Z] = k; //OC is the trace back coordinates
    // Initialize the signs for track process
    // Completed: 0; In process: 1
    COOD[X] = 1; 
    COOD[Y] = 1; 
    COOD[Z] = 1; 
    // Initialize the signs for recording if the track back hits the boundary
    // Hit the boundary: 0; Not hit the boundary: 1
    LOC[X] = 1; 
    LOC[Y] = 1; 
    LOC[Z] = 1;
    //Initialize the number of iterations
    it=1;

    // Trace back more if the any of the trace is still in process 
    while(COOD[X]==1 || COOD[Y] ==1 || COOD[Z] ==1)
    {
      it++;
      // If trace in X is in process and donot hit the boundary
      if(COOD[X]==1 && LOC[X]==1) 
      set_x_location(para, var, flagu, gx, &u0, &i, &j, &k, OL, OC, LOC, COOD); 
      // If trace in Y is in process and donot hit the boundary
      if(COOD[Y]==1 && LOC[Y]==1) 
      set_y_location(para, var, flagu, y, &v0, &i, &j, &k, OL, OC, LOC, COOD); 
      // If trace in Z is in process and donot hit the boundary
      if(COOD[Z]==1 && LOC[Z]==1) 
      set_z_location(para, var, flagu, z, &w0, &i, &j, &k, OL, OC, LOC, COOD); 
    
      //if(it>itmax)
      //{
      //printf("Error: advection.c, can not track the location for VX(%d, %d,%d)",
          //i, j, k);
      //printf("after %d iterations.\n", it);
      //}
    } // End of while() for backward tracing

    // Set the coordinates of previous location if it is as boundary
    if(u0>0 && LOC[X]==0) OC[X] -=1;
    if(v0>0 && LOC[Y]==0) OC[Y] -=1;
    if(w0>0 && LOC[Z]==0) OC[Z] -=1;

    // LOC=1 means that traced back point hits boundary.
    // When u is less than zero, the traced back point is towards west. 
    // When it hits the boundary, according to the code, 
    // the east cell of west boundary will be given as the traced back location
    // and the index (OC) is also marching on one cell. 
    // This means that the index and traced back location are on same cell. 
    // That is why the index needs to move backward for one cell 
    // in order to do interpolation.
    if(u0<0 && LOC[X]==1) OC[X] -=1; 
    if(v0<0 && LOC[Y]==1) OC[Y] -=1;
    if(w0<0 && LOC[Z]==1) OC[Z] -=1;

        // interpolate
    x_1 = (OL[X]-gx[IX(OC[X],OC[Y],OC[Z])]) 
      / (gx[IX(OC[X]+1,OC[Y],OC[Z])]-gx[IX(OC[X],OC[Y],OC[Z])]); 
    y_1 = (OL[Y]-y[IX(OC[X],OC[Y],OC[Z])])
      / (y[IX(OC[X],OC[Y]+1,OC[Z])]-y[IX(OC[X],OC[Y],OC[Z])]);
    z_1 = (OL[Z]-z[IX(OC[X],OC[Y],OC[Z])]) 
      / (z[IX(OC[X],OC[Y],OC[Z]+1)]-z[IX(OC[X],OC[Y],OC[Z])]);
      // checkpoint: adve_VX()
    d[IX(i,j,k)] = interpolation(para, d0, &x_1, &y_1, &z_1, &OC[X], &OC[Y], &OC[Z]);
  } //end of FOR_U_CELL
} //end of kernel adve_VX()

///////////////////////////////////////////////////////////////////////////////
/// kenerl of advection for VY
/// 6/16/2015
///////////////////////////////////////////////////////////////////////////////
__kernel void adve_VY(__global PARA_DATA_SIMP *para, __global REAL *var, __global  int *BINDEX) {
  int i = get_global_id(0);
  int j = get_global_id(1);
  int k = get_global_id(2);
  int it;
  int itmax = 20000; // Max number of iterations for backward tracing 
  int imax = para->geom.imax, jmax = para->geom.jmax;
  int kmax = para->geom.kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int size=(imax+2)*(jmax+2)*(kmax+2);
  REAL x_1, y_1, z_1;
  REAL dt = para->mytime.dt; 
  REAL u0, v0, w0;

  int  COOD[3], LOC[3];
  REAL OL[3];
  int  OC[3];

  // global variables
  __global REAL *x = &var[X*size], *y = &var[Y*size],  *z = &var[Z*size]; 
  __global REAL *gy = &var[GY*size]; 
  __global REAL *u = &var[VX*size], *v = &var[VY*size], *w = &var[VZ*size];
  __global REAL *flagv = &var[FLAGV*size], *flagp = &var[FLAGP*size];

  // Acustomized variables
  // In advection, VY is used to calcalute new value to update temporal TMP2
  // Then TMP2 is used in diffusion to update the VY
  // Therefore in diffusion the order of d and d0 is reversed
  // note by wei, 6/17/2015
  // checkpoint: adve_VY()

  __global REAL *d = &var[TMP2*size];
  __global REAL *d0 = &var[VY*size];

  //FOR_V_CELL
  if (i>0 && i<imax+1 && j>0 && j<jmax && k>0 && k<kmax+1 && flagv[IX(i,j,k)] < 0) {

    /*-------------------------------------------------------------------------
    | Step 1: Tracing Back
    -------------------------------------------------------------------------*/
    // Get velocities at the location of VY
    u0 = (REAL) 0.5
       * ((u[IX(i,j,k)]+u[IX(i-1,j,  k)])*(y [IX(i,j+1,k)]-gy[IX(i,j,k)])
         +(u[IX(i,j+1,k)]+u[IX(i-1,j+1,k)])*(gy[IX(i,j,  k)]-y[IX(i,j,k)]))
       / (y[IX(i,j+1,k)]-y[IX(i,j,k)]);
    v0 = v[IX(i,j,k)]; 
    w0 = (REAL) 0.5
       * ((w[IX(i,j,k)]+w[IX(i,j,k-1)])*(y[IX(i,j+1,k)]-gy[IX(i,j,k)])
         +(w[IX(i,j+1,k)]+w[IX(i,j+1,k-1)])*(gy[IX(i,j,k)]-y[IX(i,j,k)]))
       / (y[IX(i,j+1,k)]-y[IX(i,j,k)]); 
    // Find the location at previous time step
    OL[X] = x[IX(i,j,k)] - u0*dt; 
    OL[Y] = gy[IX(i,j,k)] - v0*dt;
    OL[Z] = z[IX(i,j,k)] - w0*dt;
    // Initialize the coordinates of previous step
    OC[X] = i;
    OC[Y] = j;
    OC[Z] = k;
    // Initialize the signs for track process
    // Completed: 0; In process: 1
    COOD[X] = 1; 
    COOD[Y] = 1; 
    COOD[Z] = 1;
    // Initialize the signs for recording if the track back hits the boundary
    // Hit the boundary: 0; Not hit the boundary: 1
    LOC[X] = 1; 
    LOC[Y] = 1; 
    LOC[Z] = 1;
    //Initialize the number of iterations
    it=1;

    // Trace back more if the any of the trace is still in process 
    while(COOD[X]==1 || COOD[Y] ==1 || COOD[Z] == 1) {
      it++;
      if(COOD[X]==1 && LOC[X]==1)
        set_x_location(para, var, flagv, x, &u0, &i, &j, &k, OL, OC, LOC, COOD); 
      if(COOD[Y]==1 && LOC[Y]==1)
        set_y_location(para, var, flagv, gy, &v0, &i, &j, &k, OL, OC, LOC, COOD); 
      if(COOD[Z]==1 && LOC[Z]==1)
        set_z_location(para, var, flagv, z, &w0, &i, &j, &k, OL, OC, LOC, COOD); 

      //if(it>itmax) {
        //printf("Error: advection.c can not track the location for VY(%d, %d,%d)",
                //i, j, k);
        //printf("after %d iterations.\n", it);
      //}
    } // End of while() loop

    // Set the coordinates of previous location if it is as boundary
    if(u0>=0 && LOC[X]==0) OC[X] -= 1; 
    if(v0>=0 && LOC[Y]==0) OC[Y] -= 1;
    if(w0>=0 && LOC[Z]==0) OC[Z] -= 1;

    if(u0<0 && LOC[X]==1) OC[X] -= 1;
    if(v0<0 && LOC[Y]==1) OC[Y] -= 1;
    if(w0<0 && LOC[Z]==1) OC[Z] -= 1;
            
    /*-------------------------------------------------------------------------
    | Interpolating for all variables
    -------------------------------------------------------------------------*/
    x_1 = (OL[X]-x[IX(OC[X],OC[Y],OC[Z])])
        / (x[IX(OC[X]+1,OC[Y],OC[Z])]-x[IX(OC[X],OC[Y],OC[Z])]); 
    y_1 = (OL[Y]-gy[IX(OC[X],OC[Y],OC[Z])])
        / (gy[IX(OC[X],OC[Y]+1,OC[Z])]-gy[IX(OC[X],OC[Y],OC[Z])]);
    z_1 = (OL[Z]-z[IX(OC[X],OC[Y],OC[Z])])
        / (z[IX(OC[X],OC[Y],OC[Z]+1)]-z[IX(OC[X],OC[Y],OC[Z])]);
    d[IX(i,j,k)] = interpolation(para, d0, &x_1, &y_1, &z_1, &OC[X],&OC[Y],&OC[Z]);
  } // end of FOR_V_CELL
} // End of trace_vy()

///////////////////////////////////////////////////////////////////////////////
/// kenerl of advection for VZ
/// 6/17/2015
///////////////////////////////////////////////////////////////////////////////
__kernel void adve_VZ(__global PARA_DATA_SIMP *para, __global REAL *var, __global  int *BINDEX) {
  int i = get_global_id(0);
  int j = get_global_id(1);
  int k = get_global_id(2);
  int it;
  int itmax = 20000; // Max number of iterations for backward tracing 
  int imax = para->geom.imax, jmax = para->geom.jmax;
  int kmax = para->geom.kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int size=(imax+2)*(jmax+2)*(kmax+2);
  REAL x_1, y_1, z_1;
  REAL dt = para->mytime.dt; 
  REAL u0, v0, w0;

  int  COOD[3], LOC[3];
  REAL OL[3];
  int  OC[3];

  // global variables
  __global REAL *x = &var[X*size], *y = &var[Y*size],  *z = &var[Z*size]; 
  __global REAL *gz = &var[GZ*size]; 
  __global REAL *u = &var[VX*size], *v = &var[VY*size], *w = &var[VZ*size];
  __global REAL *flagw = &var[FLAGW*size], *flagp = &var[FLAGP*size];

  // Acustomized variables
  // In advection, VZ is used to calcalute new value to update temporal TMP3
  // Then TMP3 is used in diffusion to update the VZ
  // Therefore in diffusion the order of d and d0 is reversed
  // note by wei, 6/17/2015
  // checkpoint: adve_VZ()

  __global REAL *d = &var[TMP3*size];
  __global REAL *d0 = &var[VZ*size];

  //FOR_W_CELL
  if (i>0 && i<imax+1 && j>0 && j<jmax+1 && k>0 && k<kmax && flagw[IX(i,j,k)] < 0) {
    u0 = (REAL) 0.5 
       * ((u[IX(i,j,k  )]+u[IX(i-1,j,k  )])*(z [IX(i,j,k+1)]-gz[IX(i,j,k)])
         +(u[IX(i,j,k+1)]+u[IX(i-1,j,k+1)])*(gz[IX(i,j,k  )]- z[IX(i,j,k)]))
       /  (z[IX(i,j,k+1)]-z[IX(i,j,k)]);
    v0 = (REAL) 0.5
       * ((v[IX(i,j,k  )]+v[IX(i,j-1,k  )])*(z [IX(i,j,k+1)]-gz[IX(i,j,k)])
       +(v[IX(i,j,k+1)]+v[IX(i,j-1,k+1)])*(gz[IX(i,j,k  )]-z [IX(i,j,k)]))
       /  (z[IX(i,j,k+1)]-z[IX(i,j,k)]); 
    w0 = w[IX(i,j,k)]; 
    // Find the location at previous time step
    OL[X] = x[IX(i,j,k)] - u0*dt; 
    OL[Y] = y[IX(i,j,k)] - v0*dt;
    OL[Z] = gz[IX(i,j,k)] - w0*dt;
    // Initialize the coordinates of previous step
    OC[X] = i;
    OC[Y] = j;
    OC[Z] = k;
    // Initialize the signs for track process
    // Completed: 0; In process: 1
    COOD[X] = 1;
    COOD[Y] = 1;
    COOD[Z] = 1;
    // Initialize the signs for recording if the track back hits the boundary
    // Hit the boundary: 0; Not hit the boundary: 1
    LOC[X] = 1;
    LOC[Y] = 1;
    LOC[Z] = 1;
    //Initialize the number of iterations
    it=1;

    // Trace back more if the any of the trace is still in process 
    while(COOD[X]==1 || COOD[Y] ==1 || COOD[Z] == 1) {
      it++;
      if(COOD[X]==1 && LOC[X]==1)
        set_x_location(para, var, flagw, x, &u0, &i, &j, &k, OL, OC, LOC, COOD); 
      if(COOD[Y]==1 && LOC[Y]==1)
        set_y_location(para, var, flagw, y, &v0, &i, &j, &k, OL, OC, LOC, COOD); 
      if(COOD[Z]==1 && LOC[Z]==1)
        set_z_location(para, var, flagw, gz, &w0, &i, &j, &k, OL, OC, LOC, COOD); 

      //if(it>itmax) {
        //printf("Error: advection.c can not track the location for VZ(%d, %d,%d)",
                //i, j, k);
        //printf("after %d iterations.\n", it);
      //}
    } // End of while() loop

    // Set the coordinates of previous location if it is as boundary
    if(u0>=0 && LOC[X] == 0) OC[X] -=1;
    if(v0>=0 && LOC[Y] == 0) OC[Y] -=1;
    if(w0>=0 && LOC[Z] == 0) OC[Z] -=1;

    if(u0<0 && LOC[X]==1) OC[X] -=1;
    if(v0<0 && LOC[Y]==1) OC[Y] -=1;
    if(w0<0 && LOC[Z]==1) OC[Z] -=1;

    /*-------------------------------------------------------------------------
    | Interpolating for all variables
    -------------------------------------------------------------------------*/
    x_1 = (OL[X]- x[IX(OC[X],OC[Y],OC[Z])])
        / ( x[IX(OC[X]+1,OC[Y],   OC[Z]  )]- x[IX(OC[X],OC[Y],OC[Z])]); 
    y_1 = (OL[Y]- y[IX(OC[X],OC[Y],OC[Z])])
        / ( y[IX(OC[X],  OC[Y]+1, OC[Z]  )]- y[IX(OC[X],OC[Y],OC[Z])]);
    z_1 = (OL[Z]-gz[IX(OC[X],OC[Y],OC[Z])])
        / (gz[IX(OC[X],  OC[Y],   OC[Z]+1)]-gz[IX(OC[X],OC[Y],OC[Z])]);
     d[IX(i,j,k)] = interpolation(para, d0, &x_1, &y_1, &z_1, &OC[X], &OC[Y], &OC[Z]);
  }// end of FOR_W_CELL

}// end of adve_VZ()

///////////////////////////////////////////////////////////////////////////////
/// kenerl of advection for T
/// 6/17/2015
///////////////////////////////////////////////////////////////////////////////
__kernel void adve_T(__global PARA_DATA_SIMP *para, __global REAL *var, __global  int *BINDEX) {
  int i = get_global_id(0);
  int j = get_global_id(1);
  int k = get_global_id(2);
  int it;
  int itmax = 20000; // Max number of iterations for backward tracing 
  int imax = para->geom.imax, jmax = para->geom.jmax;
  int kmax = para->geom.kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int size=(imax+2)*(jmax+2)*(kmax+2);
  REAL x_1, y_1, z_1;
  REAL dt = para->mytime.dt; 
  REAL u0, v0, w0;

  int  COOD[3], LOC[3];
  REAL OL[3];
  int  OC[3];

  // global variables
  __global REAL *x = &var[X*size], *y = &var[Y*size],  *z = &var[Z*size]; 
  __global REAL *gz = &var[GZ*size]; 
  __global REAL *u = &var[VX*size], *v = &var[VY*size], *w = &var[VZ*size];
  __global REAL *flagp = &var[FLAGP*size];

  // Acustomized variables
  // In advection, TEMP is used to calcalute new value to update temporal TMP1
  // Then TMP3 is used in diffusion to update the TEMP
  // Therefore in diffusion the order of d and d0 is reversed
  // note by wei, 6/17/2015
  // checkpoint: adve_T()

  __global REAL *d = &var[TMP1*size];
  __global REAL *d0 = &var[TEMP*size];

  //FOR_EACH_CELL
  if (i>0 && i<imax+1 && j>0 && j<jmax+1 && k>0 && k<kmax+1 && flagp[IX(i,j,k)] < 0) {
    u0 = (REAL) 0.5 * (u[IX(i,j,k)]+u[IX(i-1,j,k  )]);
    v0 = (REAL) 0.5 * (v[IX(i,j,k)]+v[IX(i,j-1,k  )]);
    w0 = (REAL) 0.5 * (w[IX(i,j,k)]+w[IX(i,j  ,k-1)]);
    // Find the location at previous time step
    OL[X] = x[IX(i,j,k)] - u0*dt; 
    OL[Y] = y[IX(i,j,k)] - v0*dt;
    OL[Z] = z[IX(i,j,k)] - w0*dt;
    // Initialize the coordinates of previous step
    OC[X] = i; 
    OC[Y] = j; 
    OC[Z] = k;  
    // Initialize the signs for tracing process
    // Completed: 0; In process: 1
    COOD[X] = 1; 
    COOD[Y] = 1; 
    COOD[Z] = 1;
    // Initialize the flags for recording if the tracing back hits the boundary
    // Hit the boundary: 0; Not hit the boundary: 1
    LOC[X] = 1; 
    LOC[Y] = 1; 
    LOC[Z] = 1;
    //Initialize the number of iterations
    it=1;

    // Trace back more if the any of the trace is still in process 
    while(COOD[X]==1 || COOD[Y] ==1 || COOD[Z] == 1) {
      it++;
      // If trace in X is in process and donot hit the boundary
      if(COOD[X]==1 && LOC[X]==1)
        set_x_location(para, var, flagp, x, &u0, &i, &j, &k, OL, OC, LOC, COOD);
      // If trace in Y is in process and donot hit the boundary
      if(COOD[Y]==1 && LOC[Y]==1)
        set_y_location(para, var, flagp, y, &v0, &i, &j, &k, OL, OC, LOC, COOD);
      // If trace in Z is in process and donot hit the boundary
      if(COOD[Z]==1 && LOC[Z]==1)
        set_z_location(para, var, flagp, z, &w0, &i, &j, &k, OL, OC, LOC, COOD); 
      //if(it>itmax) {
        //printf("trace_scalar(): Could not track the location for scalar "
          //"variable %d at cell(%d, %d,%d) after %d iterations", 
          //TEMP, i, j, k, it);
      //}
    } // End of while() for backward tracing

    // Set the coordinates of previous location if it is as boundary
    if(u0>=0 && LOC[X]==0) OC[X] -=1;
    if(v0>=0 && LOC[Y]==0) OC[Y] -=1;
    if(w0>=0 && LOC[Z]==0) OC[Z] -=1;

    if(u0<0 && LOC[X]==1) OC[X] -=1;
    if(v0<0 && LOC[Y]==1) OC[Y] -=1;
    if(w0<0 && LOC[Z]==1) OC[Z] -=1;

    //Store the local minimum and maximum values
    //var[LOCMIN][IX(i,j,k)]=check_min(para, d0, OC[X], OC[Y], OC[Z]); 
    //var[LOCMAX][IX(i,j,k)]=check_max(para, d0, OC[X], OC[Y], OC[Z]); 

    /*-------------------------------------------------------------------------
    | Interpolate
    -------------------------------------------------------------------------*/
    x_1 = (OL[X]- x[IX(OC[X],OC[Y],OC[Z])])
        / ( x[IX(OC[X]+1,OC[Y],   OC[Z]  )] - x[IX(OC[X],OC[Y],OC[Z])]); 
    y_1 = (OL[Y]- y[IX(OC[X],OC[Y],OC[Z])])
        / ( y[IX(OC[X],  OC[Y]+1, OC[Z]  )] - y[IX(OC[X],OC[Y],OC[Z])]);
    z_1 = (OL[Z]- z[IX(OC[X],OC[Y],OC[Z])])
        / ( z[IX(OC[X],  OC[Y],   OC[Z]+1)] - z[IX(OC[X],OC[Y],OC[Z])]);
    d[IX(i,j,k)] = interpolation(para, d0, &x_1, &y_1, &z_1, &OC[X], &OC[Y], &OC[Z]);

  } //end of FOR_EACH_CELL
} //end of adve_T()

__kernel void project(__global PARA_DATA_SIMP *para, __global REAL *var, __global int *BINDEX) {
  int i = get_global_id(0);
  int j = get_global_id(1);
  int k = get_global_id(2);
  int bar;
  int imax = para->geom.imax, jmax = para->geom.jmax;
  int kmax = para->geom.kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int size=(imax+2)*(jmax+2)*(kmax+2);
  REAL dt= para->mytime.dt;
  __global REAL *x = &var[X*size], *y = &var[Y*size], *z = &var[Z*size];
  __global REAL *gx = &var[GX*size], *gy = &var[GY*size], *gz = &var[GZ*size];
  __global REAL *u = &var[VX*size], *v = &var[VY*size], *w = &var[VZ*size];
  __global REAL *p = &var[IP*size], *b = &var[B*size], *ap = &var[AP*size], *ab = &var[AB*size], *af = &var[AF*size];
  __global REAL *ae = &var[AE*size], *aw = &var[AW*size], *an = &var[AN*size], *as = &var[AS*size];
  REAL dxe,dxw, dyn,dys,dzf,dzb,Dx,Dy,Dz;
  __global REAL *flagu = &var[FLAGU*size],*flagv = &var[FLAGV*size],*flagw = &var[FLAGW*size],*flagp = &var[FLAGP*size];

  /****************************************************************************
  | Calculate all coefficients
  ****************************************************************************/
  //FOR_EACH_CELL
  if (i>0 && i<imax+1 && j>0 && j<jmax+1 && k>0 && k<kmax+1) {
    dxe =  x[IX(i+1,j,  k)]   -  x[IX(i,  j,  k)];
    dxw =  x[IX(i,  j,  k)]   -  x[IX(i-1,j,  k)];
    dyn =  y[IX(i,  j+1,k)]   -  y[IX(i,  j,  k)];
    dys =  y[IX(i,  j,  k)]   -  y[IX(i,  j-1,k)];
    dzf =  z[IX(i,  j,  k+1)] -  z[IX(i,  j,  k)];
    dzb =  z[IX(i,  j,  k)]   -  z[IX(i,  j,  k-1)];
    Dx  = gx[IX(i,  j,  k)]   - gx[IX(i-1,j,  k)];
    Dy  = gy[IX(i,  j,  k)]   - gy[IX(i,  j-1,k)];
    Dz  = gz[IX(i,  j,  k)]   - gz[IX(i,  j,  k-1)];

    ae[IX(i,j,k)] = Dy*Dz/dxe;
    aw[IX(i,j,k)] = Dy*Dz/dxw;
    an[IX(i,j,k)] = Dx*Dz/dyn;
    as[IX(i,j,k)] = Dx*Dz/dys;
    af[IX(i,j,k)] = Dx*Dy/dzf;
    ab[IX(i,j,k)] = Dx*Dy/dzb;
    b[IX(i,j,k)] = Dx*Dy*Dz/dt*((u[IX(i-1,j,k)]-u[IX(i,j,k)])/Dx
                 + (v[IX(i,j-1,k)]-v[IX(i,j,k)])/Dy
                 + (w[IX(i,j,k-1)]-w[IX(i,j,k)])/Dz);

    
    ap[IX(i,j,k)] = ae[IX(i,j,k)] + aw[IX(i,j,k)] + as[IX(i,j,k)] + an[IX(i,j,k)]
                  + af[IX(i,j,k)] + ab[IX(i,j,k)];
   }// end of if
} //end of project()


///////////////////////////////////////////////////////////////////////////////
/// kenerl of projection
/// 6/17/2015
///////////////////////////////////////////////////////////////////////////////
__kernel void project_velo_corr(__global PARA_DATA_SIMP *para, __global REAL *var, __global int *BINDEX) {  
  int i = get_global_id(0);
  int j = get_global_id(1);
  int k = get_global_id(2);
  int bar;
  int imax = para->geom.imax, jmax = para->geom.jmax;
  int kmax = para->geom.kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int size=(imax+2)*(jmax+2)*(kmax+2);
  REAL dt= para->mytime.dt;
  __global REAL *x = &var[X*size], *y = &var[Y*size], *z = &var[Z*size];
  __global REAL *u = &var[VX*size], *v = &var[VY*size], *w = &var[VZ*size];
  __global REAL *p = &var[IP*size];
  __global REAL *flagu = &var[FLAGU*size],*flagv = &var[FLAGV*size],*flagw = &var[FLAGW*size],*flagp = &var[FLAGP*size];
  
  /****************************************************************************
  | Correct the velocity
  ****************************************************************************/
  //FOR_U_CELL
  if (i>0 && i<imax && j>0 && j<jmax+1 && k>0 && k<kmax+1 && flagu[IX(i,j,k)] < 0) {
    u[IX(i,j,k)] = u[IX(i,j,k)] - dt*(p[IX(i+1,j,k)]-p[IX(i,j,k)]) / (x[IX(i+1,j,k)]-x[IX(i,j,k)]);
  }

  //FOR_V_CELL
  if (i>0 && i<imax+1 && j>0 && j<jmax && k>0 && k<kmax+1 && flagv[IX(i,j,k)] < 0) {
    v[IX(i,j,k)] = v[IX(i,j,k)] - dt*(p[IX(i,j+1,k)]-p[IX(i,j,k)]) / (y[IX(i,j+1,k)]-y[IX(i,j,k)]);
  }

  //FOR_W_CELL
  if (i>0 && i<imax+1 && j>0 && j<jmax+1 && k>0 && k<kmax && flagw[IX(i,j,k)] < 0) {
    w[IX(i,j,k)] = w[IX(i,j,k)] - dt*(p[IX(i,j,k+1)]-p[IX(i,j,k)]) / (z[IX(i,j,k+1)]-z[IX(i,j,k)]);
  }
} // End of project_velo_corr( )

///////////////////////////////////////////////////////////////////////////////
/// kenerl solving linear equation set
/// 6/7/2015
///////////////////////////////////////////////////////////////////////////////
__kernel void Ax_BSolver(__global PARA_DATA_SIMP *para, __global REAL *var,__global int *index, __global int *WF)
{

  int i = get_global_id(0);
  int j = get_global_id(1);
  int k = get_global_id(2);
  int imax = para->geom.imax, jmax= para->geom.jmax;
  int kmax = para->geom.kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int size=(imax+2)*(jmax+2)*(kmax+2);
  
  __global REAL *as = &var[AS*size], *aw = &var[AW*size], *ae = &var[AE*size], *an = &var[AN*size];
  __global REAL *ap = &var[AP*size], *af = &var[AF*size], *ab = &var[AB*size], *b = &var[B*size];
  __global REAL *flag;
  __global float *x;
 
  // set the variables
  //wei tian added
  //6/16/2015
  if (index[0] == VX) {
    flag = &var[FLAGU*size];
		if (WF[0] == 0) {
		    x = &var[TMP1*size]; 
		}
		else{
        x = &var[VX*size];  		
		}
  }// end of if
  else if  (index[0] == VY) {
    flag = &var[FLAGV*size];
		if (WF[0] == 0) {
		    x = &var[TMP2*size]; 
		}
		else{
        x = &var[VY*size];  		
		} 
  }// end of else if
  else if  (index[0] == VZ) {
    flag = &var[FLAGW*size];
		if (WF[0] == 0) {
		    x = &var[TMP3*size]; 
		}
		else{
        x = &var[VZ*size]; 		
		}  
  }// end of else if
  else if  (index[0] == TEMP) {
    flag = &var[FLAGP*size];
		if (WF[0] == 0) {
		    x = &var[TMP1*size];
		}
		else{
        x = &var[TEMP*size]; 	
		}      
  }// end of else if
  
  // jacobian solver
  if (i>0 && i<imax+1 && j>0 && j<jmax+1 && k>0 && k<kmax+1 && flag[IX(i,j,k)] < 0) {
     // var[TMP4] are newly created exclusively for jacobian solver
  var[TMP4*size+IX(i,j,k)] = ( ae[IX(i,j,k)]*x[IX(i+1,j,k)]
             + aw[IX(i,j,k)]*x[IX(i-1,j,k)]
             + an[IX(i,j,k)]*x[IX(i,j+1,k)]
             + as[IX(i,j,k)]*x[IX(i,j-1,k)]
             + af[IX(i,j,k)]*x[IX(i,j,k+1)]
             + ab[IX(i,j,k)]*x[IX(i,j,k-1)]
             + b[IX(i,j,k)] ) / ap[IX(i,j,k)];

  } // end of if
} // end of kernel Ax_BSolver()

///////////////////////////////////////////////////////////////////////////////
/// update results for jacobian solver
/// 29/9/2015
/// wei tian
///////////////////////////////////////////////////////////////////////////////
__kernel void Ax_BSolver_upd(__global PARA_DATA_SIMP *para, __global REAL *var,__global int *index, __global int *WF)
{

  int i = get_global_id(0);
  int j = get_global_id(1);
  int k = get_global_id(2);
  int imax = para->geom.imax, jmax= para->geom.jmax;
  int kmax = para->geom.kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int size=(imax+2)*(jmax+2)*(kmax+2);
  
  __global REAL *as = &var[AS*size], *aw = &var[AW*size], *ae = &var[AE*size], *an = &var[AN*size];
  __global REAL *ap = &var[AP*size], *af = &var[AF*size], *ab = &var[AB*size], *b = &var[B*size];
  __global REAL *flag;
  __global float *x;
 
  // set the variables
  //wei tian added
  //6/16/2015
  if (index[0] == VX) {
    flag = &var[FLAGU*size];
		if (WF[0] == 0) {
		    x = &var[TMP1*size]; 
		}
		else{
        x = &var[VX*size];  		
		}
  }// end of if
  else if  (index[0] == VY) {
    flag = &var[FLAGV*size];
		if (WF[0] == 0) {
		    x = &var[TMP2*size]; 
		}
		else{
        x = &var[VY*size];  		
		} 
  }// end of else if
  else if  (index[0] == VZ) {
    flag = &var[FLAGW*size];
		if (WF[0] == 0) {
		    x = &var[TMP3*size]; 
		}
		else{
        x = &var[VZ*size]; 		
		}  
  }// end of else if
  else if  (index[0] == TEMP) {
    flag = &var[FLAGP*size];
		if (WF[0] == 0) {
		    x = &var[TMP1*size];
		}
		else{
        x = &var[TEMP*size]; 	
		}      
  }// end of else if
  
  // jacobian solver
  if (i>0 && i<imax+1 && j>0 && j<jmax+1 && k>0 && k<kmax+1 && flag[IX(i,j,k)] < 0) {
     // var[TMP4] are newly created exclusively for jacobian solver
   x[IX(i,j,k)] =  var[TMP4*size+IX(i,j,k)];
  } // end of if
} // end of kernel Ax_BSolver_upd()

///////////////////////////////////////////////////////////////////////////////
/// kenerl solving linear equation set especially for pressure
/// 6/17/2015
///////////////////////////////////////////////////////////////////////////////
__kernel void Ax_BSolver_P(__global PARA_DATA_SIMP *para, __global REAL *var)
{
  int i = get_global_id(0);
  int j = get_global_id(1);
  int k = get_global_id(2);
  
  //printf("Ax_Solver: now solving [%d,%d,%d]\n",i, j, k);

  int imax = para->geom.imax, jmax= para->geom.jmax;
  int kmax = para->geom.kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int size=(imax+2)*(jmax+2)*(kmax+2);
  int it=0;
  
  __global REAL *as = &var[AS*size], *aw = &var[AW*size], *ae = &var[AE*size], *an = &var[AN*size];
  __global REAL *ap = &var[AP*size], *af = &var[AF*size], *ab = &var[AB*size], *b = &var[B*size];
  __global REAL *flag = &var[FLAGP*size];
  __global float *x = &var[IP*size];
 
  // jacobian solver
  if (i>0 && i<imax+1 && j>0 && j<jmax+1 && k>0 && k<kmax+1 && flag[IX(i,j,k)] < 0) {
  //printf("%d: the cell calc is [%d, %d, %d]\n",it, i, j,k);
  var[TMP4*size+IX(i,j,k)] = ( ae[IX(i,j,k)]*x[IX(i+1,j,k)]
             + aw[IX(i,j,k)]*x[IX(i-1,j,k)]
             + an[IX(i,j,k)]*x[IX(i,j+1,k)]
             + as[IX(i,j,k)]*x[IX(i,j-1,k)]
             + af[IX(i,j,k)]*x[IX(i,j,k+1)]
             + ab[IX(i,j,k)]*x[IX(i,j,k-1)]
             + b[IX(i,j,k)] ) / ap[IX(i,j,k)];
  } // end of if
} // end of kernel Ax_BSolver_P()

///////////////////////////////////////////////////////////////////////////////
/// kenerl solving linear equation set especially for pressure
/// 6/17/2015
///////////////////////////////////////////////////////////////////////////////
__kernel void Ax_BSolver_P_upd(__global PARA_DATA_SIMP *para, __global REAL *var)
{
  int i = get_global_id(0);
  int j = get_global_id(1);
  int k = get_global_id(2);
  
  //printf("Ax_Solver: now solving [%d,%d,%d]\n",i, j, k);

  int imax = para->geom.imax, jmax= para->geom.jmax;
  int kmax = para->geom.kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int size=(imax+2)*(jmax+2)*(kmax+2);
  int it=0;
  
  __global REAL *as = &var[AS*size], *aw = &var[AW*size], *ae = &var[AE*size], *an = &var[AN*size];
  __global REAL *ap = &var[AP*size], *af = &var[AF*size], *ab = &var[AB*size], *b = &var[B*size];
  __global REAL *flag = &var[FLAGP*size];
  __global float *x = &var[IP*size];
 
  // jacobian solver
  if (i>0 && i<imax+1 && j>0 && j<jmax+1 && k>0 && k<kmax+1 && flag[IX(i,j,k)] < 0) {
     x[IX(i,j,k)] = var[TMP4*size+IX(i,j,k)];
  } // end of if
} // end of kernel Ax_BSolver_P_upd()

///////////////////////////////////////////////////////////////////////////////
/// set boundary conditions for temperature
/// 6/7/2015
/// Wei Tian
///////////////////////////////////////////////////////////////////////////////
__kernel void set_bnd_T(__global PARA_DATA_SIMP *para, __global REAL *var, __global int *WF, __global int *BINDEX) {
    int it = get_global_id(0);
  int i, j, k;
  int imax = para->geom.imax, jmax = para->geom.jmax;
  int kmax = para->geom.kmax;
  int size = (imax + 2)*(jmax + 2)*(kmax + 2);
  int IMAX = imax + 2, IJMAX = (imax + 2)*(jmax + 2);
  __global REAL *aw = &var[AW*size], *ae = &var[AE*size], *as = &var[AS*size], *an = &var[AN*size];
  __global REAL *af = &var[AF*size], *ab = &var[AB*size], *b = &var[B*size],
    *qflux = &var[QFLUX*size], *qfluxbc = &var[QFLUXBC*size];
  __global REAL *flagp = &var[FLAGP*size];
  __global REAL *psi;
  REAL axy, ayz, azx; 
  REAL h ; // fixme, other flow needs to be remodelled
  REAL rhoCp_1 = 1 / (para->prob.rho * para->prob.Cp);
  REAL D;
  int tmp1=0;
  
  // get the index
  i = BINDEX[0*size+it];
    j = BINDEX[1*size+it];
    k = BINDEX[2*size+it]; 

  // if Advection
  if (WF[0] == 0) {
    psi = &var[TMP1*size];
  }
  // else if diffusion
  else if (WF[0] == 1) {
    psi = &var[TEMP*size];  
  }
    
  axy = area_xy(para, var, &i, &j, &k);
  ayz = area_yz(para, var, &i, &j, &k);
  azx = area_zx(para, var, &i, &j, &k);

  /*-------------------------------------------------------------------------
  | Inlet boundary
  | 0: Inlet, -1: Fluid,  1: Solid Wall or Block, 2: Outlet
  -------------------------------------------------------------------------*/
  if (flagp[IX(i, j, k)] == INLET || flagp[IX(i,j,k)]==RACK_OUTLET) psi[IX(i, j, k)] = var[TEMPBC*size + IX(i, j, k)];

  /*-------------------------------------------------------------------------
  | Solid wall or block
  -------------------------------------------------------------------------*/
  if (flagp[IX(i, j, k)] == SOLID) {
    /*......................................................................
    | Constant temperature
    ......................................................................*/
    if (BINDEX[3*size + it] == 1) {
      psi[IX(i, j, k)] = var[TEMPBC*size + IX(i, j, k)];

      // West boundary wall and eastern neighbor cell is fluid
      if (i == 0) {
        if (flagp[IX(i + 1, j, k)] == FLUID) {
            tmp1 = i + 1;
          D = 0.5f * length_x(para, var, &tmp1, &j, &k);
          h = para->prob.Cp * para->prob.rho * para->prob.alpha / D;
          h = para->prob.Cp * para->prob.rho *para->prob.coeff_h;
          aw[IX(i + 1, j, k)] = h * rhoCp_1 * ayz;
          qflux[IX(i, j, k)] = h * (psi[IX(i + 1, j, k)] - psi[IX(i, j, k)]);
        }
      } // End of if(i==0)
      // East boundary wall and western neighbor cell is fluid
      else if (i == imax + 1) {
        if (flagp[IX(i - 1, j, k)] == FLUID) {
            tmp1 = i -1;
          D = 0.5f * length_x(para, var, &tmp1, &j, &k);
          h = para->prob.Cp * para->prob.rho * para->prob.alpha / D;
          h = para->prob.Cp * para->prob.rho *para->prob.coeff_h;
          ae[IX(i - 1, j, k)] = h * rhoCp_1 * ayz;
          qflux[IX(i, j, k)] = h * (psi[IX(i - 1, j, k)] - psi[IX(i, j, k)]);
        }
      } // End of else if(i==imax+1)
      // Between West and East
      else {
        // Eastern neighbor cell is fluid
        if (flagp[IX(i + 1, j, k)] == FLUID) {
            tmp1 = i + 1;
          D = 0.5f * length_x(para, var, &tmp1, &j, &k);
          h = para->prob.Cp * para->prob.rho * para->prob.alpha / D;
          h = para->prob.Cp * para->prob.rho *para->prob.coeff_h;
          aw[IX(i + 1, j, k)] = h * rhoCp_1 * ayz;
          qflux[IX(i, j, k)] = h * (psi[IX(i + 1, j, k)] - psi[IX(i, j, k)]);
        }
        // Western neighbor cell is fluid
        if (flagp[IX(i - 1, j, k)] == FLUID) {
            tmp1 = i + 1;
          D = 0.5f * length_x(para, var, &tmp1, &j, &k);
          h = para->prob.Cp * para->prob.rho * para->prob.alpha / D;
          h = para->prob.Cp * para->prob.rho *para->prob.coeff_h;
          ae[IX(i - 1, j, k)] = h * rhoCp_1 * ayz;
          qflux[IX(i, j, k)] = h * (psi[IX(i - 1, j, k)] - psi[IX(i, j, k)]);
        }
      } // End of 0<i<imax+1
      // South wall boundary and northern neighbor is fluid
      if (j == 0) {
        if (flagp[IX(i, j + 1, k)] == FLUID) {
            tmp1 = j + 1;
          D = 0.5f * length_y(para, var, &i, &tmp1, &k);
          h = para->prob.Cp * para->prob.rho * para->prob.alpha / D;
          h = para->prob.Cp * para->prob.rho *para->prob.coeff_h;
          as[IX(i, j + 1, k)] = h * rhoCp_1 * azx;
          qflux[IX(i, j, k)] = h * (psi[IX(i, j + 1, k)] - psi[IX(i, j, k)]);
        }
      }
      // North wall boundary and southern neighbor is fluid
      else if (j == jmax + 1) {
        if (flagp[IX(i, j - 1, k)] == FLUID) {
            tmp1 = j - 1;
          D = 0.5f * length_y(para, var, &i, &tmp1, &k);
          h = para->prob.Cp * para->prob.rho * para->prob.alpha / D;
          h = para->prob.Cp * para->prob.rho *para->prob.coeff_h;
          an[IX(i, j - 1, k)] = h * rhoCp_1 * azx;
          qflux[IX(i, j, k)] = h * (psi[IX(i, j - 1, k)] - psi[IX(i, j, k)]);
        }
      }
      // Between South and North
      else {
        // Southern neighbor is fluid
        if (flagp[IX(i, j - 1, k)] == FLUID) {
            tmp1 = j - 1;
          D = 0.5f * length_y(para, var, &i, &tmp1, &k);
          h = para->prob.Cp * para->prob.rho * para->prob.alpha / D;
          h = para->prob.Cp * para->prob.rho *para->prob.coeff_h;
          an[IX(i, j - 1, k)] = h * rhoCp_1 * azx;
          qflux[IX(i, j, k)] = h * (psi[IX(i, j - 1, k)] - psi[IX(i, j, k)]);
        }
        // Northern neighbor is fluid
        if (flagp[IX(i, j + 1, k)] == FLUID) {
            tmp1 = j + 1;
          D = 0.5f * length_y(para, var, &i, &tmp1, &k);
          h = para->prob.Cp * para->prob.rho * para->prob.alpha / D;
          h = para->prob.Cp * para->prob.rho *para->prob.coeff_h;
          as[IX(i, j + 1, k)] = h * rhoCp_1 * azx;
          qflux[IX(i, j, k)] = h * (psi[IX(i, j + 1, k)] - psi[IX(i, j, k)]);
        }
      }
      // Floor and ceiling neighbor is fluid
      if (k == 0) {
        if (flagp[IX(i, j, k + 1)] == FLUID) {
            tmp1 = k + 1;
          D = 0.5f * length_z(para, var, &i, &j, &tmp1);
          h = para->prob.Cp * para->prob.rho * para->prob.alpha / D;
          h = para->prob.Cp * para->prob.rho *para->prob.coeff_h;
          ab[IX(i, j, k + 1)] = h * rhoCp_1 * axy;
          qflux[IX(i, j, k)] = h * (psi[IX(i, j, k + 1)] - psi[IX(i, j, k)]);
        }
      }
      // Ceiling and floor neighbor is fluid
      else if (k == kmax + 1) {
        if (flagp[IX(i, j, k - 1)] == FLUID) {
            tmp1 = k + 1;
          D = 0.5f * length_z(para, var, &i, &j, &tmp1);
          h = para->prob.Cp * para->prob.rho * para->prob.alpha / D;
          h = para->prob.Cp * para->prob.rho *para->prob.coeff_h;
          af[IX(i, j, k - 1)] = h * rhoCp_1 * axy;
          qflux[IX(i, j, k)] = h * (psi[IX(i, j, k - 1)] - psi[IX(i, j, k)]);
        }
      }
      // Between Floor and Ceiling
      else {
        // Ceiling neighbor is fluid
        if (flagp[IX(i, j, k + 1)] == FLUID) {
            tmp1 = k + 1;
          D = 0.5f * length_z(para, var, &i, &j, &tmp1);
          h = para->prob.Cp * para->prob.rho * para->prob.alpha / D;
          h = para->prob.Cp * para->prob.rho *para->prob.coeff_h;
          ab[IX(i, j, k + 1)] = h * rhoCp_1 * axy;
          qflux[IX(i, j, k)] = h * (psi[IX(i, j, k + 1)] - psi[IX(i, j, k)]);
        }
        // Floor neighbor is fluid
        if (flagp[IX(i, j, k - 1)] == FLUID) {
            tmp1 = k - 1;
          D = 0.5f * length_z(para, var, &i, &j, &tmp1);
          h = para->prob.Cp * para->prob.rho * para->prob.alpha / D;
          h = para->prob.Cp * para->prob.rho *para->prob.coeff_h;
          af[IX(i, j, k - 1)] = h * rhoCp_1 * axy;
          qflux[IX(i, j, k)] = h * (psi[IX(i, j, k - 1)] - psi[IX(i, j, k)]);
        }
      }
    } // End of constant temperature wall
    /*.......................................................................
    | Constant heat flux
    .......................................................................*/
    if (BINDEX[3*size + it] == 0) {
      // West wall boundary and eastern neighbor is fluid
      if (i == 0) {
        if (flagp[IX(i + 1, j, k)] == FLUID) {
          aw[IX(i + 1, j, k)] = 0;
          tmp1 = i + 1;
          D = 0.5f * length_z(para, var, &tmp1, &j, &k);
          h = para->prob.Cp * para->prob.rho * para->prob.alpha / D;
          h = para->prob.Cp * para->prob.rho *para->prob.coeff_h;
          b[IX(i + 1, j, k)] += rhoCp_1 * qfluxbc[IX(i, j, k)] * ayz;
          // get the temperature of solid surface
          psi[IX(i, j, k)] = qfluxbc[IX(i, j, k)] / h + psi[IX(i + 1, j, k)];
        }
      } // End of if(i==0)
      // East wall boundary and western neighbor is fluid
      else if (i == imax + 1) {
        if (flagp[IX(i - 1, j, k)] == FLUID) {
          ae[IX(i - 1, j, k)] = 0;
          tmp1 = i - 1;
          D = 0.5f * length_z(para, var, &tmp1, &j, &k);
          h = para->prob.Cp * para->prob.rho * para->prob.alpha / D;
          h = para->prob.Cp * para->prob.rho *para->prob.coeff_h;
          b[IX(i - 1, j, k)] += rhoCp_1 * qfluxbc[IX(i, j, k)] * ayz;
          // get the temperature of solid surface
          psi[IX(i, j, k)] = qfluxbc[IX(i, j, k)] / h + psi[IX(i - 1, j, k)];
        }
      } // End of else if(i==imax+1)
      // Between West and East
      else {
        // Eastern neighbor is fluid
        if (flagp[IX(i + 1, j, k)] == FLUID) {
          aw[IX(i + 1, j, k)] = 0;
          tmp1 = i + 1;
          D = 0.5f * length_z(para, var, &tmp1, &j, &k);
          h = para->prob.Cp * para->prob.rho * para->prob.alpha / D;
          h = para->prob.Cp * para->prob.rho *para->prob.coeff_h;
          b[IX(i + 1, j, k)] += rhoCp_1 * qfluxbc[IX(i, j, k)] * ayz;
          // get the temperature of solid surface
          psi[IX(i, j, k)] = qfluxbc[IX(i, j, k)] / h + psi[IX(i + 1, j, k)];
        }
        // Western neighbor is fluid
        if (flagp[IX(i - 1, j, k)] == FLUID) {
          ae[IX(i - 1, j, k)] = 0;
          tmp1 = i - 1;
          D = 0.5f * length_z(para, var, &tmp1, &j, &k);
          h = para->prob.Cp * para->prob.rho * para->prob.alpha / D;
          h = para->prob.Cp * para->prob.rho *para->prob.coeff_h;
          b[IX(i - 1, j, k)] += rhoCp_1 * qfluxbc[IX(i, j, k)] * ayz;
          // get the temperature of solid surface
          psi[IX(i, j, k)] = qfluxbc[IX(i, j, k)] / h + psi[IX(i - 1, j, k)];
        }
      } // End of else
      // South wall boundary and northern neighbor is fluid
      if (j == 0) {
        if (flagp[IX(i, j + 1, k)] == FLUID) {
          as[IX(i, j + 1, k)] = 0;
          tmp1 = j + 1;
          D = 0.5f * length_z(para, var, &i, &tmp1, &k);
          h = para->prob.Cp * para->prob.rho * para->prob.alpha / D;
          h = para->prob.Cp * para->prob.rho *para->prob.coeff_h;
          b[IX(i, j + 1, k)] += rhoCp_1 * qfluxbc[IX(i, j, k)] * azx;
          // get the temperature of solid surface
          psi[IX(i, j, k)] = qfluxbc[IX(i, j, k)] / h + psi[IX(i, j + 1, k)];
        }
      }
      // North wall boundary and southern neighbor is fluid
      else if (j == jmax + 1) {
        if (flagp[IX(i, j - 1, k)] == FLUID) {
          an[IX(i, j - 1, k)] = 0;
          tmp1 = j - 1;
          D = 0.5f * length_z(para, var, &i, &tmp1, &k);
          h = para->prob.Cp * para->prob.rho * para->prob.alpha / D;
          h = para->prob.Cp * para->prob.rho *para->prob.coeff_h;
          b[IX(i, j - 1, k)] += rhoCp_1 * qfluxbc[IX(i, j, k)] * azx;
          // get the temperature of solid surface
          psi[IX(i, j, k)] = qfluxbc[IX(i, j, k)] / h + psi[IX(i, j - 1, k)];
        }
      }
      // Between South and North
      else {
        // Southern neighbor is fluid
        if (flagp[IX(i, j - 1, k)] == FLUID) {
          an[IX(i, j - 1, k)] = 0;
          tmp1 = j - 1;
          D = 0.5f * length_z(para, var,  &i, &tmp1, &k);
          h = para->prob.Cp * para->prob.rho * para->prob.alpha / D;
          h = para->prob.Cp * para->prob.rho *para->prob.coeff_h;
          b[IX(i, j - 1, k)] += rhoCp_1 * qfluxbc[IX(i, j, k)] * azx;
          // get the temperature of solid surface
          psi[IX(i, j, k)] = qfluxbc[IX(i, j, k)] / h + psi[IX(i, j - 1, k)];
        }
        // Northern neighbor is fluid
        if (flagp[IX(i, j + 1, k)] == FLUID) {
          as[IX(i, j + 1, k)] = 0;
          tmp1 = j + 1;
          D = 0.5f * length_z(para, var,  &i, &tmp1, &k);
          h = para->prob.Cp * para->prob.rho * para->prob.alpha / D;
          h = para->prob.Cp * para->prob.rho *para->prob.coeff_h;
          b[IX(i, j + 1, k)] += rhoCp_1 * qfluxbc[IX(i, j, k)] * azx;
          // get the temperature of solid surface
          psi[IX(i, j, k)] = qfluxbc[IX(i, j, k)] / h + psi[IX(i, j + 1, k)];
        }
      }
      // Floor boundary and ceiling neighbor is fluid
      if (k == 0) {
        if (flagp[IX(i, j, k + 1)] == FLUID) {
          ab[IX(i, j, k + 1)] = 0;
          tmp1 = k + 1;
          D = 0.5f * length_z(para, var, &i, &j, &tmp1);
          h = para->prob.Cp * para->prob.rho * para->prob.alpha / D;
          h = para->prob.Cp * para->prob.rho *para->prob.coeff_h;
          b[IX(i, j, k + 1)] += rhoCp_1 * qfluxbc[IX(i, j, k)] * axy;
          // Get the temperature on the solid surface
          psi[IX(i, j, k)] = qfluxbc[IX(i, j, k)] / h + psi[IX(i, j, k + 1)];
        }
      }
      // Ceiling boundary and floor neighbor is fluid
      else if (k == kmax + 1) {
        if (flagp[IX(i, j, k - 1)] == FLUID) {
          af[IX(i, j, k - 1)] = 0;
          tmp1 = k - 1;
          D = 0.5f * length_z(para, var, &i, &j, &tmp1);
          h = para->prob.Cp * para->prob.rho * para->prob.alpha / D;
          h = para->prob.Cp * para->prob.rho *para->prob.coeff_h;
          b[IX(i, j, k - 1)] += rhoCp_1 * qfluxbc[IX(i, j, k)] * axy;
          // Get the temperature on the solid surface
          psi[IX(i, j, k)] = qfluxbc[IX(i, j, k)] / h + psi[IX(i, j, k - 1)];
        }
      }
      // Between Floor and Ceiling
      else {
        // Ceiling neighbor is fluid
        if (flagp[IX(i, j, k + 1)] == FLUID) {
          ab[IX(i, j, k + 1)] = 0;
          tmp1 = k + 1;
          D = 0.5f * length_z(para, var, &i, &j, &tmp1);
          h = para->prob.Cp * para->prob.rho * para->prob.alpha / D;
          h = para->prob.Cp * para->prob.rho *para->prob.coeff_h;
          b[IX(i, j, k + 1)] += rhoCp_1 * qfluxbc[IX(i, j, k)] * axy;
          // Get the temperature on the solid surface
          psi[IX(i, j, k)] = qfluxbc[IX(i, j, k)] / h + psi[IX(i, j, k + 1)];
        }
        // Floor neighbor is fluid
        if (flagp[IX(i, j, k - 1)] == FLUID) {
          af[IX(i, j, k - 1)] = 0;
          tmp1 = k - 1;
          D = 0.5f * length_z(para, var, &i, &j, &tmp1);
          h = para->prob.Cp * para->prob.rho * para->prob.alpha / D;
          h = para->prob.Cp * para->prob.rho *para->prob.coeff_h;
          b[IX(i, j, k - 1)] += rhoCp_1 * qfluxbc[IX(i, j, k)] * axy;
          // Get the temperature on the solid surface
          psi[IX(i, j, k)] = qfluxbc[IX(i, j, k)] / h + psi[IX(i, j, k - 1)];
        }
      }
    } // End of constant heat flux
  } // End of wall boundary

  /*-------------------------------------------------------------------------
  | Outlet boundary
  -------------------------------------------------------------------------*/
  if (flagp[IX(i, j, k)] == OUTLET) {
    // West
    if (i == 0) {
      aw[IX(i + 1, j, k)] = 0;
      psi[IX(i, j, k)] = psi[IX(i + 1, j, k)];
    }
    // North
    if (i == imax + 1) {
      ae[IX(i - 1, j, k)] = 0;
      psi[IX(i, j, k)] = psi[IX(i - 1, j, k)];
    }
    // South
    if (j == 0) {
      as[IX(i, j + 1, k)] = 0;
      psi[IX(i, j, k)] = psi[IX(i, j + 1, k)];
    }
    // North
    if (j == jmax + 1) {
      an[IX(i, j - 1, k)] = 0;
      psi[IX(i, j, k)] = psi[IX(i, j - 1, k)];
    }
    // Floor
    if (k == 0) {
      ab[IX(i, j, k + 1)] = 0;
      psi[IX(i, j, k)] = psi[IX(i, j, k + 1)];
    }
    // Ceiling
    if (k == kmax + 1) {
      af[IX(i, j, k - 1)] = 0;
      psi[IX(i, j, k)] = psi[IX(i, j, k - 1)];
    }
  } // End of boundary for outlet

} //end of set_bnd_T()

///////////////////////////////////////////////////////////////////////////////
/// set boundary conditions for advection of temperature, solved by implicit scheme
/// 8/8/2017
/// Wei Tian
///////////////////////////////////////////////////////////////////////////////
__kernel void set_bnd_T_im(__global PARA_DATA_SIMP *para, __global REAL *var, __global int *BINDEX) {
  int it = get_global_id(0);
  int i, j, k;
  int imax = para->geom.imax, jmax = para->geom.jmax;
  int kmax = para->geom.kmax;
  int size = (imax + 2)*(jmax + 2)*(kmax + 2);
  int IMAX = imax + 2, IJMAX = (imax + 2)*(jmax + 2);
  __global REAL *aw = &var[AW*size], *ae = &var[AE*size], *as = &var[AS*size], *an = &var[AN*size];
  __global REAL *af = &var[AF*size], *ab = &var[AB*size], *b = &var[B*size];
  __global REAL *flagp = &var[FLAGP*size];
  __global REAL *psi;
  
  // get the index
  i = BINDEX[0*size+it];
  j = BINDEX[1*size+it];
  k = BINDEX[2*size+it]; 

  // if Advection, psi is assigned the pointer of tmp1
  psi = &var[TMP1*size];
  
	// start assigning the boundary conditions
	// : inlet
	if (flagp[IX(i, j, k)] == INLET || flagp[IX(i,j,k)]==RACK_OUTLET) psi[IX(i, j, k)] = var[TEMPBC*size + IX(i, j, k)];
	
	// : wall
	if (flagp[IX(i, j, k)] == SOLID) {
		//: temperature as boundary condition type, then give the BC set by user
		if (BINDEX[3*size + it] == 0) {
		  psi[IX(i, j, k)] = var[TEMPBC*size + IX(i, j, k)];
		}
		//: heat flux as boundary condition type, then give the temperature at wall from last step
		else {
		  psi[IX(i, j, k)] = var[TEMP*size + IX(i, j, k)];
		}
  } // end of solid wall

	// : outlet
	if (flagp[IX(i, j, k)] == OUTLET || flagp[IX(i, j, k)] == TILE) {
	  // :west
		if (i == 0) {
		  psi[IX(i, j, k)] = psi[IX(i + 1, j, k)];
		} // end of i==0
		// :east
		if (i == imax + 1) {
		  psi[IX(i, j, k)] = psi[IX(i - 1, j, k)];
		} // end of i == imax+1
		// :south
		if (j == 0) {
		  psi[IX(i, j, k)] = psi[IX(i, j + 1, k)];
		} // end of j==0
		// :north
		if (j == jmax + 1) {
		  psi[IX(i, j, k)] = psi[IX(i, j - 1, k)];
		} // end of j == jmax+1
		// :floor
		if (k == 0) {
		  psi[IX(i, j, k)] = psi[IX(i, j, k + 1)];
		} // end of k ==0
		//: ceiling
		if (k == kmax + 1) {
		  psi[IX(i, j, k)] = psi[IX(i, j, k - 1)];
		} // end of k ==kmax+1
	} // end of outlet

} //end of set_bnd_T_im

///////////////////////////////////////////////////////////////////////////////
/// set boundary conditions for pressure
/// 6/19/2015
/// Wei Tian
///////////////////////////////////////////////////////////////////////////////

__kernel void set_bnd_pressure(__global PARA_DATA_SIMP *para, __global REAL *var,__global int *BINDEX) {
    int it = get_global_id(0);
  int i, j, k;
  int imax = para->geom.imax, jmax = para->geom.jmax;
  int kmax = para->geom.kmax;
  int size = (imax + 2)*(jmax + 2)*(kmax + 2);
  int IMAX = imax + 2, IJMAX = (imax + 2)*(jmax + 2);
  __global REAL *aw = &var[AW*size], *ae = &var[AE*size], *as = &var[AS*size], *an = &var[AN*size];
  __global REAL *af = &var[AF*size], *ab = &var[AB*size];
  __global REAL *flagp = &var[FLAGP*size];
  __global REAL *p = &var[IP*size];

    // get the index
  i = BINDEX[0*size+it];
    j = BINDEX[1*size+it];
    k = BINDEX[2*size+it]; 

  /*-------------------------------------------------------------------------
    | For X direction
    -------------------------------------------------------------------------*/
    if(i>0) {
      if(flagp[IX(i-1,j,k)]<0) {
        p[IX(i,j,k)] = p[IX(i-1,j,k)];
        ae[IX(i-1,j,k)] = 0;
      }
    }
    if(i<imax+1) {
      if(flagp[IX(i+1,j,k )]<0) {
        p[IX(i,j,k)] = p[IX(i+1,j,k)];
        aw[IX(i+1,j,k)] = 0;
      }
    }
    /*-------------------------------------------------------------------------
    | For Y direction
    -------------------------------------------------------------------------*/
    if(j>0) {
      if(flagp[IX(i,j-1,k)]<0) {
        p[IX(i,j,k)] = p[IX(i,j-1,k )];
        an[IX(i,j-1,k)] = 0;
      }
    }
    if(j<jmax+1) {
      if(flagp[IX(i,j+1,k)]<0) {
        p[IX(i,j,k )] = p[IX(i,j+1,k )];
        as[IX(i,j+1,k)] = 0;
      }
    }
    /*-------------------------------------------------------------------------
    | For Z direction
    -------------------------------------------------------------------------*/
    if(k>0) {
      if(flagp[IX(i,j,k-1)]<0) {
        p[IX(i,j,k)] = p[IX(i,j,k-1)];
        af[IX(i,j,k-1)] = 0;
      }
    }
    if(k<kmax+1) {
      if(flagp[IX(i,j,k+1 )]<0) {
        p[IX(i,j,k)] = p[IX(i,j,k+1)];
        ab[IX(i,j,k+1)] = 0;
      }
    }
} //end of set_bnd_pressure()

///////////////////////////////////////////////////////////////////////////////
/// set boundary conditions for VX
/// 6/12/2015
/// Wei Tian
///////////////////////////////////////////////////////////////////////////////
__kernel void set_bnd_VX(__global PARA_DATA_SIMP *para, __global REAL *var, __global int *WF, __global int *BINDEX) {
    int it = get_global_id(0);
  int i, j, k;
  int imax = para->geom.imax, jmax = para->geom.jmax;
  int kmax = para->geom.kmax;
  int size = (imax + 2)*(jmax + 2)*(kmax + 2);
  int IMAX = imax + 2, IJMAX = (imax + 2)*(jmax + 2);
  __global REAL *aw = &var[AW*size], *ae = &var[AE*size], *as = &var[AS*size], *an = &var[AN*size];
  __global REAL *af = &var[AF*size], *ab = &var[AB*size], *b = &var[B*size];
  __global REAL *flagp = &var[FLAGP*size];
  __global REAL *psi;
  
  // get the index
  i = BINDEX[0*size+it];
    j = BINDEX[1*size+it];
    k = BINDEX[2*size+it]; 


  // if Advection
  if (WF[0] == 0) {
    psi = &var[TMP1*size];
  }
  // else if diffusion
  else if (WF[0] == 1) {
    psi = &var[VX*size];  
  }

  //printf("the FUNC is %d\n", WF[0]);
      // Inlet
    if (flagp[IX(i, j, k)] == 0) {
      psi[IX(i, j, k)] = var[VXBC*size+IX(i, j, k)];
      if (i != 0) psi[IX(i - 1, j, k)] = var[VXBC*size+IX(i, j, k)];
      }
    // Rack outlet
    if(flagp[IX(i,j,k)]==RACK_OUTLET) {
      psi[IX(i, j, k)] = var[VXBC*size+IX(i, j, k)];
      }	  
    // Solid wall
    // non-flip wall boundary conditions for walls
    // fixme: for lid-driven-cavity case the B.C. needs revised
    // wei tian noted on 6/12/2015
    if (flagp[IX(i, j, k)] == 1) {
        if (k==kmax+1) {
        psi[IX(i, j, k)] = 0.0;
      }
      else {
        psi[IX(i, j, k)] = 0;
        if (i != 0) psi[IX(i - 1, j, k)] = 0;
      }
    }
    // Outlet
    // Neumann boudnary conditions for outlets
    // wei tian noted on 6/12/2015
    if (flagp[IX(i, j, k)] == 2 || flagp[IX(i, j, k)] == TILE) {
			// if solving tile using presure correction method
			if (flagp[IX(i, j, k)] == TILE && para->solv.tile_flow_correct == PRESSURE_BASE) {
				// West
				if (i == 0) {
					psi[IX(i, j, k)] = var[TILE_FLOW_BC*size+IX(i, j, k)];
				}
				// East
				if (i == imax + 1) {
					psi[IX(i - 1, j, k)] = var[TILE_FLOW_BC*size+IX(i, j, k)];
				}
			}
			else {
				// West
				if (i == 0) {
					psi[IX(i, j, k)] = psi[IX(i + 1, j, k)];
					aw[IX(i + 1, j, k)] = 0;
				}
				// East
				if (i == imax + 1) {
					psi[IX(i - 1, j, k)] = psi[IX(i - 2, j, k)];
					ae[IX(i - 2, j, k)] = 0;
				}
			}
      // South
      if (j == 0) as[IX(i, j + 1, k)] = 0;
      // North
      if (j == jmax + 1) an[IX(i, j - 1, k)] = 0;
      // Floor
      if (k == 0) ab[IX(i, j, k + 1)] = 0;
      // Ceiling
      if (k == kmax + 1) af[IX(i, j, k - 1)] = 0;
    }
        // Rack inlet, which is outlet of DC room
        // Apply a fixed velocity
        if(flagp[IX(i,j,k)]==RACK_INLET){
        	// if west side is fluid
        	if (flagp[IX(i-1,j,k)]==-1) {
        		psi[IX(i-1,j,k)] = var[VXBC*size+IX(i,j,k)];
        	}
        	// if east side is fluid
        	if (flagp[IX(i+1,j,k)]==-1) {
        		psi[IX(i+1,j,k)] = var[VXBC*size+IX(i,j,k)];
        	}
        	// if south side is fluid
        	if (flagp[IX(i,j-1,k)]==-1) {
        		psi[IX(i,j-1,k)] = 0;
        	}
        	// if north side is fluid
        	if (flagp[IX(i,j+1,k)]==-1) {
        		psi[IX(i,j+1,k)] = 0;
        	}
        	// if back side is fluid
        	if (flagp[IX(i,j,k-1)]==-1) {
        		psi[IX(i,j,k-1)] = 0;
        	}
        	// if front side is fluid
        	if (flagp[IX(i,j,k+1)]==-1) {
        		psi[IX(i,j,k+1)] = 0;
        	}
        }// end of setting rack inlet	
  
}//end of set_bnd_VX()

///////////////////////////////////////////////////////////////////////////////
/// set boundary conditions for VX for advection solved by implicit scheme
/// 8/8/2017
/// Wei Tian
///////////////////////////////////////////////////////////////////////////////
__kernel void set_bnd_VX_im(__global PARA_DATA_SIMP *para, __global REAL *var, __global int *BINDEX) {
    int it = get_global_id(0);
  int i, j, k;
  int imax = para->geom.imax, jmax = para->geom.jmax;
  int kmax = para->geom.kmax;
  int size = (imax + 2)*(jmax + 2)*(kmax + 2);
  int IMAX = imax + 2, IJMAX = (imax + 2)*(jmax + 2);
  __global REAL *aw = &var[AW*size], *ae = &var[AE*size], *as = &var[AS*size], *an = &var[AN*size];
  __global REAL *af = &var[AF*size], *ab = &var[AB*size], *b = &var[B*size];
  __global REAL *flagp = &var[FLAGP*size];
  __global REAL *psi;
  int put_X = para->geom.tile_putX, put_Y =para->geom.tile_putY, put_Z = para->geom.tile_putZ;
  
  // get the index
  i = BINDEX[0*size+it];
  j = BINDEX[1*size+it];
  k = BINDEX[2*size+it];
	
	// if advection, psi is tmp1
	psi = &var[TMP1*size];
	
	// start assigning the boundary condition
	// :Inlet
  if (flagp[IX(i, j, k)] == INLET) {
	  psi[IX(i, j, k)] = var[VXBC*size+IX(i, j, k)];
    if (i != 0) psi[IX(i - 1, j, k)] = var[VXBC*size+IX(i, j, k)];
	} //end of inlet

	// :Rack outlet
  if (flagp[IX(i, j, k)] == RACK_OUTLET) {
	  psi[IX(i, j, k)] = var[VXBC*size+IX(i, j, k)];
    if (i != 0) psi[IX(i - 1, j, k)] = var[VXBC*size+IX(i, j, k)];
	} //end of rack outlet	
	
	// :Solid
	if (flagp[IX(i, j, k)] == 1) {
    psi[IX(i, j, k)] = 0;
    if (i != 0) psi[IX(i - 1, j, k)] = 0;
	} // end of solid

			// Tile
			if (flagp[IX(i, j, k)] == TILE && para->solv.tile_flow_correct == PRESSURE_BASE) {
				if (i == 0 ) {
					psi[IX(i, j, k)] = var[TILE_FLOW_BC*size+IX(i, j, k)];
				}
				if (i == imax+1) {
					psi[IX(i - 1, j, k)] = var[TILE_FLOW_BC*size+IX(i, j, k)];
				}
			}
			// Model the room and plenum together, assuming that tile is put in Z direction
			else if (flagp[IX(i, j, k)] == TILE && para->solv.dc_space == TILE_ROOM && put_X) {
				// treat tile as inlet for above fluid cell
				psi[IX(i,j,k)] = var[TILE_FLOW_BC*size+IX(i, j, k)];
				// treat the tile as outlet for below fluid tile
				ae[IX(i - 2, j, k)] = 0;
				psi[IX(i - 1, j, k)] = psi[IX(i - 2, j, k)];
			}
			// Tile
			else if (flagp[IX(i, j, k)] == TILE  && para->solv.tile_flow_correct == NS_SOURCE) {
				// West
				if (i == 0) {
					psi[IX(i, j, k)] = psi[IX(i + 1, j, k)];
				}
				// East
				if (i == imax + 1) {
					psi[IX(i - 1, j, k)] = psi[IX(i - 2, j, k)];
				}
			}	
	
	
	// :Outlet
	if (flagp[IX(i, j, k)] == 2) {
	  // :west
		if (i == 0) {
		  psi[IX(i, j, k)] = psi[IX(i + 1, j, k)];
		} // end of i==0
		// :east
		if (i == imax + 1) {
		  psi[IX(i - 1, j, k)] = psi[IX(i - 2, j, k)];
		} // end of i == imax+1
		// note for other directions, assume that the velocity is zero
	} // end of outlet

      // Rack inlet, which is outlet of DC room
      // Apply a fixed velocity
      if(flagp[IX(i,j,k)]==RACK_INLET){
      	// if west side is fluid
      	if (flagp[IX(i-1,j,k)]==-1) {
      		psi[IX(i-1,j,k)] = var[VXBC*size+IX(i,j,k)];
      	}
      	// if east side is fluid
      	if (flagp[IX(i+1,j,k)]==-1) {
      		psi[IX(i+1,j,k)] = var[VXBC*size+IX(i,j,k)];
      	}
      	// if south side is fluid
      	if (flagp[IX(i,j-1,k)]==-1) {
      		psi[IX(i,j-1,k)] = 0;
      	}
      	// if north side is fluid
      	if (flagp[IX(i,j+1,k)]==-1) {
      		psi[IX(i,j+1,k)] = 0;
      	}
      	// if back side is fluid
      	if (flagp[IX(i,j,k-1)]==-1) {
      		psi[IX(i,j,k-1)] = 0;
      	}
      	// if front side is fluid
      	if (flagp[IX(i,j,k+1)]==-1) {
      		psi[IX(i,j,k+1)] = 0;
      	}
      }// end of setting rack inlet	
	
} // end of set_bnd_VX_im


///////////////////////////////////////////////////////////////////////////////
/// set boundary conditions for VY
/// 6/12/2015
///////////////////////////////////////////////////////////////////////////////
__kernel void set_bnd_VY(__global PARA_DATA_SIMP *para, __global REAL *var, __global int *WF, __global int *BINDEX) {
    int it = get_global_id(0);
  int i, j, k;
  int imax = para->geom.imax, jmax = para->geom.jmax;
  int kmax = para->geom.kmax;
  int size = (imax + 2)*(jmax + 2)*(kmax + 2);
  int IMAX = imax + 2, IJMAX = (imax + 2)*(jmax + 2);
  __global REAL *aw = &var[AW*size], *ae = &var[AE*size], *as = &var[AS*size], *an = &var[AN*size];
  __global REAL *af = &var[AF*size], *ab = &var[AB*size], *b = &var[B*size];
  __global REAL *flagp = &var[FLAGP*size];
    __global REAL *psi;
    // get the index
    i = BINDEX[0*size+it];
    j = BINDEX[1*size+it];
    k = BINDEX[2*size+it]; 
    
    // if Advection
    if (WF[0] == 0) {
      psi = &var[TMP2*size];
    }
    // else if diffusion
    else if (WF[0] == 1) {
      psi = &var[VY*size];  
    }

    // Inlet
    if (flagp[IX(i, j, k)] == INLET) {
      psi[IX(i, j, k)] = var[VYBC*size+IX(i, j, k)];
      if (j != 0) psi[IX(i, j - 1, k)] = var[VYBC*size+IX(i, j, k)];
    }
	
    // rack outlet
    if (flagp[IX(i, j, k)] == RACK_OUTLET) {
      psi[IX(i, j, k)] = var[VYBC*size+IX(i, j, k)];
      if (j != 0) psi[IX(i, j - 1, k)] = var[VYBC*size+IX(i, j, k)];
    }	
	
    // Solid wall
    if (flagp[IX(i, j, k)] == 1) {
      psi[IX(i, j, k)] = 0;
      if (j != 0) psi[IX(i, j - 1, k)] = 0;
    }
    // Outlet
    if (flagp[IX(i, j, k)] == 2 || flagp[IX(i, j, k)] == TILE) {
      // West
      if (i == 0) aw[IX(i + 1, j, k)] = 0;
      // East
      if (i == imax + 1) ae[IX(i - 1, j, k)] = 0;
			
      // if solving tile using presure correction method
			if (flagp[IX(i, j, k)] == TILE && para->solv.tile_flow_correct == PRESSURE_BASE) {
				// South
				if (j == 0) {
					psi[IX(i, j, k)] = var[TILE_FLOW_BC*size+IX(i, j, k)];
				}
				// North
				if (j == jmax + 1) {
					psi[IX(i, j - 1, k)] = var[TILE_FLOW_BC*size+IX(i, j, k)];
				}
			}
			else {
				// South
				if (j == 0) {
					as[IX(i, j + 1, k)] = 0;
					psi[IX(i, j, k)] = psi[IX(i, j + 1, k)];
				}
				// North
				if (j == jmax + 1) {
					an[IX(i, j - 2, k)] = 0;
					psi[IX(i, j - 1, k)] = psi[IX(i, j - 2, k)];
				}
			}
      // Floor
      if (k == 0) ab[IX(i, j, k + 1)] = 0;
      if (k == kmax + 1) af[IX(i, j, k - 1)] = 0;
    }
	
        // Rack inlet, which is outlet of DC room
        // Apply a fixed velocity
        if(flagp[IX(i,j,k)]==RACK_INLET){
        	// if west side is fluid
        	if (flagp[IX(i-1,j,k)]==-1) {
        		psi[IX(i-1,j,k)] = 0;
        	}
        	// if east side is fluid
        	if (flagp[IX(i+1,j,k)]==-1) {
        		psi[IX(i+1,j,k)] = 0;
        	}
        	// if south side is fluid
        	if (flagp[IX(i,j-1,k)]==-1) {
        		psi[IX(i,j-1,k)] = var[VYBC*size+IX(i,j,k)];
        	}
        	// if north side is fluid
        	if (flagp[IX(i,j+1,k)]==-1) {
        		psi[IX(i,j+1,k)] = var[VYBC*size+IX(i,j,k)];
        	}
        	// if back side is fluid
        	if (flagp[IX(i,j,k-1)]==-1) {
        		psi[IX(i,j,k-1)] = 0;
        	}
        	// if front side is fluid
        	if (flagp[IX(i,j,k+1)]==-1) {
        		psi[IX(i,j,k+1)] = 0;
        	}
        }// end of setting rack inlet	
	
}//end of set_bnd_VY()

///////////////////////////////////////////////////////////////////////////////
/// set boundary conditions for VY for advection solved by implicit scheme
/// 8/8/2017
/// Wei Tian
///////////////////////////////////////////////////////////////////////////////
__kernel void set_bnd_VY_im(__global PARA_DATA_SIMP *para, __global REAL *var, __global int *BINDEX) {
    int it = get_global_id(0);
  int i, j, k;
  int imax = para->geom.imax, jmax = para->geom.jmax;
  int kmax = para->geom.kmax;
  int size = (imax + 2)*(jmax + 2)*(kmax + 2);
  int IMAX = imax + 2, IJMAX = (imax + 2)*(jmax + 2);
  __global REAL *aw = &var[AW*size], *ae = &var[AE*size], *as = &var[AS*size], *an = &var[AN*size];
  __global REAL *af = &var[AF*size], *ab = &var[AB*size], *b = &var[B*size];
  __global REAL *flagp = &var[FLAGP*size];
  __global REAL *psi;
  int put_X = para->geom.tile_putX, put_Y = para->geom.tile_putY, put_Z = para->geom.tile_putZ;  

  // get the index
  i = BINDEX[0*size+it];
  j = BINDEX[1*size+it];
  k = BINDEX[2*size+it]; 
    
  // if Advection, psi is tmp2
  psi = &var[TMP2*size];

  // start assigning boundary conditions
	// :Inlet
  if (flagp[IX(i, j, k)] == INLET) {
    psi[IX(i, j, k)] = var[VYBC*size+IX(i, j, k)];
    if (j != 0) psi[IX(i, j - 1, k)] = var[VYBC*size+IX(i, j, k)];
  }
  
	// :Rack outlet
  if (flagp[IX(i, j, k)] == RACK_OUTLET) {
    psi[IX(i, j, k)] = var[VYBC*size+IX(i, j, k)];
    if (j != 0) psi[IX(i, j - 1, k)] = var[VYBC*size+IX(i, j, k)];
  }  
  
	// :Solid
	if (flagp[IX(i, j, k)] == 1) {
    psi[IX(i, j, k)] = 0;
    if (j != 0) psi[IX(i, j - 1, k)] = 0;
  }
  
			// Tile
			if (flagp[IX(i, j, k)] == TILE && para->solv.tile_flow_correct == PRESSURE_BASE) {
				// South
				if (j == 0) {
					psi[IX(i, j, k)] = var[TILE_FLOW_BC*size+IX(i, j, k)];
				}
				// North
				if (j == jmax + 1) {
					psi[IX(i, j - 1, k)] = var[TILE_FLOW_BC*size+IX(i, j, k)];
				}
			}
			// Model the room and plenum together, assuming that tile is put in Z direction
			else if (flagp[IX(i, j, k)] == TILE && para->solv.dc_space == TILE_ROOM && put_Y) {
				// treat tile as inlet for above fluid cell
				psi[IX(i,j,k)] = var[TILE_FLOW_BC*size+IX(i, j, k)];
				// treat the tile as outlet for below fluid tile
				an[IX(i, j - 2, k)] = 0;
				psi[IX(i, j - 1, k)] = psi[IX(i, j - 2, k)];
			}
			// Tile
			else if (flagp[IX(i, j, k)] == TILE && para->solv.tile_flow_correct == NS_SOURCE) {
				// South
				if (j == 0) {
					psi[IX(i, j, k)] = psi[IX(i, j + 1, k)];
				}
				// North
				if (j == jmax + 1) {
					psi[IX(i, j - 1, k)] = psi[IX(i, j - 2, k)];
				}
			}  
  
	// :Outlet
	if (flagp[IX(i, j, k)] == 2) {
    // :south
    if (j == 0) {
      psi[IX(i, j, k)] = psi[IX(i, j + 1, k)];
    }
    // :north
    if (j == jmax + 1) {
      psi[IX(i, j - 1, k)] = psi[IX(i, j - 2, k)];
    }
	} // end of outlet

			// Rack inlet, which is outlet of DC room
			// Apply a fixed velocity
			if(flagp[IX(i,j,k)]==RACK_INLET){
				// if west side is fluid
				if (flagp[IX(i-1,j,k)]==-1) {
					psi[IX(i-1,j,k)] = 0;
				}
				// if east side is fluid
				if (flagp[IX(i+1,j,k)]==-1) {
					psi[IX(i+1,j,k)] = 0;
				}
				// if south side is fluid
				if (flagp[IX(i,j-1,k)]==-1) {
					psi[IX(i,j-1,k)] = var[VYBC*size+IX(i,j,k)];
				}
				// if north side is fluid
				if (flagp[IX(i,j+1,k)]==-1) {
					psi[IX(i,j+1,k)] = var[VYBC*size+IX(i,j,k)];
				}
				// if back side is fluid
				if (flagp[IX(i,j,k-1)]==-1) {
					psi[IX(i,j,k-1)] = 0;
				}
				// if front side is fluid
				if (flagp[IX(i,j,k+1)]==-1) {
					psi[IX(i,j,k+1)] = 0;
				}
			}// end of setting rack inlet	
	
} //end of void set_bnd_VY_im

///////////////////////////////////////////////////////////////////////////////
/// set boundary conditions for VZ
/// 6/12/2015
// Wei Tian
///////////////////////////////////////////////////////////////////////////////
__kernel void set_bnd_VZ(__global PARA_DATA_SIMP *para, __global REAL *var, __global int *WF, __global int *BINDEX) {
    int it = get_global_id(0);
  int i, j, k;
  int imax = para->geom.imax, jmax = para->geom.jmax;
  int kmax = para->geom.kmax;
  int size = (imax + 2)*(jmax + 2)*(kmax + 2);
  int IMAX = imax + 2, IJMAX = (imax + 2)*(jmax + 2);
  __global REAL *aw = &var[AW*size], *ae = &var[AE*size], *as = &var[AS*size], *an = &var[AN*size];
  __global REAL *af = &var[AF*size], *ab = &var[AB*size], *b = &var[B*size];
  __global REAL *flagp = &var[FLAGP*size];
  __global REAL *psi;
    // get the index
    i = BINDEX[0*size+it];
    j = BINDEX[1*size+it];
    k = BINDEX[2*size+it]; 
    // if Advection
    if (WF[0] == 0) {
      psi = &var[TMP3*size];
    }
    // else if diffusion
    else if (WF[0] == 1) {
      psi = &var[VZ*size];  
    }
        // Inlet
    if (flagp[IX(i, j, k)] == INLET) {
      psi[IX(i, j, k)] = var[VZBC*size+IX(i, j, k)];
      if (k != 0) psi[IX(i, j, k - 1)] = var[VZBC*size+IX(i, j, k)];
    }
        // Rack outlet
    if (flagp[IX(i, j, k)] == RACK_OUTLET) {
      psi[IX(i, j, k)] = var[VZBC*size+IX(i, j, k)];
      if (k != 0) psi[IX(i, j, k - 1)] = var[VZBC*size+IX(i, j, k)];
    }	
    if (flagp[IX(i, j, k)] == SOLID) {
      psi[IX(i, j, k)] = 0;
      if (k != 0) psi[IX(i, j, k - 1)] = 0;
    }
    if (flagp[IX(i, j, k)] == OUTLET || flagp[IX(i, j, k)] == TILE) {
      // West
      if (i == 0) aw[IX(i + 1, j, k)] = 0;
      // East
      if (i == imax + 1) ae[IX(i - 1, j, k)] = 0;
      //South
      if (j == 0) as[IX(i, j + 1, k)] = 0;
      // North
      if (j == jmax + 1) an[IX(i, j - 1, k)] = 0;
      
			// if solving tile using presure correction method
			if(flagp[IX(i, j, k)] == TILE && para->solv.tile_flow_correct == PRESSURE_BASE) {
				// Floor
				if (k == 0) {
					psi[IX(i, j, k)] = var[TILE_FLOW_BC*size+IX(i, j, k)];
				}
				// Ceiling
				if (k == kmax + 1) {
					psi[IX(i, j, k - 1)] = var[TILE_FLOW_BC*size+IX(i, j, k)];
				}
			}
			else {
				// Floor
				if (k == 0) {
					ab[IX(i, j, k + 1)] = 0;
					psi[IX(i, j, k)] = psi[IX(i, j, k + 1)];
				}
				// Ceiling
				if (k == kmax + 1) {
					af[IX(i, j, k - 2)] = 0;
					psi[IX(i, j, k - 1)] = psi[IX(i, j, k - 2)];
				}
			}
    }   

        // Rack inlet, which is outlet of DC room
        // Apply a fixed velocity
        if(flagp[IX(i,j,k)]==RACK_INLET){
        	// if west side is fluid
        	if (flagp[IX(i-1,j,k)]==-1) {
        		psi[IX(i-1,j,k)] = 0;
        	}
        	// if east side is fluid
        	if (flagp[IX(i+1,j,k)]==-1) {
        		psi[IX(i+1,j,k)] = 0;
        	}
        	// if south side is fluid
        	if (flagp[IX(i,j-1,k)]==-1) {
        		psi[IX(i,j-1,k)] = 0;
        	}
        	// if north side is fluid
        	if (flagp[IX(i,j+1,k)]==-1) {
        		psi[IX(i,j+1,k)] = 0;
        	}
        	// if back side is fluid
        	if (flagp[IX(i,j,k-1)]==-1) {
        		psi[IX(i,j,k-1)] = var[VZBC*size+IX(i,j,k)];
        	}
        	// if front side is fluid
        	if (flagp[IX(i,j,k+1)]==-1) {
        		psi[IX(i,j,k+1)] = var[VZBC*size+IX(i,j,k)];
        	}
        }// end of setting rack inlet	
	
}//end of set_bnd_VZ()

///////////////////////////////////////////////////////////////////////////////
/// set boundary conditions for VZ for advection solved by implicit scheme
/// 8/8/2017
/// Wei Tian
///////////////////////////////////////////////////////////////////////////////
__kernel void set_bnd_VZ_im(__global PARA_DATA_SIMP *para, __global REAL *var, __global int *BINDEX) {
    int it = get_global_id(0);
  int i, j, k;
  int imax = para->geom.imax, jmax = para->geom.jmax;
  int kmax = para->geom.kmax;
  int size = (imax + 2)*(jmax + 2)*(kmax + 2);
  int IMAX = imax + 2, IJMAX = (imax + 2)*(jmax + 2);
  __global REAL *aw = &var[AW*size], *ae = &var[AE*size], *as = &var[AS*size], *an = &var[AN*size];
  __global REAL *af = &var[AF*size], *ab = &var[AB*size], *b = &var[B*size];
  __global REAL *flagp = &var[FLAGP*size];
  __global REAL *psi;
  int put_X = para->geom.tile_putX, put_Y =para->geom.tile_putY, put_Z = para->geom.tile_putZ;

	// get the index
  i = BINDEX[0*size+it];
  j = BINDEX[1*size+it];
  k = BINDEX[2*size+it]; 
  
	// if Advection
  psi = &var[TMP3*size];

  // :Inlet
  if (flagp[IX(i, j, k)] == INLET) {
    psi[IX(i, j, k)] = var[VZBC*size+IX(i, j, k)];
    if (k != 0) psi[IX(i, j, k - 1)] = var[VZBC*size+IX(i, j, k)];
  } // end of inlet

  // :Rack outlet
  if (flagp[IX(i, j, k)] == RACK_OUTLET) {
    psi[IX(i, j, k)] = var[VZBC*size+IX(i, j, k)];
    if (k != 0) psi[IX(i, j, k - 1)] = var[VZBC*size+IX(i, j, k)];
  } // end of rack outlet  
  
	// :Solid
  if (flagp[IX(i, j, k)] == SOLID) {
    psi[IX(i, j, k)] = 0;
    if (k != 0) psi[IX(i, j, k - 1)] = 0;
  } //end of solid

			// Tile
			if (flagp[IX(i, j, k)] == TILE && para->solv.tile_flow_correct == PRESSURE_BASE) {
				// Floor
				if (k == 0) {
						psi[IX(i, j, k)] = var[TILE_FLOW_BC*size+IX(i, j, k)];
				}
				// Ceiling
				if (k == kmax + 1) {
						psi[IX(i, j, k - 1)] = var[TILE_FLOW_BC*size+IX(i, j, k)];
				}
			}
			// Model the room and plenum together, assuming that tile is put in Z direction
			else if (flagp[IX(i, j, k)] == TILE && para->solv.dc_space == TILE_ROOM && put_Z) {
				// treat tile as inlet for above fluid cell
				psi[IX(i,j,k)] = var[TILE_FLOW_BC*size+IX(i, j, k)];
				// treat the tile as outlet for below fluid tile
				af[IX(i, j, k - 2)] = 0;
				psi[IX(i, j, k - 1)] = psi[IX(i, j, k - 2)];
			}
			// Tile
			else if (flagp[IX(i, j, k)] == TILE && para->solv.tile_flow_correct == NS_SOURCE) {
				// Floor
				if (k == 0) {
					psi[IX(i, j, k)] = psi[IX(i, j, k + 1)];
				}
				// Ceiling
				if (k == kmax + 1) {
					psi[IX(i, j, k - 1)] = psi[IX(i, j, k - 2)];
				}
			}  
  
	// :Outlet
  if (flagp[IX(i, j, k)] == OUTLET) {
    // :Floor
    if (k == 0) {
      psi[IX(i, j, k)] = psi[IX(i, j, k + 1)];
    }
    // :Ceiling
    if (k == kmax + 1) {
      psi[IX(i, j, k - 1)] = psi[IX(i, j, k - 2)];
    }
  } // end of outlet

      // Rack inlet, which is outlet of DC room
      // Apply a fixed velocity
      if(flagp[IX(i,j,k)]==RACK_INLET){
      	// if west side is fluid
      	if (flagp[IX(i-1,j,k)]==-1) {
      		psi[IX(i-1,j,k)] = 0;
      	}
      	// if east side is fluid
      	if (flagp[IX(i+1,j,k)]==-1) {
      		psi[IX(i+1,j,k)] = 0;
      	}
      	// if south side is fluid
      	if (flagp[IX(i,j-1,k)]==-1) {
      		psi[IX(i,j-1,k)] = 0;
      	}
      	// if north side is fluid
      	if (flagp[IX(i,j+1,k)]==-1) {
      		psi[IX(i,j+1,k)] = 0;
      	}
      	// if back side is fluid
      	if (flagp[IX(i,j,k-1)]==-1) {
      		psi[IX(i,j,k-1)] = var[VZBC*size+IX(i,j,k)];
      	}
      	// if front side is fluid
      	if (flagp[IX(i,j,k+1)]==-1) {
      		psi[IX(i,j,k+1)] = var[VZBC*size+IX(i,j,k)];
      	}
      }// end of setting rack inlet
	  
} // end of set_bnd_VZ_im

///////////////////////////////////////////////////////////////////////////////
/// kernels to set the AP coefficients
/// 6/12/2015
// Wei Tian
///////////////////////////////////////////////////////////////////////////////
__kernel void ap_coeff(__global PARA_DATA_SIMP *para, __global REAL *var, __global int *var_type) {
  
  int i = get_global_id(0);
  int j = get_global_id(1);
  int k = get_global_id(2);
  int imax = para->geom.imax, jmax = para->geom.jmax;
  int kmax = para->geom.kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int size=(imax+2)*(jmax+2)*(kmax+2);
  __global float *aw = &var[AW*size];
  __global float *ae = &var[AE*size];
  __global float *as = &var[AS*size];
  __global float *an = &var[AN*size];
  __global float *af = &var[AF*size];
  __global float *ab = &var[AB*size]; 
  __global float *ap = &var[AP*size];
  __global float *ap0 = &var[AP0*size]; 
    // U
  if (var_type[0] == VX) {
    // FOR_U_CELL
    if (i>0 && i<imax && j>0 && j<jmax+1 && k>0 && k<kmax+1) {
          ap[IX(i,j,k)] = ap0[IX(i,j,k)] + ae[IX(i,j,k)] + aw[IX(i,j,k)]
                        + an[IX(i,j,k)]  + as[IX(i,j,k)] + af[IX(i,j,k)]
                        + ab[IX(i,j,k)];
    } //end of FOR_U_CELL

  } //end of U
  // V
  else if (var_type[0] == VY) {
    // FOR_V_CELL
    if (i>0 && i<imax+1 && j>0 && j<jmax && k>0 && k<kmax+1) {
          ap[IX(i,j,k)] = ap0[IX(i,j,k)] + ae[IX(i,j,k)] + aw[IX(i,j,k)]
                        + an[IX(i,j,k)]  + as[IX(i,j,k)] + af[IX(i,j,k)]
                        + ab[IX(i,j,k)];
    } //end of FOR_V_CELL
  } //end of V
  // W
  if (var_type[0] == VZ) {
    // FOR_W_CELL
    if (i>0 && i<imax+1 && j>0 && j<jmax+1 && k>0 && k<kmax) {
          ap[IX(i,j,k)] = ap0[IX(i,j,k)] + ae[IX(i,j,k)] + aw[IX(i,j,k)]
                        + an[IX(i,j,k)]  + as[IX(i,j,k)] + af[IX(i,j,k)]
                        + ab[IX(i,j,k)];
    } //end of FOR_W_CELL
  } //end of W
  // T
  if (var_type[0] == TEMP){
      // FOR_ALL_CELL
    if (i>0 && i<imax+1 && j>0 && j<jmax+1 && k>0 && k<kmax+1) {
          ap[IX(i,j,k)] = ap0[IX(i,j,k)] + ae[IX(i,j,k)] + aw[IX(i,j,k)]
                        + an[IX(i,j,k)]  + as[IX(i,j,k)] + af[IX(i,j,k)]
                        + ab[IX(i,j,k)];
    } //end of FOR_W_CELL
  } //end of T
  // Pressure
  if (var_type[0] == IP){
    // FOR_ALL_CELL
    if (i>0 && i<imax+1 && j>0 && j<jmax+1 && k>0 && k<kmax+1) {
       ap[IX(i,j,k)] = ae[IX(i,j,k)] + aw[IX(i,j,k)] + as[IX(i,j,k)] + an[IX(i,j,k)]
                  + af[IX(i,j,k)] + ab[IX(i,j,k)];
    }
  } //end of Pressure
} //end of __kernel void ap_coeff

///////////////////////////////////////////////////////////////////////////////
/// kernels to set the AP coefficients
/// 6/12/2015
// Wei Tian
///////////////////////////////////////////////////////////////////////////////
__kernel void ap_im_coeff(__global PARA_DATA_SIMP *para, __global REAL *var, __global int *var_type) {
  
  int i = get_global_id(0);
  int j = get_global_id(1);
  int k = get_global_id(2);
  int imax = para->geom.imax, jmax = para->geom.jmax;
  int kmax = para->geom.kmax;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int size=(imax+2)*(jmax+2)*(kmax+2);
  __global float *aw = &var[AW*size];
  __global float *ae = &var[AE*size];
  __global float *as = &var[AS*size];
  __global float *an = &var[AN*size];
  __global float *af = &var[AF*size];
  __global float *ab = &var[AB*size]; 
  __global float *ap = &var[AP*size];
  __global float *ap0 = &var[AP0*size]; 

	__global float *x = &var[X*size];
  __global float *y = &var[Y*size];
  __global float *z = &var[Z*size];
	__global float *gx = &var[GX*size];
  __global float *gy = &var[GY*size];
  __global float *gz = &var[GZ*size];
	__global float *u = &var[VX*size];
	__global float *v = &var[VY*size];
	__global float *w = &var[VZ*size];
   
  REAL dxe, dxw, dyn, dys, dzf, dzb, Dx, Dy, Dz;
	REAL uw, ue, vs, vn, wb, wf;
	REAL Fw, Fe, Fs, Fn, Fb, Ff;

    // U
  if (var_type[0] == VX) {
    // FOR_U_CELL
    if (i>0 && i<imax && j>0 && j<jmax+1 && k>0 && k<kmax+1) {
				// define the dimensions
				dxe = gx[IX(i + 1, j, k)] - gx[IX(i, j, k)];
				dxw = gx[IX(i, j, k)] - gx[IX(i - 1, j, k)];
				Dx = x[IX(i + 1, j, k)] - x[IX(i, j, k)];
				Dy = gy[IX(i, j, k)] - gy[IX(i, j - 1, k)];
				Dz = gz[IX(i, j, k)] - gz[IX(i, j, k - 1)];

				// define the velocity at the surface
				uw = 0.5*(u[IX(i - 1, j, k)] + u[IX(i, j, k)]);
				ue = 0.5*(u[IX(i, j, k)] + u[IX(i + 1, j, k)]);
				vs = 0.5*(v[IX(i, j-1, k)] + v[IX(i + 1, j-1, k)]);
				vn = 0.5*(v[IX(i, j, k)] + v[IX(i + 1, j, k)]);
				wb = 0.5*(w[IX(i, j, k-1)] + w[IX(i+1, j, k-1)]);
				wf = 0.5*(w[IX(i, j, k)] + w[IX(i + 1, j, k)]);

				// define the flow rate at the surface
				Fw = uw * Dy*Dz;
				Fe = ue * Dy*Dz;
				Fs = vs * Dx*Dz;
				Fn = vn * Dx*Dz;
				Fb = wb * Dx*Dy;
				Ff = wf * Dx*Dy;
        
				// update the ap coefficient
				ap[IX(i,j,k)] = ap0[IX(i,j,k)] + ae[IX(i,j,k)] + aw[IX(i,j,k)]
                        + an[IX(i,j,k)]  + as[IX(i,j,k)] + af[IX(i,j,k)]
                        + ab[IX(i,j,k)] + Fe - Fw + Fn - Fs + Ff - Fb;
    } //end of FOR_U_CELL

  } //end of U
  // V
  else if (var_type[0] == VY) {
    // FOR_V_CELL
    if (i>0 && i<imax+1 && j>0 && j<jmax && k>0 && k<kmax+1) {
		    // define the dimensions
				dyn = gy[IX(i, j + 1, k)] - gy[IX(i, j, k)];
				dys = gy[IX(i, j, k)] - gy[IX(i, j - 1, k)];
				Dx = gx[IX(i, j, k)] - gx[IX(i - 1, j, k)];
				Dy = y[IX(i, j + 1, k)] - y[IX(i, j, k)];
				Dz = gz[IX(i, j, k)] - gz[IX(i, j, k - 1)];

				// define the velocity at the surface
				uw = 0.5*(u[IX(i - 1, j, k)] + u[IX(i - 1, j + 1, k)]);
				ue = 0.5*(u[IX(i, j, k)] + u[IX(i, j + 1, k)]);
				vs = 0.5*(v[IX(i, j - 1, k)] + v[IX(i, j, k)]);
				vn = 0.5*(v[IX(i, j, k)] + v[IX(i, j + 1, k)]);
				wb = 0.5*(w[IX(i, j, k - 1)] + w[IX(i, j + 1, k - 1)]);
				wf = 0.5*(w[IX(i, j, k)] + w[IX(i, j +1, k)]);

				// define the flow rate at the surface
				Fw = uw * Dy*Dz;
				Fe = ue * Dy*Dz;
				Fs = vs * Dx*Dz;
				Fn = vn * Dx*Dz;
				Fb = wb * Dx*Dy;
				Ff = wf * Dx*Dy;

				// update the ap coefficient
        ap[IX(i,j,k)] = ap0[IX(i,j,k)] + ae[IX(i,j,k)] + aw[IX(i,j,k)]
                        + an[IX(i,j,k)]  + as[IX(i,j,k)] + af[IX(i,j,k)]
                        + ab[IX(i,j,k)] + Fe - Fw + Fn - Fs + Ff - Fb;
    } //end of FOR_V_CELL
  } //end of V
  // W
  if (var_type[0] == VZ) {
    // FOR_W_CELL
    if (i>0 && i<imax+1 && j>0 && j<jmax+1 && k>0 && k<kmax) {
        // define the dimensions
				dzf = gz[IX(i, j, k + 1)] - gz[IX(i, j, k)];
				dzb = gz[IX(i, j, k)] - gz[IX(i, j, k - 1)];
				Dx = gx[IX(i, j, k)] - gx[IX(i - 1, j, k)];
				Dy = gy[IX(i, j, k)] - gy[IX(i, j - 1, k)];
				Dz = z[IX(i, j, k + 1)] - z[IX(i, j, k)];

				// define the velocity at the surface
				uw = 0.5*(u[IX(i - 1, j, k)] + u[IX(i - 1, j, k + 1)]);
				ue = 0.5*(u[IX(i, j, k)] + u[IX(i, j, k + 1)]);
				vs = 0.5*(v[IX(i, j - 1, k)] + v[IX(i, j - 1, k + 1)]);
				vn = 0.5*(v[IX(i, j, k)] + v[IX(i, j, k + 1)]);
				wb = 0.5*(w[IX(i, j, k - 1)] + w[IX(i, j, k)]);
				wf = 0.5*(w[IX(i, j, k)] + w[IX(i, j, k + 1)]);


				// define the flow rate at the surface
				Fw = uw * Dy*Dz;
				Fe = ue * Dy*Dz;
				Fs = vs * Dx*Dz;
				Fn = vn * Dx*Dz;
				Fb = wb * Dx*Dy;
				Ff = wf * Dx*Dy;

				// update the ap coefficient
				ap[IX(i,j,k)] = ap0[IX(i,j,k)] + ae[IX(i,j,k)] + aw[IX(i,j,k)]
                        + an[IX(i,j,k)]  + as[IX(i,j,k)] + af[IX(i,j,k)]
                        + ab[IX(i,j,k)] + Fe - Fw + Fn - Fs + Ff - Fb;
    } //end of FOR_W_CELL
  } //end of W
  // T
  if (var_type[0] == TEMP){
      // FOR_ALL_CELL
    if (i>0 && i<imax+1 && j>0 && j<jmax+1 && k>0 && k<kmax+1) {
        // define the dimensions
				Dx = gx[IX(i, j, k)] - gx[IX(i - 1, j, k)];
				Dy = gy[IX(i, j, k)] - gy[IX(i, j - 1, k)];
				Dz = gz[IX(i, j, k)] - gz[IX(i, j, k - 1)];

				// define the velocity at the surface
				uw = u[IX(i - 1, j, k)];
				ue = u[IX(i, j, k )];
				vs = v[IX(i, j - 1, k)];
				vn = v[IX(i, j, k)];
				wb = w[IX(i, j, k - 1)];
				wf = w[IX(i, j, k)];
				
				// define the flow rate at the surface
				Fw = uw * Dy*Dz;
				Fe = ue * Dy*Dz;
				Fs = vs * Dx*Dz;
				Fn = vn * Dx*Dz;
				Fb = wb * Dx*Dy;
				Ff = wf * Dx*Dy;

				// update the ap coefficient
				// change the demoninator by removing "+ Fe - Fw + Fn - Fs + Ff - Fb;"
				// this is to avoid the incorrect cacluation when the mass balance is not
				//  strictly met
				ap[IX(i,j,k)] = ap0[IX(i,j,k)] + ae[IX(i,j,k)] + aw[IX(i,j,k)]
                        + an[IX(i,j,k)]  + as[IX(i,j,k)] + af[IX(i,j,k)]
                        + ab[IX(i,j,k)]; 
    } //end of FOR_ALL_CELL
  } //end of T

} //end of __kernel void ap_im_coeff

///////////////////////////////////////////////////////////////////////////////
/// Kernel to adjust velocity
/// 7/1/2015
// Wei Tian
///////////////////////////////////////////////////////////////////////////////
__kernel void adjust_velocity(__global PARA_DATA_SIMP *para, __global REAL *var, __global int *BINDEX) {
  int i, j, k;
  int it = get_global_id(0);;
  int imax = para->geom.imax, jmax = para->geom.jmax;
  int kmax = para->geom.kmax;
  int index= para->geom.index;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int size=(imax+2)*(jmax+2)*(kmax+2);
   __global REAL *u = &var[VX*size], *v = &var[VY*size], *w = &var[VZ*size];
  REAL dvel;
   __global REAL *flagp = &var[FLAGP*size];

  dvel = para->bc_simp.mass_corr;

  /*---------------------------------------------------------------------------
  | Adjust the outflow
  ---------------------------------------------------------------------------*/
    i = BINDEX[0*size+it];
    j = BINDEX[1*size+it];
    k = BINDEX[2*size+it];

    if(flagp[IX(i,j,k)]==2) {
      if(i==0) u[IX(i,j,k)] -= dvel;
      if(i==imax+1) u[IX(i-1,j,k)]+= dvel;
      if(j==0) v[IX(i,j,k)] -= dvel;
      if(j==jmax+1) v[IX(i,j-1,k)] += dvel;
      if(k==0) w[IX(i,j,k)] -= dvel;
      if(k==kmax+1) w[IX(i,j,k-1)] += dvel;
    }
} // End of mass_conservation()


///////////////////////////////////////////////////////////////////////////////
/// Kernel to set the force mass balance
/// Use only one core
/// 7/1/2015
// Wei Tian
///////////////////////////////////////////////////////////////////////////////
__kernel void mass_conservation(__global PARA_DATA_SIMP *para, __global REAL *var, __global int *BINDEX) {
  int i, j, k;
  int it;
  int imax = para->geom.imax, jmax = para->geom.jmax;
  int kmax = para->geom.kmax;
  int index= para->geom.index;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int size=(imax+2)*(jmax+2)*(kmax+2);
   __global REAL *u = &var[VX*size], *v = &var[VY*size], *w = &var[VZ*size];
  REAL mass_in = (REAL) 0.0, mass_out = (REAL) 0.00000001;
  REAL area_out=0;
   __global REAL *flagp = &var[FLAGP*size];
  REAL axy, ayz, azx;

  // Go through all the inlets and outlets
  for(it=0; it<index; it++) {
    i = BINDEX[0*size+it];
    j = BINDEX[1*size+it];
    k = BINDEX[2*size+it];

    axy = area_xy(para, var, &i, &j, &k);
    ayz = area_yz(para, var, &i, &j, &k);
    azx = area_zx(para, var, &i, &j, &k);
    /*-------------------------------------------------------------------------
    | Compute the total inflow
    -------------------------------------------------------------------------*/
    if(flagp[IX(i,j,k)]==0) {
      // West
      if(i==0) mass_in += u[IX(i,j,k)] * ayz;
      // East
      if(i==imax+1) mass_in += (-u[IX(i,j,k)]) * ayz;
      // South
      if(j==0) mass_in += v[IX(i,j,k)] * azx;
      // North
      if(j==jmax+1) mass_in += (-v[IX(i,j,k)]) * azx;
      // Floor
      if(k==0) mass_in += w[IX(i,j,k)] * axy;
      // Ceiling
      if(k==kmax+1) mass_in += (-w[IX(i,j,k)]) * axy;
    }
    /*-------------------------------------------------------------------------
    | Compute the total outflow
    -------------------------------------------------------------------------*/
    if(flagp[IX(i,j,k)]==2) {
      // West
      if(i==0) {
        mass_out += (-u[IX(i,j,k)]) * ayz;
        area_out +=  ayz;
      }
      // East
      if(i==imax+1) {
        mass_out += u[IX(i-1,j,k)] * ayz;
        area_out += ayz;
      }
      // South
      if(j==0) {
        mass_out += (-v[IX(i,j,k)]) * azx;
        area_out += azx;
      }
      // North
      if(j==jmax+1) {
        mass_out += v[IX(i,j-1,k)] * azx;
        area_out += azx;
      }
      // Floor
      if(k==0) {
        mass_out += (-w[IX(i,j,k)]) * axy;
        area_out += axy;
      }
      // Ceiling
      if(k==kmax+1) {
        mass_out += w[IX(i,j,k-1)] * axy;
        area_out += axy;
      }
    } // End of computing outflow
  } // End of for loop for going through all the inlets and outlets

  /*---------------------------------------------------------------------------
  | Return the adjusted velocity for mass conservation
  ---------------------------------------------------------------------------*/
  para->bc_simp.mass_corr = (mass_in-mass_out)/area_out;

  //printf("the velocity correction is %f\n",para->bc_simp.mass_corr );
} // End of adjust_velocity()

///////////////////////////////////////////////////////////////////////////////
/// Kernel to reset average data
/// 8/9/2015
// Wei Tian
///////////////////////////////////////////////////////////////////////////////
__kernel void reset_time_averaged_data(__global PARA_DATA_SIMP *para, __global REAL *var) {
  int i = get_global_id(0);
  int j = get_global_id(1);
  int k = get_global_id(2);
  int it;
  int step = para->mytime.step_mean;
  int imax = para->geom.imax, jmax = para->geom.jmax;
  int kmax = para->geom.kmax;
  int index= para->geom.index;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int size=(imax+2)*(jmax+2)*(kmax+2);
   __global REAL *u = &var[VX*size], *v = &var[VY*size], *w = &var[VZ*size], *T = &var[TEMP*size];
   __global REAL *um = &var[VXM*size], *vm = &var[VYM*size], *wm = &var[VZM*size], *Tm = &var[TEMPM*size];

   step = 0;
   um[IX(i,j,k)] = 0;
   vm[IX(i,j,k)] = 0;
   wm[IX(i,j,k)] = 0;
   Tm[IX(i,j,k)] = 0;

} //end of reset_time_averaged_data()


///////////////////////////////////////////////////////////////////////////////
/// Kernel to average data
/// 8/9/2015
// Wei Tian
///////////////////////////////////////////////////////////////////////////////
__kernel void add_time_averaged_data(__global PARA_DATA_SIMP *para, __global REAL *var) {
  int i = get_global_id(0);
  int j = get_global_id(1);
  int k = get_global_id(2);
  int it;
  int step = para->mytime.step_mean;
  int imax = para->geom.imax, jmax = para->geom.jmax;
  int kmax = para->geom.kmax;
  int index= para->geom.index;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int size=(imax+2)*(jmax+2)*(kmax+2);
   __global REAL *u = &var[VX*size], *v = &var[VY*size], *w = &var[VZ*size], *T = &var[TEMP*size];
   __global REAL *um = &var[VXM*size], *vm = &var[VYM*size], *wm = &var[VZM*size], *Tm = &var[TEMPM*size];

   um[IX(i,j,k)] = um[IX(i,j,k)] + u[IX(i,j,k)];
   vm[IX(i,j,k)] = vm[IX(i,j,k)] + v[IX(i,j,k)];
   wm[IX(i,j,k)] = wm[IX(i,j,k)] + w[IX(i,j,k)];
   Tm[IX(i,j,k)] = Tm[IX(i,j,k)] + T[IX(i,j,k)];
   //Mum[IX(i,j,k)] = Mum[IX(i,j,k)] + sqrt(u[IX(i,j,k)]*u[IX(i,j,k)]+v[IX(i,j,k)]*v[IX(i,j,k)]+w[IX(i,j,k)]*w[IX(i,j,k)]);   
   para->mytime.step_mean = para->mytime.step_mean +1;
} //end of add_time_averaged_data()

///////////////////////////////////////////////////////////////////////////////
/// Kernel to average data (void, not used at the moment)
/// 8/9/2015
// Wei Tian
///////////////////////////////////////////////////////////////////////////////
__kernel void time_averaged(__global PARA_DATA_SIMP *para, __global REAL *var) {
  int i = get_global_id(0);
  int j = get_global_id(1);
  int k = get_global_id(2);
  int step = para->mytime.step_mean;
  int imax = para->geom.imax, jmax = para->geom.jmax;
  int kmax = para->geom.kmax;
  int index= para->geom.index;
  int IMAX = imax+2, IJMAX = (imax+2)*(jmax+2);
  int size=(imax+2)*(jmax+2)*(kmax+2);
   __global REAL *u = &var[VX*size], *v = &var[VY*size], *w = &var[VZ*size], *T = &var[TEMP*size];
   __global REAL *um = &var[VXM*size], *vm = &var[VYM*size], *wm = &var[VZM*size], *Tm = &var[TEMPM*size];

   um[IX(i,j,k)] = um[IX(i,j,k)]/para->mytime.step_mean;
   vm[IX(i,j,k)] = vm[IX(i,j,k)]/para->mytime.step_mean;
   wm[IX(i,j,k)] = wm[IX(i,j,k)]/para->mytime.step_mean;
   Tm[IX(i,j,k)] = Tm[IX(i,j,k)]/para->mytime.step_mean;
   printf("the step is %d\n", para->mytime.step_mean);
} //end of time_averaged()