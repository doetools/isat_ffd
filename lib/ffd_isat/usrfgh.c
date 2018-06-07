//Note: User do not need edit this file.

#include"usrfgh.h"

extern int useNumericalDifferentiation;  // pass global variable

///////////////////////////////////////////////////////////////////////////////
//subroutine usrfgh( need, nx, x, nf, nh, iusr, rusr, fa, ga, ha )
//==========================================================================
//-------  INPUT
//   need       - integer array indicating function values to be returned
//   need(1)      = 1, return  f   
//   need(2)      = 1, return  g = df/dx
//   nx         - number of components of x  ( nx >= 1 )
//   x          - components of x
//   nf         - number of components of f  ( nf >= 1 )
//   nh         - number of components of h  ( nh >= 0 )
//   iusr       - integer array: iusr(1) = idtab, iusr(2) = ifsin
//   rusr       - pointer for ffd
//-------  OUTPUT
//   f          - f(x)
//   g          - g(nf,nx)
//   h          - not used
//!==========================================================================
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
/// Solver to evaluate f, g, and h which called by ISATAB
///
///\param All pointers which definde in ISATAB
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void USRFGH(int need[], int *nx, double x[], int *nf, int *nh, int iusr[], double rusr[], double f[], double** g, double h[]){

  if (need[1] == 1 && useNumericalDifferentiation == 1) {
    need[1] = 0; //disable g evaluation form solver
    numericalDifferentiation(g);
  }

  ffd_ISAT(need, x, f, g, (void *)iusr);

  return;
}

///////////////////////////////////////////////////////////////////////////////
/// Evaluate Numerical Differentiation of f(x)
///
///\param Pointer of g[nx][nf]
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void numericalDifferentiation (double g[][nf_SIZE]){
  // Evaluate Numerical Differentiation of f(x)

  mpc_log("numericalDifferentiation(): Numerical Differentiation has not been implemented", MPC_ERROR);
  exit(1);
}
