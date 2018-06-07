#include "io.h"


///////////////////////////////////////////////////////////////////////////////
/// Write the log file
///
///\param message Pointer the message
///\param msg_type Type ogf message
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////


void mpc_log(char *message, MPC_MSG_TYPE msg_type) {
  if(msg_type==MPC_NEW) {
    if((file_log=fopen("log.mpc","w+"))==NULL) {
        fprintf(stderr, "Error:can not open error file!\n");
        exit(1);
    }
  }
  else if((file_log=fopen("log.mpc","a+"))==NULL) {
    fprintf(stderr,"Error:can not open error file!\n");
    exit(1);
  }

  switch(msg_type) {
    case MPC_WARNING:
      fprintf(file_log, "WARNING in %s\n", message);
      break;
    case MPC_ERROR:
      fprintf(file_log, "ERROR in %s\n", message);
      break;
    // Normal log
    default:
      fprintf(file_log, "%s\n", message);
  }
  fclose(file_log);
} // End of mpc_log()


///////////////////////////////////////////////////////////////////////////////
/// Record FFD simulation data at each step
///
///\param para Pointer to FFD parameters
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////

//int InitializeFFDData(PARA_DATA *para) {
//  //Check total steps set value
//  if (totalTimeSteps != para->mytime->step_total){
//    mpc_log("InitializeFFDData(): Total simulation steps setting has different value at common.h and FFD", MPC_ERROR);
//  }
//
//  //Initilize FFDRecord to zero
//  memset(FFDRecord, 0, sizeof(FFDRecord));
//   
//}
///////////////////////////////////////////////////////////////////////////////
/// Record FFD simulation data at each step
///
///\param para Pointer to FFD parameters
///\param var Pointer to FFD simulation variables
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////


//void recordFFDData(PARA_DATA *para, REAL **var) {
//
//  int step = para->mytime->step_current;
//
//  if (step != 0) {
//
//    int imax, jmax, kmax; 
//    int IMAX, IJMAX;
//
//    imax = para->geom->imax; 
//    imax = para->geom->jmax; 
//    imax = para->geom->kmax; 
//    IMAX = imax+2;
//    IJMAX = (imax+2)*(jmax+2);
//
//
//    FFDRecord[1][step-skipSteps] = var[TEMP][IX(imax/2,jmax/2,kmax*3/4)];
//
//    //ffdOutput[1] = var[VX][IX(imax/2,jmax/2,kmax*3/4)];
//    //ffdOutput[2] = var[VY][IX(imax/2,jmax/2,kmax*3/4)];
//    //ffdOutput[3] = var[VZ][IX(imax/2,jmax/2,kmax*3/4)];
//    //ffdOutput[0] = var[TEMP][IX(10,10,14)];
//
//  }
//}
//
/////////////////////////////////////////////////////////////////////////////////
///// Caculate and write simulation result to ffdOutput for Fortran program read
/////
/////\param No parameter needed
/////
/////\return No return needed
/////////////////////////////////////////////////////////////////////////////////
//void writeFFDOutput() {
//
//  int a =1;
//
//
//
//}