#include "mpc.h"
#include <math.h>

///////////////////////////////////////////////////////////////////////////////
//recursive subroutine isatab( idtab, mode, nx, x, nf, nh, nhd, usrfgh, iusr, rusr, info, rinfo, fa, ga, ha, stats )
//-------  INPUT
//
//   idtab       - unique identifier of the table (idtab >= 0 )
//   mode        - integer determining the action to be taken 
//                 = 0 for a "regular query call" (i.e., given x, return fa)
//                 > 0 for "special calls" (described below)
//   nx          - number of components of x  ( nx >= 1 )
//   x           - components of x
//   nf          - number of components of f  ( nf >= 1 )
//   nh          - number of components of h  ( nh >= 0 )
//   nhd         - dimension of ha ( nhd >= max(1,nh) )
//   usrfgh      - name of the user-supplied subroutine that returns
//                 f(x), g(x) = df(x)/dx, and h(x) (see below)
//   iusr        - user-defined integer array passed to usrfgh
//                 Note: instead passing integer array
//                 ffd exchange struct pointer here
//   rusr        - user-defined real array passed to usrfgh
//   info        - integer array controlling ISATAB operations (see below)
//   rinfo       - real array controlling ISATAB operations (see below)
//
//-------  OUTPUT
//
//   fa          - piecewise-linear approximation to f(x)
//   ga          - piecewise-constant approximation to g(x) (for info(2)=1)
//                 (ga must be dimensioned at least nf*nx)
//   ha          - piecewise-constant approximation to h(x) (for info(3)=1)
//   stats       - statistics of ISATAB performance (see below)
//
//  Note: Depending on  mode  and other settings, some of the input may not be
//  referenced, and some of the output may not be set. (Write for ISAT Fortan lib)
///////////////////////////////////////////////////////////////////////////////

/******************************************************************************
| Declaration of ISATAB
******************************************************************************/
void ISATAB(int *idtab, int *mode, const int *nx, double x[], const int *nf, const int *nh, const int *nhd, void *usrfgh, \
             int iusr[], double rusr[], int info[], double rinfo[], double fa[], double ga[nx_SIZE][nf_SIZE], double ha[], double stats[] );

/* global variables */

/******************************************************************************
| Initialize parameter for calling ISAT or slover
******************************************************************************/

  //Please define nx_SIZE, nf_SIZE, nh_SIZE in common hearder file.
  //Note: usually only change nx_SIZE and nf_SIZE

  //#define nx_SIZE 1                   // Dimension of x
  //#define nf_SIZE 1                   // Dimension of f
  //#define nh_SIZE 1                   // Dimension of h
  
  #define ISAT_SMALL 1e-6
  // define the ita
  float ita = 5 / 10000;              // allow 5 add or grow per 10000 queries
  float VPower = 1.0;
  int idtab = 1;                      // Just set 1 if only one table
  int mode = 0;                       // mode = 0 is defualt

  int nx;                             // Calculate in main
  int nf;                             // Calculate in main

  int nh = 0;                         // If the additional function h(x) is not required, set nh=0, nhd=1.
  int nhd = 1;                        // Function h is not used so far

  double x[nx_SIZE] = {0};            // Initialize x
  double fa[nf_SIZE] = {0};           // Initialize fa
  double ga[nx_SIZE][nf_SIZE] = {0};  // Initialize ga. Note: In Fortran, g(nf,nx). Initialize in C by a reverse matrix
  double ha[nh_SIZE] = {0};           // Initialize fa

  //void *usrfgh                      // Replaced by unused pointer
  //int iusr[1] = {0};                // User defined struct here

  double rusr[nf_SIZE] = {nf_SIZE};   // Note rusr[1]=its length, so >=1. It may used to contain f in some special usage
  int info[70] = {0};                 // Initialize by dafult value 0
  double rinfo[100] = {0};            // Initialize by dafult value 0
  double stats[100] = {0};            // Initialize memory for store performance info
  double unusedPointer = 0;           // Initialize a unused pointer
  int need[3] = {1,0,0};              // Initialize need for call solver directly
  
  double time_elipse = 0; //the time before launching the isat
  ffdIO ffdStruct;                    // A struct passed to ffd contains constant
                                      // Note: only some constant could pass to ffd through here. Input and out put must in x,fa and ga.

/******************************************************************************
| Initialize variables for performance statistics
******************************************************************************/

  clock_t tStart = 0;
  clock_t tEnd = 0;
  time_t tExeStart;                   //Start of system time
  time_t tExeEnd;                     //End of system time
  double cpuISAT = 0;                 //CPU secounds on ISAT
  double cpuDE = 0;                   //CPU secounds on direct evaluation. Direct evaluation called by ISAT is not included.
  double cpuCum = 0;
  int nQuery = 0;                     //Number of queries
  int QueryLimit = 10000;             //Number of queries
  double errSum = 0;                  //Sum of f error
  double errMax = 0;                  //Max of f error
  int caltyp =0;                      //0 training
                                      //1 test
  //Define hash table
  hashStruct l, *p, *r, *tmp, *records = NULL;

/******************************************************************************
| MPC control global variable
******************************************************************************/
  //Logistic
  int useISAT = 1;                         //If use ISAT.
  int useRandomTest = 0;                   //If use Random Test.Need: xBoundary, nRanCall.
  int useTablePreparation = 1;             //If use Table Preparation. Need: xBoundary, xStep.
  int useBinarySelectedPoint = 0;          //If use binary selected training point, like binary tree.

  int useNumericalDifferentiation = 0;     //If use numerical differentiation for ga. Note: Not be implemented yet.
  int useAccuracyTest = 1;                 //If use accuracy test. Note: Perform a direct evaluation at each query. CPU time added in cpuDE
  int useErrorContraolCorrection = 1;      /*if use error contraol correction. The total error defind as vector length. Use this to crroect order of magnitudes of each component of output vector.
                                             please implement this at ffd_ISAT.c */
  int useBoundaryCenterRange=0;            //if use Boundary center range. The maximum inputs value difference equals xBoundaryCenterRange/2. Here use it to set limititon of different wall temperature.

  int useNormalDistribution=1;             /* = 0, use uniform distribution at random test
                                              = 1, use normal distribution at random test */
  int useUnboundedTest = 1;                   /* = 0, use xBoundary for random test
                                                = 1, use xBoundary2 for random test */
  int useRoundInput = 1;                   // If use rounded input, chose rounded digits after decimal

  //Value
  //double xBoundary[nx_SIZE][2]={{0,5},{2,5},{4,7}};
  double xBoundary[nx_SIZE][2]={{25.0,30.0},{25.0,30.0}};
                                           //Contains lower and upper boundary of x.
  double xBoundaryCenterRange=10;
                                           //Contains centered range of inpust in a rectangular domian. To set maximum difference of inputs at different dimensions.
  //double xStep[nx_SIZE]={1,1,1};
  double xStep[nx_SIZE]={0.1,0.1};
                                           //Step size for gerarate table.
  int nRanCall = 65;                     //Number of Random call for testing

  double xBoundary2[nx_SIZE][2]={{15.0,35.0},{15.0,35.0}};
                                           //Contains lower and upper boundary of x for "unbounderd case".
  double mu=25, sigma=10/3;  // for bounded
  double sigma2=10/3;  // for unbounded
  
//****************************************************  
  int digAftdec = 1;                               //*
												   //*
	#define RoundDigitsControl 1                   //*
	#define MANUAL_TEST                            //*
	#ifdef MANUAL_TEST   //Called by ISAT          //*
	  double *input1, *input2;                     //*
	  FILE *file_params_tmp;                       //*
	  FILE *file_params_tmp1;                      //*
	  char string[400];                            //*
	  char char_tmp[400];                          //*
	  int i_ite = 0;                               //*
	  int j_ite = 0;                               //*
	#endif                                         //*
//**************************************************** 

  // new functions to update the digits dynamically
  void update_digit() {
	  // if using binary select point then skip since it will automatically change
	  // if evenly distribution and training
	  if (!useBinarySelectedPoint && !caltyp) {
		  digAftdec = 1;
	  }
	  // else if evenly distribution and evaluation
	  else if (!useBinarySelectedPoint && caltyp) {
		  digAftdec = 3;
	  }
	  //else if binary and evaluation
	  else if (useBinarySelectedPoint && caltyp) {
		  digAftdec = 3;
	  }
	  else {
		  digAftdec = digAftdec;
	  }
  }

int main () {

	/****************************************************************************
	| Setting for ISAT table
	|
	| Note: Setting in nameList file isat_1.nml is preferable if file exists
	****************************************************************************/

	//Note: C is zero based and Fortran is one based. Use manual num-1

	//Basic setting
	info[1] = 0;                      /*info(2) if_g     = 0, no need gradient
									  = 1, use gradient*/
	info[11] = 2;                     /*info(12) isatop  = 1, write stats array in .op file for isat performance output on node 0
									  = 2, for isat performance output on all nodes*/

	info[21] = 1;                     /*info(22) istats  = 0, stats and fa is not returned
									  = 1, return stats as well as fa*/
	// may be overwritten
	rinfo[0] = 0.4;                   //rinfo(1) errtol, Note: total error
									  //Need muanally add a factor in f(x) if we want different accuray in myltiple out puts
									  //Need double check

									  //Checkpointing
	info[9] = 0;                      /*info(10) ichin    = 0, the ISAT table is created from scratch
									  = 1,Read previous table, required tab file existed and identical case condition*/

	info[10] = 2;                     /*info(11) ichout   = 1, To checkpoint the table occasionally on node 0
									  = 2, to checkpoint the table occasionally on all nodes*/

	info[19] = 1;                     /*info(20) stats0   = 0, for ichin=1, continue from  stats  read from isat_#_P.tab
									  = 1  for ichin=1, re-initialize  stats*/

	/****************************************************************************
	| Initialize log file and executive time
	****************************************************************************/
	tExeStart = time(NULL);
	sprintf(logMsg, "Executable start: %s", ctime(&tExeStart));
	mpc_log(logMsg, MPC_NEW);


	/****************************************************************************
	| overwrite the parameters by user settings
	****************************************************************************/
	if ((file_params_tmp = fopen("mpc.dat", "r")) == NULL) {
		return 1;
	}
	mpc_log("------------------------------------------------------------", MPC_NORMAL);
	while (fgets(string, 400, file_params_tmp) != NULL) {
		if (EOF == sscanf(string, "%s", char_tmp)) {
			return 0;
		}

		if (!strcmp(char_tmp, "Err_Global")) {
			sscanf(string, "%s%lf", char_tmp, &rinfo[0]);
			printf("%s is %f\n", char_tmp, rinfo[0]);
			sprintf(logMsg, "Total Error: %f", rinfo[0]);
			mpc_log(logMsg, MPC_NORMAL);
		}
		else if (!strcmp(char_tmp, "EvaSize")) {
			sscanf(string, "%s%d", char_tmp, &nRanCall);
			printf("%s is %d\n", char_tmp, nRanCall);
			sprintf(logMsg, "Total Evaluation Size: %d", nRanCall);
			mpc_log(logMsg, MPC_NORMAL);
		}
		else if (!strcmp(char_tmp, "VPower")) {
			sscanf(string, "%s%f", char_tmp, &VPower);
			printf("%s is %f\n", char_tmp, VPower);
			sprintf(logMsg, "The coefficient to velocity: %f", VPower);
			mpc_log(logMsg, MPC_NORMAL);
		}
		else if (!strcmp(char_tmp, "ita")) {
			sscanf(string, "%s%f", char_tmp, &ita);
			printf("%s is %f\n", char_tmp, ita);
			sprintf(logMsg, "The ita is: %f", ita);
			mpc_log(logMsg, MPC_NORMAL);
		}
		else if (!strcmp(char_tmp, "BinaryTreeTra")) {
			sscanf(string, "%s%d", char_tmp, &useBinarySelectedPoint);
			printf("%s is %d\n", char_tmp, useBinarySelectedPoint);
			if (useBinarySelectedPoint) {
				sprintf(logMsg, "Train method is: BinaryTree");
				mpc_log(logMsg, MPC_NORMAL);
			}
			else {
				sprintf(logMsg, "Train method is: Evenly Distribution");
				mpc_log(logMsg, MPC_NORMAL);
			}

		}
		else
			continue;
	}
	fclose(file_params_tmp);
	mpc_log("------------------------------------------------------------", MPC_NORMAL);

    #ifdef MANUAL_TEST
	/*For Test Purpose*/
	// allocate memory
	input1 = (double *)calloc(nRanCall, sizeof(double));
	input2 = (double *)calloc(nRanCall, sizeof(double));
	// read from evaluate file
	if ((file_params_tmp = fopen("input.txt", "r")) == NULL) {
		return 1;
	}
	for (i_ite = 0;i_ite < nRanCall;i_ite++) {
		fgets(string, 400, file_params_tmp);
		sscanf(string, "%lf%lf", &input1[i_ite], &input2[i_ite]);
		printf("Input 1 and input 2 is %f,%f\n", input1[i_ite], input2[i_ite]);
	}
	fclose(file_params_tmp);
	#endif
  
	//read the existing ffd result file
	if(READ_FFD_RESULT) read_existing();

  /****************************************************************************
  | Calculate nx and nf
  ****************************************************************************/
  nx = sizeof(x)/sizeof(double);
  nf = sizeof(fa)/sizeof(double);

  //-----------------------------Old Example -----------------------------------
  //intitialize struct
  // Assume this is the struct that FDD need
  // Here is a old example for input and output. Struct should used for constant.

  //ffdStruct.nx = &nx;
  //ffdStruct.nf = &nf;
  //ffdStruct.x = (double *)&x;
  //ffdStruct.f = (double *)&fa;
  ////ffdStruct.g = (double *)&ga;
  //---------------------------------------------------------------------------

  /****************************************************************************
  | Update non-rectangular input domian
  ****************************************************************************/
  if (useBoundaryCenterRange){
    //Need be fixed
  }
  /****************************************************************************
  | Table Preparation
  ****************************************************************************/
  if (RoundDigitsControl) update_digit();

  //Check useISAT to protect mistaken setting
  if (useTablePreparation && useISAT && !useBinarySelectedPoint)
    nDemArrEva(nx,xStep);

  /****************************************************************************
  | Use Binary Selected Point Training
  ****************************************************************************/
  if (useBinarySelectedPoint)
      binaryTrain ();

  /****************************************************************************
  | Random Test
  ****************************************************************************/
  if (useRandomTest){
    caltyp = 0;
	if (RoundDigitsControl) update_digit();
    //Initialize random function
    srand (time(NULL));

    randomCall(nRanCall, useNormalDistribution, useUnboundedTest);
  }
  else {
	  caltyp = 1;
	  if (RoundDigitsControl) update_digit();
	  //Initialize random function
	  srand(time(NULL));
	  for (i_ite = 0;i_ite < nRanCall;i_ite++) {
		  x[0] = input1[i_ite];
		  x[1] = input2[i_ite];
		  evaluate();
	  }//end of for
  }

  /****************************************************************************
  | Post precess of log
  ****************************************************************************/
  //Accuracy test result
  if (useAccuracyTest){
    sprintf(logMsg, "errSum: %f", errSum);
    mpc_log(logMsg, MPC_NORMAL);
    sprintf(logMsg, "errMax: %f", errMax);
    mpc_log(logMsg, MPC_NORMAL);
    sprintf(logMsg, "errAve: %f", errSum/nQuery);
    mpc_log(logMsg, MPC_NORMAL);
  }

  //Statistical data from ISAT
  if (info[21]){
    sprintf(logMsg, "stats 12345: \n%f \n%f \n%f \n%f \n%f", stats[0],stats[1],stats[2],stats[3],stats[4]);
    mpc_log(logMsg, MPC_NORMAL);
  }

  //Other statistical data
  sprintf(logMsg, "nQuery: %d", nQuery);
  mpc_log(logMsg, MPC_NORMAL);
  sprintf(logMsg, "cpuISAT: %f", cpuISAT);
  mpc_log(logMsg, MPC_NORMAL);
  sprintf(logMsg, "cpuDE: %f", cpuDE);
  mpc_log(logMsg, MPC_NORMAL);

  // Executive time and Ending
  tExeEnd = time(NULL);
  sprintf(logMsg, "\nExecutive time: %f seconds", difftime(tExeEnd,tExeStart));
  mpc_log(logMsg, MPC_NORMAL);

  sprintf(logMsg, "Executable end: %s", ctime(&tExeEnd));
  mpc_log(logMsg, MPC_NORMAL);

  //printf("End");
  //getchar();
  return 0;

} // End of Main


///////////////////////////////////////////////////////////////////////////////
/// Prepara table by calling ISATAB
///
/// Note: Table will prepared with xBoundary and xStep
///       Calling with input nx
///\param Demantion number of current recursion
///
///\return No return needed

//Fixme: Update 
///////////////////////////////////////////////////////////////////////////////
void nDemArrEva(int dimension, double xStep[]) {

  int i = dimension - 1;
  double currentX=0.0;

  if (i == -1){
    //Check hash table
    findHash();
    if (!p){
      addHash();
      evaluate();
    }
  } 
  else {

    for (currentX = xBoundary[i][0]; currentX <= xBoundary[i][1]+ISAT_SMALL; currentX += xStep[i]) {
      x[i] = currentX;
	  //printf("currentX[%d]\tbond[%d][1] is %f\t%f\n", i,i,currentX, xBoundary[i][1]);
      nDemArrEva(i,xStep);
    }

  } //End of if
} // End of iterate



///////////////////////////////////////////////////////////////////////////////
/// Perform a evaluation with current input in global variable
///
/// Note: Use solver or ISAT depends on useISAT
///
///////////////////////////////////////////////////////////////////////////////
void evaluate (){
  int i;
  //check if use rounded input
  if (useRoundInput){
    for (i=0; i<nx; i++){
      x[i] = my_round(x[i],digAftdec);
    }
  }

  if (useISAT){
    tStart = clock();
	//printf("start time in C is %.4f\n", tStart);
    ISATAB(&idtab, &mode, &nx, x, &nf, &nh, &nhd, (void *)&unusedPointer, (int *)&ffdStruct, rusr, info, rinfo, fa, ga, ha, stats);
	tEnd = clock();
	cpuCum = (double)(tEnd - tStart) / CLOCKS_PER_SEC;
	cpuISAT += cpuCum;
	//printf("ISAT VS Manual: %.4f\t%.4f\n", stats[81], cpuISAT);
	//getchar();

  } else {
    tStart = clock();
    ffd_ISAT(need, x, fa, ga, (void *)&ffdStruct);
    cpuDE += (clock() - tStart) / CLOCKS_PER_SEC;
  }

  //increase nQuery
  nQuery++;

  if (useISAT && useAccuracyTest && caltyp)
    accuracyTest();
  else
    writeRecord ();
  
  //update time_elipse
  time_elipse = cpuISAT;

}//end of evaluate

 ///////////////////////////////////////////////////////////////////////////////
 /// Perform a direct evaluation of f and record accuracy
 ///
 /// Note: This function is only called by evaluate()
 ///////////////////////////////////////////////////////////////////////////////
void accuracyTest() {
	int i;
	double fDE[nf_SIZE] = { 0 };
	double dotProduct = 0;
	double err = 0;

	//Direct evaluate
	tStart = clock();
	ffd_ISAT(need, x, fDE, ga, (void *)&ffdStruct);
	cpuDE += (clock() - tStart) / CLOCKS_PER_SEC;

	//Record accuracy

	//Use vector length to represent error: square root of the dot product
	for (i = 0; i < nf; i++)
		dotProduct += (fa[i] - fDE[i]) * (fa[i] - fDE[i]);

	err = sqrt(dotProduct);

	//Update sum and max
	errSum += err;

	if (err > errMax)
		errMax = err;
	/*
	!  stats( 1) - queries:    total number of queries
	!  stats( 2) - p_ret:	   total number of primary retrieves
	!  stats( 3) - s_ret:	   total number of secondary retrieves
	!  stats( 4) - gro:        total number of queries resulting in grows
	!  stats( 5) - add:        total number of queries resulting in adds
	!  stats( 9) - last:       last action, 2=primary retrieve, 4=grow, etc.
	*/

	sprintf(logMsg, "nQuery, Input[0], Input[1], outPut[0], Err, stats[1-5], isat_cum,cpu_cum:\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f", nQuery, x[0], x[1], fa[0], err, stats[0], stats[1], stats[2], stats[3], stats[4], stats[8], stats[81], cpuISAT);
	mpc_log(logMsg, MPC_NORMAL);
}//end of accurateTest


///////////////////////////////////////////////////////////////////////////////
/// Perform a evaluation with random input within boundary of x
///
/// Note: Use solver or ISAT depends on useISAT
///
///\param number of random calls
///
///\return No return needed
///////////////////////////////////////////////////////////////////////////////
void randomCall (int nCall, int useNormalDistribution, int useUnboundedTest){

  int i, j, r;

  ////Initialize random function
  //srand (time(NULL));

  for (i = 0; i < nCall; i++) {
    //Prepare input
    for (j = 0; j < nx; j++) {
      //r = rand() %1000;
      //x[j] = xBoundary[j][0] + (xBoundary[j][1]-xBoundary[j][0]) * r/1000;
      x[j] = getRandom(j,useNormalDistribution,useUnboundedTest);
    }
    ////Check centered range 2D
    //if (useBoundaryCenterRange && (abs(x[0]-x[1]) > xBoundaryCenterRange/2)){
    //  i--;
    //  continue;
    //}

    //Evaluate
    evaluate ();
  }

}



///////////////////////////////////////////////////////////////////////////////
/// Fixme : write description
///////////////////////////////////////////////////////////////////////////////
void writeRecord (){
	if (fabs(cpuISAT - time_elipse) > 1.000000/*threshold for the time difference*/) {
  sprintf(logMsg, "nQuery, Input[0], Input[1], outPut[0], stats[1-5], isat_cum, cpu_cum:\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f", nQuery,x[0],x[1],fa[0],stats[0],stats[1],stats[2],stats[3],stats[4],stats[8],stats[81],cpuISAT );
  mpc_log(logMsg, MPC_NORMAL);
	}
}


///////////////////////////////////////////////////////////////////////////////
/// Fixme : write description
///////////////////////////////////////////////////////////////////////////////
double getRandom (int dimension, int useNormalDistribution, int useUnboundedTest){

  //double (*xBound)[nx_SIZE][2];
  double r;

  ////Check if use bounded test
  //if (useUnboundedTest)
    //xBound = &xBoundary2;
  //else
  //  xBound = &xBoundary;

  if (!useNormalDistribution){
    r = rand() %1000;
    if (useUnboundedTest)
      return xBoundary2[dimension][0] + (xBoundary2[dimension][1]-xBoundary2[dimension][0]) * r/1000;
    else
      return xBoundary[dimension][0] + (xBoundary[dimension][1]-xBoundary[dimension][0]) * r/1000;
  } 
  else {
    if (useUnboundedTest){
      return randNormal (mu,sigma2);
    }
    else{
      r = randNormal (mu,sigma);
      while (r < xBoundary[dimension][0] || r > xBoundary[dimension][1]){
        r = randNormal (mu,sigma);
      }
      return r;

    }// end of if (useUnboundedTest)
  }// end of if (!useNormalDistribution)

}

///////////////////////////////////////////////////////////////////////////////
/// Fixme : write description
/// codes source: https://phoxis.org/2013/05/04/generating-random-numbers-from-normal-distribution-in-c/
/// Polar method: https://en.wikipedia.org/wiki/Marsaglia_polar_method
/// Box-method: https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
/// normal distribution: https://en.wikipedia.org/wiki/Normal_distribution#Generating_values_from_normal_distribution

///////////////////////////////////////////////////////////////////////////////
double randNormal (double mu, double sigma) {
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;
 
  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }
 
  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);
 
  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;
 
  call = !call;
 
  return (mu + sigma * (double) X1);
}


///////////////////////////////////////////////////////////////////////////////
// Fixme : write description
///////////////////////////////////////////////////////////////////////////////
void binaryTrain () {
  double sumOfGroAdd= 0;
  int i = 0,j = 0;
  int divide =1;
  double xStep[nx_SIZE];
  
  //initialize the digAftdec
  if (RoundDigitsControl) digAftdec = 0;

  //Round loop of evaluation
  do {
	//dynamically change the digits
	if (RoundDigitsControl) digAftdec++;

    sumOfGroAdd  = stats[3]+stats[4];

    //update xStep
    for (i=0; i<nx; i++){
      xStep[i] = (xBoundary[i][1]- xBoundary[i][0])/divide;

    }

    // evaluation with current divide
    nDemArrEva(nx,xStep);

    divide = divide*2;

  //} while (sumOfGroAdd != (stats[3]+stats[4]) /*&& nQuery <=QueryLimit*/);
  } while (sumOfGroAdd + nQuery * 5 / 1000 <= (stats[3] + stats[4]) /*&& nQuery <=QueryLimit*/);

}

///////////////////////////////////////////////////////////////////////////////
// Fixme : write description
///////////////////////////////////////////////////////////////////////////////
double my_round(double x, int digits) {
  double d,fac;

  d= digits;
  fac= pow(10, d);

  return floor(x*fac)/fac;
}

///////////////////////////////////////////////////////////////////////////////
// Fixme : write description
///////////////////////////////////////////////////////////////////////////////
void findHash() {
  int i;

  memset(&l, 0, sizeof(hashStruct));

  for (i=0;i < nx; i++)
    l.key.x[i]  = x[i];

  HASH_FIND(hh, records, &l.key, sizeof(hashKey), p);

}

///////////////////////////////////////////////////////////////////////////////
// Fixme : write description
///////////////////////////////////////////////////////////////////////////////
void addHash() {
  int i;

  r = (hashStruct*)malloc( sizeof(hashStruct) );
  memset(r, 0, sizeof(hashStruct));

  for (i=0;i < nx; i++)
    r->key.x[i]  = x[i];

  HASH_ADD(hh, records, key, sizeof(hashKey), r);

}