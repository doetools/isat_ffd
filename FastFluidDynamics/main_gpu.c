#define CL_USE_DEPRECATED_OPENCL_2_0_APIS
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL\cl.h>
#endif

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#ifdef _msc_ver
#include <windows.h>
#else
//#include <unistd.h>
#endif

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
//#include "glut.h"

#include "main_gpu.h"

#ifdef FFD_ISAT   //Called by ISAT
//Extern global variables form ISAT
extern double ffdOutput[];         // centerT, centerUx, centerUy, centerUz
#endif

//char *msg;
//using namespace std;

/* global variables */
REAL **var;
REAL *var_flat; // flattened var from two dimensions to one dimension
int  **BINDEX;
int *bindex_flat;
REAL *locmin, *locmax;
static PARA_DATA para;
static PARA_DATA_SIMP para_simp;
static GEOM_DATA geom;
static PROB_DATA prob;
static TIME_DATA mytime;
static INPU_DATA inpu;
static OUTP_DATA outp1;
static BC_DATA bc;
static SOLV_DATA solv;
static SENSOR_DATA sens;
static INIT_DATA init;

clock_t start, end;

void mypause(int x)
{
  if (x != 1){
    printf("Kernel Will Not Run... hit [ENTER] to Exit . . .");
    fflush(stdout);
    getchar();
    exit(0);
  }
  else {
    printf("Program Build Succesful To Contine with Kernel Execution hit [ENTER] . . .");
    fflush(stdout);
    getchar();
  }
}

int allocate_memory(PARA_DATA *para) {

  int nb_var, i;
  int size = (geom.imax + 2) * (geom.jmax + 2) * (geom.kmax + 2);
  int flat_size_var, flat_size_bindex;
  /****************************************************************************
  | Allocate memory for variables *var
  ****************************************************************************/
  

  /****************************************************************************
  | Allocate memory for variables **var
  ****************************************************************************/
  nb_var = C2BC + 1;
  var = (REAL **)malloc(nb_var*sizeof(REAL*));
  if (var == NULL) {
    ffd_log("allocate_memory(): Could not allocate memory for var.",
      FFD_ERROR);
    return 1;
  }

  for (i = 0; i<nb_var; i++) {
    var[i] = (REAL *)calloc(size, sizeof(REAL));
    if (var[i] == NULL) {
      sprintf(msg,
        "allocate_memory(): Could not allocate memory for var[%d]", i);
      ffd_log(msg, FFD_ERROR);
      return 1;
    }
  }

  /****************************************************************************
  | Allocate memory for variables *var
  ****************************************************************************/
  flat_size_var = nb_var*size;
  var_flat = (REAL *)malloc(flat_size_var*sizeof(REAL));
  if (var_flat == NULL) {
    ffd_log("allocate_memory(): Could not allocate memory for var_flat.",
      FFD_ERROR);
    return 1;
  }

  /****************************************************************************
  | Allocate memory for boundary cells
  | BINDEX[0]: i of global coordinate in IX(i,j,k)
  | BINDEX[1]: j of global coordinate in IX(i,j,k)
  | BINDEX[2]: k of global coordinate in IX(i,j,k)
  | BINDEX[3]: Fixed temperature or fixed heat flux
  | BINDEX[4]: Boundary ID to identify which boundary it belongs to
  ****************************************************************************/
  BINDEX = (int **)malloc(5 * sizeof(int*));
  if (BINDEX == NULL) {
    ffd_log("allocate_memory(): Could not allocate memory for BINDEX.",
      FFD_ERROR);
    return 1;
  }

  for (i = 0; i<5; i++) {
    BINDEX[i] = (int *)malloc(size*sizeof(int));
    if (BINDEX[i] == NULL) {
      sprintf(msg,
        "allocate_memory(): Could not allocate memory for BINDEX[%d]", i);
      ffd_log(msg, FFD_ERROR);
      return 1;
    }
  }

  /****************************************************************************
  | Allocate memory for variables **var
  ****************************************************************************/
  int row = 5;
  flat_size_bindex = row*size;
  bindex_flat = (int *)malloc(flat_size_bindex*sizeof(int));
  if (bindex_flat == NULL) {
    ffd_log("allocate_memory(): Could not allocate memory for bindex_flat.",
      FFD_ERROR);
    return 1;
  }
  return 0;
} // End of allocate_memory()


///////////////////////////////////////////////////////////////////////////////
/// Main routine of FFD
///
///\para coupled simulation Integer to identify the simulation type
///
///\return 0 if no error occurred
///////////////////////////////////////////////////////////////////////////////
int ffd_prep(int cosimulation) {

  // Initialize the parameters para
  para.geom = &geom;
  para.inpu = &inpu;
  para.outp = &outp1;
  para.prob = &prob;
  para.mytime = &mytime;
  para.bc = &bc;
  para.solv = &solv;
  para.sens = &sens;
  para.init = &init;

  // Initialize the parameters para

  if (initialize(&para) != 0) {
    ffd_log("ffd(): Could not initialize simulation parameters.", FFD_ERROR);
    return 1;
  }

  // Overwrite the mesh and simulation data using SCI generated file
  if (para.inpu->parameter_file_format == SCI) {
    if (read_sci_max(&para, var) != 0) {
      ffd_log("ffd(): Could not read SCI data.", FFD_ERROR);
      return 1;
    }
  }

  // Allocate memory for the variables
  if (allocate_memory(&para) != 0) {
    ffd_log("ffd(): Could not allocate memory for the simulation.", FFD_ERROR);
    return 1;
  }

  // Set the initial values for the simulation data
  if (set_initial_data(&para, var, BINDEX)) {
    ffd_log("ffd(): Could not set initial data.", FFD_ERROR);
    return 1;
  }
  ffd_log("ffd(): successfully initialize virables.", FFD_NORMAL);
  // Read previous simulation data as initial values
  //if (para.inpu->read_old_ffd_file == 1) read_ffd_data(&para, var);

  ffd_log("ffd.c: Start FFD solver.", FFD_NORMAL);
  //write_tecplot_data(&para, var, "initial");
  
  //initialize the simplified para_simp
  if (init_para_simp(&para, &para_simp) != 0) {
    ffd_log("ffd.c: fail in initializing para_simp", FFD_ERROR);
  }
    // flatten var
  flat_var(&para, var, var_flat);

  // flatten bindex
  flat_index(&para, BINDEX, bindex_flat);

  return 0;
} // End of ffd( )


int ffd()
{
  int bar = 0;
  ffd_log("Start FFD simulation on GPU", FFD_NEW);
  int size ;
  int nb_var = C2BC + 1;
  /*define the index*/
  int JACO_IT = 0; // for control for iteration
  int IT_MAX = 30; // 20 iteration for jacobian solver
  int NEXT = 1; // for control of simulation
  int VX_IND[1] = { VX };
  int VY_IND[1] = { VY };
  int VZ_IND[1] = { VZ };
  int TEMP_IND[1] = { TEMP };
  int PRE_IND[1] = { IP };
  int ADVE_IND[1] = { 0 };
  int DIFF_IND[1] = { 1 };
  /*variable for OPENCL*/
  size_t global_work_size[3];
  size_t global_work_size_bc[1];
  size_t global_work_size_mb[1];
  size_t local_group_size[3];
  size_t local_group_size_bc[1];

  /*OPENCL MEMORY ALLOCATIONS*/ 
  
  cl_kernel kernel[25] = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
    NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL , NULL}; //<--- DEFINE HERE THE NUMBER OF KERNELS THAT ARE FOUND IN THE .CL FILE

  /*CPU MEMORY ALLOCATIONS*/
  if (ffd_prep(0) != 0) {
    ffd_log("main(): Error in launching the FFD simulation", FFD_ERROR);
  }



  sprintf(msg, "main(): start getting GPU information");
  ffd_log(msg, FFD_NORMAL);

  size = (geom.imax + 2) * (geom.jmax + 2) * (geom.kmax + 2);

  /*Step1: Getting platforms and choose an available one.*/
  cl_uint numPlatforms;  //the NO. of platforms
  cl_platform_id platform = NULL;  //the chosen platform
  sprintf(msg, "main(): start getting platform information");
  ffd_log(msg, FFD_NORMAL);
  cl_int  status = clGetPlatformIDs(0, NULL, &numPlatforms);
  if (status != CL_SUCCESS)
  {
    sprintf(msg, "Error: Getting platforms!");
    ffd_log(msg, FFD_ERROR);
    return FAILURE;
  }

  /*For clarity, choose the first available platform. */
  if (numPlatforms > 0)
  {
    cl_platform_id* platforms = (cl_platform_id*)malloc(numPlatforms* sizeof(cl_platform_id));
    status = clGetPlatformIDs(numPlatforms, platforms, NULL);
    platform = platforms[0];
    free(platforms);
  }
  sprintf(msg, "main(): finishing checking avaiblable platforms");
  ffd_log(msg, FFD_NORMAL);
  /*Step 2:Query the platform and choose the first GPU device if has one.Otherwise use the CPU as device.*/
  cl_uint        numDevices = 0;
  cl_device_id        devices[1] = { NULL };
  cl_device_id        *devices_GPU = NULL;
  cl_device_id        *devices_CPU = NULL;
  status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 0, NULL, &numDevices);

  /*if (numDevices == 0)  //no GPU available.
  {
    printf("%s", "No GPU device available......");
    printf("%s", "CPU is now set as default device.\n");
    status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU, 0, NULL, &numDevices);
    devices = (cl_device_id*)malloc(numDevices * sizeof(cl_device_id));
    status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU, numDevices, devices, NULL);
  }
  else
  {
    devices = (cl_device_id*)malloc(numDevices * sizeof(cl_device_id));
    status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, numDevices, devices, NULL);
    //status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU, numDevices, devices, NULL);
    printf("%s", "GPU Obtained.....");
  }
  */

  //printf("To use GPU press [0] for CPU press [1]...\n");
  //scanf("%d",&x);

  if (numDevices==0)  //no GPU available.
  {
    //printf("%s", "CPU is now set as default device.\n");
    //status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU, 0, NULL, &numDevices);
	status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU, 0, NULL, &numDevices);
    devices_CPU = (cl_device_id*)malloc(numDevices * sizeof(cl_device_id));
    status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU, numDevices, devices_CPU, NULL);
  }
  else 
  {
    devices_GPU = (cl_device_id*)malloc((numDevices) * sizeof(cl_device_id));
    status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, numDevices, devices_GPU, NULL);
    //status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU, numDevices, devices, NULL);
    //printf("%s", "GPU Obtained.....");
  }

  /*menue*/
  // 0: old gpu
  // 1: new gpu
  // 2: cpu
  /*menu*/
  int x = 1;

  // old gpu
  if (x == 0)  
	  devices[0] = devices_GPU[0];
  // new gpu
  else if (x == 1)
	  devices[0] = devices_GPU[1];
  // new cpu
  else if (x == 2)
	  devices[0] = devices_CPU[0];
  

  sprintf(msg, "main(): finishing chosing the avaiblable platform");
  ffd_log(msg, FFD_NORMAL);

  /*Step 3: Create context.*/
  cl_context context = clCreateContext(NULL, 1, devices, NULL, NULL, NULL);
  sprintf(msg, "main(): finishing creating context for the platform");
  ffd_log(msg, FFD_NORMAL);

  /*Step 4: Creating command queue associate with the context.*/
  cl_command_queue commandQueue = clCreateCommandQueue(context, devices[0], NULL, NULL);
  sprintf(msg, "main(): finishing creating queue associated with the context");
  ffd_log(msg, FFD_NORMAL);

  /*Create kernel from source*/
  sprintf(msg, "main(): start creating the Write buffers associated with the kernel arguments");
  ffd_log(msg, FFD_NORMAL);
  //LOADING THE KERNEL FILE INTO SYSTEM MEMMORY
  
  
  // Program Setup
  
  sprintf(msg, "main(): start cl kernel file reading");
  ffd_log(msg, FFD_NORMAL);
  FILE *fp;

  const char filename[] = "./Kernels_3D.cl"; //<--- NAME OF THE .CL FILE USED INSIDE OF THE PROJECT 
  size_t source_size;

  char *source_str;
  fp = fopen(filename, "r");
  if (!fp)
    {
    printf("Failed to Load Kernel.\n");
    return FAILURE;
    }
  source_str = (char*)malloc(MAX_SOURCE_SIZE);
  source_size = fread(source_str, 1, MAX_SOURCE_SIZE, fp) ;
  fclose(fp);

  sprintf(msg, "main(): finish cl kernel file reading");
  ffd_log(msg, FFD_NORMAL);

  cl_program program = clCreateProgramWithSource(context, 1, (const char **)&source_str, (const size_t *)&source_size, NULL);

  status = clBuildProgram(program, 1, &devices[0], "-I ./", NULL, NULL);
  
  //free the allocated source_str
  free(source_str);

  //Debug Information for the kernel

  if (status != CL_SUCCESS){
    printf("\nFailed to build the program...Debug Information below:\n");

    size_t len;
    char *buffer;
    buffer = (char*)malloc(MAX_SOURCE_SIZE);
    clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_LOG, 4048*100 * sizeof(char), buffer, &len);
    printf("\nPROGRAM_BUILD_LOG\n");
    printf("%s\n", buffer);
    ffd_log(buffer, FFD_NORMAL);
    clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_STATUS, 4048 * 20 * sizeof(char), buffer, &len);
    printf("\nPROGRAM_BUILD_STATUS\n");
    printf("%s\n", buffer);
    clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_OPTIONS, 4048 * 20 * sizeof(char), buffer, &len);
    printf("\nPROGRAM_BUILD_OPTIONS\n");
    printf("%s\n", buffer);
    mypause(0);
	free(buffer);
  }
    else if (status == CL_SUCCESS) {
      //printf("\nProgram Build SUCCESS...Debug Information below:\n");
      size_t len;
      char *buffer;
      buffer = (char*)malloc(MAX_SOURCE_SIZE);
      clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_LOG, 4048 * 100 * sizeof(char), buffer, &len);
      //printf("\nCL_PROGRAM_BUILD_LOG\n");
      //printf("%s\n", buffer);
      clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_STATUS, 4048 * 20 * sizeof(char), buffer, &len);
      //printf("\nCL_PROGRAM_BUILD_STATUS\n");
      //printf("%s\n", buffer);
      clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_OPTIONS, 4048 * 20 * sizeof(char), buffer, &len);
      //printf("\nCL_PROGRAM_BUILD_OPTIONS\n");
      //printf("%s\n", buffer);
      //mypause(1);
	  free(buffer);
    }
	

  /*Create parallel OpenCL kernel*/
  sprintf(msg, "main(): start creating kernels executable");
  ffd_log(msg, FFD_NORMAL);
  status = clFlush(commandQueue);
  kernel[0] = clCreateKernel(program, "adve_VX", NULL); //<---- advection for U
  kernel[1] = clCreateKernel(program, "adve_VY", NULL); //<---- advection for V
  kernel[2] = clCreateKernel(program, "adve_VZ", NULL); //<---- advection for W
  kernel[3] = clCreateKernel(program, "adve_T", NULL); //<---- advection for T
  kernel[4] = clCreateKernel(program, "diff_VX", NULL); //<---- diffusion for U
  kernel[5] = clCreateKernel(program, "diff_VY", NULL); //<---- diffusion for V
  kernel[6] = clCreateKernel(program, "diff_VZ", NULL); //<---- diffusion for W
  kernel[7] = clCreateKernel(program, "diff_T", NULL); //<---- diffusion for T
  kernel[8] = clCreateKernel(program, "project", NULL); //<---- projection 
  kernel[9] = clCreateKernel(program, "project_velo_corr", NULL); //<---- velocity correction after projection
  kernel[10] = clCreateKernel(program, "Ax_BSolver", NULL); //<---- jacobian solver for other variables except for pressure
  kernel[11] = clCreateKernel(program, "Ax_BSolver_P", NULL); //<---- jacobian solver for pressure projection
  kernel[12] = clCreateKernel(program, "ap_coeff", NULL); //<---- set coefficients for AP
  kernel[13] = clCreateKernel(program, "set_bnd_T", NULL); //<---- T BC
  kernel[14] = clCreateKernel(program, "set_bnd_pressure", NULL); //<---- P BC
  kernel[15] = clCreateKernel(program, "set_bnd_VX", NULL); //<---- U BC
  kernel[16] = clCreateKernel(program, "set_bnd_VY", NULL); //<---- V BC
  kernel[17] = clCreateKernel(program, "set_bnd_VZ", NULL); //<---- W BC
  kernel[18] = clCreateKernel(program, "adjust_velocity", NULL); //<---- adjust velocity
  kernel[19] = clCreateKernel(program, "mass_conservation", NULL); //<---- mass balance mass_conservation
  kernel[20] = clCreateKernel(program, "reset_time_averaged_data", NULL); //<---- adjust velocity
  kernel[21] = clCreateKernel(program, "add_time_averaged_data", NULL); //<---- mass balance mass_conservation
  kernel[22] = clCreateKernel(program, "time_averaged", NULL); //<---- adjust velocity
  kernel[23] = clCreateKernel(program, "Ax_BSolver_upd", NULL); //<---- mass balance mass_conservation
  kernel[24] = clCreateKernel(program, "Ax_BSolver_P_upd", NULL); //<---- adjust velocity
  status = clFinish(commandQueue);

  /*GPU MEMORY CONTROL*/
  //Create buffer Object
  sprintf(msg, "main(): start creating the buffer objects");
  ffd_log(msg, FFD_NORMAL);
  status = clFlush(commandQueue);
  cl_mem var_mobj = clCreateBuffer(context, CL_MEM_READ_WRITE, nb_var*size*sizeof(REAL), NULL, &status);  //<--- var
  cl_mem para_mobj = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(PARA_DATA_SIMP), NULL, &status); //<--- Para
  cl_mem bindex_mobj = clCreateBuffer(context, CL_MEM_READ_WRITE, 5*size*sizeof(int), NULL, &status); //<--- BINDEX
  cl_mem vx_ind_mobj = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(int), NULL, &status); //<--- VX index
  cl_mem vy_ind_mobj = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(int), NULL, &status); //<--- VY index
  cl_mem vz_ind_mobj = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(int), NULL, &status); //<--- VZ index
  cl_mem T_ind_mobj = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(int), NULL, &status); //<--- T index
  cl_mem P_ind_mobj = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(int), NULL, &status); //<--- Pressure index
  cl_mem ADVE_ind_mobj = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(int), NULL, &status); //<--- advection function index
  cl_mem DIFF_ind_mobj = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(int), NULL, &status); //<--- diffusion fucntion index
  status = clFinish(commandQueue);

  /*Set global group size */
  global_work_size[0] = (para.geom->imax + 2)*sizeof(REAL) / sizeof(REAL);
  global_work_size[1] = (para.geom->jmax + 2)*sizeof(REAL) / sizeof(REAL);
  global_work_size[2] = (para.geom->kmax + 2)*sizeof(REAL) / sizeof(REAL);
  global_work_size_bc[0] = (para.geom->index)*sizeof(REAL) / sizeof(REAL);
  global_work_size_mb[0] = 1;

  local_group_size[0] = 4;
  local_group_size[1] = 4;
  local_group_size[2] = 4;
  local_group_size_bc[0] = 4;
  
  //Write buffer Object
  status = clFlush(commandQueue);
  status = clEnqueueWriteBuffer(commandQueue, var_mobj, CL_TRUE, 0, nb_var*size*sizeof(REAL), var_flat, 0, NULL, NULL); //<---var
  status = clEnqueueWriteBuffer(commandQueue, para_mobj, CL_TRUE, 0, sizeof(PARA_DATA_SIMP), &para_simp, 0, NULL, NULL); //<---para
  status = clEnqueueWriteBuffer(commandQueue, bindex_mobj, CL_TRUE, 0, 5 * size*sizeof(int), bindex_flat, 0, NULL, NULL); //<---BINDEX
  status = clEnqueueWriteBuffer(commandQueue, vx_ind_mobj, CL_TRUE, 0, sizeof(int), VX_IND, 0, NULL, NULL); //<---VX index
  status = clEnqueueWriteBuffer(commandQueue, vy_ind_mobj, CL_TRUE, 0, sizeof(int), VY_IND, 0, NULL, NULL); //<---VY index
  status = clEnqueueWriteBuffer(commandQueue, vz_ind_mobj, CL_TRUE, 0, sizeof(int), VZ_IND, 0, NULL, NULL); //<---VZ index
  status = clEnqueueWriteBuffer(commandQueue, T_ind_mobj, CL_TRUE, 0, sizeof(int), TEMP_IND, 0, NULL, NULL); //<---T index
  status = clEnqueueWriteBuffer(commandQueue, P_ind_mobj, CL_TRUE, 0, sizeof(int), PRE_IND, 0, NULL, NULL); //<---P index
  status = clEnqueueWriteBuffer(commandQueue, ADVE_ind_mobj, CL_TRUE, 0, sizeof(int), ADVE_IND, 0, NULL, NULL); //<---advection function index
  status = clEnqueueWriteBuffer(commandQueue, DIFF_ind_mobj, CL_TRUE, 0, sizeof(int), DIFF_IND, 0, NULL, NULL); //<---diffusion fucntion index
  status = clFinish(commandQueue);

  //for control of add-on of kernels
  int T_ON = 1;
  int P_ON = 1;
  int M_ON = 0;
  int print_control = 1;
  int ifAdv = 1;
  int ifDif = 1;
  int ifDeb = 0;
  /*start of simulation*/
  while (NEXT == 1)
  {
	  /*
	  if (para.mytime->step_current == para.mytime->step_total-1) {
		  ifDeb = 1;
	  }
	  */
	  if (ifAdv == 1) {
		  /*STEP 1.0 :Advection U [kernel 0]*/
		  //Set arguments
		  status = clSetKernelArg(kernel[0], 0, sizeof(cl_mem), (void*)&para_mobj);
		  status = clSetKernelArg(kernel[0], 1, sizeof(cl_mem), (void*)&var_mobj);
		  status = clSetKernelArg(kernel[0], 2, sizeof(cl_mem), (void*)&bindex_mobj);
		  //Run kernel
		  clEnqueueNDRangeKernel(commandQueue, kernel[0], 3, NULL, global_work_size, NULL, 0, NULL, NULL);
		  /*STEP 1.1 :Apply TMP1 BC [kernel 15]*/
		  //Set arguments
		  status = clSetKernelArg(kernel[15], 0, sizeof(cl_mem), (void*)&para_mobj);
		  status = clSetKernelArg(kernel[15], 1, sizeof(cl_mem), (void*)&var_mobj);
		  status = clSetKernelArg(kernel[15], 2, sizeof(cl_mem), (void*)&ADVE_ind_mobj);
		  status = clSetKernelArg(kernel[15], 3, sizeof(cl_mem), (void*)&bindex_mobj);
		  //Run kernel
		  clEnqueueNDRangeKernel(commandQueue, kernel[15], 1, NULL, global_work_size_bc, NULL, 0, NULL, NULL);

		  /*STEP 3.0 :Advection V [kernel 1]*/
		  //Set arguments
		  status = clSetKernelArg(kernel[1], 0, sizeof(cl_mem), (void*)&para_mobj);
		  status = clSetKernelArg(kernel[1], 1, sizeof(cl_mem), (void*)&var_mobj);
		  status = clSetKernelArg(kernel[1], 2, sizeof(cl_mem), (void*)&bindex_mobj);
		  //Run kernel
		  clEnqueueNDRangeKernel(commandQueue, kernel[1], 3, NULL, global_work_size, NULL, 0, NULL, NULL);
		  /*STEP 3.1 :Apply TMP2 BC [kernel 16]*/
		  //Set arguments
		  status = clSetKernelArg(kernel[16], 0, sizeof(cl_mem), (void*)&para_mobj);
		  status = clSetKernelArg(kernel[16], 1, sizeof(cl_mem), (void*)&var_mobj);
		  status = clSetKernelArg(kernel[16], 2, sizeof(cl_mem), (void*)&ADVE_ind_mobj);
		  status = clSetKernelArg(kernel[16], 3, sizeof(cl_mem), (void*)&bindex_mobj);
		  //Run kernel
		  clEnqueueNDRangeKernel(commandQueue, kernel[16], 1, NULL, global_work_size_bc, NULL, 0, NULL, NULL);

		  /*STEP 5.0 :Advection W [kernel 2]*/
		  //Set arguments
		  status = clSetKernelArg(kernel[2], 0, sizeof(cl_mem), (void*)&para_mobj);
		  status = clSetKernelArg(kernel[2], 1, sizeof(cl_mem), (void*)&var_mobj);
		  status = clSetKernelArg(kernel[2], 2, sizeof(cl_mem), (void*)&bindex_mobj);
		  //Run kernel
		  clEnqueueNDRangeKernel(commandQueue, kernel[2], 3, NULL, global_work_size, NULL, 0, NULL, NULL);
		  /*STEP 5.1 :Apply TMP3 BC [kernel 17]*/
		  //Set arguments
		  status = clSetKernelArg(kernel[17], 0, sizeof(cl_mem), (void*)&para_mobj);
		  status = clSetKernelArg(kernel[17], 1, sizeof(cl_mem), (void*)&var_mobj);
		  status = clSetKernelArg(kernel[17], 2, sizeof(cl_mem), (void*)&ADVE_ind_mobj);
		  status = clSetKernelArg(kernel[17], 3, sizeof(cl_mem), (void*)&bindex_mobj);
		  //Run kernel
		  clEnqueueNDRangeKernel(commandQueue, kernel[17], 1, NULL, global_work_size_bc, NULL, 0, NULL, NULL);
	  }
	  if (ifDeb == 1) {
		  status = clEnqueueReadBuffer(commandQueue, var_mobj, CL_TRUE, 0, nb_var*size*sizeof(REAL), var_flat, 0, NULL, NULL);
		  sprintf(msg, "main(): start writing results for advection");
		  unflat_var(&para, var, var_flat);
		  ffd_log(msg, FFD_NORMAL);
		  write_tecplot_data(&para, var, "result_adv");
	  }
	  if (ifDif == 1) {
		  /*STEP 2.0 :Diffusion U [kernel 4]*/
		  //Set arguments
		  status = clSetKernelArg(kernel[4], 0, sizeof(cl_mem), (void*)&para_mobj);
		  status = clSetKernelArg(kernel[4], 1, sizeof(cl_mem), (void*)&var_mobj);
		  status = clSetKernelArg(kernel[4], 2, sizeof(cl_mem), (void*)&bindex_mobj);
		  //Run kernel
		  clEnqueueNDRangeKernel(commandQueue, kernel[4], 3, NULL, global_work_size, NULL, 0, NULL, NULL);
		  /*STEP 2.1 :Apply U BC [kernel 15]*/
		  //Set arguments
		  status = clSetKernelArg(kernel[15], 0, sizeof(cl_mem), (void*)&para_mobj);
		  status = clSetKernelArg(kernel[15], 1, sizeof(cl_mem), (void*)&var_mobj);
		  status = clSetKernelArg(kernel[15], 2, sizeof(cl_mem), (void*)&DIFF_ind_mobj);
		  status = clSetKernelArg(kernel[15], 3, sizeof(cl_mem), (void*)&bindex_mobj);
		  //Run kernel
		  clEnqueueNDRangeKernel(commandQueue, kernel[15], 1, NULL, global_work_size_bc, NULL, 0, NULL, NULL);
		  /*STEP 2.2 :Diffusion U AP [kernel 12]*/
		  //Set arguments
		  status = clSetKernelArg(kernel[12], 0, sizeof(cl_mem), (void*)&para_mobj);
		  status = clSetKernelArg(kernel[12], 1, sizeof(cl_mem), (void*)&var_mobj);
		  status = clSetKernelArg(kernel[12], 2, sizeof(cl_mem), (void*)&vx_ind_mobj);
		  //Run kernel
		  clEnqueueNDRangeKernel(commandQueue, kernel[12], 3, NULL, global_work_size, NULL, 0, NULL, NULL);
		  /*STEP 2.3 :Diffusion U Solver [kernel 10]*/
		  JACO_IT = 0; // initialize index before calculation
		  while (JACO_IT < IT_MAX) {
			  // set the argument for the kenerl
			  status = clSetKernelArg(kernel[10], 0, sizeof(cl_mem), (void*)&para_mobj);
			  status = clSetKernelArg(kernel[10], 1, sizeof(cl_mem), (void*)&var_mobj);
			  status = clSetKernelArg(kernel[10], 2, sizeof(cl_mem), (void*)&vx_ind_mobj);
			  // run the kernel 
			  clEnqueueNDRangeKernel(commandQueue, kernel[10], 3, NULL, global_work_size, NULL, 0, NULL, NULL);

			  status = clSetKernelArg(kernel[23], 0, sizeof(cl_mem), (void*)&para_mobj);
			  status = clSetKernelArg(kernel[23], 1, sizeof(cl_mem), (void*)&var_mobj);
			  status = clSetKernelArg(kernel[23], 2, sizeof(cl_mem), (void*)&vx_ind_mobj);
			  clEnqueueNDRangeKernel(commandQueue, kernel[23], 3, NULL, global_work_size, NULL, 0, NULL, NULL);
			  // iteration marches on
			  JACO_IT++;
		  } //end of while loop for jacobian loop

		  /*STEP 2.4 :Apply U BC [kernel 15]*/
		  //Set arguments
		  status = clSetKernelArg(kernel[15], 0, sizeof(cl_mem), (void*)&para_mobj);
		  status = clSetKernelArg(kernel[15], 1, sizeof(cl_mem), (void*)&var_mobj);
		  status = clSetKernelArg(kernel[15], 2, sizeof(cl_mem), (void*)&DIFF_ind_mobj);
		  status = clSetKernelArg(kernel[15], 3, sizeof(cl_mem), (void*)&bindex_mobj);
		  //Run kernel
		  clEnqueueNDRangeKernel(commandQueue, kernel[15], 1, NULL, global_work_size_bc, NULL, 0, NULL, NULL);

		  /*STEP 4.0 :Diffusion V [kernel 5]*/
		  //Set arguments
		  status = clSetKernelArg(kernel[5], 0, sizeof(cl_mem), (void*)&para_mobj);
		  status = clSetKernelArg(kernel[5], 1, sizeof(cl_mem), (void*)&var_mobj);
		  status = clSetKernelArg(kernel[5], 2, sizeof(cl_mem), (void*)&bindex_mobj);
		  //Run kernel
		  clEnqueueNDRangeKernel(commandQueue, kernel[5], 3, NULL, global_work_size, NULL, 0, NULL, NULL);
		  /*STEP 4.1 :Apply V BC [kernel 16]*/
		  //Set arguments
		  status = clSetKernelArg(kernel[16], 0, sizeof(cl_mem), (void*)&para_mobj);
		  status = clSetKernelArg(kernel[16], 1, sizeof(cl_mem), (void*)&var_mobj);
		  status = clSetKernelArg(kernel[16], 2, sizeof(cl_mem), (void*)&DIFF_ind_mobj);
		  status = clSetKernelArg(kernel[16], 3, sizeof(cl_mem), (void*)&bindex_mobj);
		  //Run kernel
		  clEnqueueNDRangeKernel(commandQueue, kernel[16], 1, NULL, global_work_size_bc, NULL, 0, NULL, NULL);
		  /*STEP 4.2 :Diffusion V AP [kernel 12]*/
		  //Set arguments
		  status = clSetKernelArg(kernel[12], 0, sizeof(cl_mem), (void*)&para_mobj);
		  status = clSetKernelArg(kernel[12], 1, sizeof(cl_mem), (void*)&var_mobj);
		  status = clSetKernelArg(kernel[12], 2, sizeof(cl_mem), (void*)&vy_ind_mobj);
		  //Run kernel
		  clEnqueueNDRangeKernel(commandQueue, kernel[12], 3, NULL, global_work_size, NULL, 0, NULL, NULL);
		  /*STEP 4.3 :Diffusion V Solver [kernel 10]*/
		  JACO_IT = 0; // initialize index before calculation
		  while (JACO_IT < IT_MAX) {
			  // set the argument for the kenerl
			  status = clSetKernelArg(kernel[10], 0, sizeof(cl_mem), (void*)&para_mobj);
			  status = clSetKernelArg(kernel[10], 1, sizeof(cl_mem), (void*)&var_mobj);
			  status = clSetKernelArg(kernel[10], 2, sizeof(cl_mem), (void*)&vy_ind_mobj);
			  // run the kernel for 5 times
			  clEnqueueNDRangeKernel(commandQueue, kernel[10], 3, NULL, global_work_size, NULL, 0, NULL, NULL);

			  status = clSetKernelArg(kernel[23], 0, sizeof(cl_mem), (void*)&para_mobj);
			  status = clSetKernelArg(kernel[23], 1, sizeof(cl_mem), (void*)&var_mobj);
			  status = clSetKernelArg(kernel[23], 2, sizeof(cl_mem), (void*)&vy_ind_mobj);
			  clEnqueueNDRangeKernel(commandQueue, kernel[23], 3, NULL, global_work_size, NULL, 0, NULL, NULL);
			  // iteration marches on
			  JACO_IT++;
		  } //end of while loop for jacobian loop
		  /*STEP 4.4 :Apply V BC [kernel 16]*/
		  //Set arguments
		  status = clSetKernelArg(kernel[16], 0, sizeof(cl_mem), (void*)&para_mobj);
		  status = clSetKernelArg(kernel[16], 1, sizeof(cl_mem), (void*)&var_mobj);
		  status = clSetKernelArg(kernel[16], 2, sizeof(cl_mem), (void*)&DIFF_ind_mobj);
		  status = clSetKernelArg(kernel[16], 3, sizeof(cl_mem), (void*)&bindex_mobj);
		  //Run kernel
		  clEnqueueNDRangeKernel(commandQueue, kernel[16], 1, NULL, global_work_size_bc, NULL, 0, NULL, NULL);


		  /*STEP 6.0 :Diffusion W [kernel 6]*/
		  //Set arguments
		  status = clSetKernelArg(kernel[6], 0, sizeof(cl_mem), (void*)&para_mobj);
		  status = clSetKernelArg(kernel[6], 1, sizeof(cl_mem), (void*)&var_mobj);
		  status = clSetKernelArg(kernel[6], 2, sizeof(cl_mem), (void*)&bindex_mobj);
		  //Run kernel
		  clEnqueueNDRangeKernel(commandQueue, kernel[6], 3, NULL, global_work_size, NULL, 0, NULL, NULL);
		  /*STEP 6.1 :Apply W BC [kernel 17]*/
		  //Set arguments
		  status = clSetKernelArg(kernel[17], 0, sizeof(cl_mem), (void*)&para_mobj);
		  status = clSetKernelArg(kernel[17], 1, sizeof(cl_mem), (void*)&var_mobj);
		  status = clSetKernelArg(kernel[17], 2, sizeof(cl_mem), (void*)&DIFF_ind_mobj);
		  status = clSetKernelArg(kernel[17], 3, sizeof(cl_mem), (void*)&bindex_mobj);
		  //Run kernel
		  clEnqueueNDRangeKernel(commandQueue, kernel[17], 1, NULL, global_work_size_bc, NULL, 0, NULL, NULL);
		  /*STEP 6.2 :Diffusion W AP [kernel 12]*/
		  //Set arguments
		  status = clSetKernelArg(kernel[12], 0, sizeof(cl_mem), (void*)&para_mobj);
		  status = clSetKernelArg(kernel[12], 1, sizeof(cl_mem), (void*)&var_mobj);
		  status = clSetKernelArg(kernel[12], 2, sizeof(cl_mem), (void*)&vz_ind_mobj);
		  //Run kernel
		  clEnqueueNDRangeKernel(commandQueue, kernel[12], 3, NULL, global_work_size, NULL, 0, NULL, NULL);
		  /*STEP 6.3 :Diffusion W Solver [kernel 10]*/
		  JACO_IT = 0; // initialize index before calculation
		  while (JACO_IT < IT_MAX) {
			  // set the argument for the kenerl
			  status = clSetKernelArg(kernel[10], 0, sizeof(cl_mem), (void*)&para_mobj);
			  status = clSetKernelArg(kernel[10], 1, sizeof(cl_mem), (void*)&var_mobj);
			  status = clSetKernelArg(kernel[10], 2, sizeof(cl_mem), (void*)&vz_ind_mobj);

			  // run the kernel for 5 times
			  clEnqueueNDRangeKernel(commandQueue, kernel[10], 3, NULL, global_work_size, NULL, 0, NULL, NULL);

			  status = clSetKernelArg(kernel[23], 0, sizeof(cl_mem), (void*)&para_mobj);
			  status = clSetKernelArg(kernel[23], 1, sizeof(cl_mem), (void*)&var_mobj);
			  status = clSetKernelArg(kernel[23], 2, sizeof(cl_mem), (void*)&vz_ind_mobj);
			  clEnqueueNDRangeKernel(commandQueue, kernel[23], 3, NULL, global_work_size, NULL, 0, NULL, NULL);
			  // iteration marches on
			  JACO_IT++;
		  } //end of while loop for jacobian loop
		  /*STEP 6.4 :Apply W BC [kernel 17]*/
		  //Set arguments
		  status = clSetKernelArg(kernel[17], 0, sizeof(cl_mem), (void*)&para_mobj);
		  status = clSetKernelArg(kernel[17], 1, sizeof(cl_mem), (void*)&var_mobj);
		  status = clSetKernelArg(kernel[17], 2, sizeof(cl_mem), (void*)&DIFF_ind_mobj);
		  status = clSetKernelArg(kernel[17], 3, sizeof(cl_mem), (void*)&bindex_mobj);
		  //Run kernel
		  clEnqueueNDRangeKernel(commandQueue, kernel[17], 1, NULL, global_work_size_bc, NULL, 0, NULL, NULL);
	  }

	  if (ifDeb == 1) {
		  status = clEnqueueReadBuffer(commandQueue, var_mobj, CL_TRUE, 0, nb_var*size*sizeof(REAL), var_flat, 0, NULL, NULL);
		  unflat_var(&para, var, var_flat);
		  sprintf(msg, "main(): start writing results for diffusion");
		  ffd_log(msg, FFD_NORMAL);
		  write_tecplot_data(&para, var, "result_dif");
	  }

    if (P_ON == 1) {
      /*STEP 7.0 :Projection [kernel 8]*/
      //Set arguments
      status = clSetKernelArg(kernel[8], 0, sizeof(cl_mem), (void*)&para_mobj);
      status = clSetKernelArg(kernel[8], 1, sizeof(cl_mem), (void*)&var_mobj);
      status = clSetKernelArg(kernel[8], 2, sizeof(cl_mem), (void*)&bindex_mobj);
      //Run kernel
      clEnqueueNDRangeKernel(commandQueue, kernel[8], 3, NULL, global_work_size, NULL, 0, NULL, NULL);
      /*STEP 7.1 :Pressure BC [kernel 14]*/
      //Set arguments
      status = clSetKernelArg(kernel[14], 0, sizeof(cl_mem), (void*)&para_mobj);
      status = clSetKernelArg(kernel[14], 1, sizeof(cl_mem), (void*)&var_mobj);
      status = clSetKernelArg(kernel[14], 2, sizeof(cl_mem), (void*)&bindex_mobj);
      //Run kernel
      clEnqueueNDRangeKernel(commandQueue, kernel[14], 1, NULL, global_work_size_bc, NULL, 0, NULL, NULL);
      /*STEP 7.2 :Projection W AP [kernel 12]*/
      //Set arguments
      status = clSetKernelArg(kernel[12], 0, sizeof(cl_mem), (void*)&para_mobj);
      status = clSetKernelArg(kernel[12], 1, sizeof(cl_mem), (void*)&var_mobj);
      status = clSetKernelArg(kernel[12], 2, sizeof(cl_mem), (void*)&P_ind_mobj);
      //Run kernel
      clEnqueueNDRangeKernel(commandQueue, kernel[12], 3, NULL, global_work_size, NULL, 0, NULL, NULL);
      /*STEP 7.3 :Projection Solver [kernel 11]*/
        JACO_IT = 0; // initialize index before calculation
      while (JACO_IT < IT_MAX) {
        // set the argument for the kenerl
        status = clSetKernelArg(kernel[11], 0, sizeof(cl_mem), (void*)&para_mobj);
        status = clSetKernelArg(kernel[11], 1, sizeof(cl_mem), (void*)&var_mobj);
        // run the kernel for 5 times
        clEnqueueNDRangeKernel(commandQueue, kernel[11], 3, NULL, global_work_size, NULL, 0, NULL, NULL);

        status = clSetKernelArg(kernel[24], 0, sizeof(cl_mem), (void*)&para_mobj);
        status = clSetKernelArg(kernel[24], 1, sizeof(cl_mem), (void*)&var_mobj);
        clEnqueueNDRangeKernel(commandQueue, kernel[24], 3, NULL, global_work_size, NULL, 0, NULL, NULL);
        // iteration marches on
        JACO_IT++;
      } //end of while loop for jacobian loop
      /*STEP 7.4 :Pressure BC [kernel 14]*/
      //Set arguments
      status = clSetKernelArg(kernel[14], 0, sizeof(cl_mem), (void*)&para_mobj);
      status = clSetKernelArg(kernel[14], 1, sizeof(cl_mem), (void*)&var_mobj);
      status = clSetKernelArg(kernel[14], 2, sizeof(cl_mem), (void*)&bindex_mobj);
      //Run kernel
      clEnqueueNDRangeKernel(commandQueue, kernel[14], 1, NULL, global_work_size_bc, NULL, 0, NULL, NULL);
      /*STEP 7.5 :Velocity after Projection [kernel 9]*/
      //Set arguments
      status = clSetKernelArg(kernel[9], 0, sizeof(cl_mem), (void*)&para_mobj);
      status = clSetKernelArg(kernel[9], 1, sizeof(cl_mem), (void*)&var_mobj);
      status = clSetKernelArg(kernel[9], 2, sizeof(cl_mem), (void*)&bindex_mobj);
      //Run kernel
      clEnqueueNDRangeKernel(commandQueue, kernel[9], 3, NULL, global_work_size, NULL, 0, NULL, NULL);
    }

	if (ifDeb == 1) {
		status = clEnqueueReadBuffer(commandQueue, var_mobj, CL_TRUE, 0, nb_var*size*sizeof(REAL), var_flat, 0, NULL, NULL);
		unflat_var(&para, var, var_flat);
		sprintf(msg, "main(): start writing results for projection");
		ffd_log(msg, FFD_NORMAL);
		write_tecplot_data(&para, var, "result_pro");
	}

    /*STEP ## :Mass Balance Disabled*/
    if (para.bc->nb_outlet != 0) {
      if (M_ON == 1) {
        /* Mass balance running on one core Kernel[19]*/
        status = clSetKernelArg(kernel[19], 0, sizeof(cl_mem), (void*)&para_mobj);
        status = clSetKernelArg(kernel[19], 1, sizeof(cl_mem), (void*)&var_mobj);
        status = clSetKernelArg(kernel[19], 2, sizeof(cl_mem), (void*)&bindex_mobj);
        //Run kernel
        clEnqueueNDRangeKernel(commandQueue, kernel[19], 1, NULL, global_work_size_mb, NULL, 0, NULL, NULL);
        /* Adjust velocity Kernel[18]*/
        status = clSetKernelArg(kernel[18], 0, sizeof(cl_mem), (void*)&para_mobj);
        status = clSetKernelArg(kernel[18], 1, sizeof(cl_mem), (void*)&var_mobj);
        status = clSetKernelArg(kernel[18], 2, sizeof(cl_mem), (void*)&bindex_mobj);
        //Run kernel
        clEnqueueNDRangeKernel(commandQueue, kernel[18], 1, NULL, global_work_size_bc, NULL, 0, NULL, NULL);
      } //end of if M_ON
    }// end of mass balance

    if (T_ON == 1) {
      /*STEP 9.0 :Advection T [kernel 3]*/
      //Set arguments
      status = clSetKernelArg(kernel[3], 0, sizeof(cl_mem), (void*)&para_mobj);
      status = clSetKernelArg(kernel[3], 1, sizeof(cl_mem), (void*)&var_mobj);
      status = clSetKernelArg(kernel[3], 2, sizeof(cl_mem), (void*)&bindex_mobj);
      //Run kernel
      clEnqueueNDRangeKernel(commandQueue, kernel[3], 3, NULL, global_work_size, NULL, 0, NULL, NULL);
      /*STEP 9.1 :Apply TMP3 BC [kernel 13]*/
      //Set arguments
      status = clSetKernelArg(kernel[13], 0, sizeof(cl_mem), (void*)&para_mobj);
      status = clSetKernelArg(kernel[13], 1, sizeof(cl_mem), (void*)&var_mobj);
      status = clSetKernelArg(kernel[13], 2, sizeof(cl_mem), (void*)&ADVE_ind_mobj);
      status = clSetKernelArg(kernel[13], 3, sizeof(cl_mem), (void*)&bindex_mobj);
      //Run kernel
      clEnqueueNDRangeKernel(commandQueue, kernel[13], 1, NULL, global_work_size_bc, NULL, 0, NULL, NULL);


      /*STEP 10.0 :Diffusion T [kernel 7]*/
      //Set arguments
      status = clSetKernelArg(kernel[7], 0, sizeof(cl_mem), (void*)&para_mobj);
      status = clSetKernelArg(kernel[7], 1, sizeof(cl_mem), (void*)&var_mobj);
      status = clSetKernelArg(kernel[7], 2, sizeof(cl_mem), (void*)&bindex_mobj);
      //Run kernel
      clEnqueueNDRangeKernel(commandQueue, kernel[7], 3, NULL, global_work_size, NULL, 0, NULL, NULL);
      /*STEP 10.1 :Apply T BC [kernel 13]*/
      //Set arguments
      status = clSetKernelArg(kernel[13], 0, sizeof(cl_mem), (void*)&para_mobj);
      status = clSetKernelArg(kernel[13], 1, sizeof(cl_mem), (void*)&var_mobj);
      status = clSetKernelArg(kernel[13], 2, sizeof(cl_mem), (void*)&DIFF_ind_mobj);
      status = clSetKernelArg(kernel[13], 3, sizeof(cl_mem), (void*)&bindex_mobj);
      //Run kernel
      clEnqueueNDRangeKernel(commandQueue, kernel[13], 1, NULL, global_work_size_bc, NULL, 0, NULL, NULL);
      /*STEP 10.2 :Diffusion T AP [kernel 12]*/
      //Set arguments
      status = clSetKernelArg(kernel[12], 0, sizeof(cl_mem), (void*)&para_mobj);
      status = clSetKernelArg(kernel[12], 1, sizeof(cl_mem), (void*)&var_mobj);
      status = clSetKernelArg(kernel[12], 2, sizeof(cl_mem), (void*)&T_ind_mobj);
      //Run kernel
      clEnqueueNDRangeKernel(commandQueue, kernel[12], 3, NULL, global_work_size, NULL, 0, NULL, NULL);
      /*STEP 10.3 :Diffusion T Solver [kernel 10]*/
      JACO_IT = 0; // initialize index before calculation
      while (JACO_IT < IT_MAX) {
        // set the argument for the kenerl
        status = clSetKernelArg(kernel[10], 0, sizeof(cl_mem), (void*)&para_mobj);
        status = clSetKernelArg(kernel[10], 1, sizeof(cl_mem), (void*)&var_mobj);
        status = clSetKernelArg(kernel[10], 2, sizeof(cl_mem), (void*)&T_ind_mobj);

        // run the kernel for 5 times
        clEnqueueNDRangeKernel(commandQueue, kernel[10], 3, NULL, global_work_size, NULL, 0, NULL, NULL);

        status = clSetKernelArg(kernel[23], 0, sizeof(cl_mem), (void*)&para_mobj);
        status = clSetKernelArg(kernel[23], 1, sizeof(cl_mem), (void*)&var_mobj);
        status = clSetKernelArg(kernel[23], 2, sizeof(cl_mem), (void*)&T_ind_mobj);
        clEnqueueNDRangeKernel(commandQueue, kernel[23], 3, NULL, global_work_size, NULL, 0, NULL, NULL);
        // iteration marches on
        JACO_IT++;
      } //end of while loop for jacobian loop
      /*STEP 10.4 :Apply T BC [kernel 13]*/
      //Set arguments
      status = clSetKernelArg(kernel[13], 0, sizeof(cl_mem), (void*)&para_mobj);
      status = clSetKernelArg(kernel[13], 1, sizeof(cl_mem), (void*)&var_mobj);
      status = clSetKernelArg(kernel[13], 2, sizeof(cl_mem), (void*)&DIFF_ind_mobj);
      status = clSetKernelArg(kernel[13], 3, sizeof(cl_mem), (void*)&bindex_mobj);
      //Run kernel
      clEnqueueNDRangeKernel(commandQueue, kernel[13], 1, NULL, global_work_size_bc, NULL, 0, NULL, NULL);
    }
	if (ifDeb == 1) {
		status = clEnqueueReadBuffer(commandQueue, var_mobj, CL_TRUE, 0, nb_var*size*sizeof(REAL), var_flat, 0, NULL, NULL);
		unflat_var(&para, var, var_flat);
		sprintf(msg, "main(): start writing results for scalar");
		ffd_log(msg, FFD_NORMAL);
		write_tecplot_data(&para, var, "result_scalar");
	}

      /*STEP 11.0: Call back date to host*/
      status = clEnqueueReadBuffer(commandQueue, ADVE_ind_mobj, CL_TRUE, 0, sizeof(int), ADVE_IND, 0, NULL, NULL);
    /*STEP 11.1 :Time Marches On*/
      timing(&para);
    /*STEP 11.2 :add averaged data ?*/
      //if (para.outp->cal_mean == 1 && para.mytime->t>para.mytime->t_steady) {
      if (para.mytime->t>para.mytime->t_steady) {
        if (print_control <= 1) {
          sprintf(msg, "main(): start averaging data at", para.mytime->t);
          ffd_log(msg, FFD_NORMAL);
		  para.mytime->step_mean = 0;
          print_control += 1;
        }
        status = clSetKernelArg(kernel[21], 0, sizeof(cl_mem), (void*)&para_mobj);
        status = clSetKernelArg(kernel[21], 1, sizeof(cl_mem), (void*)&var_mobj);
        clEnqueueNDRangeKernel(commandQueue, kernel[21], 3, NULL, global_work_size, NULL, 0, NULL, NULL);
        para.mytime->step_mean++;
      }
    /*STEP 11.2 :Simulation Continues ?*/
        NEXT = para.mytime->step_current < para.mytime->step_total ? 1 : 0;
  }// end of while loop

  /*STEP 11: Call back date to host*/
  status = clEnqueueReadBuffer(commandQueue, var_mobj, CL_TRUE, 0, nb_var*size*sizeof(REAL), var_flat, 0, NULL, NULL);
  //status = clEnqueueReadBuffer(commandQueue, para_mobj, CL_TRUE, 0, sizeof(para_simp), &para_simp, 0, NULL, NULL);
  

  /*un-flatten back the var and bindex variable*/
  unflat_var(&para, var, var_flat);
  /*STEP 11: run time average*/
  average_time(&para, var);

#ifdef FFD_ISAT //Called by ISAT
  int imax = para.geom->imax, jmax = para.geom->jmax;
  int kmax = para.geom->kmax;
  int IMAX = imax + 2, IJMAX = (imax + 2)*(jmax + 2);
  ffdOutput[0] = average_volume(&para, var, var[TEMPM]);
  //printf("ffdOuput=\t%f\n", ffdOutput[0]);
  //getchar();
  //ffdOutput[0] = var[TEMPM][IX(para.geom->imax / 2, para.geom->jmax / 2, para.geom->kmax / 2+2)];
  //printf("the ffdOutput[0] is %f\n", ffdOutput[0]);
  //getchar();
#endif

  /*For saving time do not output data*//*
  //writing the tecplot results for output
  sprintf(msg, "main(): start writing results");
  ffd_log(msg, FFD_NORMAL);
  write_tecplot_data(&para, var, "result");
  //write_tecplot_all_data(&para, var, "result_all");
  */

  /*Clean Up of memory Allocations*/
  sprintf(msg, "main(): start freeing memories in GPU");
  ffd_log(msg, FFD_NORMAL);
  status = clFlush(commandQueue);
  status = clReleaseKernel(kernel[0]);
  status = clReleaseKernel(kernel[1]);
  status = clReleaseKernel(kernel[2]);
  status = clReleaseKernel(kernel[3]);
  status = clReleaseKernel(kernel[4]);
  status = clReleaseKernel(kernel[5]);
  status = clReleaseKernel(kernel[6]);
  status = clReleaseKernel(kernel[7]);
  status = clReleaseKernel(kernel[8]);
  status = clReleaseKernel(kernel[9]);
  status = clReleaseKernel(kernel[10]);
  status = clReleaseKernel(kernel[11]);
  status = clReleaseKernel(kernel[12]);
  status = clReleaseKernel(kernel[13]);
  status = clReleaseKernel(kernel[14]);
  status = clReleaseKernel(kernel[15]);
  status = clReleaseKernel(kernel[16]);
  status = clReleaseKernel(kernel[17]);
  status = clReleaseKernel(kernel[18]);
  status = clReleaseKernel(kernel[19]);
  status = clReleaseKernel(kernel[20]);
  status = clReleaseKernel(kernel[21]);
  status = clReleaseKernel(kernel[22]);
  status = clReleaseKernel(kernel[23]);
  status = clReleaseKernel(kernel[24]);
  status = clReleaseProgram(program);
  status = clReleaseMemObject(var_mobj);
  status = clReleaseMemObject(para_mobj);
  status = clReleaseMemObject(bindex_mobj);
  status = clReleaseMemObject(vx_ind_mobj);
  status = clReleaseMemObject(vy_ind_mobj);
  status = clReleaseMemObject(vz_ind_mobj);
  status = clReleaseMemObject(T_ind_mobj);
  status = clReleaseMemObject(P_ind_mobj);
  status = clReleaseMemObject(ADVE_ind_mobj);
  status = clReleaseMemObject(DIFF_ind_mobj);
  status = clReleaseContext(context);
  status = clFinish(commandQueue);
  status = clReleaseCommandQueue(commandQueue);


  /*free the memory in CPU*/
  sprintf(msg, "main(): start freeing memories in CPU");
  ffd_log(msg, FFD_NORMAL);
  free_data(var);
  free_index(BINDEX);
  free(var_flat);
  free(bindex_flat);

  return SUCCESS;
}