#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#include "allvars.h"
#include "proto.h"
#include <sys/time.h>
#include "map.h"

void print_usage_message(int);
void get_input(int nArg, void* argv[], struct st_region* region, 
	       struct st_im* image, struct st_cam* cam, struct st_flag* flag,
	       struct st_out* out); 
void determine_hsml();
void calculate_hsml(struct st_flag* flag, MyFloat pixSize);

int Verbose = 0;

/**
 * @brief Main public function to compute particle smoothing lengths and
 *        project them onto a surface density image.
 */
int findHsmlAndProjectYB_Tau(int argc, void *argv[]) {
  
  /* Structures for camera parameters, selection region, 
   * image parameters, and configuration flags */
  struct st_cam cam;
  struct st_region region;
  struct st_im image;
  struct st_flag flag;
  struct st_out out;

  /* Parse input parameters into internally-useful form */
  get_input(argc, argv, &region, &image, &cam, &flag, &out);

  /* Calculate particles' smoothing lengths */
  calculate_hsml(&flag, region.xWidth*2/image.xPix);

  if(flag.hsml_only)
    return 0;

  /* Create actual image */
  make_image(&cam, &out, &flag);

  return 0;
}

/**
 * @brief Check that input is sensible and parse arguments if so.
 * 
 * If the input is inconsistent, a usage message is printed instead.
 * 
 * @param nArg Number of detected arguments.
 * @param argv[] The arguments passed to the main function.
 * @param region: Pointer to #st_region structure, for particle limits.
 * @param image: Pointer to #st_im structure, for image parameters.
 * @param cam: Pointer to #st_cam structure, for camera parameters.
 * @param flag: Pointer to #st_flag structure, for flags.
 *
 * Internally accessed global variables:
 *   -- NumPart, P, Hsml, Mass, Quantity, DesDensNgb, Hmax,
 *      Value, ValueQuantity
 */
void get_input(int nArg, void* argv[], struct st_region* region, 
	       struct st_im* image, struct st_cam* cam, struct st_flag* flag,
	       struct st_out* out) 
{
  
  int kk;

  /* Check that input is as expected and abort with instructions if not */
  if(nArg != 35) {
    print_usage_message(nArg);
    exit(0);
  }
  
  /* Number of particles in input data */
  NumPart = *(long *) argv[0]; 

  /* Relevant properties of particles */
  P = (struct particle_data *) argv[1];
  Hsml = (float *) argv[2];
  Mass = (float *) argv[3];
  Quantity = (float *) argv[4];
  PVel = (struct particle_vel_data *) argv[32];

  /* Limits of particles to be considered. Note that x/y limits imported 
   * here are only used in orthogonal projection. */
  region->xWidth = *(float *) argv[5];
  region->yWidth = *(float *) argv[6];
  region->zMin = *(float *) argv[7];
  region->zMax = *(float *) argv[8];
  region->vzMin = *(float *) argv[30];
  region->vzMax = *(float*) argv[31];

  /* Image parameters */
  image->xPix = *(int *) argv[9];
  image->yPix = *(int *) argv[10];
  image->zPix = *(int *) argv[11];
  image->tau = *(float *) argv[12];
  image->nPix = (long)image->xPix * image->yPix * image->zPix;
  image->zPix_out = 2;
  
  /* Smoothing parameters */
  DesDensNgb =  *(int *) argv[13];
  Hmax = *(float *) argv[14];
  
  /* Positioning of the camera */
  MyFloat* mtemp = (MyFloat*) argv[15];
  for (kk = 0; kk < 3; kk++)
    cam->pos[kk] = mtemp[kk];
  float* temp = (float*) argv[16];
  for (kk = 0; kk < 3; kk++)
    cam->dir[kk] = temp[kk];
  temp = (float*) argv[17];
  for (kk = 0; kk < 3; kk++)
    cam->angle[kk] = temp[kk];
  temp = (float*) argv[18];
  for (kk = 0; kk < 2; kk++)
    cam->fov[kk] = temp[kk];
  cam->roll_mode = *(int*) argv[19];

  /* For convenience, store pointers to image and region structs */
  cam->image = image;
  cam->region = region;

  /* Output arrays */
  Value =         (float *) argv[20];
  ValueQuantity = (float *) argv[21];
  Cube = (float *) argv[33];
  CubeQuantity = (float *) argv[34];

  /* Flags to modify the behaviour of the code */
  flag->use_hsml = *(int *) argv[22];
  flag->hsml_only = *(int *) argv[23];
  flag->tree_alloc_fac = *(float *) argv[24];
  Verbose = *(int*) argv[25];
  flag->losvd = *(int *) argv[29];

  /* Output arrays for camera basis vectors */
  out->x = (float*) argv[26];
  out->y = (float*) argv[27];
  out->z = (float*) argv[28];
  
  
  if (!Verbose)
    return;
  
  printf("\n");
  printf("========================================================\n");
  printf("           Image generator -- settings                  \n");
  printf("--------------------------------------------------------\n");
  printf("NumPart = %ld\n", NumPart);
  printf("Selection region: xWidth=%g, yWidth=%g,\n"
	 "                  zMin=%g, zMax=%g, \n"
         "                  vzMin=%g, vzMax=%g\n",
	 region->xWidth, region->yWidth, region->zMin, region->zMax,
	 region->vzMin, region->vzMax);
  printf("Image parameters: xPix=%d, yPix=%d, zPix=%d, tau=%.3e\n",
	 image->xPix, image->yPix, image->zPix, image->tau);
  printf("Basis vectors:\n");
  printf("  x --> [%g, %g, %g]\n", out->x[0], out->x[1], out->x[2]);
  printf("  y --> [%g, %g, %g]\n", out->y[0], out->y[1], out->y[2]);
  printf("  z --> [%g, %g, %g]\n", out->z[0], out->z[1], out->z[2]);
  
  printf("CamDir = [%g, %g, %g]\n", cam->dir[0], cam->dir[1], cam->dir[2]);
  printf("CamPos = [%g, %g, %g]\n", cam->pos[0], cam->pos[1], cam->pos[2]);
  printf("CamAngle = [%g, %g, %g]\n", 
	 cam->angle[0], cam->angle[1], cam->angle[2]);
  printf("CamFOV = [%g, %g]\n", cam->fov[0], cam->fov[1]);
  printf("Camera roll mode: %d\n", cam->roll_mode);
  printf("Flags: use_hsml=%d, hsml_only=%d, tree_alloc_fac=%f, losvd=%d\n",
	 flag->use_hsml, flag->hsml_only, flag->tree_alloc_fac,
	 flag->losvd);
  printf("=======================================================\n");  
  printf("\n");
  
  return;
}




/**
 * @brief Print message outlining required arguments in case their number
 *        did not match the expected number.
 */
void print_usage_message(int nArg) {
  fprintf(stderr, "\n\nWrong number of arguments! (found %d) \n\n", nArg);
  fprintf(stderr, "Need to supply the following arguments:\n\n");
  
  return;
}


void determine_hsml(void)
{
  int i;
  double h = Hmax/100;  /* Just a semi-reasonable, semi-conservative starting guess */

  /*  omp_set_num_threads(1); */
  /*#pragma omp parallel for firstprivate(h, R2list) */
  for(i = 0; i < NumPart; i++)
    {
      /*int signal = 0; */

      if((i % (NumPart / 100)) == 0)
	{
	  printf("x");
	  fflush(stdout);
	}

      /*
      if(i > (signal / 100.0) * NumPart)
        {
          printf("x");
          fflush(stdout);
          signal++;
        }
      */

      /*printf("Particle %d...\n", i); */
      if(Hsml[i] == 0)
	{
	  Hsml[i] = h = ngb_treefind(P[i].Pos, DesDensNgb, h * 1.1);
	  /*if(Hsml[i] == 0)
	    printf("Warning! Particle %d has hsml = 0!\n", i); */
	}

    }

  printf("\n");
}

/**
 * @brief Calculate smoothing lenghts of all particles.
 * 
 * These can either be calculated internally from scratch (mode==1)
 * or taken from the input Hsml array, clipped to Hmin and Hmax (mode==0)
 *
 * @param flag Pointer to flag structure (#st_flag)
 * @param pixSize The pixel size to be used for softening calculation
 *                (used internally in tree building).
 *
 * Internally accessed global variables:
 *   -- Hmax, Softening
 */
void calculate_hsml(struct st_flag* flag, MyFloat pixSize) {
  
  /* If we use input values, we only have to max-clip them */
  if (flag->use_hsml) {
    long ii;
    for (ii = 0; ii < NumPart; ii++)
      if (Hsml[ii] > Hmax)
	Hsml[ii] = Hmax;
  }
  
  /* More fun option: we need to calculate smoothing lengths from scratch */
  else {

    double time_ph_start = get_wall_time();
    double time_tree_start, time_hsml_start, time_restore_start;

    /* Bring particles in Peano-Hilbert order for tree-building */
    printf("Bring particles into Peano-Hilbert order...");
    peano_hilbert_order();
    printf(" done (%g sec.).\n", get_wall_time()-time_ph_start);
    
    /* Allocate tree nodes, use supplied allocation factor */
    tree_treeallocate(flag->tree_alloc_fac * NumPart, NumPart);
    
    /* Build tree */
    Softening = pixSize/100;   /* Used internally by tree routines */
    time_tree_start = get_wall_time();
    printf("Build tree...\n");
    tree_treebuild();
    printf("done (%g sec.)\n", get_wall_time()-time_tree_start);

    /* Actually find each particle's smoothing length */
    time_hsml_start = get_wall_time();
    printf("Calculating smoothing lengths...\n");
    determine_hsml();
    printf("done (%g sec.)\n", get_wall_time()-time_hsml_start);

    /* Free memory used for tree */
    tree_treefree();

    /* Important: restore particles back to original order. */
    time_restore_start = get_wall_time();
    printf("Restoring original particle order...");
    restore_original_order();
    printf(" done (%g sec.)\n", get_wall_time()-time_restore_start);
    
  } /* Done with need-to-calculate-hsml section */

  return;
}


double get_wall_time() {

  struct timeval time;
  if(gettimeofday(&time, NULL))
    return 0;
  
  return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
