/**
 * This file contains the 'map-making' routines.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <sys/time.h>

#include "allvars.h"
#include "map.h"
#include "map_private.h"

extern int Verbose;

/* ---------------------------------------------------------------- *
 * -------------- Externally-used function: make_image() ---------- *
 * ---------------------------------------------------------------- */

/**
 * @brief Create surface density images from individual particles.
 *
 * Actually, four maps are created. A mass and mass-weighted-quantity
 * map uses standard 2D projection, integrating along the (camera-)z*
 * axis. A second pair of maps uses simplified ray-tracing to give more
 * weight to structures close to the camera.
 * 
 * @param cam The camera parameters (as #st_cam struct pointer)
 * @param out Pointer to #st_out struct holding output camera bases.
 * @param flag Pointer to #st_flag struct holding general flags.
 * @param mass3D_out 3D mass map provided from outside.
 * @param quant3D_out 3D quant map provided from outside.
 */
void make_image(struct st_cam* cam, struct st_out* out,
		struct st_flag* flag) {

  /* Temporary 3D maps */
  float* mass3D = NULL;
  float* quant3D = NULL;

  /* Need to point mass3D and quant3D to provided output maps, if 
   * we want to store the full 3D distribution */
  if (flag->losvd > 0) {
    mass3D = Cube;
    quant3D = CubeQuantity;
  }

  /* Set up camera frame */
  get_basis_vectors(cam, out);

  /* Set up 'general-purpose' projection parameters: */
  setup_projection(cam);

  /* Sum up particle contributions to temporary (3D) maps */
  /* (NB: quant is left as quant*mass here!) */
  calculate_3D_map(&mass3D, &quant3D, cam, flag);

  /* Project 3D map to 2D: simple sum and ray-tracing */
  /* (only at the end of this is quant normalized by mass) */
  project_3D_map(mass3D, quant3D, cam->image, flag, cam->flag_perspectivic,
		   cam->region->zMax - cam->region->zMin);


  printf("\nDone generating SPH image.\n");
  return;

} /* ends make_image() */


/* ---------- get_basis_vectors()  ----------------------------------- */

/**
 * @brief Compute the basis vectors of the camera frame (x*, y*, z*)
 *
 * @param cam Pointer to camera property structure (#st_cam)
 * @param out Pointer to #st_out struct holding output camera bases.
 *
 * Convention: 
 *  -- The z* vector points along the inverse camera viewing direction
 *     (i.e. outwards towards the unit sphere, rather than inwards).
 *  -- x* and y* are both in the tangential plane to the sphere, at the 
 *     point where the camera sits (looking inwards). The direction of 
 *     x* depends on the camera roll angle (rho) and the basis definition
 *     (see below), y* is 90 degrees counter-clockwise from x*.
 *  -- The cross-product of x* (x) y* is z* (as usual). 
 *
 * There are three implemented definitions of setting the 'zero-roll' x*'
 * vector, controlled by the BasisType parameter:
 *  -- 0: standard geographic definition (x*' --> East, y*' --> North);
 *        degenerate at poles (defined as x*'=x/y*'=y at N, x*'=-x/y*'=y at S)
 *  -- 1: Parallel transport of the standard (x/y/z) bases from the North
 *        Pole to the camera point. This definition is degenerate at the 
 *        South Pole, where we define x*'=-x, y*'=y.
 *  -- 2: `y-flow' bases, where y'* is the (normalized) projection of y
 *        onto the sphere at the camera point. This corresponds to 
 *        x/y/z at the North Pole, and -x/y/z at the South Pole (like the
 *        standard geographic definition, BasisType==0). It is degenerate at
 *        the `y-entry/exit points' (0, +-1, 0), with definition of 
 *        x*'=x, y*'=z (0, -1, 0) and x*'=x, y*'=-z, which keeps the axes
 *        aligned during a rotation about the y-axis.
 */
void get_basis_vectors(struct st_cam* cam, struct st_out* out) {

  int kk;  /* Counter for loops over dimensions */

  /* First thing: check if basis vectors were supplied. In this case,
     we're done already... */
  const double len_x = fabs(out->x[0]) + fabs(out->x[1]) + fabs(out->x[2]);
  const double len_y = fabs(out->y[0]) + fabs(out->y[1]) + fabs(out->y[2]);
  const double len_z = fabs(out->z[0]) + fabs(out->z[1]) + fabs(out->z[2]);
  if (len_x > 1e-3 && len_y > 1e-3 && len_z > 1e-3) {
    printf("Input basis vectors are valid, using those...\n");
    normalize_vector_float(out->x);
    normalize_vector_float(out->y);
    normalize_vector_float(out->z);
    
    for (kk = 0; kk < 3; kk++) {
      cam->base.x[kk] = out->x[kk];
      cam->base.y[kk] = out->y[kk];
      cam->base.z[kk] = out->z[kk];
      cam->dir[kk] = -out->z[kk]; 
    }

  } else {

    /* Struct to hold sin/cos values of angles. */
    struct st_trig trig;
    
    /* Structs to hold geographic and camera basis vectors */
    struct st_base base_geo;
    
    /* Check if the cam->dir variable is set to meaningful value */
    /* If it is ~(0,0,0), initialise from cam->angle */
    /* Otherwise, normalise and compute angles. */
    if(fabs(cam->dir[0]) < 1e-3 && fabs(cam->dir[1]) < 1e-3 
       && fabs(cam->dir[2]) < 1e-3)
      get_cam_dir_from_angles(cam, &trig);
    else {
      normalize_vector_float(cam->dir);
      
      /* Also need to set up the corresponding angles, to set up 
       * x' and y' basis vectors (in camera frame) */
      get_angles_from_cam_dir(&trig, cam);
    }
    
    /* z* vector is always the *negative* camera direction */
    for (kk = 0; kk < 3; kk++)
      cam->base.z[kk] = -(cam->dir[kk]);
    
    /* Definition of x* and y* vectors depends on settings */
    if (cam->roll_mode == 2)   /* 'y-flow' basis vectors (`old-style') */ 
      form_camera_bases_yparallel(cam);
    else {
      /* Form `geographic' basis vectors at camera position (-CamDir) *
       * (see function definition for more details)                   */
      form_geographic_bases(&base_geo, &trig, cam);
      
      if (cam->roll_mode == 0)  /* Use standard geographic basis vectors */
	form_geographic_camera_bases(&base_geo, cam);
      else
	form_camera_bases(&base_geo, cam);
    }
  } /* ends section if bases had to be calculated */
  
  printf("Camera basis vectors determined as\n"
	   "[%g, %g, %g],\n[%g, %g, %g],\n[%g, %g, %g]\n" , 
	   cam->base.x[0], cam->base.x[1], cam->base.x[2], 
	   cam->base.y[0], cam->base.y[1], cam->base.y[2],
	   cam->base.z[0], cam->base.z[1], cam->base.z[2]);
    
  const double xy = cam->base.x[0]*cam->base.y[0] + 
    cam->base.x[1]*cam->base.y[1] + cam->base.x[2]*cam->base.y[2];
  const double xz = cam->base.x[0]*cam->base.z[0] + 
    cam->base.x[1]*cam->base.z[1] + cam->base.x[2]*cam->base.z[2];
  const double yz = cam->base.y[0]*cam->base.z[0] + 
    cam->base.y[1]*cam->base.z[1] + cam->base.y[2]*cam->base.z[2];
  
  printf("   Test: dot products = %g (xy), %g (xz), %g (yz)\n", xy,xz,yz);

  /* Store camera bases in output structure, to pass back to library caller */
  for (kk = 0; kk < 3; kk++) {
    out->x[kk] = cam->base.x[kk];
    out->y[kk] = cam->base.y[kk];
    out->z[kk] = cam->base.z[kk];
  }

  return;
} /* ends get_basis_vectors() */ 

 
/**
 * @brief Compute camera direction vector from supplied angles.
 *
 * Conventions:
 *  -- Angles specify the point on the (unit) sphere from which the camera 
 *     looks inwards, while CamDir points along the camera viewing axis.
 *  -- [0, 0, 1] is the North Pole, [0, 0, -1] the South Pole
 *  -- Theta is measured from the positive z axis: the North Pole has 
 *     theta=0, the equator theta=pi/2, the South Pole theta=pi.
 *  -- Phi is measured counter-clockwise from the x-axis.
 *
 * This implies that e.g. the default view from the North Pole corresponds
 * to angles (0, 0, 0) and direction (0, 0, -1).
 *
 * @param cam Pointer to camera parameter structure (#st_cam).
 * @param trig Pointer to #st_trig structure, to hold cos/sin of angles.
 */
void get_cam_dir_from_angles(struct st_cam* cam, struct st_trig* trig) {
  
  printf("... initialize line of sight from CamAngle...\n");
  
  trig->cosTheta = cos(cam->angle[0]);
  trig->cosPhi = cos(cam->angle[1]);
  trig->sinTheta = sin(cam->angle[0]);
  trig->sinPhi = sin(cam->angle[1]);
  
  /* Minus sign is there to get camera looking 'inwards'. */
  cam->dir[0] = -1*trig->sinTheta*trig->cosPhi;
  cam->dir[1] = -1*trig->sinTheta*trig->sinPhi;
  cam->dir[2] = -1*trig->cosTheta;

  printf("    CamDir recalculated as [%g, %g, %g]\n",  
	 cam->dir[0], cam->dir[1], cam->dir[2]);
  return;
}

/**
 * @brief Compute relevant angles from CamDir.
 *
 * See description of get_cam_dir_from_angles for angle/direction conventions.
 *
 * @param trig Pointer to #st_trig structure, to hold cos/sin of angles.
 * @param cam Pointer to camera parameter structure (#st_cam)
 */
void get_angles_from_cam_dir(struct st_trig* trig, struct st_cam* cam) {

  /* Minus signs because camera looks inwards from specified point on *
   * the unit sphere */

  trig->cosTheta = -(cam->dir[2]);
  trig->sinTheta = sqrt(1-(cam->dir[2])*(cam->dir[2]));
  trig->cosPhi = -cam->dir[0] / trig->sinTheta; 
  trig->sinPhi = -cam->dir[1] / trig->sinTheta; 
  cam->angle[1] = atan2(-(cam->dir[1]), -(cam->dir[0]));
  
  return;
}

/**
 * @brief Calculate the 'geographic basis vectors' at the camera position.
 * 
 * The camera position here is the point on the (infinitesimal) unit sphere 
 * around the `actual' camera coordinates on which the camera sits, looking 
 * inwards. The geographic basis vectors lie in the plane tangential to the 
 * sphere at this point, and point `south' (x': towards [0, 0, -1]) and 
 * `east' (y': towards increasing phi at constant theta).
 *
 * Note that this is *not* the standard 'East-is-right, North-is-up' 
 * configuration, because this makes it easier to specify rotations about
 * the y-axis (i.e. left/right from the North Pole).
 * 
 * The above definition breaks down at the North and South poles. There, 
 * we use the ad-hoc definition of x'=x, y'=y (North) and x'=-x, y'=y (South).
 * This keeps the bases continuous for rotations about the y axis.
 *
 * @param geo Pointer to structure holding the bases to be computed.
 * @param trig Pointer to structure holding the cos/sin of camera angles.
 * @param cam Pointer to camera parameter structure (#st_cam).
 */
void form_geographic_bases(struct st_base* geo, struct st_trig* trig,
			   struct st_cam* cam) {

  if (cam->dir[2] < -0.9999) {  /*  North pole */
    geo->x[0] = 1.0;
    geo->x[1] = 0.0;
    geo->x[2] = 0.0;
    geo->y[0] = 0.0;
    geo->y[1] = 1.0;
    geo->y[2] = 0.0;
  } else if (cam->dir[2] > 0.9999) {  /* South pole */
    geo->x[0] = -1.0;
    geo->x[1] = 0.0;
    geo->x[2] = 0.0;
    geo->y[0] = 0.0;
    geo->y[1] = 1.0;
    geo->y[2] = 0.0;
  } else {  /* Rest of the globe */
    geo->x[0] = trig->cosTheta * trig->cosPhi;
    geo->x[1] = trig->cosTheta * trig->sinPhi;
    geo->x[2] = -1*trig->sinTheta;
    geo->y[0] = -1*trig->sinPhi;
    geo->y[1] = trig->cosPhi;
    geo->y[2] = 0.0;
  }

  if (!Verbose)
    return;

  printf("(Modified) geographic vectors determined as\n"
	 "   [%g, %g, %g], [%g, %g, %g]\n", 
	 geo->x[0], geo->x[1], geo->x[2], geo->y[0], geo->y[1], geo->y[2]);
  
  const double dot = geo->x[0]*geo->y[0] + 
    geo->x[1]*geo->y[1] +
    geo->x[2]*geo->y[2]; 

  const double dot_xz = geo->x[0]*cam->base.z[0] + 
    geo->x[1]*cam->base.z[1] + geo->x[2]*cam->base.z[2];
  const double dot_yz = geo->y[0]*cam->base.z[0] + 
    geo->y[1]*cam->base.z[1] + geo->y[2]*cam->base.z[2];

  printf("   Test: dot-products = %g (xy), %g (xz), %g (yz)\n", 
	 dot, dot_xz, dot_yz);
  	 
  return;
}

/**
 * @brief Calculate the actual camera basis vectors.
 * 
 * This corresponds to parallel transport of the standard x/y/z vectors
 * from the North Pole to the camera point.
 *
 * By construction, this gives continuous frame orientations for all 
 * rotations crossing the North Pole. However, it does *not* yield
 * continuous frame orientations for other rotations, nor give the 
 * standard Cartesian orientations at all 6 poles (+-x, +-y, +-z):
 * it *does* work at +z (xy), -y (xz), and -x (zy). 
 *
 * @param geo Pointer to the geographic basis vectors (#st_base struct).
 * @param cam Pointer to the to be calculated camera basis vectors.
 */

void form_camera_bases(struct st_base* geo, struct st_cam* cam) {
  
  int kk;

  /* Tau is the (counter-clockwise) rotation relative to the input
   * geographic basis vectors */
  double tau;

  /* Since the geo-bases are set up specially at the poles, we 
   * need to treat them specially here too --> set up to not need 
   * additional modification by phi. */
  if (cam->dir[2] < -0.9999 || cam->dir[2] > 0.9999)
    tau = cam->angle[2];
  else
    tau = cam->angle[2] - cam->angle[1];

  /* Now construct final bases in frame spanned by geo-bases */
  for (kk = 0; kk < 3; kk++) {
    cam->base.x[kk] = cos(tau) * geo->x[kk] + sin(tau) * geo->y[kk];
    cam->base.y[kk] = -sin(tau) * geo->x[kk] + cos(tau) * geo->y[kk];
  }

  return;
}


/**
 * @brief Calculate the camera bases in standard geographical orientation
 *        (North is up, East is right).
 *
 * @param geo Pointer to the geographic basis vectors (#st_base struct).
 * @param cam Pointer to the to be calculated camera basis vectors.
 */

void form_geographic_camera_bases(struct st_base* geo, struct st_cam* cam) {
  
  int kk;
  double tau;  /* Rotation relative to modified geographic bases */

  /* Since the geo-bases are set up specially at the poles, we 
   * need to treat them specially here too --> set up to not need 
   * additional modification by phi. */
  if (cam->dir[2] < -0.9999 || cam->dir[2] > 0.9999)
    tau = cam->angle[2];
  else
    tau = cam->angle[2]+M_PI/2;

  /* Now construct final bases in frame spanned by geo-bases */
  for (kk = 0; kk < 3; kk++) {
    cam->base.x[kk] = cos(tau) * geo->x[kk] + sin(tau) * geo->y[kk];
    cam->base.y[kk] = -sin(tau) * geo->x[kk] + cos(tau) * geo->y[kk];
  }

  return;
}


/**
 * @brief Calculate actual camera basis vectors using alternative definition.
 *
 * Here, the y*' vector is taken as the projection of the (real) y
 * basis vector onto the sphere, with x*' pointing 90 degrees clockwise on
 * the sphere surface. The camera basis vectors y* and x* are then these 
 * vectors rotated counter-clockwise by the specified angle (rho).
 *
 * The degenerate points (0, -1, 0) and (0, 1, 0) are treated specially,
 * with x*'=x, y*'=z and x*'=x, y*'=-z, respectively; this keeps the bases 
 * consistent during a rotation about the x-axis.
 *
 */
void form_camera_bases_yparallel(struct st_cam* cam) {

  int kk;

  /* Set up preliminary basis vectors (without intentional camera roll) */
  struct st_base camDash;

  const double sinRho = sin(cam->angle[2]);
  const double cosRho = cos(cam->angle[2]);
  
  /* Deal with degenerate points at 'entrance and exit' of y-axis: */
  if(cam->dir[1] < -0.9999) {
    camDash.x[0] = 1.;
    camDash.x[1] = 0.;
    camDash.x[2] = 0.;
    camDash.y[0] = 0.;
    camDash.y[1] = 0.;
    camDash.y[2] = 1.;
  } else if(cam->dir[1] > 0.9999) {
    camDash.x[0] = 1.;
    camDash.x[1] = 0.;
    camDash.x[2] = 0.;
    camDash.y[0] = 0.;
    camDash.y[1] = 0.;
    camDash.y[2] = -1.;
  } else {  /* Normal, non-degenerate points */
    
    /* Set up y*' as projection of y onto sphere at camera point */
    /* (at least I think this is what this is doing...)         */
    camDash.y[0] = -(cam->dir[1]*cam->dir[0]);
    camDash.y[1] = 1.0 - (cam->dir[1]*cam->dir[1]);
    camDash.y[2] = -(cam->dir[1]*cam->dir[2]);
    normalize_vector_double(camDash.y);
    
    /* Form x*' as cross product of CamDir and y*' */
    camDash.x[0] = cam->dir[1]*camDash.y[2] - cam->dir[2]*camDash.y[1];
    camDash.x[1] = cam->dir[2]*camDash.y[0] - cam->dir[0]*camDash.y[2];
    camDash.x[2] = cam->dir[0]*camDash.y[1] - cam->dir[1]*camDash.y[0];
    verify_unit_length(camDash.x);
  }

  printf("Preliminary basis vectors determined as\n"
	 "[%g, %g, %g], [%g, %g, %g]\n", 
	 camDash.x[0], camDash.x[1], camDash.x[2],
	 camDash.y[0], camDash.y[1], camDash.y[2]);

  /* Form final basis vector by applying rho rotation  *
   * (Construct in frame spanned by preliminary bases) */
  for(kk = 0; kk < 3; kk++) {
    cam->base.x[kk] = cosRho*camDash.x[kk] + sinRho*camDash.y[kk];
    cam->base.y[kk] = cosRho*camDash.y[kk] - sinRho*camDash.x[kk];
  }
  return;
}


/**
 * @brief Normalize a (double) input vector to unit length.
 */
void normalize_vector_double(double* vec) {
  int kk;  /* Loop counter over dimensions */
  
  const double dLength = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
  for (kk = 0; kk < 3; kk++)
    vec[kk] /= dLength;
  return;
}

/**
 * @brief Normalize a (float) input vector to unit length.
 */
void normalize_vector_float(float* vec) {
  int kk;  /* Loop counter over dimensions */  

  const double dLength = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
  for (kk = 0; kk < 3; kk++)
    vec[kk] /= dLength;
  return;
}


/**
 * @brief Verify that an input vector is unit length (quit if not).
 *
 * @param vec Input vector (double).
 */
void verify_unit_length(double* vec) {
  
  int kk;

  double dLTemp = 0;
  for (kk = 0; kk < 3; kk++) 
    dLTemp += vec[kk]*vec[kk];

  if(dLTemp > 1.01 || dLTemp < 0.99) {
    printf("Unexpected value of vector=%g. Terminating.\n", dLTemp);
    exit(101);
  }
  return;
}


/* ---------- setup_projection() -------------------------------- */

/**
 * @brief Set up general projection parameters.
 *
 * This function determines whether we use perspectivic or orthogonal
 * projection and calculates the pixel dimensions (which must be square).
 * In perspectivic mode, it also calculates the image scale at unit distance.
 *
 * @param param Parameter structure pointer, of type #st_param
 * 
 * Internally accessed global variables:
 *  -- CamFOV, Xpixels, Ypixels, Xmax, Xmin
 */
void setup_projection(struct st_cam* cam) {
 
  struct st_im* image = cam->image;
  struct st_region* region = cam->region;

  /* Check if we are performing orthogonal or perspectivic projection: */
  if (cam->fov[0] > 1e-3) {
    cam->flag_perspectivic = 1;
    printf("Perspectivic projection (FOV=%g/%g), nX=%d, nY=%d\n", 
	   cam->fov[0], cam->fov[1], image->xPix, image->yPix);
        
    /* Re-calculate frame extent, at unit distance */
    region->xWidth = 1.0 * tan(cam->fov[0]);
    region->yWidth = 1.0 * tan(cam->fov[1]);
    
  } else {
    cam->flag_perspectivic = 0;
  }

  /* Set pixel length and size (at unit distance in perspectivic mode) */
  cam->pix_length = (region->xWidth)*2 / image->xPix;
  cam->pix_area = (cam->pix_length) * cam->pix_length;
  printf("Pixel area = %f\n", cam->pix_area);

  /* Sanity check to make sure pixels are square */
  {
    double pixRatio = region->yWidth*2 / image->yPix / cam->pix_length;
    if (pixRatio > 1.001 || pixRatio < 0.999) {
      printf("Pixels are non-square, which is not allowed:\n"
	     "x=%g, y=%g\n",
	     cam->pix_length, region->yWidth*2/image->yPix);
      exit(130619);
    }
  }
     
  return;
}


/* ---------- calculate_3D_map() ------------------------------------- */

/**
 * @brief Main interpolation routine to calculate 3D distributions. 
 * 
 * @param massMap (Output) Pointer-to-pointer of 3D mass map
 * @param quantMap (Output) Pointer-to-pointer (3D weighted quantity map).
 * @param cam Pointer to #st_cam camera parameter structure.
 */
void calculate_3D_map(float** p_massMap, float** p_quantMap, 
		      struct st_cam* cam, struct st_flag* flag) {

  long ii;  /* Loop counter over particles */
  double time_start = get_wall_time();

  printf("\nProjecting to 3D grid with %ld elements...\n",
	 cam->image->nPix);

  /* Allocate memory for temporary maps and initialize */
  if (flag->losvd == 0)  /* 3D maps provided as input otherwise */
    allocate_3D_maps(cam->image->nPix, p_massMap, p_quantMap);

  /* Memory now allocated, set up pointer to outputs */
  float* massMap = *p_massMap;
  float* quantMap = *p_quantMap;
  initialize_maps(massMap, quantMap, cam->image);

  /* Main loop over particles is embarassingly parallel --> OpenMP */
#pragma omp parallel for schedule(static, 1)
  for(ii = 0; ii < NumPart; ii++) {

    /*if (ii == 5068165) */
    /*printf("Particle %ld...\n", ii); */

    /* Particle position in camera frame, in pixel units */
    MyFloat pos[2];

    /* Particle position along camera axis, in physical units */
    MyFloat zeta;
    float zeta_vel;

    /* Particle's smoothing length */
    float hsml = Hsml[ii];

    /* Print progress signal if appropriate (not too often) */
    print_signal(ii);

    /*if (ii == 5068165) */
    /*printf("   ... getting coordinates...\n"); */
    
    /* Calculate particle's coordinates */
    if (get_particle_coordinates(ii, cam, flag, pos, &zeta, &zeta_vel, 
				 &hsml) < 0)
      continue;

    /*if (ii == 5068165) */
    /*printf("   ... checking FoV...\n"); */

    /* Check whether particle is within camera field of view */
    if(pos[0] + hsml < 0 || pos[0] - hsml > cam->image->xPix
       || pos[1] + hsml < 0 || pos[1] - hsml > cam->image->yPix)
      continue;

    /*if (ii == 5068165) */
    /*printf("   ... compute z layer...\n"); */
    
    /* Add particle to respective layer of 3D map */
    int iz = compute_z_level(zeta, zeta_vel, cam->region, cam->image,
			     flag);

    /*if (ii == 5068165) */
    /*printf("   ... find pixel range...\n"); */
      
    /* Calculate range of pixels covered by this particle */
    struct st_range range;
    find_pixel_range(pos, hsml, &range, cam->image);

    /*    if (ii == 5068165) */
    /*printf("   ... find kernel normalization...\n"); */
        
    /* Determine kernel normalizaton - stop if it is insignificant*/
    double wkNorm = get_kernel_normalization(pos, hsml, &range);
    if (wkNorm < 1.0e-10) 
      continue;

    /* Finish normalization by multiplying with pixel area
     * --> gives final maps in surface density */
    wkNorm *= cam->pix_area;
    if (cam->flag_perspectivic)
      wkNorm *= (zeta*zeta);

    /*printf(" ... add to map...\n"); */

    /* Actually add particle contributions to map pixels */
    add_particle_to_map(pos, Mass[ii], Quantity[ii], massMap, quantMap, 
			&range, cam->image, hsml, iz, wkNorm);    

  } /* ends loop through individual particles */
  
  printf("\n");
  printf("done (took %.2g sec.)\n\n", get_wall_time()-time_start);
    
  return;
} /* ends calculate_3d_map() */


/**
 * @brief Allocate memory for the temporary per-particle maps.
 *
 * @param nPixels Total number of pixels in 3D maps
 * @param pp_map Pointer-to-pointer for mass array
 * @param pp_quant Pointer-to-pointer for quantity array  
 */
void allocate_3D_maps(long nPixels, 
		      float** pp_map, float** pp_quant) {

  *pp_map = (float *) malloc(nPixels * sizeof(float));
  *pp_quant = (float *) malloc(nPixels * sizeof(float));
  if (!*pp_map) {
    printf("Could not allocate memory for TempMap!\n");
    exit(2203191335);
  }
  if (!*pp_quant) {
    printf("Count not allocate memory for TempQuantity!\n");
    exit(2203191336);
  }
  return;
}

/**
 * @brief Initialize temporary and output maps to zero.
 * 
 * @param mass3D The 3D mass map
 * @param quant3D The 3D mass-weighted-quantity map
 * @param im Pointer to image parameter structure (#st_im).
 */
void initialize_maps(float* mass3D, float* quant3D, struct st_im* im) {

  /* Loop counter over pixels */
  long ii, kk;
  
#pragma omp parallel for
  /* Initialise 3D maps to zero */
  for(ii = 0; ii < im->nPix; ii++) {
    mass3D[ii] = 0;
    quant3D[ii] = 0;
  }
  
#pragma omp parallel for
  /* Initialise output map to zero */
  for(ii = 0; ii < (long)im->xPix * im->yPix; ii++)
    for (kk = 0; kk < 2; kk++) {
      Value[kk + ii] = 0;
      ValueQuantity[kk + ii] = 0;
    }
  return;
}


/**
 * @brief Print signal (.) evenly throughout loop progress.
 *
 * @param n Current iteration index.
 */
static inline void print_signal(long n) {
  
  int sigNumPart = 100;
  if (NumPart < 100)
    sigNumPart = NumPart;
  
  if((n % (NumPart / sigNumPart)) == 0) {
    printf(".");
    fflush(stdout);
  }
  return;
}

/** 
 * @brief Calculate the coordinates of a particle in the camera frame 
 *
 * @param ii Index of particle to process.
 * @param cam Pointer to camera parameters (type #st_cam).
 * @param pos (Output): pointer to particle (pixel) coordinates.
 * @param zeta (Output): z coordinate in physical units
 * @param hsml (Output): (Clipped) smoothing length of the particle (pix units)
 *
 * Internally used global variable:
 *  -- P [Particle structure]
 */
static inline int get_particle_coordinates(long ii, struct st_cam* cam,
					   struct st_flag* flag,
					   MyFloat* pos, MyFloat* zeta,
					   float* zeta_vel,
					   float *hsml) { 
  
  int kk;               /* Dimension loop counter */
  double pos_sim[3];    /* Particle position in sim. frame (shifted) */

  /* Pixel size, and image width */
  MyFloat pix_length = cam->pix_length;
  MyFloat im_width[2] = {cam->region->xWidth, cam->region->yWidth};

  /* Extract camera basis vectors */
  const struct st_base base = cam->base;

  /*printf("      compute zeta...\n"); */
  *zeta = 0;
  *zeta_vel = 0;
  /* Shift particle position relative to camera centre, and find zeta
   * (physical distance along camera viewing axis) 
   * 
   * NB: Same projection is done for zeta_vel, so moving *along* the 
   *     camera viewing direction gives positive zeta_vel, moving
   *     *against* the camera viewing direction negative.  */
  for (kk = 0; kk < 3; kk++) {
    pos_sim[kk] = P[ii].Pos[kk] - cam->pos[kk];
    *zeta += pos_sim[kk] * cam->dir[kk];
  }
  
  /* Only project velocity along the line if it is required: */
  if (flag->losvd == 1)
    for (kk = 0; kk < 3; kk++)
      *zeta_vel += PVel[ii].Vel[kk] * cam->dir[kk];

  /* Keep input-z velocity if desired: */
  if (flag->losvd == 2)
    *zeta_vel = PVel[ii].Vel[2];

  /*printf("      check whether particle is beyond boundaries...\n"); */
  /* Can stop if particle is outside desired range in z-direction */
  if(*zeta < cam->region->zMin || *zeta > cam->region->zMax)
    return -1;

  if (flag->losvd == 1 || flag->losvd == 2)
    if (*zeta_vel < cam->region->vzMin || *zeta_vel > cam->region->vzMax)
      return -1;

  /*printf("      ...convert to actual distances...\n"); */
  /* In perspectivic projection, need to now convert at-unit-distance
   * quantities to actual lengths for this particle */
  if(cam->flag_perspectivic) {
    pix_length *= *zeta;
    for (kk = 0; kk < 2; kk++)
      im_width[kk] *= *zeta;
  }

  /*printf("      ...get minimum smoothing length...\n"); */
  /* Get minimum smoothing length of the particle */
  const float hMin = pix_length * 1.001/2;
  if (*hsml < hMin)
    *hsml = hMin;
  *hsml /= pix_length;  /* Convert to pixel units */

  /*printf("      ...project onto camera frame...\n"); */
  
  /* Project coordinates onto camera frame */
  pos[0] = 0;
  pos[1] = 0;
  for (kk = 0; kk < 3; kk++) {
    pos[0] += pos_sim[kk]*base.x[kk];
    pos[1] += pos_sim[kk]*base.y[kk];
  }

  /*
  printf("pos=(%g/%g)\n", pos[0], pos[1]);
  if (pos[1] > 0.005)
    exit(777);
  */

  /* And convert physical to pixel coordinates */
  for (kk = 0; kk < 2; kk++) 
    pos[kk] = (pos[kk] + im_width[kk])/pix_length;

  /*
  printf("pix_length=%g, im_width=(%g/%g)\n",
	 pix_length, im_width[0], im_width[1]);
  printf("pos=(%g/%g)\n", pos[0], pos[1]);
  if (pos[1] > 1240)
    exit(888);
  if (pos[1] < 1160)
    exit(999);
  */

  return 1;
} /* ends get_particle_coordinates() */
  
/**
 * @brief Compute the z level to which a particle belongs.
 * 
 * @param zeta Distance from camera along z* axis.
 * @param reg Pointer to #st_region structure for image limits.
 * @param image Pointer to #st_im structure for image parameters.
 * 
 * @Return iz Level (=pixel coordinate, int) of particle, increasing 
 *            *towards* the camera (so that ray tracing can be done in 
 *            order of increasing iz).
 */
static inline int compute_z_level(MyFloat zeta, float zeta_vel,
				  struct st_region* reg,
				  struct st_im* image,
				  struct st_flag* flag) {
  int zPix = image->zPix;

  int iz;
  if (flag->losvd == 0 || flag->losvd == 3)
    iz = (int)((zeta-reg->zMin)/(reg->zMax - reg->zMin) * ((float)zPix));
  else
    iz = (int)((zeta_vel-reg->vzMin)/(reg->vzMax - reg->vzMin)*((float)zPix));

  /* Clip level to within actual pixel range */
  if (iz < 0) {
    printf("WARNING: iz = %d", iz);
    iz = 0;
  }
  else if (iz >= zPix && zPix > 0)
    iz = zPix - 1;
      
  /* Invert z-order so we can build image from 0 upwards */
  /* No, don't do that, because it gives counter-intuitive 3D cubes... */
  /* iz = zPix - iz - 1; */ 

  return iz;
}

/**
 * @brief Determine range of pixels covered by a particle. 
 * 
 * @param pos Position of particle, in pixel units.
 * @param h Smoothing length of particle, in pixel units.
 * @param range (Output) pixel range, as #st_range struct.
 */
static inline void find_pixel_range(MyFloat* pos, float h, 
				    struct st_range* range,
				    struct st_im* im) {
  
  range->xMin = (int) floor(pos[0]-h+0.5); 
  range->xMax = (int) floor(pos[0]+h+0.5); 
  range->yMin = (int) floor(pos[1]-h+0.5); 
  range->yMax = (int) floor(pos[1]+h+0.5); 
  
  range->xMinFull = range->xMin < 0 ? 0 : range->xMin;
  range->xMaxFull = range->xMax >= im->xPix ? im->xPix-1 : range->xMax;
  range->yMinFull = range->yMin < 0 ? 0 : range->yMin;
  range->yMaxFull = range->yMax >= im->yPix ? im->yPix-1 : range->yMax;
  
  return;
}

/**
 * @brief Compute the normalization of the smoothing kernel.
 * 
 * @param pos Position of particle (in pixel units).
 * @param h Particle's smoothing length, in pixel units.
 * @param range Range of pixels covered by this particle.
 * 
 * @Return Kernel normalization (sum of kernel weights over all pixels).
 */
static inline double get_kernel_normalization(MyFloat* pos,
					      float h,
					      struct st_range* range) {
  
  /* Kernel normalization, start from zero */
  double sum = 0;
  const double h2 = h*h;

  int ix, iy;  /* Counters for looping over pixels */

  /*
  printf("xRange=(%d/%d), yRange=(%d/%d)\n",
	 range->xMin, range->xMax, range->yMin, range->yMax);
  */
  
  /* Loop through all pixels covered by particle, INCLUDING
   * ones off the edge - just to compute normalization! */
  for(ix = range->xMin; ix <= range->xMax; ix++)
    for(iy = range->yMin; iy <= range->yMax; iy++) {
      
      /* Pixel's distance from particle position */
      const double dx = ix + 0.5 - pos[0];
      const double dy = iy + 0.5 - pos[1];
      const double r2 = dx * dx + dy * dy;

      /*printf("ix=%d, iy=%d, r2=%g\n", ix, iy, r2); */
      
      /* Only consider particles within h from particle */
      if(r2 < h2) {
	const double r = sqrt(r2);
	const double u = r / h;

	/* Compute kernel weight at this position */
	double wk;
	if(u < 0.5)
	  wk = (2.546479089470 + 15.278874536822 * (u - 1) * u * u);
	else
	  wk = 5.092958178941 * (1.0 - u) * (1.0 - u) * (1.0 - u);
	
	sum += wk;
      } /* ends section for near-particle pixels */	
    } /* ends loop through pixels covered by this particle */
  
  return sum;
}

/**
 * @brief Actually add particle's contribution to pixels in 3D maps.
 * 
 * The map is built in units of *surface density*.
 * 
 * @param pos The particle's 2D position, in pixel units.
 * @param mass The particle's mass
 * @param quant The particle's qantity
 * @param tempMap The 3D mass (surface density) map to be incremented
 * @param tempQuant The 3D mass-weighted-quantity map to be incremented
 * @param range Pixel range covered by particle, as #st_range struct.
 * @param im Pointer to the image parameter structure (#st_im)
 * @param h Particle's smoothing length, in pixel units.
 * @param iz z layer to which to add the particle.
 * @param sumWeight The kernel normalization factor, including the 
 *                  factor for pixel area.
 */
static inline void add_particle_to_map(MyFloat* pos, float mass, float quant,
				       float* tempMap, float* tempQuant,
				       struct st_range* range, 
				       struct st_im* im, float h, int iz,
				       double sumWeight) {

  int ix, iy;   /* Pixel indices in x, y coordinates. */
  int yPix = im->yPix;
  int zPix = im->zPix;
  const double h2 = h*h;

  /* Loop through pixels to update output map -- now only 'real' ones! */
  for(ix = range->xMinFull; ix <= range->xMaxFull; ix++) {
    for(iy = range->yMinFull; iy <= range->yMaxFull; iy++) {

      const double dx = ix + 0.5 - pos[0];
      const double dy = iy + 0.5 - pos[1];
      const double r2 = dx*dx + dy*dy;

      /* Only care about pixel if it is within h from particle */
      if(r2 < h2) {
	const double r = sqrt(r2);
	const double u = r / h;
        
	/* Compute kernel weight for this pixel */
	double wk;
	if(u < 0.5)
	  wk = (2.546479089470 + 15.278874536822 * (u - 1) * u * u);
	else
	  wk = 5.092958178941 * (1.0 - u) * (1.0 - u) * (1.0 - u);

	wk /= sumWeight;

	/* Construct 1D pixel index */
	const long iPix = (long)iz + iy*zPix + ix*zPix*yPix;

	/* Actually increment maps */

	tempMap[iPix] += (float)(mass*wk);
	tempQuant[iPix] += (float)(mass*quant*wk);
      
      } /* ends section only for 'covered' pixels. */
    } /* ends loop through y-pixels */
  } /* ends loop through x-pixels */
  
  return;
} /* ends add_particle_to_map() */


/* ---------- project_3d_map() --------------------------------------- */

/**
 * @brief Project the 3D map to 2D, using sum and ray-tracing.
 * 
 * The ray-tracing image is stored in layer 1 of the output maps,
 * the simple sum image in layer 0.
 * 
 * @param tempMass The 3D mass map.
 * @param tempQuant The 3D mass-weighted quantity map.
 * @param im Pointer to #st_im structure with image parameters.
 * @param flag_perspectivic 1 if using perspec. projection, 0 otherwise.
 * @param zRange Distance between furthest and closest point from camera
 *               (only relevant if perspectivic projection is selected).
 */
void project_3D_map(float* tempMass, float* tempQuant, struct st_im* im,
		    struct st_flag* flag,
		    int flag_perspectivic, MyFloat zRange) {

  const double time_start = get_wall_time();
  int ix, iy, iz;  /* Counters for looping over pixels */
  const int zPix_out = im->zPix_out;
  const int yPix = im->yPix;
  const int zPix = im->zPix;
  const float tau = im->tau;  /* Optical depth parameter */

  printf("Produce 2D maps from 3D grid...");

#pragma omp parallel for private(iy, iz)
  for(ix = 0; ix < im->xPix; ix++) {
    for(iy = 0; iy < yPix; iy++) {

      /* 1D pixel index for 2D maps */
      long ip2 = (long)iy*zPix_out + ix*zPix_out*yPix;

      /* Start sums for final values of this 2D pixel */
      double massPix = 0, massPixAll = 0, quantPix = 0, quantPixAll = 0;

      /* Loop through individual z-pixels of current 'column' and 
       * compute total values 
       * NB: Loop backwards through z pixels, because we don't invert
       *     them anymore (so high iz == furthest away). */
      for(iz = zPix-1; iz >= 0; iz--) {

	/* 1D pixel index for 3D maps */
	long ip3 = (long)iz + iy*zPix + ix*zPix*yPix;

	const float currWeight = tempMass[ip3];
	const float currQuant = tempQuant[ip3];

	/* Surface density of current pixel. In orthogonal projection,
	 * this is simply mass/area. In perspectivic projection, we have
	 * to convert to area-at-unit-distance, or within a unit 'ray'
	 * that gets larger at larger distances */
	float mass_per_ray_sq;
	
	/* Simple integration through 3D box is straightforward */
	massPixAll += (double) currWeight;
	quantPixAll += (double) currQuant;
	
	/* Now ray tracing -- first need to find relevant surface density */
	if(flag_perspectivic) {
	  /* Distance from camera at centre of current z-pixel */ 
	  const MyFloat level_zeta = (iz+0.5) * zRange / zPix;
	  mass_per_ray_sq = currWeight * level_zeta * level_zeta;  
	} else
	  mass_per_ray_sq = currWeight;
        
	/* Compute transparency factor for pixels behind current one */
	const float alpha = exp(-mass_per_ray_sq/tau);

	/* Update ray-traced sum: 'through-shining' fraction of old + new */
	massPix = massPix*alpha + currWeight;
	quantPix = quantPix*alpha + currQuant;

	/* Also need to normalize NOW the 3D quant map, if it is 
	 * returned as output (in velocity-split mode) */
	if (flag->losvd > 0)
	  if (currWeight > 0)
	    tempQuant[ip3] /= currWeight;
		
      } /* ends loop through z pixels */

      /* Normalize quant values in current pixel (take out mass factor) */
      if(massPix > 0)
	quantPix /= massPix;
      if(massPixAll > 0)
	quantPixAll /= massPixAll;

      /* And store in output array */
      Value[ip2] = (float) massPix;
      Value[ip2+1] = (float) massPixAll;   /* `z' is fastest varying! */

      ValueQuantity[ip2] = (float) quantPix;
      ValueQuantity[ip2+1] = (float) quantPixAll;
   
    } /* ends loop through y... */
  } /* ... and x pixels */
  
  printf(" done (took %.2g sec.)\n", get_wall_time() - time_start);

  print_max_value(im->xPix * im->yPix * im->zPix_out);
  
  if (flag->losvd == 0) {
    /* Free memory used for temporary 3D maps */
    free((void*)tempMass);
    free((void*)tempQuant);
  } 

  return;

} /* ends project_3d_map() */


void print_basis_vectors(const struct st_base* const base) {

  printf("Basis vectors determined as\n"
	 "[%g, %g, %g],\n[%g, %g, %g],\n[%g, %g, %g]\n" , 
	 base->x[0], base->x[1], base->x[2], 
	 base->y[0], base->y[1], base->y[2],
	 base->z[0], base->z[1], base->z[2]);
  return;
}


void print_max_value(long nPix) {

  long ii;
  float valMax = 0;
  for (ii = 0; ii < nPix; ii++) {
    if (Value[ii] > valMax)
      valMax = Value[ii];
  }
  printf("max Value = %f\n", valMax);
  return;
}
