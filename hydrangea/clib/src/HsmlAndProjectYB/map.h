#ifndef _MAP_H
#define _MAP_H

/**
 * @brief Struct to hold sin/cos values of angles.
 */
struct st_trig {
  double cosTheta;
  double cosPhi;
  double sinTheta;
  double sinPhi;
};

/** 
 * @brief Struct to hold a set of three basis vectors. 
 */
struct st_base {
  double x[3];
  double y[3];
  double z[3];
};


/**
 * @brief Struct to hold output arrays (camera bases)
 */
struct st_out {
  float* x;
  float* y;
  float* z;
};


/**
 * @brief Image parameters.
 */
struct st_im {
  /* Number of (total!) pixels in x*, y*, z* directions */ 
  int xPix;
  int yPix;
  int zPix;
  int zPix_out;  /* z dimension of (2D) output map. */
  long nPix;  /* Total number of (3D) image pixels */

  /* Optical depth parameter for combining z-slices (-1 to leave in 3D) */
  float tau;
};

/**
 * @brief Particle selection region. 
 */
struct st_region {

  /* x and y specify the box (half-)width in camera frame (orthogonal)
   * or half-width at unit distance from camera (perspectivic): */
  float xWidth;
  float yWidth;

  /* For z, store separate min and max for perspectivic projections */
  float zMin;
  float zMax;

  /* Store min and max velocity, for velocity-splitting mode */
  float vzMin;
  float vzMax;

};


/**
 * @brief Holds parameters for camera position and orientation.
 */
struct st_cam {
  /* Position of camera in simulation frame (actual coordinates) */
  MyFloat pos[3];

  /* Viewing direction of the camera, in simulation frame coordinates. */
  float dir[3];

  /* Viewing direction of the camera specified in angles. The first two 
   * elements (theta, phi) give the position on a (hypothetical, 
   * infinitesimal) unit sphere around the camera position, at which the 
   * camera sits looking inwards. Theta is measured downwards from the +z axis,
   * phi counter-clockwise from +x, both in radians.
   * E.g. for a view down the z-axis, dir = (0, 0, -1) and angle = (0, 0, 0); 
   * along the x axis corresponds to dir = (1, 0, 0) and 
   * angle = (pi/2, pi, 0). The third angle, rho, specifies the camera roll
   * angle, measured counter-clockwise from the standard configuration 
   * at the chosen position (see roll_mode).
   */ 
  float angle[3];
 
  /* Camera field of view, in radians (in x* and y*, respectively) */
  float fov[2];

  /* Flag to indicate orthogonal (0) or perspectivic (1) projection */
  int flag_perspectivic;

  /* Set of camera basis vectors x*, y*, and z*. */
  struct st_base base;

  /* Flag to specify zero-roll orientation of the camera:
   *   0: standard geographic orientation (N is up, E is right).
   *   1: parallel transport from the North Pole [0, 0, 1].
   *   2: y* parallel to projection of y onto sphere at camera angles.
   */
  int roll_mode;

  /* Pointer to corresponding image parameters, for convenience. */
  struct st_im* image;

  /* Pointer to corresponding image region, for convenience. */
  struct st_region* region;

  /* Area of one pixel (at unit distance in perspectivic mode) */
  float pix_area;

  /* Side length of one pixel (which are enforced to be square) */
  float pix_length;
};

struct st_range {
  int xMin;
  int xMax;
  int xMinFull;
  int xMaxFull;
  int yMin;
  int yMax;
  int yMinFull;
  int yMaxFull;
};


/**
 * @brief Struct to hold configuration flags.
 */
struct st_flag {
  int use_hsml;   /* 0: compute hsml from scratch, 1: take supplied input */
  int hsml_only;  /* 0: compute image, 1: only compute hsml. */

  /* Over-allocation factor for tree-nodes relative to number of particles. 
   * It generally seems to work to have 5 if N < 1e6 and 2.5 otherwise. */
  float tree_alloc_fac; 

  /* Flag to control velocity-splitting mode.
   * 0 --> no velocity split (slice by cam-z coordinate instead).
   * 1 --> velocity split along input z direction
   * 2 --> velocity split along cam-z direction
   */
  int losvd;
}; 

/* =============== Function declarations ======================== */

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
 */
void make_image(struct st_cam* cam, struct st_out* out,
		struct st_flag* flag);

/**
 * @brief Compute the basis vectors of the camera frame (x*, y*, z*)
 *
 * @param cam Pointer to camera property structure (#st_cam)
 * @param out Pointer to #st_out struct holding output camera bases.
 * @param flag Pointer to #st_flag struct holding general flags.
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
void get_basis_vectors(struct st_cam* cam, struct st_out* out);

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
void get_cam_dir_from_angles(struct st_cam* cam, struct st_trig* trig);

/**
 * @brief Compute relevant angles from CamDir.
 *
 * See description of get_cam_dir_from_angles for angle/direction conventions.
 *
 * @param trig Pointer to #st_trig structure, to hold cos/sin of angles.
 * @param cam Pointer to camera parameter structure (#st_cam)
 */
void get_angles_from_cam_dir(struct st_trig* trig, struct st_cam* cam);

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
			   struct st_cam* cam);

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

void form_camera_bases(struct st_base* geo, struct st_cam* cam);

/**
 * @brief Calculate the camera bases in standard geographical orientation
 *        (North is up, East is right).
 *
 * @param geo Pointer to the geographic basis vectors (#st_base struct).
 * @param cam Pointer to the camera parameter struct (#st_cam).
 */

void form_geographic_camera_bases(struct st_base* geo, struct st_cam* cam);

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
void form_camera_bases_yparallel(struct st_cam* cam);

/**
 * @brief Normalize a (double) input vector to unit length.
 */
void normalize_vector_double(double* vec);

/**
 * @brief Normalize a (float) input vector to unit length.
 */
void normalize_vector_float(float* vec);

/**
 * @brief Verify that an input vector is unit length (quit if not).
 *
 * @param vec Input vector (double).
 */
void verify_unit_length(double* vec);

/**
 * @brief Set up general projection parameters.
 *
 * @param param Parameter structure pointer, of type #st_param
 * 
 * Internally accessed global variables:
 *  -- CamFOV, Xpixels, Ypixels, Xmax, Xmin
 */
void setup_projection(struct st_cam* cam);

/**
 * @brief Main interpolation routine to calculate 3D distributions. 
 * 
 * @param massMap (Output) 3D mass map
 * @param quantMap (Output) 3D weighted quantity map.
 * @param cam Pointer to #st_cam camera parameter structure.
 */
void calculate_3D_map(float** massMap, float** quantMap, struct st_cam* cam,
		      struct st_flag* flag);

/**
 * @brief Allocate memory for the temporary per-particle maps.
 *
 * @param nPixels Total number of pixels in 3D maps
 * @param pp_map Pointer-to-pointer for mass array
 * @param pp_quant Pointer-to-pointer for quantity array  
 */
void allocate_3D_maps(long nPixels, 
		      float** pp_map, float** pp_quant);

/**
 * @brief Initialize temporary and output maps to zero.
 * 
 * @param mass3D The 3D mass map
 * @param quant3D The 3D mass-weighted-quantity map
 * @param im Pointer to image parameter structure (#st_im).
 */
void initialize_maps(float* mass3D, float* quant3D, struct st_im* im);


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
		    int flag_perspectivic, MyFloat zRange);


/**
 * @brief Return the wall-clock time in seconds.
 */
double get_wall_time();



#endif   /* _MAP_H */




void print_basis_vectors(const struct st_base* const base);
void print_max_value(long nPix);
