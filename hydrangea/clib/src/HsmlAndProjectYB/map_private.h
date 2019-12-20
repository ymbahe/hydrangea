/* Internal-only function declarations */

#include "map.h"

/**
 * @brief Print signal (.) evenly throughout loop progress.
 *
 * @param n Current iteration index.
 */
static inline void print_signal(long n);

/** 
 * @brief Calculate the coordinates of a particle in the camera frame 
 *
 * @param ii Index of particle to process.
 * @param cam Pointer to camera parameters (type #st_cam).
 * @param pos (Output): pointer to particle (pixel) coordinates.
 * @param zeta (Output): z coordinate in physical units
 * @param hsml (Output): (Clipped) smoothing length of the particle (pix units)
 *
 * Internally used global variables: 
 *   -- P, Zmin, Zmax, Zpixels, CamDir, CamPos
 */
static inline int get_particle_coordinates(long ii, struct st_cam* cam,
					   struct st_flag* flag,
					   MyFloat* pos, MyFloat* zeta,
					   float* zeta_vel,
					   float *hsml);

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
				  struct st_flag* flag);

/**
 * @brief Determine range of pixels covered by a particle. 
 * 
 * @param pos Position of particle, in pixel units.
 * @param h Smoothing length of particle, in pixel units.
 * @param range (Output) pixel range, as #st_range struct.
 */
static inline void find_pixel_range(MyFloat* pos, float h, 
				    struct st_range* range, struct st_im* im);

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
				struct st_range* range);

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
			 double sumWeight);
