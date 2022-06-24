#ifndef PARTICLE_DIST_H
#define PARTICLE_DIST_H

#include <gsl/gsl_vector.h>
#include <stdbool.h>

/**
 * A structure to hold the per-particle properties of the system, including
 * positions, velocities etc.
 */
typedef struct
{
	int N;
    double box_width;
	gsl_vector *r;     // Position
	gsl_vector *v;     // Velocity
	gsl_vector *m;     // Mass
	gsl_vector *u;     // Internal energy
	gsl_vector *rho;   // Density
	gsl_vector *p;     // Pressure
    gsl_vector *h;     // Smoothing length
    gsl_vector *omega; // Adjustment to kernel derivative due to variable smoothing length
} ParticleDist;

typedef struct
{
    int index;
	double r, v, m, u, rho, p, h, omega;
} Particle;

typedef void (*InitFunc)(ParticleDist *dist, int i, bool verbose);

/**
 * Create a distribution with N particles in a 1D box of length L, all with the same
 * provided mass, internal energy and smoothing length.
 */
ParticleDist *particle_dist_create(int N, double L, double m, double u, double h);

/**
 * Free the memory used by the particle distribution, including all the arrays and
 * such that are pointed to.
 */
void particle_dist_free(ParticleDist *dist);

/**
 * Initialise the particles in the distribution using the provided function. The
 * function will be supplied each particle index individually.
 */
void particle_dist_init(ParticleDist *dist, InitFunc init, bool verbose);

/**
 * Returns the properties of the particle with index i.
 * i can range from -N to 2N-1, where values outide of the range 0 to N-1
 * represent ghost particles on the left and right sides of the box.
 * Changing the properties of the retuned particle has no effect.
 */
void particle_dist_getparticle(ParticleDist *dist, int i, Particle *p);

/**
 * Set the properties of the particle with index i.
 */
void particle_dist_setparticle(ParticleDist *dist, int i, Particle *p);

#endif // PARTICLE_DIST_H
