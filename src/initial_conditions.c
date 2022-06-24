#include "initial_conditions.h"
#include "particle_dist.h"
#include "constants.h"
#include "eos.h"

#include <string.h>
#include <math.h>
#include <stdbool.h>

void init_func_from_string(const char *str, EOSFunc eos, InitFunc *func)
{
	if (strcmp(str, "uniform") == 0)
	{
		*func = init_uniform;
	}
	else if (strcmp(str, "colliding_streams") == 0)
	{
		if (eos == EOS_isothermal) *func = init_colliding_streams_isothermal;
		else if (eos == EOS_ideal) *func = init_colliding_streams_ideal;
	}
	else if (strcmp(str, "shock_tube") == 0)
	{
		*func = init_shock_tube;
	}
	else
	{
		*func = NULL;
	}
}

/**
 * Initialise the box with particles equally spaced all the way along
 * and with zero velocity.
 */
void init_uniform(ParticleDist *dist, int i, bool verbose)
{
	if (verbose) printf("Initialising simulation with static uniform distribution\n");
	double particle_spacing = dist->box_width / (double)dist->N;
	double offset = particle_spacing / 2.0 - dist->box_width / 2.0;
	double r = (double)i * particle_spacing + offset;
	gsl_vector_set(dist->r, i, r);
	gsl_vector_set(dist->v, i, 0.0);
}

/**
 * Initialise with particles uniformly distributed along the box
 * and velocities towards its centre.
 */
void init_colliding_streams_isothermal(ParticleDist *dist, int i, bool verbose)
{
	if (verbose) printf("Initialising simulation for colliding streams test\n");
	double particle_spacing = dist->box_width / (double)dist->N;
	double offset = particle_spacing / 2.0 - dist->box_width / 2.0;
	double r = (double)i * particle_spacing + offset;
	double v = r < 0.0 ? C_S : -C_S;
	gsl_vector_set(dist->r, i, r);
	gsl_vector_set(dist->v, i, v);
}

/**
 * Initialise with particles uniformly distributed along the box
 * and velocities towards its centre. Sets u = 1/(gamma - 1) and
 * calculates the speed of sound as sqrt(gamma).
 */
void init_colliding_streams_ideal(ParticleDist *dist, int i, bool verbose)
{
	if (verbose) printf("Initialising simulation for colliding streams test\n");
	double particle_spacing = dist->box_width / (double)dist->N;
	double offset = particle_spacing / 2.0 - dist->box_width / 2.0;
	double r = (double)i * particle_spacing + offset;
	double c_s = sqrt(IDEAL_GAMMA);
	double v = r < 0.0 ? c_s : -c_s;
	gsl_vector_set(dist->r, i, r);
	gsl_vector_set(dist->v, i, v);
	gsl_vector_set(dist->u, i, 1.0 / (IDEAL_GAMMA - 1.0));
}

/**
 * Initialise with high pressure on the LHS and low pressure on
 * the RHS of the box. Particles have zero velocity.
 */
void init_shock_tube(ParticleDist *dist, int i, bool verbose)
{
	// N_l + N_r = N
	// N_l / N_r = DENSITY_RATIO
	int numberOnLeft = (float) dist->N / (1.0 + 1.0 / DENSITY_RATIO);
	if (verbose)
	{
		printf("Initialising simulation for shock tube test\n");
		// Warn if the density ratio is not achievable with the number of particles
		if (fmod(dist->N, 1.0 + 1.0 / DENSITY_RATIO) != 0.0)
			printf("[WARNING] Cannot achieve density ratio of %g:1 with %d particles\n", DENSITY_RATIO, dist->N);
	}
	// LHS
	if (i < numberOnLeft)
	{
		double spacing = (dist->box_width / 2.0) / (float)numberOnLeft;
		double offset = spacing / 2.0 - dist->box_width / 2.0;
		gsl_vector_set(dist->r, i, (float)i * spacing + offset);
		gsl_vector_set(dist->u, i, 1.0 / (IDEAL_GAMMA - 1.0));
	}
	// RHS
	else
	{
		double spacing = (dist->box_width / 2.0) / (dist->N - (float)numberOnLeft);
		double offset = spacing / 2.0;
		gsl_vector_set(dist->r, i, (float)(i - numberOnLeft) * spacing + offset);
		gsl_vector_set(dist->u, i, 0.8 / (IDEAL_GAMMA - 1.0));
	}
}
