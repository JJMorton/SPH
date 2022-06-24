#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <stdbool.h>

#include "particle_dist.h"
#include "density_estimator.h"
#include "eos.h"

typedef void (*IntegratorMethod)(ParticleDist *dist, DensityEstimator *density_est, double dt);
/**
 * Given a string, returns the integration method to be used, or NULL if the string's invalid
 */
void integrator_method_from_string(const char *str, IntegratorMethod *method);

/**
 * One timestep of numerical integration using the Euler method.
 */
void integrator_method_euler(ParticleDist *dist, DensityEstimator *density_est, double dt);

/**
 * One timestep of numerical integration using the second order Runge-Kutta method.
 */
void integrator_method_rk2(ParticleDist *dist, DensityEstimator *density_est, double dt);

/**
 * A structure to hold parameters for the integrator
 */
struct integrator
{
	// Interval at which to dump to the file
	double dump_interval;

	double last_dump;

	int iterations;

	// Current world time
	double time;

	// World time to stop at
	double duration;

	// Whether the simulation should continue or not
	bool proceed;

	// The file which the state should be dumped to
	FILE *logfile;

	// The integration method (Euler, Runge-Kutta, etc.)
	IntegratorMethod method;

};
typedef struct integrator Integrator;

/**
 * Step the integrator forward in time by dt, dumping to the log file if necessary
 */
void integrator_step(Integrator *integrator, ParticleDist *dist, DensityEstimator *density_est, double dt);

/**
 * Returns the required timestep of the simulation, based off maximum accelerations and the CFL condition
 */
double integrator_calc_timestep(Integrator *integrator, ParticleDist *dist, DensityEstimator *density_est);

#endif // INTEGRATOR_H
