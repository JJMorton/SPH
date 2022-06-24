#ifndef INITIAL_CONDITIONS_H
#define INITIAL_CONDITIONS_H

#include <stdbool.h>

#include "particle_dist.h"
#include "eos.h"

/**
 * A function representing an initial condition setter
 */
typedef void (*InitFunc)(ParticleDist *dist, int i, bool verbose);

/**
 * Get the initialisation function corresponding with its name as a string
 */
void init_func_from_string(const char *str, EOSFunc eos, InitFunc *func);

// Initialisation functions for various different systems
void init_uniform(ParticleDist *dist, int i, bool verbose);
void init_colliding_streams_isothermal(ParticleDist *dist, int i, bool verbose);
void init_colliding_streams_ideal(ParticleDist *dist, int i, bool verbose);
void init_shock_tube(ParticleDist *dist, int i, bool verbose);

#endif // INITIAL_CONDITIONS_H
