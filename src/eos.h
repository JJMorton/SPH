#ifndef EOS_H
#define EOS_H

#include "particle_dist.h"

/**
 * A function representing an equation of state, returning pressure
 */
typedef double (*EOSFunc)(Particle *p);

/**
 * Get the EOS function corresponding with its name as a string
 */
void EOS_from_string(const char *str, EOSFunc *func);

double EOS_isothermal(Particle *p);
double EOS_ideal(Particle *p);

#endif // EOS_H
