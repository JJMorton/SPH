#ifndef DENSITY_ESTIMATOR_H
#define DENSITY_ESTIMATOR_H

#include <stdbool.h>

#include "particle_dist.h"
#include "eos.h"

/**
 * Type definition for the density estimator kernel function (and its derivatives)
 */
typedef double (*KernelFunc)(double x, double h);

/**
 * A structure to hold parameters necessary for density & pressure estimation
 */
typedef struct
{
    KernelFunc kernel;
    KernelFunc kernel_diffx;
    KernelFunc kernel_diffh;
    bool var_smoothlength; // Whether to use variable smoothing lengths
    EOSFunc EOS;
} DensityEstimator;

/**
 * Which particle's smoothing length to use in the sum.
 */
typedef enum
{
    SMOOTHLENGTH_A,
    SMOOTHLENGTH_B,
    SMOOTHLENGTH_AVG
} SmoothType;

/**
 * Populates the kernel and kernel_diff fields in the provided density estimator struct
 */
void density_kernel_from_string(const char *str, DensityEstimator *density_est);

/**
 * Provide a function f(a, b) and a density estimator kernel k(x, h).
 * Will evaluate SUM_b [ f(a, b) * k(x_a - x_b, h) ].
 * The function will be given a particle and one of its neighbours, and should return
 * this neighbour's contribution to the quantity being calculated at position a.
 * By specifying a SmoothType, you can choose which particle's smoothing length
 * to use.
 */
double smooth_around_particle(double (*func)(Particle *, Particle *), KernelFunc, SmoothType, ParticleDist *, int i);

// The density estimator kernels, take physical position x and a smoothing length h
double density_kernel(double x, double h);
double density_kernel_diffx(double x, double h);
double density_kernel_diffh(double x, double h);

void update_particle_props(DensityEstimator *density_est, ParticleDist *dist);

#endif // DENSITY_ESTIMATOR_H
