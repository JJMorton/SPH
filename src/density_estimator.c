#include <math.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <string.h>
#include <stdbool.h>

#include "constants.h"
#include "density_estimator.h"
#include "particle_dist.h"

/*
 * Will evaluate SUM_b [ f(a, b) * k(x_a - x_b, h) ].
 * See header file for more information.
 */
double smooth_around_particle(double (*func)(Particle *, Particle *), KernelFunc kernel, SmoothType s, ParticleDist *dist, int i)
{
	Particle pa;
	particle_dist_getparticle(dist, i, &pa);

	double sum = func(&pa, &pa) * kernel(0.0, pa.h);

	// Iterate outwards from particle a in both directions
	// Stop once the kernel function drops to zero
	for (int step = -1; step <= 1; step += 2)
	{
		for (int b = i + step; b < 2 * dist->N && b >= -dist->N; b += step)
		{
			Particle pb;
			particle_dist_getparticle(dist, b, &pb);

			// Get the smoothing length to use
			double h;
			switch (s)
			{
				case SMOOTHLENGTH_A: h = pa.h; break;
				case SMOOTHLENGTH_B: h = pb.h; break;
				case SMOOTHLENGTH_AVG: h = 0.5 * (pa.h + pb.h); break;
			}

			// Evaluate the kernel and the supplied function
			double k = kernel(pa.r - pb.r, h);
			if (k == 0.0) break;

			// Add their product to the sum
			sum += func(&pa, &pb) * k;
		}
	}

	return sum;
}

// https://www.desmos.com/calculator/a88sro5gfd
/**
 * The kernel used in density estimation.
 * This is a cubic spline function, dropping to zero at x/h == 2.
 * It is normalised such that k(x=0, h) = 1
 */
double density_kernel(double x, double h)
{
	double q = fabs(x) / h;
	double sigma = 2.0 / 3.0;
	if (q < 1.0)
	{
		return sigma / h * (1.0 - 1.5*q*q + 0.75*q*q*q);
	}
	else if (q < 2.0)
	{
		double z = 2.0 - q;
		return sigma / h * 0.25 * z*z*z;
	}
	return 0.0;
}

/**
 * d/dx of the density estimator kernel
 */
double density_kernel_diffx(double x, double h)
{
	double q = fabs(x) / h;
	double sigma = 2.0 / 3.0;
	double dqdx = (x < 0.0 ? -1.0 : 1.0) / h;
	if (q < 1.0)
	{
		return dqdx * sigma / h * (-2.0*1.5*q + 3.0*0.75*q*q);
	}
	else if (q < 2.0)
	{
		double z = 2.0 - q;
		return -dqdx * sigma / h * 3.0 * 0.25 * z*z;
	}
	return 0.0;
}

/**
 * d/dh of the density estimator kernel
 */
double density_kernel_diffh(double x, double h)
{
	double q = fabs(x) / h;
	double sigma = 2.0 / 3.0;
	if (q < 1.0)
	{
		return sigma / (h*h) * (-1.0 + 9.0/2.0 * q*q - 3.0 * q*q*q);
	}
	else if (q < 2.0)
	{
		return sigma / (h*h) * (-2.0 + 6.0 * q - 9.0/2.0 * q*q + q*q*q);
	}
	return 0.0;
}

/**
 * Return the density estimating kernel to use, given a string describing it
 */
void density_kernel_from_string(const char *str, DensityEstimator *density_est)
{
	if (strcmp(str, "cubic_spline") == 0)
	{
		density_est->kernel = density_kernel;
		density_est->kernel_diffx = density_kernel_diffx;
		density_est->kernel_diffh = density_kernel_diffh;
	}
	else
	{
		density_est->kernel = density_est->kernel_diffx = density_est->kernel_diffh = NULL;
	}
}

/**
 * Particle b's contribution to the density at particle a, to be smoothed by the density kernel
 */
double weighting_density(Particle *a, Particle *b)
{
	return b->m;
}

/**
 * Calculate the pressure and density at the location of each particle
 * and save in the `p` and `rho` arrays of the particle distribution.
 */
void update_particle_props(DensityEstimator *density_est, ParticleDist *dist)
{
	// Store the densities that we calculate at the position of each particle, so we only calculate them once
	for (int i = 0; i < dist->N; i++)
	{
		Particle p;
		particle_dist_getparticle(dist, i, &p);

		if (density_est->var_smoothlength)
		{
			// eta specifies smoothing length in units of mean particle spacing
			double eta = 1.3;

			double epsilon = 0.0001 * dist->box_width / dist->N;
			int max_iter = 30;
			int iter = 0;

			/*
			 * Newton-Raphson method to simultaneously find density and smoothing length.
			 * x_n+1 = x_n - f(x_n) / f'(x_n).
			 * In this case, x == h and f == rho
			 * 
			 * We are solving the equation rho(h) - eta * m/h = 0
			 */
			for (iter = 0; iter < max_iter; iter++)
			{
				// Calculate x_n, f(x) and f'(x)
				p.rho = smooth_around_particle(weighting_density, density_est->kernel, SMOOTHLENGTH_A, dist, i);
				double rho_dh = smooth_around_particle(weighting_density, density_est->kernel_diffh, SMOOTHLENGTH_A, dist, i);
				double f = p.rho - eta * p.m / p.h;
				double f_diff = rho_dh + eta * p.m / (p.h * p.h);

				// The change in x
				double change = -f / f_diff;

				// Recalculate density with new smoothing length before we go into next iteration
				p.h = p.h + change;
				particle_dist_setparticle(dist, i, &p);

				// Finish root finding if we have converged
				if (fabs(change) <= epsilon) break;
			}

			// Print warnings if it seems like things went wrong
			if (iter == max_iter - 1) fprintf(stderr, "[WARNING] Reached maximum iterations (%d) for smoothing length estimation\n", max_iter);
			if (p.h <= 0 || fabs(p.h) > dist->box_width) printf("[WARNING] Smoothing length diverged from expected values (h = %g)\n", p.h);

			// Omega is an adjustment to d/dx terms in the equations of motion, best to just calculate it now
			p.omega = 1.0 + p.h / p.rho * smooth_around_particle(weighting_density, density_est->kernel_diffh, SMOOTHLENGTH_A, dist, i);
		}
		else
		{
			p.omega = 1.0;
			p.rho = smooth_around_particle(weighting_density, density_est->kernel, SMOOTHLENGTH_A, dist, i);
		}

		// Calculate the pressure and save in the arrays
		p.p = density_est->EOS(&p);
		particle_dist_setparticle(dist, i, &p);
	}
}
