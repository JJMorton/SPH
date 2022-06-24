#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

#include "integrator.h"
#include "constants.h"
#include "particle_dist.h"
#include "density_estimator.h"
#include "eos.h"


/*
 * The following weighting functions are to be used with the smooth_around_particle()
 * function in density_estimator.h
 */

/**
 * The P_a / Omega_a rho_a^2 term in the expression for particle a's acceleration
 * in equation 30 of Price (2010).
 * Smooth with smoothing length of particle a.
 */
double weighting_dvdt_a(Particle *pa, Particle *pb)
{
	double dvdt = -pb->m * pa->p / (pa->rho * pa->rho * pa->omega);
	return dvdt;
}

/**
 * The P_b / Omega_b rho_b^2 term in the expression for particle a's acceleration
 * in equation 30 of Price (2010).
 * Smooth with smoothing length of particle b.
 */
double weighting_dvdt_b(Particle *pa, Particle *pb)
{
	double dvdt = -pb->m * pb->p / (pb->rho * pb->rho * pb->omega);
	return dvdt;
}

/**
 * The artificial viscosity acting on particle a, when an ideal gas equation
 * of state is being used.
 * Smooth with smoothing length of particles a and b
 */
double weighting_dvdt_viscous_ideal(Particle *pa, Particle *pb)
{
	double r_ab = pa->r - pb->r;
	double v_ab = pa->v - pb->v;
	double rho_bar = 0.5 * (pa->rho + pb->rho);
	double r_unit = r_ab >= 0.0 ? 1.0 : -1.0;
	double c_sa = sqrt(IDEAL_GAMMA * pa->p / pa->rho);
	double c_sb = sqrt(IDEAL_GAMMA * pb->p / pb->rho);
	double v_sig = (r_ab * v_ab) > 0.0 ? 0.0 : 0.5 * (c_sa + c_sb - VISCOUS_BETA * v_ab * r_unit);
	double dvdt = 0.5 * pb->m * VISCOUS_ALPHA * v_sig * v_ab * r_unit / rho_bar;
	return dvdt;
}

/**
 * The artificial viscosity acting on particle a, when an isothermal equation
 * of state is being used.
 * Smooth with smoothing length of particles a and b
 */
double weighting_dvdt_viscous_isothermal(Particle *pa, Particle *pb)
{
	double r_ab = pa->r - pb->r;
	double v_ab = pa->v - pb->v;
	double rho_bar = 0.5 * (pa->rho + pb->rho);
	double r_unit = r_ab >= 0.0 ? 1.0 : -1.0;
	double v_sig = (r_ab * v_ab) > 0.0 ? 0.0 : (C_S - 0.5 * VISCOUS_BETA * v_ab * r_unit);
	double dvdt = 0.5 * pb->m * VISCOUS_ALPHA * v_sig * v_ab * r_unit / rho_bar;
	return dvdt;
}

/**
 * Half the change in internal energy of particle a.
 * Smooth with smoothing length of particle a
 */
double weighting_dudt_ideal_a(Particle *pa, Particle *pb)
{
	double r_ab = pa->r - pb->r;
	double v_ab = pa->v - pb->v;
	double dudt = pa->p / (pa->rho * pa->rho * pa->omega) * pb->m * v_ab;

	// Viscous term
	double rho_bar = 0.5 * (pa->rho + pb->rho);
	double r_unit = r_ab >= 0.0 ? 1.0 : -1.0;
	double c_sa = sqrt(IDEAL_GAMMA * pa->p / pa->rho);
	double c_sb = sqrt(IDEAL_GAMMA * pb->p / pb->rho);
	double v_sig = (r_ab * v_ab) > 0.0 ? 0.0 : 0.5 * (c_sa + c_sb - VISCOUS_BETA * v_ab * r_unit);
	dudt += -0.25 * pb->m / rho_bar * VISCOUS_ALPHA * v_sig * v_ab*v_ab * r_unit;

	// Thermal dissipation term
	double v_sig_u = sqrt(fabs(pa->p - pb->p) / rho_bar);
	// double v_sig_u = fabs(v_ab);
	dudt += 0.5 * pb->m / rho_bar * THERMAL_COND * v_sig_u * (pa->u - pb->u) * r_unit;

	return dudt;
}

/**
 * Half the change in internal energy of particle a.
 * Smooth with smoothing length of particle b
 */
double weighting_dudt_ideal_b(Particle *pa, Particle *pb)
{
	double r_ab = pa->r - pb->r;
	double v_ab = pa->v - pb->v;
	double dudt = 0.0;

	// Viscous term
	double rho_bar = 0.5 * (pa->rho + pb->rho);
	double r_unit = r_ab >= 0.0 ? 1.0 : -1.0;
	double c_sa = sqrt(IDEAL_GAMMA * pa->p / pa->rho);
	double c_sb = sqrt(IDEAL_GAMMA * pb->p / pb->rho);
	double v_sig = (r_ab * v_ab) > 0.0 ? 0.0 : 0.5 * (c_sa + c_sb - VISCOUS_BETA * v_ab * r_unit);
	dudt += -0.25 * pb->m / rho_bar * VISCOUS_ALPHA * v_sig * v_ab*v_ab * r_unit;

	// Thermal dissipation term
	double v_sig_u = sqrt(fabs(pa->p - pb->p) / rho_bar);
	// double v_sig_u = fabs(v_ab);
	dudt += 0.5 * pb->m / rho_bar * THERMAL_COND * v_sig_u * (pa->u - pb->u) * r_unit;

	return dudt;
}


/**
 * Calculates the acceleration and velocity of each particle, saves the
 * result in the given vectors `drdt` and `dvdt`
 */
void EOM(ParticleDist *dist, DensityEstimator *density_est, gsl_vector *drdt, gsl_vector *dvdt, gsl_vector *dudt)
{
	int N = dist->N;
	if (dvdt->size != N || drdt->size != N || dudt->size != N)
	{
		fprintf(stderr, "Size mismatch between particle number and length of EOM vectors\n");
		return;
	}

	/*
	 * For each particle we will sum up all the contributions from the equations of motion,
	 * after smoothing them across the neighbour particles with smooth_around_particle().
	 */

	for (int i = 0; i < N; i++)
	{
		// Set drdt as velocity
		Particle p;
		particle_dist_getparticle(dist, i, &p);
		gsl_vector_set(drdt, i, p.v);

		KernelFunc kernel = density_est->kernel_diffx;

		// Calculate dvdt without viscosity
		double dvdt_i =
			smooth_around_particle(weighting_dvdt_a, kernel, SMOOTHLENGTH_A, dist, i) +
			smooth_around_particle(weighting_dvdt_b, kernel, SMOOTHLENGTH_B, dist, i);

		// Add the viscosity contribtion according to the EOS used
		if (density_est->EOS == EOS_ideal)
		{
			dvdt_i +=
				smooth_around_particle(weighting_dvdt_viscous_ideal, kernel, SMOOTHLENGTH_A, dist, i) +
				smooth_around_particle(weighting_dvdt_viscous_ideal, kernel, SMOOTHLENGTH_B, dist, i);
		}
		else
		{
			dvdt_i +=
				smooth_around_particle(weighting_dvdt_viscous_isothermal, kernel, SMOOTHLENGTH_A, dist, i) +
				smooth_around_particle(weighting_dvdt_viscous_isothermal, kernel, SMOOTHLENGTH_B, dist, i);
		}
		gsl_vector_set(dvdt, i, dvdt_i);

		// Calculate dudt if using ideal gas EOS
		if (density_est->EOS == EOS_ideal)
		{
			double dudt_i =
				smooth_around_particle(weighting_dudt_ideal_a, kernel, SMOOTHLENGTH_A, dist, i) +
				smooth_around_particle(weighting_dudt_ideal_b, kernel, SMOOTHLENGTH_B, dist, i);
			gsl_vector_set(dudt, i, dudt_i);
		}
	}
}


/**
 * Given a string, returns the integration method to be used, or NULL if the string's invalid
 */
void integrator_method_from_string(const char *str, IntegratorMethod *method)
{
	if (strcmp(str, "euler") == 0)
	{
		*method = integrator_method_euler;
	}
	else if (strcmp(str, "rk2") == 0)
	{
		*method = integrator_method_rk2;
	}
	else
	{
		*method = NULL;
	}
}

/**
 * First-order numerical integration with the Euler method:
 * f(t + h) = f(t) + h*df(t)/dt
 */
void integrator_method_euler(ParticleDist *dist, DensityEstimator *density_est, double dt)
{
	gsl_vector *drdt = gsl_vector_alloc(dist->N);
	gsl_vector *dvdt = gsl_vector_alloc(dist->N);
	gsl_vector *dudt = gsl_vector_alloc(dist->N);
	EOM(dist, density_est, drdt, dvdt, dudt);
	gsl_blas_daxpy(dt, drdt, dist->r);
	gsl_blas_daxpy(dt, dvdt, dist->v);
	gsl_blas_daxpy(dt, dudt, dist->u);
	gsl_vector_free(drdt);
	gsl_vector_free(dvdt);
	gsl_vector_free(dudt);
	update_particle_props(density_est, dist);
}

/**
 * Second-order Runge-Kutta numerical integration:
 * k_1 = h * df(t)/dt
 * k_2 = h * df(t + h/2)/dt
 * f(t + h) = f(t) + k_2
 */
void integrator_method_rk2(ParticleDist *dist, DensityEstimator *density_est, double dt)
{

	// The first order term, Euler's method
	// Evaluate the acceleration and velocity at the current time
	gsl_vector *delta_r1 = gsl_vector_alloc(dist->N);
	gsl_vector *delta_v1 = gsl_vector_alloc(dist->N);
	gsl_vector *delta_u1 = gsl_vector_alloc(dist->N);
	EOM(dist, density_est, delta_r1, delta_v1, delta_u1);


	// The second order term, evaluate EOM halfway through timestep
	// The function gsl_blas_daxpy(a, x, y) sets y = ax + y

	// Evaluate velocity and position at t + h/2
	gsl_blas_daxpy(0.5 * dt, delta_r1, dist->r);
	gsl_blas_daxpy(0.5 * dt, delta_v1, dist->v);
	gsl_blas_daxpy(0.5 * dt, delta_u1, dist->u);

	// Update density, internal energy, pressure, smoothing length etc. of the particles
	update_particle_props(density_est, dist);

	// Calculate the changes to position and velocity at this new time
	gsl_vector *delta_r2 = gsl_vector_alloc(dist->N);
	gsl_vector *delta_v2 = gsl_vector_alloc(dist->N);
	gsl_vector *delta_u2 = gsl_vector_alloc(dist->N);
	EOM(dist, density_est, delta_r2, delta_v2, delta_u2);


	// Step back by h/2 to get back to time t
	gsl_blas_daxpy(-0.5 * dt, delta_r1, dist->r);
	gsl_blas_daxpy(-0.5 * dt, delta_v1, dist->v);
	gsl_blas_daxpy(-0.5 * dt, delta_u1, dist->u);

	// Add the second order terms to the positions and velocities
	gsl_blas_daxpy(dt, delta_r2, dist->r);
	gsl_blas_daxpy(dt, delta_v2, dist->v);
	gsl_blas_daxpy(dt, delta_u2, dist->u);

	// Update density, internal energy, pressure, smoothing length etc. of the particles
	update_particle_props(density_est, dist);

	gsl_vector_free(delta_r1);
	gsl_vector_free(delta_v1);
	gsl_vector_free(delta_u1);
	gsl_vector_free(delta_r2);
	gsl_vector_free(delta_v2);
	gsl_vector_free(delta_u2);
}


/**
 * Returns the required timestep of the simulation, based off maximum accelerations and the CFL condition
 */
double integrator_calc_timestep(Integrator *integrator, ParticleDist *dist, DensityEstimator *density_est)
{
	gsl_vector *drdt = gsl_vector_alloc(dist->N);
	gsl_vector *dvdt = gsl_vector_alloc(dist->N);
	gsl_vector *dudt = gsl_vector_alloc(dist->N);
	EOM(dist, density_est, drdt, dvdt, dudt);
	gsl_vector_free(drdt);
	gsl_vector_free(dudt);

	// Find the particle with the smallest required timestep and use that
	double dt = INFINITY;
	for (int a = 0; a < dist->N; a++)
	{
		Particle p;
		particle_dist_getparticle(dist, a, &p);

		// Makes timestep smaller when accelerations get large
		double dt_f = sqrt(p.h / fabs(gsl_vector_get(dvdt, a)));

		// CFL condition, ensures that the timestep is smaller than the time it takes for
		// a particle to cross 'one point of resolution' (in this case one smoothing length)
		double c_a = density_est->EOS == EOS_ideal ? (sqrt(IDEAL_GAMMA * p.p / p.rho)) : C_S;
		double dt_c = p.h / c_a;

		// Take the smaller of the two timestep conditions
		double dt_min = fmin(dt_f, dt_c);
		dt = fmin(dt, dt_min);
	}

	gsl_vector_free(dvdt);

	return 0.25 * dt;
}

/**
 * Step the system forward in time by one timestep of length `dt`
 */
void integrator_step(Integrator *integrator, ParticleDist *dist, DensityEstimator *density_est, double dt)
{
	// Dump to file if reached the specified interval
	if (integrator->last_dump + integrator->dump_interval - integrator->time <= 0)
	{
		// First check that everything is in order, with no NANs or INFs or anything
		for (int i = 0; i < dist->N; i++)
		{
			double rho = gsl_vector_get(dist->rho, i);
			double r = gsl_vector_get(dist->r, i);
			double v = gsl_vector_get(dist->v, i);
			double u = gsl_vector_get(dist->u, i);
			double h = gsl_vector_get(dist->h, i);
			if (isnan(rho) || isnan(r) || isnan(v) || isnan(u) || isnan(h) || isinf(rho) || isinf(r) || isinf(v) || isinf(u) || isinf(h))
			{
				printf("[ERROR] Detected non-finite quantity, aborting.\n");
				integrator->proceed = false;
				return;
			}
			if (r > dist->box_width / 2.0 || r < -dist->box_width / 2.0)
			{
				printf("[WARNING] Detected particle(s) outside of defined box\n");
				break;
			}
		}

		// Print the relevant particle properties to the log file
		printf("Dumping at t=%lf, dt=%lf\n", integrator->time, dt);
		gsl_vector * const arrs[6] = { dist->r, dist->rho, dist->v, dist->u, dist->h, dist-> p };
		for (int arr = 0; arr < 6; arr++)
		{
			for (int i = 0; i < dist->N; i++)
			{
				fprintf(integrator->logfile, "%E\t", gsl_vector_get(arrs[arr], i));
			}
			fprintf(integrator->logfile, "\n");
		}

		integrator->last_dump = integrator->time;
	}

	// Check for the simulation end
	if (integrator->time >= integrator->duration)
	{
		printf("Reached simulation end.\n");
		integrator->proceed = false;
		return;
	}

	// Integrate
	double until_next_dump = integrator->last_dump + integrator->dump_interval - integrator->time;
	double until_end = integrator->duration - integrator->time;
	double timestep = fmin(fmin(dt, until_next_dump), until_end);
	if (timestep == 0.0)
	{
		printf("[WARNING] Attempted to integrate a timestep of dt=0 (until_next_dump=%.3E, until_end=%.3E)\n", until_next_dump, until_end);
		return;
	}
	integrator->method(dist, density_est, timestep);
	integrator->time += timestep;
	integrator->iterations++;
}
