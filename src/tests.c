#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <getopt.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "constants.h"
#include "particle_dist.h"
#include "density_estimator.h"
#include "initial_conditions.h"
#include "integrator.h"
#include "eos.h"


void test_sod_density(FILE *logfile, bool var_smoothlength)
{
	printf("Calculating density estimation for Sod shock tube...\n");
	printf("Using ");
	if (var_smoothlength) printf("variable"); else printf("fixed");
	printf(" smoothing lengths\n");

	int N = 90;
	double box_width = 1.0;
	double initial_h = box_width / (double)N;
	ParticleDist *dist = particle_dist_create(N, box_width, 1.0, 1.0 / (IDEAL_GAMMA - 1.0), initial_h);
	DensityEstimator density_est =
	{
		.var_smoothlength = var_smoothlength,
		.EOS = EOS_ideal,
		.kernel = density_kernel,
		.kernel_diffx = density_kernel_diffx,
		.kernel_diffh = density_kernel_diffh
	};
	particle_dist_init(dist, init_shock_tube, false);
	update_particle_props(&density_est, dist);
	fprintf(logfile, "position,density\n");
	for (int i = 0; i < N; i++)
	{
		Particle p;
		particle_dist_getparticle(dist, i, &p);
		fprintf(logfile, "%E,%E\n", p.r, p.rho);
	}
}

void test_density_big_o(FILE *logfile, bool var_smoothlength)
{
	printf("Using ");
	if (var_smoothlength) printf("variable"); else printf("fixed");
	printf(" smoothing lengths\n");

	for (int N = 100; N <= 2000; N += 100)
	{
		printf("Timing density estimation for %d particles...\n", N);

		double box_width = 1.0;
		double initial_h = 1.3 * box_width / (double)N;
		ParticleDist *dist = particle_dist_create(N, box_width, 1.0, 1.0 / (IDEAL_GAMMA - 1.0), initial_h);
		DensityEstimator density_est =
		{
			.var_smoothlength = var_smoothlength,
			.EOS = EOS_ideal,
			.kernel = density_kernel,
			.kernel_diffx = density_kernel_diffx,
			.kernel_diffh = density_kernel_diffh
		};
		particle_dist_init(dist, init_uniform, false);

		int iterations = 0;
		double total_time = 0.0;
		while (total_time < 4.0 && iterations < 2000)
		{
			particle_dist_init(dist, init_uniform, false);
			// Set the smoothing length back to what it was originally
			gsl_vector_set_all(dist->h, initial_h);
			clock_t tic = clock();
			update_particle_props(&density_est, dist);
			clock_t toc = clock();
			total_time += (double)(toc - tic) / CLOCKS_PER_SEC;
			iterations++;
		}
		double duration = total_time / (double)iterations;
		printf("took %g seconds (averaged over %d iterations).\n", duration, iterations);
		particle_dist_free(dist);

		fprintf(logfile, "%d,%E\n", N, duration);
	}
}

void wait_for_file_copy(FILE *logfile)
{
	fflush(logfile);
	printf("Copy the output file to keep it and press enter to continue");
	getchar();
	freopen(NULL, "w", logfile);
}

void tests_run_all(FILE *logfile)
{
	printf("Running performance tests...\n");
	test_density_big_o(logfile, false);
	wait_for_file_copy(logfile);
	test_density_big_o(logfile, true);
	wait_for_file_copy(logfile);
	test_sod_density(logfile, false);
	wait_for_file_copy(logfile);
	test_sod_density(logfile, true);
}
