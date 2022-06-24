/*
 * 
 * This program is an implementation of Smoothed Particle Hydrodynamics,
 * with variable smoothing lengths, adaptive timestepping and aritifcial
 * viscosity and conductivity.
 * 
 */

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
#include "tests.h"


void print_usage(const char *progname)
{
	fprintf(stderr,
		"Usage: %s -f log_file\n"
		"  [-d dump_interval]: the interval at which the simulation state is output to the log file (default 0.01)\n"
		"  [-N num_particles]: the number of particles to create (default 200)\n"
		"  [-s smooth_length]: the width of the smoothing kernel (variable if not specified)\n"
		"  [-h timestep]: the integration timestep to use (adaptive if not specified)\n"
		"  [-t duration]: the length of the simulation (default 1.0)\n"
		"  [-w width]: the width of the box to simulate (default 1.0)\n"
		// "  [-k kernel]: the density estimation kernel to use (only 'cubic_spline' is available)\n"
		"  [-m method]: the integration method to use, either 'euler' or 'rk2' (default 'rk2')\n"
		"  [-i initial_conditions]: the initial conditions to set up, either 'uniform', 'colliding streams' or 'shock_tube' (default 'shock_tube')\n"
		"  [-e EOS]: the equation of state to use, either 'isothermal' or 'ideal' (default 'isothermal')\n"
		"  [-T]: run tests, none of the other parameters have an efffect\n"
		"  [-?]: show this usage message\n"
	, progname);
}

int main(int argc, char *argv[])
{

	/*
	 * Parse command line options and open relevant files
	 */

	FILE  *logfile = NULL;
	int    N = 200;
	double smoothlength = 0.0;
	double timestep = 0.0;
	double duration = 1;
	double dump_interval = 0.01;
	char  *density_kernel = "cubic_spline";
	char  *integration_method = "rk2";
	char  *initial_conds = "shock_tube";
	char  *eos = "isothermal";
	double box_width = 1.0;
	bool   run_tests = false;
	bool   validargs = true;
	char   opt;

	// Options with a colon ':' after them take an argument
	// While there are still options to be read, getopt stores the option letter in 'opt'
	while ((opt = getopt(argc, argv, "f:N:s:h:t:d:k:m:i:e:w:T?")) != -1)
	{
		// After calling getopt(), the argument provided is stored in optarg
		switch (opt)
		{
			case 'f':
				// User provided the input file, attempt to open it for reading
				logfile = fopen(optarg, "w");
				if (logfile == NULL)
				{
					fprintf(stderr, "Failed to open or create file '%s'\n", optarg);
				}
				break;
			case 'N':
				// User specified the number of particles to create
				sscanf(optarg, " %i\n", &N);
				if (N <= 0)
				{
					fprintf(stderr, "Number of particles cannot be <= 0\n");
					validargs = false;
				}
				break;
			case 's':
				// User specified smoothing length
				sscanf(optarg, " %lf\n", &smoothlength);
				if (smoothlength <= 0.0)
				{
					fprintf(stderr, "Smoothing length cannot be <= 0\n");
					validargs = false;
				}
				break;
			case 'h':
				// User specified integration timestep
				sscanf(optarg, " %lf\n", &timestep);
				if (timestep <= 0.0)
				{
					fprintf(stderr, "Timestep cannot be <= 0\n");
					validargs = false;
				}
				break;
			case 't':
				// User specified simulation duration
				sscanf(optarg, " %lf\n", &duration);
				if (duration < 0.0)
				{
					fprintf(stderr, "Duration cannot be < 0\n");
					validargs = false;
				}
				break;
			case 'd':
				// User specified output interval
				sscanf(optarg, " %lf\n", &dump_interval);
				if (dump_interval <= 0.0)
				{
					fprintf(stderr, "Dump interval cannot be <= 0\n");
					validargs = false;
				}
				break;
			case 'k':
				// User specified density estimator kernel
				density_kernel = optarg;
				break;
			case 'm':
				// User specified integration method
				integration_method = optarg;
				break;
			case 'i':
				// User specified initial conditions
				initial_conds = optarg;
				break;
			case 'e':
				// User specified equation of state
				eos = optarg;
				break;
			case 'w':
				sscanf(optarg, " %lf\n", &box_width);
				if (box_width <= 0.0)
				{
					fprintf(stderr, "Box width must be a positive float");
					validargs = false;
				}
				break;
			case 'T':
				// Run the performance tests
				run_tests = true;
				break;
			case '?':
				print_usage(argv[0]);
				return 0;
			default:
				validargs = false;
		}
	}

	if (logfile == NULL)
	{
		printf("Required argument: provide a file to log to with the -f option.\n");
		validargs = false;
	}

	// Only continue if we have all the needed arguments, and they are all valid
	if (!validargs)
	{
		if (logfile != NULL) fclose(logfile);
		printf("Run '%s -?' to print the usage page.\n", argv[0]);
		exit(1);
	}

	if (run_tests)
	{
		tests_run_all(logfile);
		fclose(logfile);
		return 0;
	}

	// If we have variable smoothing, initialise the smoothing lengths to average
	// particle spacing as a sane initial guess
	bool var_smoothlength = smoothlength == 0.0;
	smoothlength = var_smoothlength ? (box_width / (double)N) : smoothlength;
	// If the user didn't specify a timestep, use adaptive timestepping
	bool adap_timestep = timestep == 0.0;

	// Create the distribution of particles
	ParticleDist *dist = particle_dist_create(N, box_width, 1.0, 1.0 / (IDEAL_GAMMA - 1.0), smoothlength);

	// Create a structure to hold integration bits and pieces
	Integrator integ =
	{
		.duration = duration,
		.dump_interval = dump_interval,
		.last_dump = -dump_interval,
		.iterations = 0,
		.time = 0.0,
		.logfile = logfile,
		.method = NULL,
		.proceed = true
	};

	// Create a structure of the parameters and methods for density estimation
	DensityEstimator density_est =
	{
		.var_smoothlength = var_smoothlength,
		.EOS = NULL,
		.kernel = NULL,
		.kernel_diffx = NULL,
		.kernel_diffh = NULL
	};

	// Set density estimator kernel to use
	density_kernel_from_string(density_kernel, &density_est);
	if (density_est.kernel == NULL)
	{
		fprintf(stderr, "'%s' is not a valid density estimation kernel\n", density_kernel);
		return 1;
	}

	// Set integration method
	integrator_method_from_string(integration_method, &integ.method);
	if (integ.method == NULL)
	{
		fprintf(stderr, "'%s' is not a valid integration method\n", integration_method);
		return 1;
	}

	// Set equation of state
	EOS_from_string(eos, &density_est.EOS);
	if (density_est.EOS == NULL)
	{
		fprintf(stderr, "'%s' is not a valid equation of state\n", eos);
		return 1;
	}

	// Set initial conditions
	dist->box_width = box_width;
	InitFunc init_func;
	init_func_from_string(initial_conds, density_est.EOS, &init_func);
	if (init_func == NULL)
	{
		fprintf(stderr, "'%s' is not a valid initial condition\n", initial_conds);
		return 1;
	}
	particle_dist_init(dist, init_func, true);

	// Save the simulation parameters to the file header
	fprintf(
		logfile,
		"N=%d\tvar_smoothlength=%d\tadap_timestep=%d\tduration=%E\tbox_width=%E\tdump_interval=%E\tkernel=%s\tmethod=%s\tinit=%s\tEOS=%s",
		N, var_smoothlength, adap_timestep, duration, box_width, dump_interval, density_kernel, integration_method, initial_conds, eos
	);
	if (!var_smoothlength) fprintf(logfile, "\tsmoothlength=%E", smoothlength);
	if (!adap_timestep) fprintf(logfile, "\ttimestep=%E", timestep);
	if (density_est.EOS == EOS_ideal) fprintf(logfile, "\tgamma=%E", IDEAL_GAMMA);
	fprintf(logfile, "\n");

	// Begin integration
	printf("Integrating...\n");
	// Calculate the initial densities
	update_particle_props(&density_est, dist);
	while (integ.proceed)
	{
		// Calculate the required timestep for the current state of the simulation
		double dt = adap_timestep ? integrator_calc_timestep(&integ, dist, &density_est) : timestep;
		// Step forward in time by one timestep
		integrator_step(&integ, dist, &density_est, dt);
	}

	// Free the memory we used
	particle_dist_free(dist);

	fclose(logfile);

	return 0;
}
