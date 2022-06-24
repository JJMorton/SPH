#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <stdbool.h>
#include <stdarg.h>
#include <math.h>

#include "particle_dist.h"
#include "initial_conditions.h"
#include "density_estimator.h"
#include "constants.h"

/**
 * Test if any of the given pointers are NULL
 */
bool containsNULL(int argc, ...)
{
	va_list va;
	va_start(va, argc);
	bool anyNULL = false;
	for (int i = 0; i < argc; i++)
	{
		void *ptr = va_arg(va, void *);
		anyNULL |= ptr == NULL;
	}
	va_end(va);
	return anyNULL;
}

/**
 * Create a particle distribution, allocating all the memory for the vectors.
 */
ParticleDist *particle_dist_create(int N, double L, double m, double u, double h)
{
	ParticleDist *dist = malloc(sizeof(*dist));
	if (dist == NULL)
	{
		fprintf(stderr, "Failed to allocate memory for particle distribution\n");
		return NULL;
	}
	dist->N = N;
	dist->box_width = 1.0;
	dist->r = gsl_vector_alloc(N);
	dist->v = gsl_vector_alloc(N);
	dist->m = gsl_vector_alloc(N);
	dist->u = gsl_vector_alloc(N);
	dist->rho = gsl_vector_alloc(N);
	dist->p = gsl_vector_alloc(N);
	dist->h = gsl_vector_alloc(N);
	dist->omega = gsl_vector_alloc(N);
	if (containsNULL(6, dist->r, dist->v, dist->u, dist->rho, dist->p, dist->h, dist->omega))
	{
		fprintf(stderr, "Failed to allocate memory for particle distribution\n");
		particle_dist_free(dist);
		abort();
	}
	gsl_vector_set_all(dist->r, 0.0);
	gsl_vector_set_all(dist->v, 0.0);
	gsl_vector_set_all(dist->m, m);
	gsl_vector_set_all(dist->u, u);
	gsl_vector_set_all(dist->rho, 0.0);
	gsl_vector_set_all(dist->p, 0.0);
	gsl_vector_set_all(dist->h, h);
	gsl_vector_set_all(dist->omega, 0.0);
	return dist;
}

/**
 * Free all the memory used by the arrays
 */
void particle_dist_free(ParticleDist *dist)
{
	if (dist->r != NULL) gsl_vector_free(dist->r);
	if (dist->v != NULL) gsl_vector_free(dist->v);
	if (dist->m != NULL) gsl_vector_free(dist->m);
	if (dist->u != NULL) gsl_vector_free(dist->u);
	if (dist->rho != NULL) gsl_vector_free(dist->rho);
	if (dist->p != NULL) gsl_vector_free(dist->p);
	if (dist->p != NULL) gsl_vector_free(dist->h);
	if (dist->p != NULL) gsl_vector_free(dist->omega);
	free(dist);
}

/**
 * Initialise the particle distribution using the provided function.
 */
void particle_dist_init(ParticleDist *dist, InitFunc init, bool verbose)
{
	for (int i = 0; i < dist->N; i++)
	{
		// Only print out details for the first particle
		init(dist, i, verbose && i == 0);
	}
}

void particle_dist_getparticle(ParticleDist *dist, int i, Particle *p)
{
	int index;
	if (i >= 0 && i < dist->N)
	{
		// Requesting a real particle
		index = i;
		p->r = gsl_vector_get(dist->r, index);
		p->v = gsl_vector_get(dist->v, index);
	}
	else if (i >= -dist->N && i < 0)
	{
		// Requesting a ghost particle on the LHS of the box
		index = -i - 1;
		p->r = -gsl_vector_get(dist->r, index) - dist->box_width;
		p->v = -gsl_vector_get(dist->v, index); // Reflective
		// p->v = 0.0;
	}
	else if (i >= dist->N && i < 2 * dist->N)
	{
		// Requesting a ghost particle on the RHS of the box
		index = 2 * dist->N - i - 1;
		p->r = -gsl_vector_get(dist->r, index) + dist->box_width;
		p->v = -gsl_vector_get(dist->v, index); // Reflective
		// p->v = 0.0;
	}
	else
	{
		printf("WARNING: Invalid particle index passed to particle_dist_getparticle()\n");
		return;
	}

	p->index = i;
	p->p = gsl_vector_get(dist->p, index);
	p->m = gsl_vector_get(dist->m, index);
	p->u = gsl_vector_get(dist->u, index);
	p->rho = gsl_vector_get(dist->rho, index);
	p->h = gsl_vector_get(dist->h, index);
	p->omega = gsl_vector_get(dist->omega, index);
}

void particle_dist_setparticle(ParticleDist *dist, int i, Particle *p)
{
	if (i >= 0 && i < dist->N)
	{
		gsl_vector_set(dist->r, i, p->r);
		gsl_vector_set(dist->v, i, p->v);
		gsl_vector_set(dist->m, i, p->m);
		gsl_vector_set(dist->u, i, p->u);
		gsl_vector_set(dist->rho, i, p->rho);
		gsl_vector_set(dist->p, i, p->p);
		gsl_vector_set(dist->h, i, p->h);
		gsl_vector_set(dist->omega, i, p->omega);
	}
	else
	{
		printf("WARNING: Invalid particle index passed to particle_dist_setparticle()\n");
		return;
	}
}
