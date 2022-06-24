#include "eos.h"
#include "constants.h"
#include "particle_dist.h"

#include <string.h>

void EOS_from_string(const char *str, EOSFunc *func)
{
	if (strcmp(str, "isothermal") == 0)
	{
		*func = EOS_isothermal;
	}
	else if (strcmp(str, "ideal") == 0)
	{
		*func = EOS_ideal;
	}
	else
	{
		*func = NULL;
	}
}

double EOS_isothermal(Particle *p)
{
	return C_S * C_S * p->rho;
}

double EOS_ideal(Particle *p)
{
	return (IDEAL_GAMMA - 1.0) * p->u * p->rho;
}
