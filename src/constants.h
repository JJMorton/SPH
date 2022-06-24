#ifndef CONSTANTS_H
#define CONSTANTS_H

#define VISCOUS_ALPHA 1
#define VISCOUS_BETA 2
#define THERMAL_COND 1.0 // For thermal dissipation, the nabla^2 u term
#define DENSITY_RATIO 8.0 // Between both sides of box (shock tube), LHS/RHS
#define C_S 1.0 // Speed of sound for isothermal systems
#define IDEAL_GAMMA 5.0/3.0 // Polytropic index of ideal gas

#endif // CONSTANTS_H
