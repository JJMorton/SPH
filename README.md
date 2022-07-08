## SPH
This project is an implementation of Smoothed Particle Hydrodynamics (SPH) in one dimension, with adaptive timestepping and per-particle smoothing lengths.
The model includes artificial viscosity as well as conductivity, as described in [1].
There are Euler and second-order Runge-Kutta integration methods implemented and available for use via the command line arguments.
Boundary conditions are enforced with ghost particles that mirror the positions of the real particles across the two ends of the box.

### References

1. D. J. Price, “Smoothed particle hydrodynamics and magnetohydrodynamics,” Journal of
Computational Physics, vol. 231, pp. 759–794, Feb. 2012.
2. J. J. Monaghan, “Smoothed particle hydrodynamics.,” Annual Rev. Astron. Astrophys.,
vol. 30, pp. 543–574, Jan. 1992.
3. S. Rosswog, “Astrophysical smooth particle hydrodynamics,” New Astronomy Reviews,
vol. 53, pp. 78–104, Apr. 2009.
4. M. R. Bate, “The role of accretion in binary star formation.”
