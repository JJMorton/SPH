#include <stdio.h>
#ifndef TESTS_H
#define TESTS_H

#include <stdlib.h>
#include <stdbool.h>

/**
 * A couple of tests to produce plots with.
 */

void test_sod_density(FILE *logfile, bool var_smoothlength);
void test_density_big_o(FILE *logfile, bool var_smoothlength);

void tests_run_all(FILE *logfile);

#endif // TESTS_H
