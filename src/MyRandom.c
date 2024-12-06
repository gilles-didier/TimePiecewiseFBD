#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <gsl_rng.h>

void my_rand_init() {
	gsl_rng_env_setup();
}

void *my_rand_get_data() {
	return (void*) gsl_rng_alloc(gsl_rng_default);
}

void my_rand_free_data(void *data) {
	gsl_rng_free((gsl_rng*) data);
}

double my_rand_unif_unit(void *data) {
	return gsl_rng_uniform((gsl_rng*) data);
}

int my_rand_bern(void *data, double prob) {
	return  my_rand_unif_unit(data)<=prob;
}

int my_rand_unif_range(void *data, int max) {
	return gsl_rng_uniform_int((gsl_rng*) data, max+1);
}
