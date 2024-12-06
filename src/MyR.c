#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include "MyR.h"


void error(const char *message, ...) {
	va_list args;

	va_start(args, message);
	vprintf(message, args);
	va_end(args);
	exit(1);
}

void warning(const char *message, ...) {
	va_list args;

	va_start(args, message);
	vprintf(message, args);
	va_end(args);
}

void Rprintf(const char *message, ...) {
	va_list args;

	va_start(args, message);
	vprintf(message, args);
	va_end(args);
}

double lgammafn(double x) {
	return lgamma(x);
}

double unif_rand() {
	double res;
printf("unif_rand\n");
	exit(0);
	#pragma omp critical
	{
		res = (((double)rand())/((double)RAND_MAX));
	}
	return res;
//	return (((double)rand())/((double)RAND_MAX));
}

int bernoulli(double prob) {
	int res;
printf("bernoulli\n");
	exit(0);
//	#pragma omp critical
	{
		res = unif_rand()<=prob;
	}
	return res;
//	return unif_rand()<=prob;
}

int unif_int_rand(int max) {
	int res;
printf("unif_int_rand\n");
	exit(0);
//	#pragma omp critical
	{
		res = (int)floor(unif_rand()*(max+1)-0.00000000000001);
		res = (res<=max)?res:max;
	}
	return res;
	//int res = (int)floor(unif_rand()*(max+1)-0.00000000000001);
	//return (res<=max)?res:max;
}

void GetRNGstate() {
	;
}

void PutRNGstate() {
	;
}
