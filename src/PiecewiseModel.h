#ifndef PiecewiseModelF
#define PiecewiseModelF

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include "Tree.h"
#include "Fossil.h"

#define RINFTY 1E99

typedef struct MODEL_PARAM {
	double birth, death, fossil, sampling;
} TypeModelParam;



/*between times startTime[i] and startTime[i+1] parameters are param[i]*/

typedef struct PIECEWISE_MODEL_PARAM {
	int size;
	double *startTime;
	TypeModelParam *param;
} TypePiecewiseModelParam;


typedef struct SAMPLING_SCHEME {
	TypePiecewiseModelParam param;
	int sizeParamTime, sizeParamBirth, sizeParamDeath, sizeParamFossil, sizeParamSampling, *indexTime, *schemeTime, *schemeBirth, *schemeDeath, *schemeFossil, *schemeSampling;
} TypeSamplingScheme;


typedef struct SAMPLING_CURRENT {
	double *currentTime, *currentBirth, *currentDeath, *currentFossil, *currentSampling;
} TypeSamplingCurrent;


#ifdef __cplusplus
extern "C" {
#endif

TypePiecewiseModelParam simple2piecewise(TypeModelParam *param, double startTime, double endTime);
int getPieceIndex(double v, TypePiecewiseModelParam *param);
void printPiecewiseModel(FILE *f, TypePiecewiseModelParam *param);
TypePiecewiseModelParam readPiecewiseModelParam(FILE *f);
void printSamplingScheme(FILE *f, TypeSamplingScheme sc);
TypeSamplingScheme readSamplingScheme(FILE *f);
TypeSamplingCurrent getSamplingCurrent(TypeSamplingScheme sampScheme);
void freeSamplingScheme(TypeSamplingScheme *sampScheme);
void freeSamplingCurrent(TypeSamplingCurrent *sampCurrent);
void setParamFromSamplingCurrent(TypePiecewiseModelParam *param, TypeSamplingCurrent sampCurrent, TypeSamplingScheme sampScheme);
TypePiecewiseModelParam getPiecewiseParamFromSamplingScheme(TypeSamplingScheme sampScheme);
void freePiecewiseParam(TypePiecewiseModelParam *param);
#ifdef __cplusplus
}
#endif

#endif
