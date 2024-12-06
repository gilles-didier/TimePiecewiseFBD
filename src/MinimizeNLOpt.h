#ifndef MinimizeNLOptF
#define MinimizeNLOptF

#include "PiecewiseModel.h"
#include "FossilInt.h"
#include "Tree.h"

typedef struct NLOPT_OPTION {
    double infSpe, supSpe, infExt, supExt, infFos, supFos, tolOptim;
    int trials, maxIter;
} TypeNLOptOption;

typedef struct ESTIMATION {
    TypePiecewiseModelParam param;
    double logLikelihood;
} TypeEstimation;

#ifdef __cplusplus
extern "C" {
#endif


typedef double TypeLikelihoodSetTreeFosFunction(TypeTree **, int nTree, TypeFossilFeature *, TypePiecewiseModelParam *);

void fprintNLoptOption(FILE *f, TypeNLOptOption *option);
void sprintNLoptOption(char *buffer, TypeNLOptOption *option);
void fprintNLoptOptionTag(FILE *f, TypeNLOptOption *option);
void fscanNLoptOptionTag(FILE *f, TypeNLOptOption *option);
int minimizePiecewiseParamFromSetTreeFossil(TypeLikelihoodSetTreeFosFunction *f, TypeTree **tree, int nTree, TypeFossilFeature *fos,TypeSamplingCurrent sampCurrentInit,  TypeSamplingScheme sampScheme, TypeNLOptOption *option, TypeEstimation *estim);
int minimizePiecewiseParamFromSetTreeFossilInt(TypeLikelihoodSetTreeFosFunction *f, TypeTree **tree, int nTree, TypeFossilIntFeature *fint, TypeFossilFeature *foinit, TypeSamplingCurrent sampCurrentInit,  TypeSamplingScheme sampScheme, TypeNLOptOption *option, TypeEstimation *estim);
#ifdef __cplusplus
}
#endif

#endif
