#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <nlopt.h>

#include "Utils.h"
#include "MyR.h"

#include "MinimizeNLOpt.h"


//#define NLOPT_ALGO NLOPT_GN_ISRES
//#define NLOPT_ALGO NLOPT_GN_ESCH
//#define NLOPT_ALGO NLOPT_LN_BOBYQA
//#define NLOPT_ALGO NLOPT_LN_COBYLA
//#define NLOPT_ALGO NLOPT_AUGLAG
//#define NLOPT_ALGO NLOPT_GN_DIRECT_L
//#define NLOPT_ALGO NLOPT_LN_SBPLX
//#define NLOPT_ALGO NLOPT_GN_CRS2_LM
//#define NLOPT_ALGO NLOPT_LN_NEWUOA
//#define NLOPT_ALGO NLOPT_G_MLSL_LDS
//#define NLOPT_ALGO NLOPT_LN_PRAXIS

/**/
#define NLOPT_ALGO NLOPT_LN_NELDERMEAD

//#define  NLOPT_ALGO NLOPT_LD_CCSAQ

//#define  NLOPT_ALGO NLOPT_LD_MMA

//#define  NLOPT_ALGO NLOPT_LD_SLSQP

//#define  NLOPT_ALGO NLOPT_LD_LBFGS


//#define  NLOPT_ALGO NLOPT_LD_TNEWTON_PRECOND_RESTART

//#define  NLOPT_ALGO NLOPT_LD_VAR2

#define MINVAL 0.01
#define INFTY 1E99
#define RINFTY 1E99
#define DEF 10
#define MIN_VAL 0.000001
#define TOLERANCE_CONSTRAINT 0.000000001
#define TOLERANCE_OPTIM 0.001


typedef struct MINIMIZATION_SET_TREE_FOS_DATA {
	TypeTree **tree;
	TypeFossilFeature *fos;
	int nTree;
	TypePiecewiseModelParam param;
	TypeSamplingCurrent sampCurrent;
	TypeSamplingScheme sampScheme;
	TypeLikelihoodSetTreeFosFunction *likelihood;
} TypeMinimizationSetTreeFossilData;


typedef struct MINIMIZATION_SET_TREE_FOS_INT_DATA {
	TypeTree **tree;
	TypeFossilFeature *fos;
	TypeFossilIntFeature *fint;
	int nTree, *index, *revert, sizeFossil, nParam;
	TypePiecewiseModelParam param;
	TypeSamplingCurrent sampCurrent;
	TypeSamplingScheme sampScheme;
	TypeLikelihoodSetTreeFosFunction *likelihood;
} TypeMinimizationSetTreeFossilIntData;

#define TAG_SPE "SPE"
#define TAG_EXT "EXT"
#define TAG_FOS "FOS"
#define TAG_TRI "TRI"
#define TAG_TOL "TOL"
#define TAG_ITE "ITE"
#define SIZE_TAG 20
#define SIZE_VAL 100

void fprintNLoptOptionTag(FILE *f, TypeNLOptOption *option) {
    fprintf(f, ":%s [%lE;%lE]\n", TAG_SPE, option->infSpe, option->supSpe);
    fprintf(f, ":%s [%lE;%lE]\n", TAG_EXT, option->infExt, option->supExt);
    fprintf(f, ":%s [%lE;%lE]\n", TAG_FOS, option->infFos, option->supFos);
    fprintf(f, ":%s %d\n", TAG_TRI, option->trials);
    fprintf(f, ":%s %lE\n", TAG_TOL, option->tolOptim);
    fprintf(f, ":%s %d\n", TAG_ITE, option->maxIter);
}

void fscanNLoptOptionTag(FILE *f, TypeNLOptOption *option) {
    char c, tag[SIZE_TAG+1], val[SIZE_VAL+1];
    for(c=fgetc(f); c!=EOF && isspace(c); c=fgetc(f));
    while(c == ':') {
        int i;
        c=fgetc(f);
        for(i=0; c!=EOF && !isspace(c) && i<SIZE_TAG; c=fgetc(f))
            tag[i++] = c;
        tag[i] = '\0';
        if(i>=SIZE_TAG) {
            fprintf(stderr, "Error when reading an optimizer options file - Tag too long:\n%s...\n", tag);
            exit(1);
        }
        for(; c!=EOF && isspace(c); c=fgetc(f));
        for(i=0; c!=EOF && !isspace(c) && i<SIZE_VAL; c=fgetc(f))
            val[i++] = c;
        val[i] = '\0';
        if(i>=SIZE_VAL) {
            fprintf(stderr, "Error when reading an optimizer options file - value too long:\n%s...\n", val);
            exit(1);
        }
        if(strcmp(tag, TAG_SPE) == 0)
            toInterval(val, &(option->infSpe), &(option->supSpe));
        if(strcmp(tag, TAG_EXT) == 0)
            toInterval(val, &(option->infExt), &(option->supExt));
        if(strcmp(tag, TAG_FOS) == 0)
            toInterval(val, &(option->infFos), &(option->supFos));
        if(strcmp(tag, TAG_TRI) == 0)
            option->trials = atoi(val);
        if(strcmp(tag, TAG_TOL) == 0)
            option->tolOptim = atof(val);
        if(strcmp(tag, TAG_ITE) == 0)
            option->maxIter = atoi(val);
        for(; c!=EOF && isspace(c); c=fgetc(f));
    }
}

void fprintNLoptOption(FILE *f, TypeNLOptOption *option) {
    fprintf(f, "Speciation rates are sampled in [%.2lE:%.2lE]\n", option->infSpe, option->supSpe);
    fprintf(f, "Extinction rates are sampled in [%.2lE:%.2lE]\n", option->infExt, option->supExt);
    fprintf(f, "Fossil rates are sampled in [%.2lE:%.2lE]\n", option->infFos, option->supFos);
    fprintf(f, "Optimizer runs %d trials and stops with tolerance %.lE or after more than %d iterations.\n", option->trials, option->tolOptim, option->maxIter);
}

void sprintNLoptOption(char *buffer, TypeNLOptOption *option) {
    buffer += sprintf(buffer, "Speciation rates are sampled in [%.2lE:%.2lE]\n", option->infSpe, option->supSpe);
    buffer += sprintf(buffer, "Extinction rates are sampled in [%.2lE:%.2lE]\n", option->infExt, option->supExt);
    buffer += sprintf(buffer, "Fossil rates are sampled in [%.2lE:%.2lE]\n", option->infFos, option->supFos);
    buffer += sprintf(buffer, "Optimizer runs %d trials and stops with tolerance %.lE or after more than %d iterations.\n", option->trials, option->tolOptim, option->maxIter);
}



void setSamplingCurrentFromVector(TypeSamplingCurrent *sampCurrent, double *x, int start, TypeSamplingScheme sampScheme) {
	int i, ind = start;
	for(i=0; i<sampScheme.sizeParamTime; i++)
		sampCurrent->currentTime[i] = x[ind++];
	for(i=0; i<sampScheme.sizeParamBirth; i++)
		sampCurrent->currentBirth[i] = x[ind++];
	for(i=0; i<sampScheme.sizeParamDeath; i++)
		sampCurrent->currentDeath[i] = x[ind++];
	for(i=0; i<sampScheme.sizeParamFossil; i++)
		sampCurrent->currentFossil[i] = x[ind++];
	for(i=0; i<sampScheme.sizeParamSampling; i++)
		sampCurrent->currentSampling[i] = x[ind++];
		
	//for(i=0; i<sampScheme.sizeParamTime; i++) {
		//printf("time %d %d %.2lf\n", i, ind, x[ind]);
		//sampCurrent->currentTime[i] = x[ind++];
	//}
	//for(i=0; i<sampScheme.sizeParamBirth; i++){
		//printf("birth %d %d %.2lf\n", i, ind, x[ind]);
		//sampCurrent->currentBirth[i] = x[ind++];
	//}
	//for(i=0; i<sampScheme.sizeParamDeath; i++){
		//printf("death %d %d %.2lf\n", i, ind, x[ind]);
		//sampCurrent->currentDeath[i] = x[ind++];
	//}
	//for(i=0; i<sampScheme.sizeParamFossil; i++){
		//printf("fos %d %d %.2lf\n", i, ind, x[ind]);
		//sampCurrent->currentFossil[i] = x[ind++];
	//}
	//for(i=0; i<sampScheme.sizeParamSampling; i++)
		//sampCurrent->currentSampling[i] = x[ind++];
}

void setVectorFromSamplingCurrent(double *x, int start, TypeSamplingCurrent sampCurrent, TypeSamplingScheme sampScheme) {
	int i, ind = start;
	for(i=0; i<sampScheme.sizeParamBirth; i++)
		x[ind++] = sampCurrent.currentBirth[i];
	for(i=0; i<sampScheme.sizeParamDeath; i++)
		x[ind++] = sampCurrent.currentDeath[i];
	for(i=0; i<sampScheme.sizeParamFossil; i++)
		x[ind++] = sampCurrent.currentFossil[i];
	for(i=0; i<sampScheme.sizeParamSampling; i++)
		x[ind++] = sampCurrent.currentSampling[i];
}

void setSamplingCurrentFossilFromVector(TypeTree **tree, TypeFossilIntFeature *fint, TypeFossilFeature *fos, int *revert, int sizeFossil, int sizeNode, TypeSamplingCurrent *sampCurrent, double *x, TypeSamplingScheme sampScheme) {
	int i, ind = 0, *tmp, n, f;
	size_t *index_tmp;
	for(i=0; i<sizeFossil; i++) {
		fos->fossilList[i].time = x[revert[ind++]];
		fos->fossilList[i].prec = fint->fossilIntList[i].prec;
	}
    for(i=0; i<sizeNode; i++)
        fos->fossil[i] = fint->fossilInt[i];
	index_tmp = qsortindex(fos->fossilList, fos->size, sizeof(TypeFossilList), compareFossilList);
    for(n=0; n<sizeNode; n++)
        if(fos->fossil[n]>=0)
            fos->fossil[n] = index_tmp[fos->fossil[n]];
    for(f=0; f<fos->size; f++)
        if(fos->fossilList[f].prec>=0)
            fos->fossilList[f].prec = index_tmp[fos->fossilList[f].prec];
    free((void*) index_tmp);
    tmp = (int*) malloc((fos->size+1)*sizeof(int));
    for(n=0; n<sizeNode; n++) {
        int ind = 0;
        for(f=fos->fossil[n]; f>=0; f=fos->fossilList[f].prec)
            tmp[ind++] = f;
        if(ind>0) {
            qsort(tmp, ind, sizeof(int), compareInt);
            fos->fossilList[tmp[0]].prec = -1;
            for(f=1; f<ind; f++)
                fos->fossilList[tmp[f]].prec = tmp[f-1];
            fos->fossil[n] = tmp[ind-1];
        }
    }
    free((void*) tmp);
//printf("x[0] %lf\n", x[0]);
//fprintFossilFeature(stdout, fos, tree[0]->name, tree[0]->size);
//for(i=0; i<sizeFossil; i++)
	//printf("x[%d] %lf\n", i, x[i]);
//for(i=sizeFossil; i<sizeFossil+3; i++)
	//printf("\tx[%d] %lf\n", i, x[i]);
//printf("\n\n");
	//for(i=0; i<sampScheme.sizeParamTime; i++) {
		//printf("time %d %d %.2lf\n", i, ind, x[ind]);
		//sampCurrent->currentTime[i] = x[ind++];
	//}
	//for(i=0; i<sampScheme.sizeParamBirth; i++){
		//printf("birth %d %d %.2lf\n", i, ind, x[ind]);
		//sampCurrent->currentBirth[i] = x[ind++];
	//}
	//for(i=0; i<sampScheme.sizeParamDeath; i++){
		//printf("death %d %d %.2lf\n", i, ind, x[ind]);
		//sampCurrent->currentDeath[i] = x[ind++];
	//}
	//for(i=0; i<sampScheme.sizeParamFossil; i++){
		//printf("fos %d %d %.2lf\n", i, ind, x[ind]);
		//sampCurrent->currentFossil[i] = x[ind++];
	//}
	//for(i=0; i<sampScheme.sizeParamSampling; i++)
		//sampCurrent->currentSampling[i] = x[ind++];

	for(i=0; i<sampScheme.sizeParamTime; i++)
		sampCurrent->currentTime[i] = x[ind++];
	for(i=0; i<sampScheme.sizeParamBirth; i++)
		sampCurrent->currentBirth[i] = x[ind++];
	for(i=0; i<sampScheme.sizeParamDeath; i++)
		sampCurrent->currentDeath[i] = x[ind++];
	for(i=0; i<sampScheme.sizeParamFossil; i++)
		sampCurrent->currentFossil[i] = x[ind++];
	for(i=0; i<sampScheme.sizeParamSampling; i++)
		sampCurrent->currentSampling[i] = x[ind++];
}

void setVectorFromSamplingCurrentFossil(double *x, TypeFossilFeature *fos, int *index, int sizeFossil, TypeSamplingCurrent sampCurrent, TypeSamplingScheme sampScheme) {
	int i, ind = 0;
	for(i=0; i<sizeFossil; i++)
		x[ind++] = fos->fossilList[index[i]].time;
	for(i=0; i<sampScheme.sizeParamTime; i++)
		x[ind++] = sampCurrent.currentTime[i];
	for(i=0; i<sampScheme.sizeParamBirth; i++)
		x[ind++] = sampCurrent.currentBirth[i];
	for(i=0; i<sampScheme.sizeParamDeath; i++)
		x[ind++] = sampCurrent.currentDeath[i];
	for(i=0; i<sampScheme.sizeParamFossil; i++)
		x[ind++] = sampCurrent.currentFossil[i];
	for(i=0; i<sampScheme.sizeParamSampling; i++)
		x[ind++] = sampCurrent.currentSampling[i];
}

double toMaximizeSetTreeFossil(unsigned n, const double *x, double *grad, void *data) {
	double tmp = 0.;
////	old = ((TypeMinimizationSetTreeFossilData*)data)->likelihood(((TypeMinimizationSetTreeFossilData*)data)->tree, ((TypeMinimizationSetTreeFossilData*)data)->nTree, ((TypeMinimizationSetTreeFossilData*)data)->fos, &(((TypeMinimizationSetTreeFossilData*)data)->param));
	//setSamplingCurrentFromVector(&(((TypeMinimizationSetTreeFossilData*)data)->sampCurrent), x, ((TypeMinimizationSetTreeFossilData*)data)->sampScheme);
	//setParamFromSamplingCurrent(&(((TypeMinimizationSetTreeFossilData*)data)->param), ((TypeMinimizationSetTreeFossilData*)data)->sampCurrent, ((TypeMinimizationSetTreeFossilData*)data)->sampScheme);
	//tmp = ((TypeMinimizationSetTreeFossilData*)data)->likelihood(((TypeMinimizationSetTreeFossilData*)data)->tree, ((TypeMinimizationSetTreeFossilData*)data)->nTree, ((TypeMinimizationSetTreeFossilData*)data)->fos, &(((TypeMinimizationSetTreeFossilData*)data)->param));
////	printf("like %.4le -> %.4le\n", old, tmp);
	return tmp;
}

double toMaximizeSetTreeFossilInt(unsigned n, const double *x, double *grad, void *data) {
	double tmp, old;
	int i, m, f;
	setSamplingCurrentFossilFromVector(((TypeMinimizationSetTreeFossilIntData*)data)->tree, ((TypeMinimizationSetTreeFossilIntData*)data)->fint, ((TypeMinimizationSetTreeFossilIntData*)data)->fos, ((TypeMinimizationSetTreeFossilIntData*)data)->revert, ((TypeMinimizationSetTreeFossilIntData*)data)->sizeFossil, ((TypeMinimizationSetTreeFossilIntData*)data)->tree[0]->size, &(((TypeMinimizationSetTreeFossilIntData*)data)->sampCurrent), x, ((TypeMinimizationSetTreeFossilIntData*)data)->sampScheme);
	for(f=1; f<((TypeMinimizationSetTreeFossilIntData*)data)->fos->size; f++)
		if(((TypeMinimizationSetTreeFossilIntData*)data)->fos->fossilList[f].time == ((TypeMinimizationSetTreeFossilIntData*)data)->fos->fossilList[f-1].time)
			return -DBL_MAX;

//	setParamFromSamplingCurrent(&(((TypeMinimizationSetTreeFossilIntData*)data)->param), ((TypeMinimizationSetTreeFossilIntData*)data)->sampCurrent, ((TypeMinimizationSetTreeFossilIntData*)data)->sampScheme);
//printPiecewiseModel(stdout, &(((TypeMinimizationSetTreeFossilIntData*)data)->param));
////exit(0);
//printf("\n comp\n");
		//fprintFossilFeature(stdout, ((TypeMinimizationSetTreeFossilIntData*)data)->fos, ((TypeMinimizationSetTreeFossilIntData*)data)->tree[0]->name, ((TypeMinimizationSetTreeFossilIntData*)data)->tree[0]->size);
//printf("\n");
	for(i=0; i<((TypeMinimizationSetTreeFossilIntData*)data)->nTree; i++) {
		for(m=0; m<((TypeMinimizationSetTreeFossilIntData*)data)->tree[i]->size; m++) {
			if(((TypeMinimizationSetTreeFossilIntData*)data)->tree[i]->node[m].child == NOSUCH) {
				switch(((TypeMinimizationSetTreeFossilIntData*)data)->fint->status[m]) {
					case contempNodeStatus:
						((TypeMinimizationSetTreeFossilIntData*)data)->tree[i]->time[m] = ((TypeMinimizationSetTreeFossilIntData*)data)->tree[i]->maxTime;
					break;
					case unknownNodeStatus:
//						error("Node %d (%s) has unknown status\n", n, tree[i]->name[m]);
					break;
					case extinctNodeStatus:
						((TypeMinimizationSetTreeFossilIntData*)data)->tree[i]->time[m] = ((TypeMinimizationSetTreeFossilIntData*)data)->fos->fossilList[((TypeMinimizationSetTreeFossilIntData*)data)->fos->fossil[m]].time;
					break;
					default:
						error("Node %d has no status\n", n);
				}
			} else
				((TypeMinimizationSetTreeFossilIntData*)data)->tree[i]->time[m] = NO_TIME;
		}
	}
	for(i=0; i<=((TypeMinimizationSetTreeFossilIntData*)data)->param.size; i++) {
		if(((TypeMinimizationSetTreeFossilIntData*)data)->sampScheme.schemeTime[i] >= 0)
			((TypeMinimizationSetTreeFossilIntData*)data)->param.startTime[i] = ((TypeMinimizationSetTreeFossilIntData*)data)->sampCurrent.currentTime[((TypeMinimizationSetTreeFossilIntData*)data)->sampScheme.schemeTime[i]];
		else
			((TypeMinimizationSetTreeFossilIntData*)data)->param.startTime[i] = ((TypeMinimizationSetTreeFossilIntData*)data)->sampScheme.param.startTime[i];
	}

	for(i=0; i<((TypeMinimizationSetTreeFossilIntData*)data)->param.size; i++) {
		if(((TypeMinimizationSetTreeFossilIntData*)data)->sampScheme.schemeBirth[i] >= 0)
			((TypeMinimizationSetTreeFossilIntData*)data)->param.param[i].birth = ((TypeMinimizationSetTreeFossilIntData*)data)->sampCurrent.currentBirth[((TypeMinimizationSetTreeFossilIntData*)data)->sampScheme.schemeBirth[i]];
		else
			((TypeMinimizationSetTreeFossilIntData*)data)->param.param[i].birth = ((TypeMinimizationSetTreeFossilIntData*)data)->sampScheme.param.param[i].birth;
		if(((TypeMinimizationSetTreeFossilIntData*)data)->sampScheme.schemeDeath[i] >= 0)
			((TypeMinimizationSetTreeFossilIntData*)data)->param.param[i].death = ((TypeMinimizationSetTreeFossilIntData*)data)->sampCurrent.currentDeath[((TypeMinimizationSetTreeFossilIntData*)data)->sampScheme.schemeDeath[i]];
		else
			((TypeMinimizationSetTreeFossilIntData*)data)->param.param[i].death = ((TypeMinimizationSetTreeFossilIntData*)data)->sampScheme.param.param[i].death;
		if(((TypeMinimizationSetTreeFossilIntData*)data)->sampScheme.schemeFossil[i] >= 0)
			((TypeMinimizationSetTreeFossilIntData*)data)->param.param[i].fossil = ((TypeMinimizationSetTreeFossilIntData*)data)->sampCurrent.currentFossil[((TypeMinimizationSetTreeFossilIntData*)data)->sampScheme.schemeFossil[i]];
		else
			((TypeMinimizationSetTreeFossilIntData*)data)->param.param[i].fossil = ((TypeMinimizationSetTreeFossilIntData*)data)->sampScheme.param.param[i].fossil;
		if(((TypeMinimizationSetTreeFossilIntData*)data)->sampScheme.schemeSampling[i] >= 0)
			((TypeMinimizationSetTreeFossilIntData*)data)->param.param[i].sampling = ((TypeMinimizationSetTreeFossilIntData*)data)->sampCurrent.currentSampling[((TypeMinimizationSetTreeFossilIntData*)data)->sampScheme.schemeSampling[i]];
		else
			((TypeMinimizationSetTreeFossilIntData*)data)->param.param[i].sampling = ((TypeMinimizationSetTreeFossilIntData*)data)->sampScheme.param.param[i].sampling;
	}
	tmp = ((TypeMinimizationSetTreeFossilIntData*)data)->likelihood(((TypeMinimizationSetTreeFossilIntData*)data)->tree, ((TypeMinimizationSetTreeFossilIntData*)data)->nTree, ((TypeMinimizationSetTreeFossilIntData*)data)->fos, &(((TypeMinimizationSetTreeFossilIntData*)data)->param));
//printf("\nlike %le\n", tmp);
//printPiecewiseModel(stdout, &(((TypeMinimizationSetTreeFossilIntData*)data)->param));
	return tmp;
}
int fillFossilIndex(int *index, int *revert, TypeTree **tree, TypeFossilIntFeature *fint) {
	int n, tot=0;
	for(n=0; n<tree[0]->size; n++)
		if(fint->fossilInt[n]>=0) {
			int f;
			for(f=fint->fossilInt[n]; f>=0; f=fint->fossilIntList[f].prec) {
				index[tot] = f;
				revert[f] = tot;
				tot++;
			}
		}
	return tot;
}

int minimizePiecewiseParamFromSetTreeFossilInt(TypeLikelihoodSetTreeFosFunction *f, TypeTree **tree, int nTree, TypeFossilIntFeature *fint, TypeFossilFeature *foinit, TypeSamplingCurrent sampCurrentInit,  TypeSamplingScheme sampScheme, TypeNLOptOption *option, TypeEstimation *estim) {
	double *lb, *lu, *x, minLikelihood;
	nlopt_opt opt;
	TypeMinimizationSetTreeFossilIntData data;
	TypeFossilFeature *fos;
	int result, i, ind, nParam;
	//fprintNLoptOptionTag(stdout, option);
	//printf("\n");
	data.index = (int*) malloc(fint->sizeFossil*sizeof(int));
	data.revert = (int*) malloc(fint->sizeFossil*sizeof(int));
	data.sizeFossil = fillFossilIndex(data.index, data.revert, tree, fint);
	nParam = data.sizeFossil+sampScheme.sizeParamTime+sampScheme.sizeParamBirth+sampScheme.sizeParamDeath+sampScheme.sizeParamFossil+sampScheme.sizeParamSampling;
	lb = (double*) malloc(nParam*sizeof(double));
	lu = (double*) malloc(nParam*sizeof(double));
	x = (double*) malloc(nParam*sizeof(double));
	ind = 0;
	for(i=0; i<data.sizeFossil; i++) {
		lb[ind] = fint->fossilIntList[data.index[i]].fossilInt.inf;
		lu[ind] = fint->fossilIntList[data.index[i]].fossilInt.sup;
//printf("fos %d %.2lf %.2lf\n", i, lb[ind], lu[ind]);
		ind++;
	}	
	for(i=0; i<sampScheme.sizeParamTime; i++) {
		lb[ind] = sampScheme.param.startTime[sampScheme.indexTime[i]-1]+0.0000001;
		lu[ind] = sampScheme.param.startTime[sampScheme.indexTime[i]+1]-0.0000001;
//printf("time %d %.2lf %.2lf\n", i, lb[ind], lu[ind]);
		ind++;
	}
	for(i=0; i<sampScheme.sizeParamBirth; i++) {
		lb[ind] = 0.000001;
		lu[ind] = option->supSpe;
//printf("birth %d %.2lf %.2lf\n", i, lb[ind], lu[ind]);
		ind++;
	}
	for(i=0; i<sampScheme.sizeParamDeath; i++) {
		lb[ind] = 0.000001;
		lu[ind] = option->supExt;
//printf("death %d %.2lf %.2lf\n", i, lb[ind], lu[ind]);
		ind++;
	}
	for(i=0; i<sampScheme.sizeParamFossil; i++) {
		lb[ind] = 0.000001;
		lu[ind] = option->supFos;
//printf("foss %d %.2lf %.2lf\n", i, lb[ind], lu[ind]);
		ind++;
	}
	for(i=0; i<sampScheme.sizeParamSampling; i++) {
		lb[ind] = 0.000001;
		lu[ind] = 1.;
		ind++;
	}
	data.tree = tree;
	data.fint = fint;
	data.nTree = nTree;
	data.nParam = nParam;
//	data.fos = sampleFossilIntSpecial(fint, tree[0]->size);
	data.fos = foinit;
	data.param = getPiecewiseParamFromSamplingScheme(sampScheme);
	data.sampScheme = sampScheme;
//	data.sampCurrent = getSamplingCurrent(sampScheme);
	data.sampCurrent = getSamplingCurrent(sampScheme);
	data.likelihood = f;
	opt = nlopt_create(NLOPT_ALGO, nParam); /* algorithm and dimensionality */
	nlopt_set_lower_bounds(opt, lb);
	nlopt_set_upper_bounds(opt, lu);
	nlopt_set_max_objective(opt, toMaximizeSetTreeFossilInt, &data);
	nlopt_set_xtol_abs1(opt, option->tolOptim);
	nlopt_set_maxeval(opt, option->maxIter);
	setVectorFromSamplingCurrentFossil(x, data.fos, data.index, data.sizeFossil, sampCurrentInit, sampScheme);
//for(i=0; i<data.sizeFossil; i++)
	//printf("%.2lf\t%.2lf\t%.2lf\n", lb[i], x[i], lu[i]);

	estim->logLikelihood = -1.e20;
	setSamplingCurrentFromVector(&(data.sampCurrent), x, data.sizeFossil, data.sampScheme);
	//printf("nlopt %u arguments\n", nlopt_get_dimension(opt));
	//printf("nlopt %s\n%u arguments\n", nlopt_algorithm_name(nlopt_get_algorithm(opt)),nlopt_get_dimension(opt));
	//printf("nlopt %u arguments (%d)\n", nlopt_get_dimension(opt), nParam);
	if(((result = nlopt_optimize(opt, x, &minLikelihood)) >= 0)) {
		estim->logLikelihood = minLikelihood;
		setSamplingCurrentFossilFromVector(data.tree, data.fint, data.fos, data.revert, data.sizeFossil, data.tree[0]->size, &(data.sampCurrent), x, data.sampScheme);

//		setSamplingCurrentFromVector(&(data.sampCurrent), x, data.sizeFossil, data.sampScheme);
		setParamFromSamplingCurrent(&(estim->param), data.sampCurrent, data.sampScheme);
		//printf("result (%d) %.4le\n", result, minLikelihood);
	} else {
		switch(result) {
			case NLOPT_SUCCESS:
				printf("Success\n");
				break;
			case NLOPT_STOPVAL_REACHED:
				printf("Stop val reached\n");
				break;
			case NLOPT_FTOL_REACHED:
				printf("Ftol reached\n");
				break;
			case NLOPT_XTOL_REACHED:
				printf("Xtol reached\n");
				break;
			case NLOPT_MAXEVAL_REACHED:
				printf("Max eval reached\n");
				break;
			case NLOPT_MAXTIME_REACHED:
				printf("Max time reached\n");
				break;
			case NLOPT_FAILURE:
				printf("General failure\n");
				break;
			case NLOPT_INVALID_ARGS:
				printf("Invalid args\n");
				break;
			case NLOPT_OUT_OF_MEMORY:
				printf("Out of memory\n");
				break;
			case NLOPT_ROUNDOFF_LIMITED:
				printf("Roundoff limited\n");
				break;
			case NLOPT_FORCED_STOP:
				printf("Forced stop\n");
				break;
			default:
				printf("failure %d\n", result);
		}
		printf("result (%d)\n", result);
		printf("like %le\n", estim->logLikelihood);
	}
//	estim->logLikelihood = estim->logLikelihood;
	free((void*)lb);
	free((void*)lu);
	free((void*)x);
	freeSamplingCurrent(&(data.sampCurrent));
	freePiecewiseParam(&(data.param));
	nlopt_destroy(opt);
	return result;
}
