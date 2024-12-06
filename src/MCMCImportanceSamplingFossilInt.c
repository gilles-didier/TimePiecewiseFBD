#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <omp.h>

#include "Utils.h"
#include "MyR.h"
#include "MyRandom.h"
#include "FBDDensity.h"
#include "MCMCImportanceSamplingFossilInt.h"
int totProp = 0, rejProp = 0, totPropT = 0, rejPropT = 0, totPropF = 0, rejPropF = 0, totPropB = 0, rejPropB = 0, totPropD = 0, rejPropD = 0, totPropG = 0, rejPropG = 0, totPropS = 0, rejPropS = 0;

/*Entries of fossilList are assumed sorted in increasing order with regard to the times while the chain prec have to go in decreasing order again with time*/
typedef struct TYPE_VARIABLE_DATA {
	char type[5];
    int sizeVar[5], size;
    double prob[5];
} TypeTypeVariableData;

/*Entries of fossilList are assumed sorted in increasing order with regard to the times while the chain prec have to go in decreasing order again with time*/
typedef struct EXTENDED_FOSSIL_FEATURE {
	TypeFossilFeature *ff;
    int *first, *foll, **pre, **fol, *sizePre, *sizeFol;
    double *boundInf, *boundSup;
} TypeExtendedFossilFeature;

typedef struct MINIMIZATION_PARAM_DATA {
    TypeTree **tree;
    TypeFossilFeature **fos;
    int nTree;
} TypeMinimizationParamData;

static void fillBranchNode(TypeExtendedFossilFeature *eff, int size, int *branch);
static void getProposalFossil(TypeExtendedFossilFeature *eff, TypeFossilIntFeature *fi, int *branch, double al, int *move, double *newTime, void *rand_data);
static double getLogDensityOne(TypeTree *tree, TypeFossilFeature *fp, TypePiecewiseModelParam *param);
static TypeExtendedFossilFeature *extendFossilFeature(TypeFossilIntFeature *fos, TypeTree **tree, int size, void *rand_data);
static void freeExtendedFossilFeature(TypeExtendedFossilFeature *eff, int size);
static void getSchemeTypeAndVariable(char *type, int *var, TypeTypeVariableData tvd, void *rand_data);
static double updateSchemeSample(TypeTree **tree, int nTree, TypeExtendedFossilFeature *eff, TypeFossilIntFeature *fi, int *branch, double prob, double al, double propParam, TypePiecewiseModelParam *param, TypeModelParam *windSize, double probSpe, double probExt, double probFos, TypeTypeVariableData tvd, TypeSamplingScheme sampScheme, TypeSamplingCurrent sampCurrent, void *rand_data);

#define INC_SIZEBUF_SC 10
#define MAX_NAME_SIZE 100

/*return log(a+b) in an accurate way*/
double logSum(double a,double b) {
	double max, min;
	if(a>b) {
		max = a;
		min = b;
	} else {
		min = a;
		max = b;
	}
	return log(max) + log1p(min/max);
}

/*return log(exp(a)+exp(b)) in an accurate way*/
double logSumLog(double a, double b) {
	double max, min;
	if(a == NEG_INFTY)
		return b;
	if(b == NEG_INFTY)
		return a;
	if(a>b) {
		max = a;
		min = b;
	} else {
		min = a;
		max = b;
	}
	return max + log1p(exp(min-max));
}

void sampleFossilIntToExtendedFossilFeature(int n, int *sampled, TypeFossilIntFeature *fos, TypeExtendedFossilFeature *eff, void *rand_data) {
	int f, k, size, *tmpInd;
	size_t *index;
	double *tmpTime;

	if(sampled[n] || fos->fossilInt[n] == NOSUCH)
		return;
	for(k=0; k<eff->sizePre[n]; k++)
		if(!sampled[eff->pre[n][k]])
			sampleFossilIntToExtendedFossilFeature(eff->pre[n][k], sampled, fos, eff, rand_data);
	for(k=0; k<eff->sizePre[n]; k++)
		if(eff->ff->fossilList[eff->ff->fossil[eff->pre[n][k]]].time > eff->boundInf[n])
			eff->boundInf[n] = eff->ff->fossilList[eff->ff->fossil[eff->pre[n][k]]].time;
	size = 0;
	for(f=fos->fossilInt[n]; f!=NOSUCH; f=fos->fossilIntList[f].prec)
		size++;
	tmpInd = (int*) malloc(size*sizeof(int));
	tmpTime = (double*) malloc(size*sizeof(double));
	size = 0;

	for(f=fos->fossilInt[n]; f>=0; f=fos->fossilIntList[f].prec) {
		double bi, bs;
		tmpInd[size] = f;
		bi = utils_MAX(fos->fossilIntList[f].fossilInt.inf, eff->boundInf[n]);
		bs = utils_MIN(fos->fossilIntList[f].fossilInt.sup, eff->boundSup[n]);
		tmpTime[size] = bi+my_rand_unif_unit(rand_data)*(bs-bi);
		eff->ff->fossilList[f].time = tmpTime[size];
		size++;
	}
	index = getTableIndex(tmpTime, size, sizeof(double), compareDouble);
	eff->first[n] = tmpInd[index[0]];
	eff->ff->fossil[n] = tmpInd[index[size-1]];
 	eff->ff->fossilList[tmpInd[index[0]]].prec = NOSUCH;
  	eff->ff->fossilList[tmpInd[index[0]]].time = tmpTime[index[0]];
	eff->foll[tmpInd[index[size-1]]] = NOSUCH;
	for(f=1; f<size; f++) {
		eff->ff->fossilList[tmpInd[index[f]]].prec = tmpInd[index[f-1]];
		eff->ff->fossilList[tmpInd[index[f]]].time = tmpTime[index[f]];
		eff->foll[tmpInd[index[f-1]]] = tmpInd[index[f]];
	}
	for(k=0; k<eff->sizePre[n]; k++)
		if(eff->ff->fossilList[eff->first[n]].time < eff->boundSup[eff->pre[n][k]])
			eff->boundSup[eff->pre[n][k]] = eff->ff->fossilList[eff->first[n]].time;
	sampled[n] = 1;
    free((void*) tmpInd);
    free((void*) tmpTime);
    free((void*) index);
}

void fillPreFol(int n, int depth, int *ancestor, TypeTree *tree, double *min, double *max, double *minMax, TypeExtendedFossilFeature *eff) {
	int c, a;
	for(a=0; a<depth; a++)
		if(max[ancestor[a]] > min[n]) {
			int k;
			for(k=0; k<eff->sizePre[n] && eff->pre[n][k] != ancestor[a]; k++)
				;
			if(k == eff->sizePre[n])
				eff->pre[n][eff->sizePre[n]++] = ancestor[a];
			for(k=0; k<eff->sizeFol[ancestor[a]] && eff->fol[ancestor[a]][k] != n; k++)
				;
			if(k == eff->sizeFol[ancestor[a]]) {
				eff->fol[ancestor[a]][eff->sizeFol[ancestor[a]]++] = n;
				 if(eff->boundSup[ancestor[a]]>minMax[n])
					eff->boundSup[ancestor[a]] = minMax[n];
			}
		}
	ancestor[depth] = n;
	for(c=tree->node[n].child; c!=NOSUCH; c=tree->node[c].sibling)
		fillPreFol(c, depth+1, ancestor, tree, min, max, minMax, eff);
}

void freeExtendedFossilFeature(TypeExtendedFossilFeature *eff, int size) {
	int n;
	freeFossilFeature(eff->ff);
	free((void*) eff->foll);
	free((void*) eff->first);
	free((void*) eff->sizePre);
	free((void*) eff->sizeFol);
	free((void*) eff->boundInf);
	free((void*) eff->boundSup);
	for(n=0; n<size; n++) {
		if(eff->pre[n] != NULL)
			free((void*) eff->pre[n]);
		if(eff->fol[n] != NULL)
			free((void*) eff->fol[n]);
	}
	free((void*) eff->pre);
	free((void*) eff->fol);
	free((void*) eff);
}

TypeExtendedFossilFeature *extendFossilFeature(TypeFossilIntFeature *fos, TypeTree **tree, int size, void *rand_data) {
	TypeExtendedFossilFeature *eff;
	int i, n, *ancestor, *sampled;
	double *min, *max, *minMax;
	eff = (TypeExtendedFossilFeature*) malloc(sizeof(TypeExtendedFossilFeature));
	eff->ff = (TypeFossilFeature*) malloc(sizeof(TypeFossilFeature));
	eff->ff->size = fos->sizeFossil;
	eff->ff->sizeBuf = fos->sizeFossil;
	eff->ff->fossilList = (TypeFossilList*) malloc(eff->ff->size*sizeof(TypeFossilList));
	eff->ff->fossil = (int*) malloc(tree[0]->size*sizeof(int));
	eff->ff->status = (TypeNodeStatus*) malloc(tree[0]->size*sizeof(TypeNodeStatus));
	eff->foll = (int*) malloc(eff->ff->size*sizeof(int));
	eff->first = (int*) malloc(tree[0]->size*sizeof(int));
	eff->sizePre = (int*) malloc(tree[0]->size*sizeof(int));
	eff->sizeFol = (int*) malloc(tree[0]->size*sizeof(int));
	eff->pre = (int**) malloc(tree[0]->size*sizeof(int*));
	eff->fol = (int**) malloc(tree[0]->size*sizeof(int*));
	eff->boundInf = (double*) malloc(tree[0]->size*sizeof(double));
	eff->boundSup = (double*) malloc(tree[0]->size*sizeof(double));
	sampled = (int*) malloc(tree[0]->size*sizeof(int));
	ancestor = (int*) malloc(tree[0]->size*sizeof(int));
	min = (double*) malloc(tree[0]->size*sizeof(double));
	max = (double*) malloc(tree[0]->size*sizeof(double));
	minMax = (double*) malloc(tree[0]->size*sizeof(double));
	for(n=0; n<tree[0]->size; n++) {
		int k;
		min[n] = DBL_MAX;
		max[n] = -DBL_MAX;
		minMax[n] = DBL_MAX;
		eff->boundInf[n] = -DBL_MAX;
		eff->boundSup[n] = DBL_MAX;
		eff->ff->status[n] = fos->status[n];
		eff->ff->fossil[n] = NOSUCH;
		eff->first[n] = NOSUCH;
		for(k=fos->fossilInt[n]; k!=NOSUCH; k=fos->fossilIntList[k].prec) {
			if(fos->fossilIntList[k].fossilInt.inf<min[n])
				min[n] = fos->fossilIntList[k].fossilInt.inf;
			if(fos->fossilIntList[k].fossilInt.sup>max[n])
				max[n] = fos->fossilIntList[k].fossilInt.sup;
			if(fos->fossilIntList[k].fossilInt.sup<minMax[n])
				minMax[n] = fos->fossilIntList[k].fossilInt.sup;
		}
		eff->sizePre[n] = 0;
		eff->sizeFol[n] = 0;
		if(fos->fossilInt[n] != NOSUCH) {
			eff->pre[n] = (int*) malloc(tree[0]->size*sizeof(int));
			eff->fol[n] = (int*) malloc(tree[0]->size*sizeof(int));
		} else {
			eff->pre[n] = NULL;
			eff->fol[n] = NULL;
		}
	}
	for(i=0; i<size; i++)
		fillPreFol(tree[i]->root, 0, ancestor, tree[i], min, max, minMax, eff);
	for(n=0; n<tree[0]->size; n++) {
		if((fos->fossilInt[n] != NOSUCH) && (eff->sizePre[n] == 0)) {
			free((void*) eff->pre[n]);
			eff->pre[n] = NULL;
		} else
			eff->pre[n] = (int*) realloc((void*)eff->pre[n], eff->sizePre[n]*sizeof(int));
		if((fos->fossilInt[n] != NOSUCH) && (eff->sizeFol[n] == 0)) {
			free((void*) eff->fol[n]);
			eff->fol[n] = NULL;
		} else
			eff->fol[n] = (int*) realloc((void*)eff->fol[n], eff->sizeFol[n]*sizeof(int));
	}
	for(n=0; n<tree[0]->size; n++)
		sampled[n] = 0;
	
	for(n=0; n<tree[0]->size; n++)
		sampleFossilIntToExtendedFossilFeature(n, sampled, fos, eff, rand_data);
	free((void*) sampled);
	free((void*) ancestor);
	free((void*) min);
	free((void*) max);
	free((void*) minMax);
	return eff;
}

/*fol <- reciprocal of prec*/
void fillFolNode(TypeFossilFeature *f, int size, int *fol) {
	int n;
	for(n=0; n<size; n++) {
		int k, tmp;
		tmp = NOSUCH;
		for(k=f->fossil[n]; k!=NOSUCH; k=f->fossilList[k].prec) {
			fol[k] = tmp;
			tmp = k;
		}
	}
}

/*branch[k] <- branch bearing k*/
void fillBranchNode(TypeExtendedFossilFeature *eff, int size, int *branch) {
	int n;
	for(n=0; n<size; n++) {
		int k;
		for(k=eff->first[n]; k!=NOSUCH; k=eff->foll[k])
			branch[k] = n;
	}
}

void getProposalFossil(TypeExtendedFossilFeature *eff, TypeFossilIntFeature *fi, int *branch, double al, int *move, double *newTime, void *rand_data) {
	double inf, sup, length;
    *move = RANGE_RAND(eff->ff->size);
    inf = utils_MAX(fi->fossilIntList[*move].fossilInt.inf, eff->boundInf[branch[*move]]);
    sup = utils_MIN(fi->fossilIntList[*move].fossilInt.sup, eff->boundSup[branch[*move]]);
    length = (sup-inf)*al/2.;
 //   *newTime = my_rand_unif_unit(rand_data)*(2.*length)+eff->ff->fossilList[*move].time-length;
    *newTime = my_rand_unif_unit(rand_data)*((sup-inf))+inf;

    if(*newTime<inf)
		*newTime = 2.*inf - *newTime;
    if(*newTime>sup)
		*newTime = 2.*sup - *newTime;
}

void updateExtendedFossilFeature(TypeTree **tree, int nTree, TypeExtendedFossilFeature *eff, int *branch, int move, double newTime) {
	int xB, i, k;
	double oldInf, oldSup;
	oldInf = eff->ff->fossilList[eff->first[branch[move]]].time;
	oldSup = eff->ff->fossilList[eff->ff->fossil[branch[move]]].time;
	if(newTime > eff->ff->fossilList[move].time) {
		for(xB=eff->foll[move]; xB!=NOSUCH && newTime>eff->ff->fossilList[xB].time; xB=eff->foll[xB])
			;
		if(xB != eff->foll[move]) {
			if(eff->ff->fossilList[move].prec != NOSUCH)
				eff->foll[eff->ff->fossilList[move].prec] = eff->foll[move];
			else
				eff->first[branch[move]] = eff->foll[move];
			eff->ff->fossilList[eff->foll[move]].prec = eff->ff->fossilList[move].prec;
			eff->foll[move] = xB;
			if(xB != NOSUCH) {
				eff->ff->fossilList[move].prec = eff->ff->fossilList[xB].prec;
				eff->foll[eff->ff->fossilList[xB].prec] = move;
				eff->ff->fossilList[xB].prec = move;
			} else {
				eff->ff->fossilList[move].prec = eff->ff->fossil[branch[move]];
				eff->foll[eff->ff->fossil[branch[move]]] = move;
				eff->ff->fossil[branch[move]] = move;
			}
		}
	} else {
		for(xB=eff->ff->fossilList[move].prec; xB!=NOSUCH && newTime<eff->ff->fossilList[xB].time; xB=eff->ff->fossilList[xB].prec)
			;
		if(xB != eff->ff->fossilList[move].prec) {
			if(eff->foll[move] != NOSUCH)
				eff->ff->fossilList[eff->foll[move]].prec = eff->ff->fossilList[move].prec;
			else
				eff->ff->fossil[branch[move]] = eff->ff->fossilList[move].prec;
			eff->foll[eff->ff->fossilList[move].prec] = eff->foll[move];
			eff->ff->fossilList[move].prec = xB;
			if(xB != NOSUCH) {
				eff->ff->fossilList[eff->foll[xB]].prec = move;
				eff->foll[move] = eff->foll[xB];
				eff->foll[xB] = move;
			} else {
				eff->foll[move] = eff->first[branch[move]];
				eff->ff->fossilList[eff->first[branch[move]]].prec = move;
				eff->first[branch[move]] = move;
			}
		}
	}
	eff->ff->fossilList[move].time = newTime;
	if(oldInf != eff->ff->fossilList[eff->first[branch[move]]].time)
		for(k=0; k<eff->sizePre[branch[move]]; k++) {
			int l;
			eff->boundSup[eff->pre[branch[move]][k]] = DBL_MAX;
			for(l=0; l<eff->sizeFol[eff->pre[branch[move]][k]]; l++)
				if(eff->ff->fossilList[eff->first[eff->fol[eff->pre[branch[move]][k]][l]]].time < eff->boundSup[eff->pre[branch[move]][k]])
					eff->boundSup[eff->pre[branch[move]][k]] = eff->ff->fossilList[eff->first[eff->fol[eff->pre[branch[move]][k]][l]]].time;
		}
	if(oldSup != eff->ff->fossilList[eff->ff->fossil[branch[move]]].time) {
		for(k=0; k<eff->sizeFol[branch[move]]; k++) {
			int l;
			eff->boundInf[eff->fol[branch[move]][k]] = -DBL_MAX;
			for(l=0; l<eff->sizePre[eff->fol[branch[move]][k]]; l++)
				if(eff->ff->fossilList[eff->ff->fossil[eff->pre[eff->fol[branch[move]][k]][l]]].time > eff->boundInf[eff->fol[branch[move]][k]])
					eff->boundInf[eff->fol[branch[move]][k]] = eff->ff->fossilList[eff->ff->fossil[eff->pre[eff->fol[branch[move]][k]][l]]].time;
		}
		for(i=0; i<nTree; i++)
			if(tree[i]->node[branch[move]].child == NOSUCH && eff->ff->status[branch[move]] == extinctNodeStatus)
				tree[i]->time[branch[move]] = eff->ff->fossilList[eff->ff->fossil[branch[move]]].time;
	}
}

TypeTypeVariableData getTypeVariableData(double probTime, double probSpe, double probExt, double probFos, TypeSamplingScheme sampScheme) {
	TypeTypeVariableData tvd;
	double totProb = 0., probSampling = 1.-probTime-probSpe-probExt-probFos;
	int i;
	tvd.size = 0;
	if(sampScheme.sizeParamTime>0) {
		tvd.type[tvd.size] = 't';
		tvd.sizeVar[tvd.size] = sampScheme.sizeParamTime;
		totProb += probTime;
		tvd.prob[tvd.size] = totProb;
		tvd.size++;
	}
	if(sampScheme.sizeParamBirth>0) {
		tvd.type[tvd.size] = 'b';
		tvd.sizeVar[tvd.size] = sampScheme.sizeParamBirth;
		totProb += probSpe;
		tvd.prob[tvd.size] = totProb;
		tvd.size++;
	}
	if(sampScheme.sizeParamDeath>0) {
		tvd.type[tvd.size] = 'd';
		tvd.sizeVar[tvd.size] = sampScheme.sizeParamDeath;
		totProb += probExt;
		tvd.prob[tvd.size] = totProb;
		tvd.size++;
	}
	if(sampScheme.sizeParamFossil>0) {
		tvd.type[tvd.size] = 'f';
		tvd.sizeVar[tvd.size] = sampScheme.sizeParamFossil;
		totProb += probFos;
		tvd.prob[tvd.size] = totProb;
		tvd.size++;
	}
	if(sampScheme.sizeParamSampling>0) {
		tvd.type[tvd.size] = 's';
		tvd.sizeVar[tvd.size] = sampScheme.sizeParamSampling;
		totProb += probSampling;
		tvd.prob[tvd.size] = totProb;
		tvd.size++;
	}
	for(i=0; i<tvd.size; i++)
		tvd.prob[i] /= tvd.prob[tvd.size-1];
	return tvd;
}

void getSchemeTypeAndVariable(char *type, int *var, TypeTypeVariableData tvd, void *rand_data) {
	double p = my_rand_unif_unit(rand_data);
	int i;
	for(i=0; i<tvd.size && p>tvd.prob[i]; i++)
		;
	*type = tvd.type[i];
	*var = my_rand_unif_range(rand_data, tvd.sizeVar[i]-1);
}

double updateSchemeSample(TypeTree **tree, int nTree, TypeExtendedFossilFeature *eff, TypeFossilIntFeature *fi, int *branch, double prob, double al, double propParam, TypePiecewiseModelParam *param, TypeModelParam *windSize, double probSpe, double probExt, double probFos, TypeTypeVariableData tvd, TypeSamplingScheme sampScheme, TypeSamplingCurrent sampCurrent, void *rand_data) {
	int move;
	double newTime, oldTime, newProb, windTime = 30.;
totProp++;
	if(my_rand_unif_unit(rand_data) < propParam) {
		char type;
		int var, i;
		double old;
		getSchemeTypeAndVariable(&type, &var, tvd, rand_data);
		switch(type) {
			case 't':
totPropT++;
				old = sampCurrent.currentTime[var];
				sampCurrent.currentTime[var] = my_rand_unif_unit(rand_data)*windTime+sampCurrent.currentTime[var]-0.5*windTime;
				if(sampCurrent.currentTime[var] < param->startTime[sampScheme.indexTime[var]-1])
					sampCurrent.currentTime[var] = 2.*param->startTime[sampScheme.indexTime[var]-1]-sampCurrent.currentTime[var];
				if(sampCurrent.currentTime[var] > param->startTime[sampScheme.indexTime[var]+1])
					sampCurrent.currentTime[var] = 2.*param->startTime[sampScheme.indexTime[var]+1]-sampCurrent.currentTime[var];
				if(sampCurrent.currentTime[var] < param->startTime[sampScheme.indexTime[var]-1])
					sampCurrent.currentTime[var] = (param->startTime[sampScheme.indexTime[var]-1]+param->startTime[sampScheme.indexTime[var]+1])/2.;
				for(i=0; i<=param->size; i++)
					if(sampScheme.schemeTime[i] == var)
						param->startTime[i] = sampCurrent.currentTime[var];
			break;
			case 'b':
totPropB++;
				old = sampCurrent.currentBirth[var];
				sampCurrent.currentBirth[var] = my_rand_unif_unit(rand_data)*windSize->birth+sampCurrent.currentBirth[var]-0.5*windSize->birth;
				if(sampCurrent.currentBirth[var] < 0.)
					sampCurrent.currentBirth[var] = -sampCurrent.currentBirth[var];
				for(i=0; i<param->size; i++)
					if(sampScheme.schemeBirth[i] == var)
						param->param[i].birth = sampCurrent.currentBirth[var];
			break;
			case 'd':
totPropD++;
				old = sampCurrent.currentDeath[var];
				sampCurrent.currentDeath[var] = my_rand_unif_unit(rand_data)*windSize->death+sampCurrent.currentDeath[var]-0.5*windSize->death;
				if(sampCurrent.currentDeath[var] < 0.)
					sampCurrent.currentDeath[var] = -sampCurrent.currentDeath[var];
				for(i=0; i<param->size; i++)
					if(sampScheme.schemeDeath[i] == var)
						param->param[i].death = sampCurrent.currentDeath[var];
			break;
			case 'f':
totPropG++;
				old = sampCurrent.currentFossil[var];
				sampCurrent.currentFossil[var] = my_rand_unif_unit(rand_data)*windSize->fossil+sampCurrent.currentFossil[var]-0.5*windSize->fossil;
				if(sampCurrent.currentFossil[var] < 0.)
					sampCurrent.currentFossil[var] = -sampCurrent.currentFossil[var];
				for(i=0; i<param->size; i++)
					if(sampScheme.schemeFossil[i] == var)
						param->param[i].fossil = sampCurrent.currentFossil[var];
			break;
			case 's':
totPropS++;
				old = sampCurrent.currentSampling[var];
				sampCurrent.currentSampling[var] = my_rand_unif_unit(rand_data)*windSize->sampling+sampCurrent.currentSampling[var]-0.5*windSize->sampling;
				if(sampCurrent.currentSampling[var] < 0.)
					sampCurrent.currentSampling[var] = -sampCurrent.currentSampling[var];
				if(sampCurrent.currentSampling[var] > 1.)
					sampCurrent.currentSampling[var] = 2.-sampCurrent.currentSampling[var];
				for(i=0; i<param->size; i++)
					if(sampScheme.schemeSampling[i] == var)
						param->param[i].sampling = sampCurrent.currentSampling[var];
			break;
			default:
				error("Unknown type '%c' while scheme sampling\n", type);
		}
		if(nTree == 1)
			newProb = getLogDensityOne(*tree, eff->ff, param);
		else
			newProb = getLogDensitySum(tree, nTree, eff->ff, param);
		if(isnan(newProb))
			error("Param %d %.2le %.2le %.2le\n", i, param->param[i].birth, param->param[i].death, param->param[i].fossil);
//double ru = my_rand_unif_unit(rand_data);
//printf("%lf\t%lf %le (%le %le)\n", ru, exp(newProb-prob), newProb-prob, newProb, prob);
		//if(ru < exp(newProb-prob)) {
		if(my_rand_unif_unit(rand_data) < exp(newProb-prob)) {
			return newProb;
		} else {
rejProp++;
			switch(type) {
				case 't':
rejPropT++;
					sampCurrent.currentTime[var] = old;
					for(i=0; i<param->size; i++)
						if(sampScheme.schemeTime[i] == var)
							param->startTime[i] = sampCurrent.currentTime[var];
				break;
				case 'b':
rejPropB++;
					sampCurrent.currentBirth[var] = old;
					for(i=0; i<param->size; i++)
						if(sampScheme.schemeBirth[i] == var)
							param->param[i].birth = sampCurrent.currentBirth[var];
				break;
				case 'd':
rejPropD++;
					sampCurrent.currentDeath[var] = old;
					for(i=0; i<param->size; i++)
						if(sampScheme.schemeDeath[i] == var)
							param->param[i].death = sampCurrent.currentDeath[var];
				break;
				case 'f':
rejPropG++;
					sampCurrent.currentFossil[var] = old;
					for(i=0; i<param->size; i++)
						if(sampScheme.schemeFossil[i] == var)
							param->param[i].fossil = sampCurrent.currentFossil[var];
				break;
				case 's':
rejPropS++;
					sampCurrent.currentSampling[var] = old;
					for(i=0; i<param->size; i++)
						if(sampScheme.schemeSampling[i] == var)
							param->param[i].sampling = sampCurrent.currentSampling[var];
				break;
				default:
					error("Unknown type '%c' while scheme sampling\n", type);
			}
			return prob;
		}
	} else { 	
totPropF++;
		getProposalFossil(eff, fi, branch, al, &move, &newTime, rand_data);
		oldTime = eff->ff->fossilList[move].time;
		updateExtendedFossilFeature(tree, nTree, eff, branch, move, newTime);
		if(nTree == 1)
			newProb = getLogDensityOne(*tree, eff->ff, param);
		else
			newProb = getLogDensitySum(tree, nTree, eff->ff, param);
		if(isnan(newProb) || isinf(newProb))
			error("NaN issue (old prob %.2le) : fossil %d %.5lf -> %.5lf (%d %s) time %.5lf\n", prob, move, oldTime, newTime, branch[move], tree[0]->name[branch[move]], tree[0]->time[branch[move]]);
		if(my_rand_unif_unit(rand_data) < exp(newProb-prob)) {
			return newProb;
		} else {
rejProp++;
rejPropF++;
			updateExtendedFossilFeature(tree, nTree, eff, branch, move, oldTime);
			return prob;
		}
	}
	return prob;
}

double getLogDensityOne(TypeTree *tree, TypeFossilFeature *fp, TypePiecewiseModelParam *param) {
	return logLikelihoodTreeFossil(tree, fp, NULL, param);
}

double getLogDensitySum(TypeTree **tree, int nTree, TypeFossilFeature *ff, TypePiecewiseModelParam *param) {
	double sum = NEG_INFTY;
#pragma omp parallel shared(sum)
	{
		int i;
#pragma omp for
		for(i=0; i<nTree; i++) {
			double ll;
			ll = logLikelihoodTreeFossil(tree[i], ff, NULL, param);
#pragma omp critical
			{
				sum = logSumLog(sum, ll);
			}
		}
	}
	return sum;
//	return sum-log((double)nTree);
}

/*return the max likelihood observed during the MCMC*/
/*it is assumed that all the branch with identifiers in the trees have the same index, the tree nodes with fossils must have names*/
double MCMCGetMaxLogLikelihoodSC(TypeTree **tree, int nTree, TypeFossilIntFeature *fi, TypeSamplingScheme sampScheme, double al, int iter, double prop, TypeModelParam *windSize, TypeModelParam *init, double probSpe, double probExt, double probFos, void *rand_data) {
	int *branch, i, s, n;
	double prob, max = 0.;
	TypeExtendedFossilFeature *eff;
	TypePiecewiseModelParam param;
	TypeTypeVariableData tvd;
	TypeSamplingCurrent sampCurrent;
	//0. -> probTime
	tvd = getTypeVariableData(0., probSpe, probExt, probFos, sampScheme);
	
	param.size =  sampScheme.param.size;
	param.startTime =  sampScheme.param.startTime;
	param.param =  (TypeModelParam*) malloc(param.size*sizeof(TypeModelParam));
	if(sampScheme.sizeParamBirth>0) {
		int i;
		sampCurrent.currentBirth = (double*) malloc(sampScheme.sizeParamBirth*sizeof(double));
		for(i=0; i<sampScheme.sizeParamBirth; i++)
			sampCurrent.currentBirth[i] = my_rand_unif_unit(rand_data)*init->birth;
	} else
		sampCurrent.currentBirth = NULL;
	if(sampScheme.sizeParamDeath>0) {
		int i;
		sampCurrent.currentDeath = (double*) malloc(sampScheme.sizeParamDeath*sizeof(double));
		for(i=0; i<sampScheme.sizeParamDeath; i++)
			sampCurrent.currentDeath[i] = my_rand_unif_unit(rand_data)*init->death;
	} else
		sampCurrent.currentDeath = NULL;
	if(sampScheme.sizeParamFossil>0) {
		int i;
		sampCurrent.currentFossil = (double*) malloc(sampScheme.sizeParamFossil*sizeof(double));
		for(i=0; i<sampScheme.sizeParamFossil; i++)
			sampCurrent.currentFossil[i] = my_rand_unif_unit(rand_data)*init->fossil;
	} else
		sampCurrent.currentFossil = NULL;
	if(sampScheme.sizeParamSampling>0) {
		int i;
		sampCurrent.currentSampling = (double*) malloc(sampScheme.sizeParamSampling*sizeof(double));
		for(i=0; i<sampScheme.sizeParamSampling; i++)
			sampCurrent.currentSampling[i] = my_rand_unif_unit(rand_data);
	} else
		sampCurrent.currentSampling = NULL;
	for(i=0; i<param.size; i++) {
		if(sampScheme.schemeBirth[i] >= 0)
			param.param[i].birth = sampCurrent.currentBirth[sampScheme.schemeBirth[i]];
		else
			param.param[i].birth = sampScheme.param.param[i].birth;
		if(sampScheme.schemeDeath[i] >= 0)
			param.param[i].death = sampCurrent.currentDeath[sampScheme.schemeDeath[i]];
		else
			param.param[i].death = sampScheme.param.param[i].death;
		if(sampScheme.schemeFossil[i] >= 0)
			param.param[i].fossil = sampCurrent.currentFossil[sampScheme.schemeFossil[i]];
		else
			param.param[i].fossil = sampScheme.param.param[i].fossil;
		if(sampScheme.schemeSampling[i] >= 0)
			param.param[i].sampling = sampCurrent.currentSampling[sampScheme.schemeSampling[i]];
		else
			param.param[i].sampling = sampScheme.param.param[i].sampling;
	}
	eff = extendFossilFeature(fi, tree, nTree, rand_data);
	for(i=0; i<nTree; i++) {
		for(n=0; n<tree[i]->size; n++) {
			if(tree[i]->node[n].child == NOSUCH) {
				switch(fi->status[n]) {
					case contempNodeStatus:
						tree[i]->time[n] = tree[i]->maxTime;
					break;
					case unknownNodeStatus:
//						error("Node %d (%s) has unknown status\n", n, tree[i]->name[n]);
					break;
					case extinctNodeStatus:
						tree[i]->time[n] = eff->ff->fossilList[eff->ff->fossil[n]].time;
					break;
					default:
						error("Node %d has no status\n", n);
				}
			} else
				tree[i]->time[n] = NO_TIME;
		}
	}
	branch = (int*) malloc(eff->ff->size*sizeof(int));
	fillBranchNode(eff, tree[0]->size, branch);
	
	if(nTree == 1)
		prob = getLogDensityOne(*tree, eff->ff, &param);
	else
		prob = getLogDensitySum(tree, nTree, eff->ff, &param);
	max = prob;
	for(s=0; s<iter; s++) {
		prob = updateSchemeSample(tree, nTree, eff, fi, branch, prob, al, prop, &param, windSize, probSpe, probExt, probFos, tvd, sampScheme, sampCurrent, rand_data);
		if(prob>max)
			max = prob;
//if(s%10 == 0) {
	//fprintf(stderr, "it %-7d/%d %le %le\r", s, iter, max, prob); fflush(stderr);
//}
	}
	if(sampCurrent.currentBirth != NULL)
		free((void*) sampCurrent.currentBirth);
	if(sampCurrent.currentDeath != NULL)
		free((void*) sampCurrent.currentDeath);
	if(sampCurrent.currentFossil != NULL)
		free((void*) sampCurrent.currentFossil);
	if(sampCurrent.currentSampling != NULL)
		free((void*) sampCurrent.currentSampling);
	free((void*) param.param);
	free((void*)branch);
	freeExtendedFossilFeature(eff, tree[0]->size);
	return max;
}


/*return the max likelihood observed during the MCMC*/
/*it is assumed that all the branch with identifiers in the trees have the same index, the tree nodes with fossils must have names*/
double MCMCGetMaxLogLikelihoodSCBis(TypeTree **tree, int nTree, TypeFossilIntFeature *fi, TypeSamplingScheme sampScheme, double al, int trial, int iter, double prop, TypeModelParam *windSize, TypeModelParam *init, double probSpe, double probExt, double probFos, void *rand_data) {
	int *branch, i, s, n, t;
	double prob, max = 0.;
	TypeExtendedFossilFeature *eff;
	TypePiecewiseModelParam param;
	TypeTypeVariableData tvd;
	TypeSamplingCurrent sampCurrent;
	
	//0. -> probTime
	tvd = getTypeVariableData(0., probSpe, probExt, probFos, sampScheme);
	
	param.size =  sampScheme.param.size;
	param.startTime =  sampScheme.param.startTime;
	param.param =  (TypeModelParam*) malloc(param.size*sizeof(TypeModelParam));
	if(sampScheme.sizeParamBirth>0)
		sampCurrent.currentBirth = (double*) malloc(sampScheme.sizeParamBirth*sizeof(double));
	else
		sampCurrent.currentBirth = NULL;
	if(sampScheme.sizeParamDeath>0)
		sampCurrent.currentDeath = (double*) malloc(sampScheme.sizeParamDeath*sizeof(double));
	else
		sampCurrent.currentDeath = NULL;
	if(sampScheme.sizeParamFossil>0)
		sampCurrent.currentFossil = (double*) malloc(sampScheme.sizeParamFossil*sizeof(double));
	else
		sampCurrent.currentFossil = NULL;
	if(sampScheme.sizeParamSampling>0)
		sampCurrent.currentSampling = (double*) malloc(sampScheme.sizeParamSampling*sizeof(double));
	else
		sampCurrent.currentSampling = NULL;
	for(i=0; i<param.size; i++) {
		if(sampScheme.schemeBirth[i] >= 0)
			param.param[i].birth = sampCurrent.currentBirth[sampScheme.schemeBirth[i]];
		else
			param.param[i].birth = sampScheme.param.param[i].birth;
		if(sampScheme.schemeDeath[i] >= 0)
			param.param[i].death = sampCurrent.currentDeath[sampScheme.schemeDeath[i]];
		else
			param.param[i].death = sampScheme.param.param[i].death;
		if(sampScheme.schemeFossil[i] >= 0)
			param.param[i].fossil = sampCurrent.currentFossil[sampScheme.schemeFossil[i]];
		else
			param.param[i].fossil = sampScheme.param.param[i].fossil;
		if(sampScheme.schemeSampling[i] >= 0)
			param.param[i].sampling = sampCurrent.currentSampling[sampScheme.schemeSampling[i]];
		else
			param.param[i].sampling = sampScheme.param.param[i].sampling;
	}
	eff = extendFossilFeature(fi, tree, nTree, rand_data);
	for(i=0; i<nTree; i++) {
		for(n=0; n<tree[i]->size; n++) {
			if(tree[i]->node[n].child == NOSUCH) {
				switch(fi->status[n]) {
					case contempNodeStatus:
						tree[i]->time[n] = tree[i]->maxTime;
					break;
					case unknownNodeStatus:
//						error("Node %d (%s) has unknown status\n", n, tree[i]->name[n]);
					break;
					case extinctNodeStatus:
						tree[i]->time[n] = eff->ff->fossilList[eff->ff->fossil[n]].time;
					break;
					default:
						error("Node %d has no status\n", n);
				}
			} else
				tree[i]->time[n] = NO_TIME;
		}
	}
	branch = (int*) malloc(eff->ff->size*sizeof(int));
	fillBranchNode(eff, tree[0]->size, branch);
	
	for(t=0; t<trial; t++) {
		int i;
		if(sampScheme.sizeParamBirth>0)
			for(i=0; i<sampScheme.sizeParamBirth; i++)
				sampCurrent.currentBirth[i] = my_rand_unif_unit(rand_data)*init->birth;
		else
			sampCurrent.currentBirth = NULL;
		if(sampScheme.sizeParamDeath>0)
			for(i=0; i<sampScheme.sizeParamDeath; i++)
				sampCurrent.currentDeath[i] = my_rand_unif_unit(rand_data)*init->death;
		else
			sampCurrent.currentDeath = NULL;
		if(sampScheme.sizeParamFossil>0)
			for(i=0; i<sampScheme.sizeParamFossil; i++)
				sampCurrent.currentFossil[i] = my_rand_unif_unit(rand_data)*init->fossil;
		else
			sampCurrent.currentFossil = NULL;
		if(sampScheme.sizeParamSampling>0)
			for(i=0; i<sampScheme.sizeParamSampling; i++)
				sampCurrent.currentSampling[i] = my_rand_unif_unit(rand_data);
		else
			sampCurrent.currentSampling = NULL;
		for(i=0; i<param.size; i++) {
			if(sampScheme.schemeBirth[i] >= 0)
				param.param[i].birth = sampCurrent.currentBirth[sampScheme.schemeBirth[i]];
			else
				param.param[i].birth = sampScheme.param.param[i].birth;
			if(sampScheme.schemeDeath[i] >= 0)
				param.param[i].death = sampCurrent.currentDeath[sampScheme.schemeDeath[i]];
			else
				param.param[i].death = sampScheme.param.param[i].death;
			if(sampScheme.schemeFossil[i] >= 0)
				param.param[i].fossil = sampCurrent.currentFossil[sampScheme.schemeFossil[i]];
			else
				param.param[i].fossil = sampScheme.param.param[i].fossil;
			if(sampScheme.schemeSampling[i] >= 0)
				param.param[i].sampling = sampCurrent.currentSampling[sampScheme.schemeSampling[i]];
			else
				param.param[i].sampling = sampScheme.param.param[i].sampling;
		}
		if(nTree == 1)
			prob = getLogDensityOne(*tree, eff->ff, &param);
		else
			prob = getLogDensitySum(tree, nTree, eff->ff, &param);
		max = prob;
		for(s=0; s<iter; s++) {
			prob = updateSchemeSample(tree, nTree, eff, fi, branch, prob, al, prop, &param, windSize, probSpe, probExt, probFos, tvd, sampScheme, sampCurrent, rand_data);
			if(prob>max)
				max = prob;
	//if(s%10 == 0) {
		//fprintf(stderr, "it %-7d/%d %le %le\r", s, iter, max, prob); fflush(stderr);
	//}
		}
	}
	if(sampCurrent.currentBirth != NULL)
		free((void*) sampCurrent.currentBirth);
	if(sampCurrent.currentDeath != NULL)
		free((void*) sampCurrent.currentDeath);
	if(sampCurrent.currentFossil != NULL)
		free((void*) sampCurrent.currentFossil);
	if(sampCurrent.currentSampling != NULL)
		free((void*) sampCurrent.currentSampling);
	free((void*) param.param);
	free((void*)branch);
	freeExtendedFossilFeature(eff, tree[0]->size);
	return max;
}


/*return the max likelihood observed during the MCMC and put the correstponding parameters in maxParam*/
/*it is assumed that all the branch with identifiers in the trees have the same index, the tree nodes with fossils must have names*/
double MCMCFillMaxParamSampleSC(TypeTree **tree, int nTree, TypeFossilIntFeature *fi, TypeSamplingScheme sampScheme, double al, int iter, double prop, TypeModelParam *windSize, TypeModelParam *init,double probTime,  double probSpe, double probExt, double probFos, TypePiecewiseModelParam *maxParam, TypeFossilFeature **foss, TypeSamplingCurrent *saCu, void *rand_data) {
	int *branch, i, s, n;
	double prob, max = 0.;
	TypeExtendedFossilFeature *eff;
	TypePiecewiseModelParam param;
	TypeTypeVariableData tvd;
	TypeSamplingCurrent sampCurrent;
	TypeFossilFeature *fsave;

	tvd = getTypeVariableData(probTime, probSpe, probExt, probFos, sampScheme);
	param.size =  sampScheme.param.size;
	param.startTime =  sampScheme.param.startTime;
	param.param =  (TypeModelParam*) malloc(param.size*sizeof(TypeModelParam));
	if(sampScheme.sizeParamTime>0)
		saCu->currentTime = (double*) malloc(sampScheme.sizeParamTime*sizeof(double));
	else
		saCu->currentTime = NULL;
	if(sampScheme.sizeParamBirth>0)
		saCu->currentBirth = (double*) malloc(sampScheme.sizeParamBirth*sizeof(double));
	else
		saCu->currentBirth = NULL;
	if(sampScheme.sizeParamDeath>0)
		saCu->currentDeath = (double*) malloc(sampScheme.sizeParamDeath*sizeof(double));
	else
		saCu->currentDeath = NULL;
	if(sampScheme.sizeParamFossil>0)
		saCu->currentFossil = (double*) malloc(sampScheme.sizeParamFossil*sizeof(double));
	else
		saCu->currentFossil = NULL;
	if(sampScheme.sizeParamSampling>0)
		saCu->currentSampling = (double*) malloc(sampScheme.sizeParamSampling*sizeof(double));
	else
		saCu->currentSampling = NULL;
	if(sampScheme.sizeParamTime>0) {
		int i;
		sampCurrent.currentTime = (double*) malloc(sampScheme.sizeParamTime*sizeof(double));
		for(i=0; i<sampScheme.sizeParamTime; i++) {
			sampCurrent.currentTime[i] = param.startTime[sampScheme.indexTime[i]-1]+(param.startTime[sampScheme.indexTime[i]+1]-param.startTime[sampScheme.indexTime[i]-1])*my_rand_unif_unit(rand_data);
			sampCurrent.currentTime[i] = param.startTime[sampScheme.indexTime[i]-1]+(param.startTime[sampScheme.indexTime[i]+1]-param.startTime[sampScheme.indexTime[i]-1])*my_rand_unif_unit(rand_data);
//			sampCurrent.currentTime[i] = -262;
		}
	} else
		sampCurrent.currentTime = NULL;
	if(sampScheme.sizeParamBirth>0) {
		int i;
		sampCurrent.currentBirth = (double*) malloc(sampScheme.sizeParamBirth*sizeof(double));
		for(i=0; i<sampScheme.sizeParamBirth; i++)
			sampCurrent.currentBirth[i] = my_rand_unif_unit(rand_data)*init->birth;
	} else
		sampCurrent.currentBirth = NULL;
	if(sampScheme.sizeParamDeath>0) {
		int i;
		sampCurrent.currentDeath = (double*) malloc(sampScheme.sizeParamDeath*sizeof(double));
		for(i=0; i<sampScheme.sizeParamDeath; i++)
			sampCurrent.currentDeath[i] = my_rand_unif_unit(rand_data)*init->death;
	} else
		sampCurrent.currentDeath = NULL;
	if(sampScheme.sizeParamFossil>0) {
		int i;
		sampCurrent.currentFossil = (double*) malloc(sampScheme.sizeParamFossil*sizeof(double));
		for(i=0; i<sampScheme.sizeParamFossil; i++)
			sampCurrent.currentFossil[i] = my_rand_unif_unit(rand_data)*init->fossil;
	} else
		sampCurrent.currentFossil = NULL;
	if(sampScheme.sizeParamSampling>0) {
		int i;
		sampCurrent.currentSampling = (double*) malloc(sampScheme.sizeParamSampling*sizeof(double));
		for(i=0; i<sampScheme.sizeParamSampling; i++)
			sampCurrent.currentSampling[i] = my_rand_unif_unit(rand_data);
	} else
		sampCurrent.currentSampling = NULL;
	for(i=0; i<=param.size; i++) {
		if(sampScheme.schemeTime[i] >= 0)
			param.startTime[i] = sampCurrent.currentTime[sampScheme.schemeTime[i]];
		else
			param.startTime[i] = sampScheme.param.startTime[i];
	}
	for(i=0; i<param.size; i++) {
		if(sampScheme.schemeBirth[i] >= 0)
			param.param[i].birth = sampCurrent.currentBirth[sampScheme.schemeBirth[i]];
		else
			param.param[i].birth = sampScheme.param.param[i].birth;
		if(sampScheme.schemeDeath[i] >= 0)
			param.param[i].death = sampCurrent.currentDeath[sampScheme.schemeDeath[i]];
		else
			param.param[i].death = sampScheme.param.param[i].death;
		if(sampScheme.schemeFossil[i] >= 0)
			param.param[i].fossil = sampCurrent.currentFossil[sampScheme.schemeFossil[i]];
		else
			param.param[i].fossil = sampScheme.param.param[i].fossil;
		if(sampScheme.schemeSampling[i] >= 0)
			param.param[i].sampling = sampCurrent.currentSampling[sampScheme.schemeSampling[i]];
		else
			param.param[i].sampling = sampScheme.param.param[i].sampling;
	}
	eff = extendFossilFeature(fi, tree, nTree, rand_data);
	for(i=0; i<nTree; i++) {
		for(n=0; n<tree[i]->size; n++) {
			if(tree[i]->node[n].child == NOSUCH) {
				switch(fi->status[n]) {
					case contempNodeStatus:
						tree[i]->time[n] = tree[i]->maxTime;
					break;
					case unknownNodeStatus:
//						error("Node %d (%s) has unknown status\n", n, tree[i]->name[n]);
					break;
					case extinctNodeStatus:
						tree[i]->time[n] = eff->ff->fossilList[eff->ff->fossil[n]].time;
					break;
					default:
						error("Node %d has no status\n", n);
				}
			} else
				tree[i]->time[n] = NO_TIME;
		}
	}
	branch = (int*) malloc(eff->ff->size*sizeof(int));
	fillBranchNode(eff, tree[0]->size, branch);
	if(nTree == 1)
		prob = getLogDensityOne(*tree, eff->ff, &param);
	else
		prob = getLogDensitySum(tree, nTree, eff->ff, &param);
	max = prob;
	fsave = cpyFossilFeature(eff->ff, tree[0]->size);
	maxParam->size = param.size;
	for(i=0; i<param.size; i++)
		maxParam->param[i] = param.param[i];
	for(s=0; s<iter; s++) {
		prob = updateSchemeSample(tree, nTree, eff, fi, branch, prob, al, prop, &param, windSize, probSpe, probExt, probFos, tvd, sampScheme, sampCurrent, rand_data);
		if(prob>max) {
			int k;
			max = prob;
			for(k=0; k<param.size; k++)
				maxParam->param[k] = param.param[k];
			for(k=0; k<=param.size; k++)
				maxParam->startTime[k] = param.startTime[k];
			if(fsave != NULL)
				freeFossilFeature(fsave);
			fsave = cpyFossilFeature(eff->ff, tree[0]->size);
			if(sampScheme.sizeParamTime>0) {
				int i;
				for(i=0; i<sampScheme.sizeParamTime; i++)
					saCu->currentTime[i] = sampCurrent.currentTime[i];
			}
			if(sampScheme.sizeParamBirth>0) {
				int i;
				for(i=0; i<sampScheme.sizeParamBirth; i++)
					saCu->currentBirth[i] = sampCurrent.currentBirth[i];
			}
			if(sampScheme.sizeParamDeath>0) {
				int i;
				for(i=0; i<sampScheme.sizeParamDeath; i++)
					saCu->currentDeath[i] = sampCurrent.currentDeath[i];
			}
			if(sampScheme.sizeParamFossil>0) {
				int i;
				for(i=0; i<sampScheme.sizeParamFossil; i++)
					saCu->currentFossil[i] = sampCurrent.currentFossil[i];
			}
			if(sampScheme.sizeParamSampling>0) {
				int i;
				for(i=0; i<sampScheme.sizeParamSampling; i++)
					saCu->currentSampling[i] = sampCurrent.currentSampling[i];
			}		
		}
if(s%10 == 0) {
	fprintf(stderr, "it %-7d/%d %le / %le (%.2lf/%.2lf)\r", s, iter, max, prob, maxParam->startTime[1], param.startTime[1]); fflush(stderr);
}
	}
	FILE *fo;
	if((fo = fopen("FossiML.csv", "w"))) {
		fprintFossilFeature(fo, fsave, tree[0]->name, tree[0]->size);
		fclose(fo);
	} else
		error("Can't open %s\n", "FossiML.csv");
	if(sampCurrent.currentBirth != NULL)
		free((void*) sampCurrent.currentBirth);
	if(sampCurrent.currentDeath != NULL)
		free((void*) sampCurrent.currentDeath);
	if(sampCurrent.currentFossil != NULL)
		free((void*) sampCurrent.currentFossil);
	if(sampCurrent.currentSampling != NULL)
		free((void*) sampCurrent.currentSampling);
	free((void*) param.param);
	free((void*)branch);
	freeExtendedFossilFeature(eff, tree[0]->size);
//	freeFossilFeature(fsave);
	*foss = fsave;
	return max;
}

/*compute the posterior distribution of the free parameters of the sampScheme (not that of the fossils*/
/*it is assumed that all the branch with identifiers in the trees have the same index, the tree nodes with fossils must have names*/
void MCMCSamplingSamplePosteriorParametersOnly(FILE *fout, FILE *find, TypeTree **tree, int nTree, TypeFossilIntFeature *fi, TypeSamplingScheme sampScheme, double al, int burn, int gap, int iter, double prop, TypeModelParam *windSize, TypeModelParam *init, double probTime, double probSpe, double probExt, double probFos, void *rand_data) {
	int *branch, i, s, n;
	double prob;
	TypeExtendedFossilFeature *eff;
	TypePiecewiseModelParam param;

	TypeTypeVariableData tvd;
	TypeSamplingCurrent sampCurrent, *sampCurrentSave;
	tvd = getTypeVariableData(probTime, probSpe, probExt, probFos, sampScheme);
	sampCurrentSave = (TypeSamplingCurrent*) malloc(iter*sizeof(TypeSamplingCurrent));
	param.size =  sampScheme.param.size;
	param.startTime =  sampScheme.param.startTime;
	param.param =  (TypeModelParam*) malloc(param.size*sizeof(TypeModelParam));
	if(sampScheme.sizeParamTime>0) {
		int i;
		sampCurrent.currentTime = (double*) malloc(sampScheme.sizeParamTime*sizeof(double));
		for(i=0; i<sampScheme.sizeParamTime; i++)
			sampCurrent.currentTime[i] = param.startTime[sampScheme.indexTime[i]-1]+(param.startTime[sampScheme.indexTime[i]+1]-param.startTime[sampScheme.indexTime[i]-1])*my_rand_unif_unit(rand_data);
			sampCurrent.currentTime[i] = -300;
		for(i=0; i<iter; i++)
			sampCurrentSave[i].currentTime = (double*) malloc(sampScheme.sizeParamTime*sizeof(double));
	} else {
		sampCurrent.currentTime = NULL;
		for(i=0; i<iter; i++)
			sampCurrentSave[i].currentTime = NULL;
	}
	if(sampScheme.sizeParamBirth>0) {
		int i;
		sampCurrent.currentBirth = (double*) malloc(sampScheme.sizeParamBirth*sizeof(double));
		for(i=0; i<sampScheme.sizeParamBirth; i++)
			sampCurrent.currentBirth[i] = my_rand_unif_unit(rand_data)*init->birth;
		for(i=0; i<iter; i++)
			sampCurrentSave[i].currentBirth = (double*) malloc(sampScheme.sizeParamBirth*sizeof(double));
	} else {
		sampCurrent.currentBirth = NULL;
		for(i=0; i<iter; i++)
			sampCurrentSave[i].currentBirth = NULL;
	}
	if(sampScheme.sizeParamDeath>0) {
		int i;
		sampCurrent.currentDeath = (double*) malloc(sampScheme.sizeParamDeath*sizeof(double));
		for(i=0; i<sampScheme.sizeParamDeath; i++)
			sampCurrent.currentDeath[i] = my_rand_unif_unit(rand_data)*init->death;
		for(i=0; i<iter; i++)
			sampCurrentSave[i].currentDeath = (double*) malloc(sampScheme.sizeParamDeath*sizeof(double));
	} else {
		sampCurrent.currentDeath = NULL;
		for(i=0; i<iter; i++)
			sampCurrentSave[i].currentDeath = NULL;
	}
	if(sampScheme.sizeParamFossil>0) {
		int i;
		sampCurrent.currentFossil = (double*) malloc(sampScheme.sizeParamFossil*sizeof(double));
		for(i=0; i<sampScheme.sizeParamFossil; i++)
			sampCurrent.currentFossil[i] = my_rand_unif_unit(rand_data)*init->fossil;
		for(i=0; i<iter; i++)
			sampCurrentSave[i].currentFossil = (double*) malloc(sampScheme.sizeParamFossil*sizeof(double));
	} else {
		sampCurrent.currentFossil = NULL;
		for(i=0; i<iter; i++)
			sampCurrentSave[i].currentFossil = NULL;
	}
	if(sampScheme.sizeParamSampling>0) {
		int i;
		sampCurrent.currentSampling = (double*) malloc(sampScheme.sizeParamSampling*sizeof(double));
		for(i=0; i<sampScheme.sizeParamSampling; i++)
			sampCurrent.currentSampling[i] = my_rand_unif_unit(rand_data);
		for(i=0; i<iter; i++)
			sampCurrentSave[i].currentSampling = (double*) malloc(sampScheme.sizeParamSampling*sizeof(double));
	} else {
		sampCurrent.currentSampling = NULL;
		for(i=0; i<iter; i++)
			sampCurrentSave[i].currentSampling = NULL;
	}
	for(i=0; i<=param.size; i++) {
		if(sampScheme.schemeTime[i] >= 0)
			param.startTime[i] = sampCurrent.currentTime[sampScheme.schemeTime[i]];
		else
			param.startTime[i] = sampScheme.param.startTime[i];
	}
	for(i=0; i<param.size; i++) {
		if(sampScheme.schemeBirth[i] >= 0)
			param.param[i].birth = sampCurrent.currentBirth[sampScheme.schemeBirth[i]];
		else
			param.param[i].birth = sampScheme.param.param[i].birth;
		if(sampScheme.schemeDeath[i] >= 0)
			param.param[i].death = sampCurrent.currentDeath[sampScheme.schemeDeath[i]];
		else
			param.param[i].death = sampScheme.param.param[i].death;
		if(sampScheme.schemeFossil[i] >= 0)
			param.param[i].fossil = sampCurrent.currentFossil[sampScheme.schemeFossil[i]];
		else
			param.param[i].fossil = sampScheme.param.param[i].fossil;
		if(sampScheme.schemeSampling[i] >= 0)
			param.param[i].sampling = sampCurrent.currentSampling[sampScheme.schemeSampling[i]];
		else
			param.param[i].sampling = sampScheme.param.param[i].sampling;
	}
	eff = extendFossilFeature(fi, tree, nTree, rand_data);
	for(i=0; i<nTree; i++) {
		for(n=0; n<tree[i]->size; n++) {
			if(tree[i]->node[n].child == NOSUCH) {
				switch(fi->status[n]) {
					case contempNodeStatus:
						tree[i]->time[n] = tree[i]->maxTime;
					break;
					case unknownNodeStatus:
					break;
					case extinctNodeStatus:
						tree[i]->time[n] = eff->ff->fossilList[eff->ff->fossil[n]].time;
					break;
					default:
						error("Node %d has no status\n", n);
				}
			} else
				tree[i]->time[n] = NO_TIME;
		}
	}
	branch = (int*) malloc(eff->ff->size*sizeof(int));
	fillBranchNode(eff, tree[0]->size, branch);

	if(nTree == 1)
		prob = getLogDensityOne(*tree, eff->ff, &param);
	else
		prob = getLogDensitySum(tree, nTree, eff->ff, &param);
	for(i=0; i<burn; i++) {
		prob = updateSchemeSample(tree, nTree, eff, fi, branch, prob, al, prop, &param, windSize, probSpe, probExt, probFos, tvd, sampScheme, sampCurrent, rand_data);
if(i%10 == 0) {
//	fprintf(stderr, "it %-7d/%d %le (%.2lf)\ttot %d\trej %d\trt %le\r", i, burn, prob, param.startTime[1], totPropT, rejPropT, ((double)rejPropT)/((double)totPropT)); fflush(stderr);
//	fprintf(stderr, "it %-7d/%d %le (%.2lf)\r", i, burn, prob, param.startTime[1]); fflush(stderr);
}
	}
	FILE *fo;
	if((fo = fopen("log_like.csv", "w"))) {
	for(s=0; s<iter; s++) {
		int j;
		fprintf(fo, "%le\n", prob);
		for(i=0; i<sampScheme.sizeParamTime; i++)
			sampCurrentSave[s].currentTime[i] = sampCurrent.currentTime[i];
		for(i=0; i<sampScheme.sizeParamBirth; i++)
			sampCurrentSave[s].currentBirth[i] = sampCurrent.currentBirth[i];
		for(i=0; i<sampScheme.sizeParamDeath; i++)
			sampCurrentSave[s].currentDeath[i] = sampCurrent.currentDeath[i];
		for(i=0; i<sampScheme.sizeParamFossil; i++)
			sampCurrentSave[s].currentFossil[i] = sampCurrent.currentFossil[i];
		for(i=0; i<sampScheme.sizeParamSampling; i++)
			sampCurrentSave[s].currentSampling[i] = sampCurrent.currentSampling[i];
		for(j=0; j<gap; j++) {
			prob = updateSchemeSample(tree, nTree, eff, fi, branch, prob, al, prop, &param, windSize, probSpe, probExt, probFos, tvd, sampScheme, sampCurrent, rand_data);
		}
if(s%10 == 0) {
//	fprintf(stderr, "it %-7d/%d %le (%.2lf)\ttot %d\trej %d\trt %le\r", s, iter, prob, param.startTime[1], totProp, rejProp, ((double)rejProp)/((double)totProp)); fflush(stderr);
	fprintf(stderr, "it %-7d/%d %le (%.2lf)\r", s, iter, prob, param.startTime[1]); fflush(stderr);
}
	}
	fclose(fo);
	}
printf("\n\ntot %d\trej %d\trt %le\n", totProp, rejProp, ((double)rejProp)/((double)totProp));
printf("\n\nTTtot %d\trej %d\trt %le\n", totPropT, rejPropT, ((double)rejPropT)/((double)totPropT));
printf("\n\nBBtot %d\trej %d\trt %le\n", totPropB, rejPropB, ((double)rejPropB)/((double)totPropB));
printf("\n\nDDtot %d\trej %d\trt %le\n", totPropD, rejPropD, ((double)rejPropD)/((double)totPropD));
printf("\n\nGGtot %d\trej %d\trt %le\n", totPropG, rejPropG, ((double)rejPropG)/((double)totPropG));
printf("\n\nSStot %d\trej %d\trt %le\n", totPropS, rejPropS, ((double)rejPropS)/((double)totPropS));
printf("\n\nFFtot %d\trej %d\trt %le\n", totPropF, rejPropF, ((double)rejPropF)/((double)totPropF));
printf("\n\nPPtot %d\trej %d\trt %le\n", totProp-totPropT-totPropF, rejProp-rejPropT-rejPropF, ((double)rejProp-rejPropT-rejPropF)/((double)totProp-totPropT-totPropF));

	for(i=0; i<sampScheme.sizeParamTime; i++)
		for(s=0; s<iter; s++)
			fprintf(fout, "%d\t%lf\n",  s+burn+1, sampCurrentSave[s].currentTime[i]);
	for(i=0; i<sampScheme.sizeParamBirth; i++)
		for(s=0; s<iter; s++)
			fprintf(fout, "%d\t%lf\n",  s+burn+1, sampCurrentSave[s].currentBirth[i]);
	for(i=0; i<sampScheme.sizeParamDeath; i++)
		for(s=0; s<iter; s++)
			fprintf(fout, "%d\t%lf\n",  s+burn+1, sampCurrentSave[s].currentDeath[i]);
	for(i=0; i<sampScheme.sizeParamFossil; i++)
		for(s=0; s<iter; s++)
			fprintf(fout, "%d\t%lf\n",  s+burn+1, sampCurrentSave[s].currentFossil[i]);
	for(i=0; i<sampScheme.sizeParamSampling; i++)
		for(s=0; s<iter; s++)
			fprintf(fout, "%d\t%lf\n",  s+burn+1, sampCurrentSave[s].currentSampling[i]);
	free((void*)param.param);
	if(sampCurrent.currentTime != NULL) {
		free((void*) sampCurrent.currentTime);
		for(s=0; s<iter; s++)
			free((void*)sampCurrentSave[s].currentTime);	
	}
	if(sampCurrent.currentBirth != NULL) {
		free((void*) sampCurrent.currentBirth);
		for(s=0; s<iter; s++)
			free((void*)sampCurrentSave[s].currentBirth);	
	}
	if(sampCurrent.currentDeath != NULL) {
		free((void*) sampCurrent.currentDeath);
		for(s=0; s<iter; s++)
			free((void*)sampCurrentSave[s].currentDeath);	
	}
	if(sampCurrent.currentFossil != NULL) {
		free((void*) sampCurrent.currentFossil);
		for(s=0; s<iter; s++)
			free((void*)sampCurrentSave[s].currentFossil);	
	}
	if(sampCurrent.currentSampling != NULL) {
		free((void*) sampCurrent.currentSampling);
		for(s=0; s<iter; s++)
			free((void*)sampCurrentSave[s].currentSampling);	
	}
	int off = 0;
	for(i=0; i<sampScheme.sizeParamTime; i++) {
		fprintf(find, "time_%d\t%d\t%d\n", i, off*iter+1, (off+1)*iter);
		off++;
	}
	for(i=0; i<sampScheme.sizeParamBirth; i++) {
		fprintf(find, "birth_%d\t%d\t%d\n", i, off*iter+1, (off+1)*iter);
		off++;
	}
	for(i=0; i<sampScheme.sizeParamDeath; i++) {
		fprintf(find, "death_%d\t%d\t%d\n", i, off*iter+1, (off+1)*iter);
		off++;
	}
	for(i=0; i<sampScheme.sizeParamFossil; i++) {
		fprintf(find, "fossil_%d\t%d\t%d\n", i, off*iter+1, (off+1)*iter);
		off++;
	}
	for(i=0; i<sampScheme.sizeParamSampling; i++) {
		fprintf(find, "sampling_%d\t%d\t%d\n", i, off*iter+1, (off+1)*iter);
		off++;
	}
	free((void*)branch);
	freeExtendedFossilFeature(eff, tree[0]->size);
}



/*compute the posterior distribution of the free parameters of the sampScheme (not that of the fossils*/
/*it is assumed that all the branch with identifiers in the trees have the same index, the tree nodes with fossils must have names*/
void MCMCSamplingSamplePosteriorParametersOnlyInit(FILE *fout, FILE *find, TypeTree **tree, int nTree, TypeFossilIntFeature *fi, TypeFossilIntFeature *fo, TypeSamplingScheme sampScheme, double al, int burn, int gap, int iter, double prop, TypeModelParam *windSize, TypeModelParam *init, double probTime, double probSpe, double probExt, double probFos, void *rand_data) {
	int *branch, i, s, n;
	double prob;
	TypeExtendedFossilFeature *eff;
	TypePiecewiseModelParam param;

	TypeTypeVariableData tvd;
	TypeSamplingCurrent sampCurrent, *sampCurrentSave;
	tvd = getTypeVariableData(probTime, probSpe, probExt, probFos, sampScheme);
	sampCurrentSave = (TypeSamplingCurrent*) malloc(iter*sizeof(TypeSamplingCurrent));
	param.size =  sampScheme.param.size;
	param.startTime =  sampScheme.param.startTime;
	param.param =  (TypeModelParam*) malloc(param.size*sizeof(TypeModelParam));
	if(sampScheme.sizeParamTime>0) {
		int i;
		sampCurrent.currentTime = (double*) malloc(sampScheme.sizeParamTime*sizeof(double));
		for(i=0; i<sampScheme.sizeParamTime; i++)
			sampCurrent.currentTime[i] = param.startTime[sampScheme.indexTime[i]-1]+(param.startTime[sampScheme.indexTime[i]+1]-param.startTime[sampScheme.indexTime[i]-1])*my_rand_unif_unit(rand_data);
		for(i=0; i<iter; i++)
			sampCurrentSave[i].currentTime = (double*) malloc(sampScheme.sizeParamTime*sizeof(double));
	} else {
		sampCurrent.currentTime = NULL;
		for(i=0; i<iter; i++)
			sampCurrentSave[i].currentTime = NULL;
	}
	if(sampScheme.sizeParamBirth>0) {
		int i;
		sampCurrent.currentBirth = (double*) malloc(sampScheme.sizeParamBirth*sizeof(double));
		for(i=0; i<sampScheme.sizeParamBirth; i++)
			sampCurrent.currentBirth[i] = my_rand_unif_unit(rand_data)*init->birth;
		for(i=0; i<iter; i++)
			sampCurrentSave[i].currentBirth = (double*) malloc(sampScheme.sizeParamBirth*sizeof(double));
	} else {
		sampCurrent.currentBirth = NULL;
		for(i=0; i<iter; i++)
			sampCurrentSave[i].currentBirth = NULL;
	}
	if(sampScheme.sizeParamDeath>0) {
		int i;
		sampCurrent.currentDeath = (double*) malloc(sampScheme.sizeParamDeath*sizeof(double));
		for(i=0; i<sampScheme.sizeParamDeath; i++)
			sampCurrent.currentDeath[i] = my_rand_unif_unit(rand_data)*init->death;
		for(i=0; i<iter; i++)
			sampCurrentSave[i].currentDeath = (double*) malloc(sampScheme.sizeParamDeath*sizeof(double));
	} else {
		sampCurrent.currentDeath = NULL;
		for(i=0; i<iter; i++)
			sampCurrentSave[i].currentDeath = NULL;
	}
	if(sampScheme.sizeParamFossil>0) {
		int i;
		sampCurrent.currentFossil = (double*) malloc(sampScheme.sizeParamFossil*sizeof(double));
		for(i=0; i<sampScheme.sizeParamFossil; i++)
			sampCurrent.currentFossil[i] = my_rand_unif_unit(rand_data)*init->fossil;
		for(i=0; i<iter; i++)
			sampCurrentSave[i].currentFossil = (double*) malloc(sampScheme.sizeParamFossil*sizeof(double));
	} else {
		sampCurrent.currentFossil = NULL;
		for(i=0; i<iter; i++)
			sampCurrentSave[i].currentFossil = NULL;
	}
	if(sampScheme.sizeParamSampling>0) {
		int i;
		sampCurrent.currentSampling = (double*) malloc(sampScheme.sizeParamSampling*sizeof(double));
		for(i=0; i<sampScheme.sizeParamSampling; i++)
			sampCurrent.currentSampling[i] = my_rand_unif_unit(rand_data);
		for(i=0; i<iter; i++)
			sampCurrentSave[i].currentSampling = (double*) malloc(sampScheme.sizeParamSampling*sizeof(double));
	} else {
		sampCurrent.currentSampling = NULL;
		for(i=0; i<iter; i++)
			sampCurrentSave[i].currentSampling = NULL;
	}
	for(i=0; i<=param.size; i++) {
		if(sampScheme.schemeTime[i] >= 0)
			param.startTime[i] = sampCurrent.currentTime[sampScheme.schemeTime[i]];
		else
			param.startTime[i] = sampScheme.param.startTime[i];
	}
	for(i=0; i<param.size; i++) {
		if(sampScheme.schemeBirth[i] >= 0)
			param.param[i].birth = sampCurrent.currentBirth[sampScheme.schemeBirth[i]];
		else
			param.param[i].birth = sampScheme.param.param[i].birth;
		if(sampScheme.schemeDeath[i] >= 0)
			param.param[i].death = sampCurrent.currentDeath[sampScheme.schemeDeath[i]];
		else
			param.param[i].death = sampScheme.param.param[i].death;
		if(sampScheme.schemeFossil[i] >= 0)
			param.param[i].fossil = sampCurrent.currentFossil[sampScheme.schemeFossil[i]];
		else
			param.param[i].fossil = sampScheme.param.param[i].fossil;
		if(sampScheme.schemeSampling[i] >= 0)
			param.param[i].sampling = sampCurrent.currentSampling[sampScheme.schemeSampling[i]];
		else
			param.param[i].sampling = sampScheme.param.param[i].sampling;
	}
	eff = extendFossilFeature(fi, tree, nTree, rand_data);
	for(i=0; i<nTree; i++) {
		for(n=0; n<tree[i]->size; n++) {
			if(tree[i]->node[n].child == NOSUCH) {
				switch(fi->status[n]) {
					case contempNodeStatus:
						tree[i]->time[n] = tree[i]->maxTime;
					break;
					case unknownNodeStatus:
					break;
					case extinctNodeStatus:
						tree[i]->time[n] = eff->ff->fossilList[eff->ff->fossil[n]].time;
					break;
					default:
						error("Node %d has no status\n", n);
				}
			} else
				tree[i]->time[n] = NO_TIME;
		}
	}
	branch = (int*) malloc(eff->ff->size*sizeof(int));
	fillBranchNode(eff, tree[0]->size, branch);

	if(nTree == 1)
		prob = getLogDensityOne(*tree, eff->ff, &param);
	else
		prob = getLogDensitySum(tree, nTree, eff->ff, &param);
	for(i=0; i<burn; i++) {
		prob = updateSchemeSample(tree, nTree, eff, fi, branch, prob, al, prop, &param, windSize, probSpe, probExt, probFos, tvd, sampScheme, sampCurrent, rand_data);
if(i%10 == 0) {
	fprintf(stderr, "it %-7d/%d %le (%.2lf)\r", i, burn, prob, param.startTime[1]); fflush(stderr);
}
	}
	for(s=0; s<iter; s++) {
		int j;
		for(i=0; i<sampScheme.sizeParamTime; i++)
			sampCurrentSave[s].currentTime[i] = sampCurrent.currentTime[i];
		for(i=0; i<sampScheme.sizeParamBirth; i++)
			sampCurrentSave[s].currentBirth[i] = sampCurrent.currentBirth[i];
		for(i=0; i<sampScheme.sizeParamDeath; i++)
			sampCurrentSave[s].currentDeath[i] = sampCurrent.currentDeath[i];
		for(i=0; i<sampScheme.sizeParamFossil; i++)
			sampCurrentSave[s].currentFossil[i] = sampCurrent.currentFossil[i];
		for(i=0; i<sampScheme.sizeParamSampling; i++)
			sampCurrentSave[s].currentSampling[i] = sampCurrent.currentSampling[i];
		for(j=0; j<gap; j++) {
			prob = updateSchemeSample(tree, nTree, eff, fi, branch, prob, al, prop, &param, windSize, probSpe, probExt, probFos, tvd, sampScheme, sampCurrent, rand_data);
		}
if(s%10 == 0) {
	fprintf(stderr, "it %-7d/%d %le (%.2lf)\r", s, iter, prob, param.startTime[1]); fflush(stderr);
}
	}

	for(i=0; i<sampScheme.sizeParamTime; i++)
		for(s=0; s<iter; s++)
			fprintf(fout, "%d\t%lf\n",  s+burn+1, sampCurrentSave[s].currentTime[i]);
	for(i=0; i<sampScheme.sizeParamBirth; i++)
		for(s=0; s<iter; s++)
			fprintf(fout, "%d\t%lf\n",  s+burn+1, sampCurrentSave[s].currentBirth[i]);
	for(i=0; i<sampScheme.sizeParamDeath; i++)
		for(s=0; s<iter; s++)
			fprintf(fout, "%d\t%lf\n",  s+burn+1, sampCurrentSave[s].currentDeath[i]);
	for(i=0; i<sampScheme.sizeParamFossil; i++)
		for(s=0; s<iter; s++)
			fprintf(fout, "%d\t%lf\n",  s+burn+1, sampCurrentSave[s].currentFossil[i]);
	for(i=0; i<sampScheme.sizeParamSampling; i++)
		for(s=0; s<iter; s++)
			fprintf(fout, "%d\t%lf\n",  s+burn+1, sampCurrentSave[s].currentSampling[i]);
	free((void*)param.param);
	if(sampCurrent.currentTime != NULL) {
		free((void*) sampCurrent.currentTime);
		for(s=0; s<iter; s++)
			free((void*)sampCurrentSave[s].currentTime);	
	}
	if(sampCurrent.currentBirth != NULL) {
		free((void*) sampCurrent.currentBirth);
		for(s=0; s<iter; s++)
			free((void*)sampCurrentSave[s].currentBirth);	
	}
	if(sampCurrent.currentDeath != NULL) {
		free((void*) sampCurrent.currentDeath);
		for(s=0; s<iter; s++)
			free((void*)sampCurrentSave[s].currentDeath);	
	}
	if(sampCurrent.currentFossil != NULL) {
		free((void*) sampCurrent.currentFossil);
		for(s=0; s<iter; s++)
			free((void*)sampCurrentSave[s].currentFossil);	
	}
	if(sampCurrent.currentSampling != NULL) {
		free((void*) sampCurrent.currentSampling);
		for(s=0; s<iter; s++)
			free((void*)sampCurrentSave[s].currentSampling);	
	}
	int off = 0;
	for(i=0; i<sampScheme.sizeParamTime; i++) {
		fprintf(find, "time_%d\t%d\t%d\n", i, off*iter+1, (off+1)*iter);
		off++;
	}
	for(i=0; i<sampScheme.sizeParamBirth; i++) {
		fprintf(find, "birth_%d\t%d\t%d\n", i, off*iter+1, (off+1)*iter);
		off++;
	}
	for(i=0; i<sampScheme.sizeParamDeath; i++) {
		fprintf(find, "death_%d\t%d\t%d\n", i, off*iter+1, (off+1)*iter);
		off++;
	}
	for(i=0; i<sampScheme.sizeParamFossil; i++) {
		fprintf(find, "fossil_%d\t%d\t%d\n", i, off*iter+1, (off+1)*iter);
		off++;
	}
	for(i=0; i<sampScheme.sizeParamSampling; i++) {
		fprintf(find, "sampling_%d\t%d\t%d\n", i, off*iter+1, (off+1)*iter);
		off++;
	}
	free((void*)branch);
	freeExtendedFossilFeature(eff, tree[0]->size);
}


/*compute the posterior distribution of the free parameters of the sampScheme (not that of the fossils*/
/*it is assumed that all the branch with identifiers in the trees have the same index, the tree nodes with fossils must have names*/
void MCMCSamplingSamplePosteriorParametersFossils(FILE *fout, FILE *find, TypeTree **tree, int nTree, TypeFossilIntFeature *fi, TypeSamplingScheme sampScheme, double al, int burn, int gap, int iter, double prop, TypeModelParam *windSize, TypeModelParam *init, double probSpe, double probExt, double probFos, void *rand_data) {
	int *branch, i, f, s, n;
	double prob, **save;
	TypeExtendedFossilFeature *eff;
	TypePiecewiseModelParam param;

	TypeTypeVariableData tvd;
	TypeSamplingCurrent sampCurrent, *sampCurrentSave;

	tvd = getTypeVariableData(0., probSpe, probExt, probFos, sampScheme);
	param.size =  sampScheme.param.size;
	param.startTime =  sampScheme.param.startTime;
	param.param =  (TypeModelParam*) malloc(param.size*sizeof(TypeModelParam));
	sampCurrentSave = (TypeSamplingCurrent*) malloc(iter*sizeof(TypeSamplingCurrent));
	if(sampScheme.sizeParamBirth>0) {
		int i;
		sampCurrent.currentBirth = (double*) malloc(sampScheme.sizeParamBirth*sizeof(double));
		for(i=0; i<sampScheme.sizeParamBirth; i++)
			sampCurrent.currentBirth[i] = my_rand_unif_unit(rand_data)*init->birth;
		for(i=0; i<iter; i++)
			sampCurrentSave[i].currentBirth = (double*) malloc(sampScheme.sizeParamBirth*sizeof(double));
	} else {
		sampCurrent.currentBirth = NULL;
		for(i=0; i<iter; i++)
			sampCurrentSave[i].currentBirth = NULL;
	}
	if(sampScheme.sizeParamDeath>0) {
		int i;
		sampCurrent.currentDeath = (double*) malloc(sampScheme.sizeParamDeath*sizeof(double));
		for(i=0; i<sampScheme.sizeParamDeath; i++)
			sampCurrent.currentDeath[i] = my_rand_unif_unit(rand_data)*init->death;
		for(i=0; i<iter; i++)
			sampCurrentSave[i].currentDeath = (double*) malloc(sampScheme.sizeParamDeath*sizeof(double));
	} else {
		sampCurrent.currentDeath = NULL;
		for(i=0; i<iter; i++)
			sampCurrentSave[i].currentDeath = NULL;
	}
	if(sampScheme.sizeParamFossil>0) {
		int i;
		sampCurrent.currentFossil = (double*) malloc(sampScheme.sizeParamFossil*sizeof(double));
		for(i=0; i<sampScheme.sizeParamFossil; i++)
			sampCurrent.currentFossil[i] = my_rand_unif_unit(rand_data)*init->fossil;
		for(i=0; i<iter; i++)
			sampCurrentSave[i].currentFossil = (double*) malloc(sampScheme.sizeParamFossil*sizeof(double));
	} else {
		sampCurrent.currentFossil = NULL;
		for(i=0; i<iter; i++)
			sampCurrentSave[i].currentFossil = NULL;
	}
	if(sampScheme.sizeParamSampling>0) {
		int i;
		sampCurrent.currentSampling = (double*) malloc(sampScheme.sizeParamSampling*sizeof(double));
		for(i=0; i<sampScheme.sizeParamSampling; i++)
			sampCurrent.currentSampling[i] = my_rand_unif_unit(rand_data);
		for(i=0; i<iter; i++)
			sampCurrentSave[i].currentSampling = (double*) malloc(sampScheme.sizeParamSampling*sizeof(double));
	} else {
		sampCurrent.currentSampling = NULL;
		for(i=0; i<iter; i++)
			sampCurrentSave[i].currentSampling = NULL;
	}
	for(i=0; i<param.size; i++) {
		if(sampScheme.schemeBirth[i] >= 0)
			param.param[i].birth = sampCurrent.currentBirth[sampScheme.schemeBirth[i]];
		else
			param.param[i].birth = sampScheme.param.param[i].birth;
		if(sampScheme.schemeDeath[i] >= 0)
			param.param[i].death = sampCurrent.currentDeath[sampScheme.schemeDeath[i]];
		else
			param.param[i].death = sampScheme.param.param[i].death;
		if(sampScheme.schemeFossil[i] >= 0)
			param.param[i].fossil = sampCurrent.currentFossil[sampScheme.schemeFossil[i]];
		else
			param.param[i].fossil = sampScheme.param.param[i].fossil;
		if(sampScheme.schemeSampling[i] >= 0)
			param.param[i].sampling = sampCurrent.currentSampling[sampScheme.schemeSampling[i]];
		else
			param.param[i].sampling = sampScheme.param.param[i].sampling;
	}
	eff = extendFossilFeature(fi, tree, nTree, rand_data);
	for(i=0; i<nTree; i++) {
		for(n=0; n<tree[i]->size; n++) {
			if(tree[i]->node[n].child == NOSUCH) {
				switch(fi->status[n]) {
					case contempNodeStatus:
						tree[i]->time[n] = tree[i]->maxTime;
					break;
					case unknownNodeStatus:
					break;
					case extinctNodeStatus:
						tree[i]->time[n] = eff->ff->fossilList[eff->ff->fossil[n]].time;
					break;
					default:
						error("Node %d has no status\n", n);
				}
			} else
				tree[i]->time[n] = NO_TIME;
		}
	}
	branch = (int*) malloc(eff->ff->size*sizeof(int));
	fillBranchNode(eff, tree[0]->size, branch);
	save = (double**) malloc(eff->ff->size*sizeof(double*));
	for(f=0; f<eff->ff->size; f++)
		save[f] = (double*) malloc(iter*sizeof(double));

	if(nTree == 1)
		prob = getLogDensityOne(*tree, eff->ff, &param);
	else
		prob = getLogDensitySum(tree, nTree, eff->ff, &param);
	for(i=0; i<burn; i++) {
		prob = updateSchemeSample(tree, nTree, eff, fi, branch, prob, al, prop, &param, windSize, probSpe, probExt, probFos, tvd, sampScheme, sampCurrent, rand_data);
if(i%10 == 0) {
	fprintf(stderr, "it %-7d\r", i); fflush(stderr);
}
	}
	for(s=0; s<iter; s++) {
		int j;
		for(i=0; i<sampScheme.sizeParamBirth; i++)
			sampCurrentSave[s].currentBirth[i] = sampCurrent.currentBirth[i];
		for(i=0; i<sampScheme.sizeParamDeath; i++)
			sampCurrentSave[s].currentDeath[i] = sampCurrent.currentDeath[i];
		for(i=0; i<sampScheme.sizeParamFossil; i++)
			sampCurrentSave[s].currentFossil[i] = sampCurrent.currentFossil[i];
		for(i=0; i<sampScheme.sizeParamSampling; i++)
			sampCurrentSave[s].currentSampling[i] = sampCurrent.currentSampling[i];
		for(f=0; f<eff->ff->size; f++)
			save[f][s] = eff->ff->fossilList[f].time;
		for(j=0; j<gap; j++) {
			prob = updateSchemeSample(tree, nTree, eff, fi, branch, prob, al, prop, &param, windSize, probSpe, probExt, probFos, tvd, sampScheme, sampCurrent, rand_data);
		}

if(i%10 == 0) {
	fprintf(stderr, "it %-7d\r", i); fflush(stderr);
}
	}

	for(i=0; i<sampScheme.sizeParamBirth; i++)
		for(s=0; s<iter; s++)
			fprintf(fout, "%d\t%lf\n",  s+burn+1, sampCurrentSave[s].currentBirth[i]);
	for(i=0; i<sampScheme.sizeParamDeath; i++)
		for(s=0; s<iter; s++)
			fprintf(fout, "%d\t%lf\n",  s+burn+1, sampCurrentSave[s].currentDeath[i]);
	for(i=0; i<sampScheme.sizeParamFossil; i++)
		for(s=0; s<iter; s++)
			fprintf(fout, "%d\t%lf\n",  s+burn+1, sampCurrentSave[s].currentFossil[i]);
	for(i=0; i<sampScheme.sizeParamSampling; i++)
		for(s=0; s<iter; s++)
			fprintf(fout, "%d\t%lf\n",  s+burn+1, sampCurrentSave[s].currentSampling[i]);
	free((void*)param.param);
	for(f=0; f<eff->ff->size; f++)
		for(s=0; s<iter; s++)
			fprintf(fout, "%d\t%lf\n",  s+burn+1, save[f][s]);
	for(f=0; f<eff->ff->size; f++)
		free((void*)save[f]);
	free((void*)save);
	if(sampCurrent.currentBirth != NULL) {
		free((void*) sampCurrent.currentBirth);
		for(s=0; s<iter; s++)
			free((void*)sampCurrentSave[s].currentBirth);	
	}
	if(sampCurrent.currentDeath != NULL) {
		free((void*) sampCurrent.currentDeath);
		for(s=0; s<iter; s++)
			free((void*)sampCurrentSave[s].currentDeath);	
	}
	if(sampCurrent.currentFossil != NULL) {
		free((void*) sampCurrent.currentFossil);
		for(s=0; s<iter; s++)
			free((void*)sampCurrentSave[s].currentFossil);	
	}
	if(sampCurrent.currentSampling != NULL) {
		free((void*) sampCurrent.currentSampling);
		for(s=0; s<iter; s++)
			free((void*)sampCurrentSave[s].currentSampling);	
	}
	int off = 0;
	for(i=0; i<sampScheme.sizeParamBirth; i++) {
		fprintf(find, "birth_%d\t%d\t%d\n", i, off*iter+1, (off+1)*iter);
		off++;
	}
	for(i=0; i<sampScheme.sizeParamDeath; i++) {
		fprintf(find, "death_%d\t%d\t%d\n", i, off*iter+1, (off+1)*iter);
		off++;
	}
	for(i=0; i<sampScheme.sizeParamFossil; i++) {
		fprintf(find, "fossil_%d\t%d\t%d\n", i, off*iter+1, (off+1)*iter);
		off++;
	}
	for(i=0; i<sampScheme.sizeParamSampling; i++) {
		fprintf(find, "sampling_%d\t%d\t%d\n", i, off*iter+1, (off+1)*iter);
		off++;
	}
	for(f=0; f<eff->ff->size; f++) {
		fprintf(find, "%d", f+1);
		if(branch[f] != NOSUCH && tree[0]->name[branch[f]] != NULL)
			fprintf(find, "_%s_%d_%d", tree[0]->name[branch[f]], (int) round(fi->fossilIntList[f].fossilInt.inf), (int) round(fi->fossilIntList[f].fossilInt.sup));
		fprintf(find, "\t%d\t%d\n", off*iter+1, (off+1)*iter); off++;
	}
	free((void*)branch);
	freeExtendedFossilFeature(eff, tree[0]->size);
}
