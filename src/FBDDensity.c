#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>


#include "Utils.h"
#include "MyR.h"
//#include "TreeExtras.h"
#include "FBDDensity.h"



#define FBD_DBL_PREC 0.000000000001

typedef enum {
	FBD_none_type = 0, /*div time can be anytime*/
	FBD_constraint_type, /*div time is anterior and/or posterior*/
	FBD_fossil_type,   /*tip time is a fossil dated at time*/
	FBD_stop_type,     /*tip time is a lineage alive at time but with fate unknown after*/
	FBD_end_type       /*tip time goes until the end of diversification, i.e., the present time*/
} TypeFBDInfoType;

typedef struct FBD_COEFF_ROW {
    double c, alpha, beta, omega, bma;
} TypeFBDCoeffRow;

typedef struct FBD_COEFF {
    TypePiecewiseModelParam *model;
    TypeFBDCoeffRow *row;
    double *logProbNotObs;
} TypeFBDCoeff;

typedef struct FBD_INFO {
	TypeFBDInfoType type; /*if the corresponding node is internal may be FBD_none_type or FBD_constraint_type and may be FBD_fossil_type, FBD_stop_type or FBD_end_type if it is a tip see above*/
	int nTips; /*number tips descending from the node*/
	double min, max; /*range of possible datation of the node*/
} TypeFBDInfo;

typedef struct FBD_ANNEXE {
	double *W; /*table of quantities required to compute the density for numbre of lineages from min to min+size-1*/
	int min, size, nStops, nFossils; /*min: min number of outing lineages; size: min+size-1 is the max number of outing lineages; nStops: number of lineages with unknown fate after a given time; nFossils: number of lineages ending with a fossil*/
} TypeFBDAnnexe;


static double lfactorial(int n);
/*return log(exp(a)+exp(b)) in an accurate way*/
static double logSumLog(double a, double b);
/*return log(exp(a)-exp(b)) in an accurate way*/
static double logDiffLog(double a,double b);
static double logProbNotObs(double t, int cur, TypeFBDCoeff *coeff);
static double logProbObs(double t, int cur, TypeFBDCoeff *coeff);
static double probObs(double t, int cur, TypeFBDCoeff *coeff);
static TypeFBDAnnexe *fillFBDAnnexe(TypeTree *tree, TypeFBDInfo *info, double *time, int nTimes, TypeFBDCoeff *coeff, TypeFBDAnnexe *annexeA, TypeFBDAnnexe *annexeB);
static void splitTreeFossilFBD(TypeTree *tree, TypeFBDInfo *info, TypeFossilFeature *fos, TypeTree ***treeList, int *size);
static void fillFBDInfo(int n, double curMin, TypeTree *tree, TypeFBDInfo *info);
static double logProbStandard(int k, int l, double tstart, double tend, int cur, TypeFBDCoeff *coeff, double lpno);
static double logProbFossil(int k, int l, double tstart, double tend, int cur, TypeFBDCoeff *coeff, double lpno);
static TypeFBDCoeff getCoeff(TypePiecewiseModelParam *model);
static int fillBasicTreeFBD(TypeTree *basic, int n, TypeTree *tree, TypeFossilTab *ftab, int toCompute, int *newIndex);
static double logProbExtinctCond(double tA, int curA, double tB, int curB, double maxTime, TypeFBDCoeff *coeff);
static TypeTree *getBasicTreeDivergenceFBD(int n, TypeTree *tree, TypeFossilTab *ftab, int *newIndex);
static double logLikelihoodBasicTreeFossil(TypeTree *tree, double *stop, TypeFBDCoeff *coeff);



typedef struct INTEGRATION_FBD_DATA {
	int k, l;
    double tend;
    TypeFBDCoeff *coeff;
} TypeIntegrationFBDData;


double getMinTimeTree(int n, TypeTree *tree, TypeFossilTab *ftab) {
	if(ftab!= NULL && ftab[n].size>0)
		return ftab[n].time[0];
	if(tree->time[n] != NO_TIME)
		return tree->time[n];
	if(tree->node[n].child == NOSUCH) {
		error("Execution error leaf %d (%s) with no time\n", n, (tree->name && tree->name[n])?tree->name[n]:"??");
	}
	return utils_MIN(getMinTimeTree(tree->node[n].child, tree, ftab), getMinTimeTree(tree->node[tree->node[n].child].sibling, tree, ftab));
}

#define PRECISION_LOG_LIKE -5*log(10.)

void fillExtinctionLogDistributionFBD(TypeDistribution *logD, int n, TypeTree *tree, TypeFossilFeature *fos, TypePiecewiseModelParam *model) {
	if(tree->node[n].child == NOSUCH && tree->time[n] != NO_TIME) {
		TypeFBDCoeff coeff = getCoeff(model);
		int cur, i;
		cur = getPieceIndex(tree->time[n], model);
		for(i=0; i<logD->size; i++)
			if(logD->item[i].val >= tree->time[n])
				logD->item[i].dens = logProbExtinctCond(tree->time[n], cur, logD->item[i].val, getPieceIndex(logD->item[i].val, model), tree->maxTime, &coeff);
			else
				logD->item[i].dens = NEG_INFTY;				
	}
}

void fillDivergenceLogDistributionFBD(TypeDistribution *logD, double *logCond, int n, TypeTree *tree, TypeFossilFeature *fos, TypePiecewiseModelParam *model) {
	int index, i;
	TypeTree *basicU;
	double *stop, diff, max;
	TypeFBDCoeff coeff = getCoeff(model);
	
	TypeFossilTab *ftab = listToFossilTab(fos, tree->size);
	*logCond = logLikelihoodTreeFossil(tree, fos, NULL, model);
	basicU = getBasicTreeDivergenceFBD(n, tree, ftab, &index);
	stop = (double*) malloc((3*basicU->size+fos->size+model->size+5)*sizeof(double));
	diff = *logCond-logLikelihoodBasicTreeFossil(basicU, stop, &coeff);
	max = getMinTimeTree(index, basicU, NULL);
	for(i=0; i<logD->size && logD->item[i].val<basicU->minTime; i++)
		logD->item[i].dens = NEG_INFTY;
	int a=i, b=logD->size-1, m;
	while(b>a+1) {
		m = (a+b)/2;
		if(logD->item[m].val<max)
			a = m;
		else
			b = m;
	}
	b=a; a=i;
	while(b>a+1) {
		double tmp;
		m = (a+b)/2;
		((TypeFBDInfo*)basicU->info)[index].max = logD->item[m].val;
		tmp = logLikelihoodBasicTreeFossil(basicU, stop, &coeff)+diff;
		if(tmp<PRECISION_LOG_LIKE)
			a = m;
		else
			b = m;
	}
	for(; i<a; i++)
		logD->item[i].dens = NEG_INFTY;
	for(; i<logD->size && logD->item[i].val<max; i++) {
		((TypeFBDInfo*)basicU->info)[index].max = logD->item[i].val;
		logD->item[i].dens = logLikelihoodBasicTreeFossil(basicU, stop, &coeff)+diff;
	}
	for(; i<logD->size; i++)
		logD->item[i].dens = 0.;
	free((void*)stop);
	for(i=0; i<tree->size; i++)
		if(ftab[i].time != NULL)
			free((void*)ftab[i].time);
	free((void*)ftab);
	if(basicU->info!=NULL)
		free((void*)basicU->info);
	freeTree(basicU);
}

double lfactorial(int n) {
	return  lgammafn((double)n+1.);
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

/*return log(exp(a)-exp(b)) in an accurate way*/
double logDiffLog(double a,double b) {
	return a + log1p(-exp(b-a));
}

/*return the probability that a lineage alive at last tA (e.g. last fossil age) goes extinct before time tB without  leaving any fossil between startTime and  t, conditionned on the that lineage goes extinct before maxTime*/
double logProbExtinctCond(double tA, int curA, double tB, int curB, double maxTime, TypeFBDCoeff *coeff) {
	if(curA == curB)
		return log(coeff->row[curA].alpha)+log(coeff->row[curA].beta)+log1p(-exp(coeff->row[curA].omega*(tB-tA)))-log(coeff->row[curA].beta-coeff->row[curA].alpha*exp(coeff->row[curA].omega*(tB-tA)));
	else {
		int cur;
		double d;
		d = 1.-coeff->model->param[curB-1].sampling*(1-coeff->row[curB].alpha*coeff->row[curB].beta*(1.-exp(coeff->row[curB].omega*(tB-coeff->model->startTime[curB])))/(coeff->row[curB].beta-coeff->row[curB].alpha*exp(coeff->row[curB].omega*(tB-coeff->model->startTime[curB]))));
		for(cur=curB-1; cur>curA; cur--)
			d = 1.-coeff->model->param[cur-1].sampling*(1-(coeff->row[cur].alpha*(coeff->row[cur].beta-d) -coeff->row[cur].beta*(coeff->row[cur].alpha-d)*(1.-exp(coeff->row[cur].omega*(coeff->model->startTime[cur+1]-coeff->model->startTime[cur]))))/(coeff->row[cur].beta-d-(coeff->row[cur].alpha-d)*exp(coeff->row[cur].omega*(coeff->model->startTime[cur+1]-coeff->model->startTime[cur]))));
		return log((coeff->row[curA].alpha*(coeff->row[curA].beta-d) -coeff->row[curA].beta*(coeff->row[curA].alpha-d)*(1.-exp(coeff->row[curA].omega*(coeff->model->startTime[curA+1]-tA))))/(coeff->row[curA].beta-d-(coeff->row[curA].alpha-d)*exp(coeff->row[curA].omega*(coeff->model->startTime[curA+1]-tA))));
	}
}

double logProbNotObs(double t, int cur, TypeFBDCoeff *coeff) {
	return log(coeff->row[cur].alpha*(coeff->row[cur].beta-coeff->row[cur].c)-coeff->row[cur].beta*(coeff->row[cur].alpha-coeff->row[cur].c)*exp(coeff->row[cur].omega*(coeff->model->startTime[cur+1]-t)))-log(coeff->row[cur].beta-coeff->row[cur].c-(coeff->row[cur].alpha-coeff->row[cur].c)*exp(coeff->row[cur].omega*(coeff->model->startTime[cur+1]-t)));
}

double logProbObs(double t, int cur, TypeFBDCoeff *coeff) {
	return log1p(-exp(logProbNotObs(t, cur, coeff)));
}			

double probObs(double t, int cur, TypeFBDCoeff *coeff) {
	return exp(logProbObs(t, cur, coeff));
}
			
double logProbStandard(int k, int l, double tstart, double tend, int cur, TypeFBDCoeff *coeff, double lpno) {
		return
			((double)l)*log1p(-exp(lpno))
			+2.*log(coeff->row[cur].bma)
			+coeff->row[cur].omega*(tend-tstart)
			+((double)k-1)*log1p(-exp(coeff->row[cur].omega*(tend-tstart)))
			-((double)k+1)*log(coeff->row[cur].beta-exp(lpno)-(coeff->row[cur].alpha-exp(lpno))*exp(coeff->row[cur].omega*(tend-tstart)));
}

double logProbFossil(int k, int l, double tstart, double tend, int cur, TypeFBDCoeff *coeff, double lpno) {
	return
		((double)l)*log1p(-exp(lpno))
		+log((double)k)
		+log(coeff->model->param[cur].fossil)
		+2.*log(coeff->row[cur].bma)
		+coeff->row[cur].omega*(tend-tstart)
		+((double)k-1)*log1p(-exp(coeff->row[cur].omega*(tend-tstart)))
		-((double)k+1)*log(coeff->row[cur].beta-exp(lpno)-(coeff->row[cur].alpha-exp(lpno))*exp(coeff->row[cur].omega*(tend-tstart)));
}

TypeFBDCoeff getCoeff(TypePiecewiseModelParam *model) {
	TypeFBDCoeff res;
	int i;
	res.model = model;
	res.row = (TypeFBDCoeffRow*) malloc(res.model->size*sizeof(TypeFBDCoeffRow));
	for(i=0; i<res.model->size; i++) {
		res.row[i].alpha = (res.model->param[i].birth+res.model->param[i].death+res.model->param[i].fossil-sqrt(pow(res.model->param[i].birth+res.model->param[i].death+res.model->param[i].fossil, 2.)-4.*res.model->param[i].birth*res.model->param[i].death))/(2*res.model->param[i].birth);
		res.row[i].beta = (res.model->param[i].birth+res.model->param[i].death+res.model->param[i].fossil+sqrt(pow(res.model->param[i].birth+res.model->param[i].death+res.model->param[i].fossil, 2.)-4.*res.model->param[i].birth*res.model->param[i].death))/(2*res.model->param[i].birth);
		res.row[i].bma = res.row[i].beta-res.row[i].alpha;
		res.row[i].omega = -res.model->param[i].birth*res.row[i].bma;
	}
	res.row[res.model->size-1].c = 1.-res.model->param[res.model->size-1].sampling;
	for(i=res.model->size-2; i>=0; i--)
		res.row[i].c = 1.-res.model->param[i].sampling*probObs(res.model->startTime[i+1],i+1, &res);
	return res;
}



/*Fill min and max times and number of tips of the subtree originating at node n and its descendants*/
void fillFBDInfo(int n, double curMin, TypeTree *tree, TypeFBDInfo *info) {
	if(tree->time[n] != NO_TIME) {
		info[n].min = tree->time[n];
	} else {
		if(info[n].min != NO_TIME)
			info[n].min = utils_MAX(info[n].min, curMin);
		else
			info[n].min = curMin;
	}
	if(tree->node[n].child == NOSUCH) {
		info[n].nTips = 1;
		if(tree->time[n] != NO_TIME)
			info[n].max = tree->time[n];
		else
			error("Error in FillFBDInfo function: tip with no time!\n");
	} else {
		int c;
		for(c=tree->node[n].child; c != NOSUCH; c=tree->node[c].sibling)
			fillFBDInfo(c, info[n].min, tree, info);
		info[n].nTips = info[tree->node[n].child].nTips;
		for(c=tree->node[tree->node[n].child].sibling; c!=NOSUCH; c = tree->node[c].sibling)
			info[n].nTips += info[c].nTips;
		if(tree->time[n] != NO_TIME) {
			info[n].max = tree->time[n];
		} else {
			if(info[n].max != NO_TIME)
				info[n].max = utils_MIN(info[n].max, info[tree->node[n].child].max);
			else
				info[n].max = info[tree->node[n].child].max;
			for(c=tree->node[tree->node[n].child].sibling; c != NOSUCH; c = tree->node[c].sibling)
				info[n].max = utils_MIN(info[n].max, info[c].max);
		}
	}
}



/*For all nodes n of tree, info[n].indexMostRecent = tip with the smallest time in subtree n, 
 * ensures that if one of the smallest times is a fossil, info[n].indexMostRecent point on it*/
double logProbLineage(double tstart, double tend, int n, TypeTree *tree, TypeFBDInfo *info, int cur, TypeFBDCoeff *coeff, TypeFBDAnnexe *annexe, double lpno) {
	if(tree->time[n] == tstart)
		return 0.;
	if(cur == coeff->model->size-1 && tend == coeff->model->startTime[coeff->model->size]) { /*end time case*/
		return logProbStandard(info[n].nTips, info[n].nTips, tstart, tend, cur, coeff, lpno);
	} else {
		if(tree->node[n].child == NOSUCH && tend == tree->time[n]) { /*tree with a single node ending before the present time*/
			if(annexe[n].nFossils > 0)
				return logProbFossil(1, 0, tstart, tend, cur, coeff, lpno);
			else
				return logProbStandard(1, 1, tstart, tend, cur, coeff, lpno);
		} else {
			int k;
			double res = NEG_INFTY;
			if(tend == coeff->model->startTime[cur+1] && coeff->model->param[cur].sampling < 1.) {
				if(annexe[n].nFossils > 0)
					for(k=0; k<annexe[n].size; k++)
						res = logSumLog(res, annexe[n].W[k]+logProbFossil(k+annexe[n].min, annexe[n].nStops, tstart, tend, cur, coeff, lpno)-log((double)k+annexe[n].min))+((double)k+annexe[n].min-1.)*log(coeff->model->param[cur].sampling);
				else
					for(k=0; k<annexe[n].size; k++)
						res = logSumLog(res, annexe[n].W[k]+logProbStandard(k+annexe[n].min, annexe[n].nStops, tstart, tend, cur, coeff, lpno)+((double)k+annexe[n].min)*log(coeff->model->param[cur].sampling));
			} else {
				if(annexe[n].nFossils > 0)
					for(k=0; k<annexe[n].size; k++)
						res = logSumLog(res, annexe[n].W[k]+logProbFossil(k+annexe[n].min, annexe[n].nStops, tstart, tend, cur, coeff, lpno)-log((double)k+annexe[n].min));
				else
					for(k=0; k<annexe[n].size; k++)
						res = logSumLog(res, annexe[n].W[k]+logProbStandard(k+annexe[n].min, annexe[n].nStops, tstart, tend, cur, coeff, lpno));
			}
			return res-lfactorial(info[n].nTips);
		}
	}
	return 0.;
}

TypeFBDAnnexe *fillFBDAnnexe(TypeTree *tree, TypeFBDInfo *info, double *time, int nTimes, TypeFBDCoeff *coeff, TypeFBDAnnexe *annexeA, TypeFBDAnnexe *annexeB) {
	int s, n, cur;
	cur = getPieceIndex((time[nTimes-1]+time[nTimes-2])/2., coeff->model); /*cur is the current index in the piecewise model, initalized between the two last stopping times*/
	if(cur <0) {
		for(s=0; s<nTimes; s++)
			printf("Time %d = %lf\n", s, time[s]);
		printTreeDebug(stdout, tree->root, tree, NULL);
		exit(0);
	}
	for(s=nTimes-1; s>0; s--) {
		TypeFBDAnnexe *annexeTmp;
		if(time[s] < coeff->model->startTime[cur]) /*test if we have to change of slice*/
			cur--;
		for(n=tree->size-1; n>=0; n--) {
			if(tree->node[n].child == NOSUCH) {
				if(tree->time[n] > time[s]) {
					annexeA[n].min = 1;
					annexeA[n].size = 1;
					annexeA[n].nFossils = 0;
					annexeA[n].nStops = 0;
					annexeA[n].W[0] = logProbLineage(time[s], time[s+1], n, tree, info, cur, coeff, annexeB, logProbNotObs(time[s+1], cur, coeff));
				} else {
					if(tree->time[n] == time[s]) {
						annexeA[n].min = 1;
						annexeA[n].size = 1;
						annexeA[n].W[0] = 0.;
						switch(info[n].type) {
							case FBD_fossil_type:
								annexeA[n].nFossils = 1;
								annexeA[n].nStops = 0;
								break;
							case FBD_stop_type:
								annexeA[n].nFossils = 0;
								annexeA[n].nStops = 1;
								break;
							default:
								annexeA[n].nFossils = 0;
								annexeA[n].nStops = 0;
						}
					}
				}
			} else {
				if(info[n].max > time[s-1]) {
					if(info[n].min >= time[s]) {
						annexeA[n].size = 1;
						annexeA[n].min = 1;
						annexeA[n].W[0] = logProbLineage(time[s], time[s+1], n, tree, info, cur, coeff, annexeB, logProbNotObs(time[s+1], cur, coeff))+lfactorial(info[n].nTips);
					} else {
						int i, j, k;
						if(info[n].max > time[s]) {
							annexeA[n].min = 1;
							annexeA[n].W[0] = logProbLineage(time[s], time[s+1], n, tree, info, cur, coeff, annexeB, logProbNotObs(time[s+1], cur, coeff))+lfactorial(info[n].nTips);
						} else
							annexeA[n].min = annexeA[tree->node[n].child].min+annexeA[tree->node[tree->node[n].child].sibling].min;
						annexeA[n].nFossils = annexeA[tree->node[n].child].nFossils+annexeA[tree->node[tree->node[n].child].sibling].nFossils;
						annexeA[n].nStops = annexeA[tree->node[n].child].nStops+annexeA[tree->node[tree->node[n].child].sibling].nStops;
						annexeA[n].size = annexeA[tree->node[n].child].min+annexeA[tree->node[n].child].size+annexeA[tree->node[tree->node[n].child].sibling].min+annexeA[tree->node[tree->node[n].child].sibling].size-annexeA[n].min-1;
						for(k=utils_MAX(2, annexeA[n].min); k<annexeA[n].min+annexeA[n].size; k++)
							annexeA[n].W[k-annexeA[n].min] = NEG_INFTY;
						for(i=0; i<annexeA[tree->node[n].child].size; i++)
							for(j=0; j<annexeA[tree->node[tree->node[n].child].sibling].size; j++)
								annexeA[n].W[i+annexeA[tree->node[n].child].min+j+annexeA[tree->node[tree->node[n].child].sibling].min-annexeA[n].min] = 
									logSumLog(annexeA[n].W[i+annexeA[tree->node[n].child].min+j+annexeA[tree->node[tree->node[n].child].sibling].min-annexeA[n].min], 
										annexeA[tree->node[n].child].W[i]+annexeA[tree->node[tree->node[n].child].sibling].W[j]);
						for(k=utils_MAX(2, annexeA[n].min); k<annexeA[n].min+annexeA[n].size; k++)
							annexeA[n].W[k-annexeA[n].min] += log(2.)- log((double) k-1.);
					}
				}
			}
		}
		annexeTmp = annexeA;
		annexeA = annexeB;
		annexeB = annexeTmp;
	}
	return annexeB;
}


/*get the likelihood of a tree ending at each fossil or contemporary time*/
double logProbSplitted(TypeTree *tree, TypeFBDInfo *info, double *time, int nTimes, TypeFBDCoeff *coeff) {
    int  i, n, cur;
	TypeFBDAnnexe *annexe[2], *annexeF;
	double res;
	if(tree->size == 0) {
		cur = getPieceIndex(tree->minTime, coeff->model);
        if(coeff->model->param[cur].death>0)
            return logProbNotObs(tree->minTime, cur, coeff);
        else
            return NEG_INFTY;
    }
	for(i=0; i<2; i++) {
		annexe[i] = (TypeFBDAnnexe*) malloc(tree->size*sizeof(TypeFBDAnnexe));
		for(n=0; n<tree->size; n++)
			annexe[i][n].W = (double*) malloc(info[n].nTips*sizeof(double));
	}
	annexeF = fillFBDAnnexe(tree, info, time, nTimes, coeff, annexe[0], annexe[1]);
	cur = getPieceIndex(tree->minTime, coeff->model);
	res = logProbLineage(tree->minTime, time[1], tree->root, tree, info, cur, coeff, annexeF, logProbNotObs(time[1], cur, coeff));
	for(i=0; i<2; i++) {
		for(n=0; n<tree->size; n++)
			free((void*)annexe[i][n].W);
		free((void*)annexe[i]);
	}
    return res;
}

/*get the likelihood of a tree ending at each fossil or contemporary time*/
double logProbSplittedX(TypeTree *tree, TypeFBDInfo *info, double *time, int nTimes, TypeFBDCoeff *coeff) {
    int  i, n, cur;
	TypeFBDAnnexe *annexe[2], *annexeF;
	double res;
	if(tree->size == 0) {
		cur = getPieceIndex(tree->minTime, coeff->model);
        if(coeff->model->param[cur].death>0)
            return logProbNotObs(tree->minTime, cur, coeff);
        else
            return NEG_INFTY;
    }
	for(i=0; i<2; i++) {
		annexe[i] = (TypeFBDAnnexe*) malloc(tree->size*sizeof(TypeFBDAnnexe));
		for(n=0; n<tree->size; n++)
			annexe[i][n].W = (double*) malloc(info[n].nTips*sizeof(double));
	}
	annexeF = fillFBDAnnexe(tree, info, time, nTimes, coeff, annexe[0], annexe[1]);
	cur = getPieceIndex(tree->minTime, coeff->model);
printPiecewiseModel(stdout, coeff->model);
	res = logProbLineage(tree->minTime, time[1], tree->root, tree, info, cur, coeff, annexeF, logProbNotObs(time[1], cur, coeff));
printf("res %le\t no %le (%.2lf\t%.2lf %d)\n", res, logProbNotObs(time[1], cur, coeff), time[1], coeff->model->startTime[cur+1], cur);
	for(i=0; i<2; i++) {
		for(n=0; n<tree->size; n++)
			free((void*)annexe[i][n].W);
		free((void*)annexe[i]);
	}
    return res;
}
double logLikelihoodBasicTreeFossil(TypeTree *tree, double *stop, TypeFBDCoeff *coeff) {
	int i, j, k, n, nstop;
	double timeMin, timeMax;
	timeMin = tree->minTime;
	timeMax = timeMin;
	nstop=0;
	stop[nstop++] = timeMin;
	for(n=0; n<tree->size; n++) { /*add all the times from tree[l] to the stop*/
		if(tree->time[n] != NO_TIME) {
			stop[nstop++] = tree->time[n];
			if(tree->time[n] > timeMax)
				timeMax = tree->time[n];
			if(tree->time[n] < timeMin)
				error("Time inconsistency (min %lf tip %lf)\n", timeMin, tree->time[n]);
		}
		if(((TypeFBDInfo*)tree->info)[n].min != NO_TIME)
			stop[nstop++] = ((TypeFBDInfo*)tree->info)[n].min;
		if(((TypeFBDInfo*)tree->info)[n].max != NO_TIME) {
			stop[nstop++] = ((TypeFBDInfo*)tree->info)[n].max;
			if(((TypeFBDInfo*)tree->info)[n].max > timeMax)
				timeMax = ((TypeFBDInfo*)tree->info)[n].max;
		}
	}
//for(i=0; i<nstop; i++)
	//printf("stop %d = %lf\n", i, stop[i]);
	for(i=0; i<=coeff->model->size; i++) /*add the times of the model between timeMin and timeMax to the stops*/
		if(coeff->model->startTime[i]>timeMin && coeff->model->startTime[i]<timeMax)
			stop[nstop++] = coeff->model->startTime[i];
//for(i=0; i<nstop; i++)
	//printf("stop %d = %lf\n", i, stop[i]);
//printf("final number of stops %d\n", nstop);
	if(nstop>0) { /*sort and suppress doublons among the stops*/
		double tmp;
		qsort(stop, nstop, sizeof(double), compareDouble);
		tmp = stop[0];
		j = 1;
		for(k=1; k<nstop; k++)
			if(stop[k] != tmp) {
				stop[j++] = stop[k];
				tmp = stop[k];
			}
		nstop = j;
	}
	fillFBDInfo(tree->root, tree->minTime, tree, (TypeFBDInfo*)tree->info); /*Fill min and max times and number of tips of tree if not empty*/
	return logProbSplitted(tree, (TypeFBDInfo*)tree->info, stop, nstop, coeff);
}

/*get the whole likelihood*/
/* by convention if a tip l is such that
 * - tree->time[l] == NO_TIME => it is extinct (it must have fossil(s)),
 * - tree->time[l] == tree->maxTime => it is contemporary,
 * - otherwise => it is unknown after but observable until tree->time[l]
 */
double logLikelihoodTreeFossil(TypeTree *tree, TypeFossilFeature *fos, TypeFBDTimeConstraint *constraint, TypePiecewiseModelParam *model) {
	TypeTree **treeList;
	int i, n, l, size;
	double *stop, res = 0.;
	TypeFBDCoeff coeff = getCoeff(model);
	TypeFBDInfo *info;
	info = (TypeFBDInfo*) malloc(tree->size*sizeof(TypeFBDInfo));
{
int f;
for(f=0; f<fos->size; f++)
	if(isnan(fos->fossilList[f].time)) {
		fprintFossilFeature(stdout, fos, tree->name, tree->size);
		error("Problem Fossil time %d is NaN\n", f);
	}
}
	for(n=0; n<tree->size; n++) {
if(tree->node[n].child == NOSUCH) {
	if(tree->time[n] > model->startTime[model->size]) {
		warning("node %d time %le > %le (set to %le)\n", n, tree->time[n], model->startTime[model->size], model->startTime[model->size]);
		tree->time[n] = model->startTime[model->size];
	}
	if(tree->time[n] < model->startTime[0]) {
		warning("node %d (%d) time %le < %le (set to %le)\n", n, tree->node[n].child, tree->time[n], model->startTime[0], model->startTime[0]);
		tree->time[n] = model->startTime[model->size];
	}
} else
	tree->time[n] = NO_TIME;

		if(tree->node[n].child != NOSUCH && tree->time[n] != NO_TIME) {
			warning("Internal node %d has time %.2lf (time discarded)\n", n, tree->time[n]);
			tree->time[n] = NO_TIME;
		}
		info[n].min = NO_TIME;
		info[n].max = NO_TIME;
	}
	if(constraint != NULL)
		for(i=0; constraint[i].time != NO_TIME; i++) {
			if(constraint[i].type == FBD_negative_constraint_type)
				info[constraint[i].node].max = constraint[i].time;
			if(constraint[i].type == FBD_positive_constraint_type)
				info[constraint[i].node].min = constraint[i].time;
		}
//printf("\nTree complet min time %lf\n", tree->minTime);
//printTreeDebug(stdout, tree->root, tree, NULL);
//{
	//int n;
	//for(n=0; n<tree->size; n++)
//if(tree->node[n].child != NOSUCH && tree->time[n] != NO_TIME && tree->time[n] > 0)
	//exit(0);
//}

	splitTreeFossilFBD(tree, info, fos, &treeList, &size);
	free((void*)info);
	stop = (double*) malloc((3*tree->size+fos->size+model->size+2)*sizeof(double));
	for(l=0; l<size; l++) {
{int n;
for(n=0; n<treeList[l]->size; n++)
	if(treeList[l]->node[n].child == NOSUCH && treeList[l]->time[n]<treeList[l]->minTime) {
		printf("Problem tree %d node %d time %lf %lf\n", l, n, treeList[l]->time[n], treeList[l]->minTime);
	FILE *ft;
	if((ft = fopen("Problem_Basic_TreeXX.newick", "w"))) {
		if(treeList[l]->name != NULL)
			free((void*)treeList[l]->name);
		treeList[l]->name = nameAll("Tip_", "Int_",treeList[l]);
		if(treeList[l]->comment != NULL)
			free((void*)treeList[l]->comment);
		treeList[l]->comment = NULL;
		fprintSubtreeNewick(ft, treeList[l]->root, treeList[l]);
		fprintf(ft, "\n");
		fclose(ft);
	}
	printf("min Time %lf\n",treeList[l]->minTime);
	error("Problem tree %d\n", l);
	exit(0);
}
}	
//printf("\nTree %d\n", l);
//printTreeDebug(stdout, treeList[l]->root, treeList[l], NULL);
		if(treeList[l]->size == 0) {
			res += logProbSplitted(treeList[l], (TypeFBDInfo*)treeList[l]->info, NULL, 0, &coeff);
		} else {
			int j, k, nstop;
			double timeMin, timeMax;
			timeMin = treeList[l]->minTime;
			timeMax = timeMin;
			nstop=0;
			stop[nstop++] = timeMin;
			for(n=0; n<treeList[l]->size; n++) { /*add all the times from tree[l] to the stop*/
				if(treeList[l]->time[n] != NO_TIME) {
					stop[nstop++] = treeList[l]->time[n];
//printf("stop tree time %d = %lf\n", n, treeList[l]->time[n]);
					if(treeList[l]->time[n] > timeMax)
						timeMax = treeList[l]->time[n];
				}
//printf("stop tree time %d min max %lf\n", l, treeList[l]->time[n);
				if(((TypeFBDInfo*)treeList[l]->info)[n].min != NO_TIME)
					stop[nstop++] = ((TypeFBDInfo*)treeList[l]->info)[n].min;
				if(((TypeFBDInfo*)treeList[l]->info)[n].max != NO_TIME) {
					stop[nstop++] = ((TypeFBDInfo*)treeList[l]->info)[n].max;
					if(((TypeFBDInfo*)treeList[l]->info)[n].max > timeMax)
						timeMax = ((TypeFBDInfo*)treeList[l]->info)[n].max;
				}
			}
			for(i=0; i<=model->size; i++) /*add the times of the model between timeMin and timeMax to the stops*/
				if(model->startTime[i]>timeMin && model->startTime[i]<timeMax)
					stop[nstop++] = model->startTime[i];
//for(i=0; i<nstop; i++)
	//printf("stop %d = %lf\n", i, stop[i]);
//printf("final number of stops %d\n", nstop);

			if(nstop>0) { /*sort and suppress doublons among the stops*/
				double tmp;
				qsort(stop, nstop, sizeof(double), compareDouble);
				tmp = stop[0];
				j = 1;
				for(k=1; k<nstop; k++)
					if(stop[k] != tmp) {
						stop[j++] = stop[k];
						tmp = stop[k];
					}
				nstop = j;
			}
//for(i=0; i<nstop; i++)
	//printf("stop %d = %lf\n", i, stop[i]);
//printf("final number of stops %d\n", nstop);
//printPiecewiseModel(stdout, model);
//exit(0);
			fillFBDInfo(treeList[l]->root, treeList[l]->minTime, treeList[l], (TypeFBDInfo*)treeList[l]->info); /*Fill min and max times and number of tips of treeList[l] if not empty*/
//{
		//if(treeList[l]->name != NULL)
			//free((void*)treeList[l]->name);
		//treeList[l]->name = nameAll("Tip_", "Int_",treeList[l]);
		//if(treeList[l]->comment != NULL)
			//free((void*)treeList[l]->comment);
		//treeList[l]->comment = NULL;
		//fprintSubtreeNewick(stdout, treeList[l]->root, treeList[l]);
		//fprintf(stdout, "\n");
//}
if(isinf(logProbSplitted(treeList[l], (TypeFBDInfo*)treeList[l]->info, stop, nstop, &coeff)) || isnan(logProbSplitted(treeList[l], (TypeFBDInfo*)treeList[l]->info, stop, nstop, &coeff))) {
	FILE *ft;
	if((ft = fopen("Problem_Basic_Tree.newick", "w"))) {
		if(treeList[l]->name != NULL)
			free((void*)treeList[l]->name);
		treeList[l]->name = nameAll("Tip_", "Int_",treeList[l]);
		if(treeList[l]->comment != NULL)
			free((void*)treeList[l]->comment);
		treeList[l]->comment = NULL;
		fprintSubtreeNewick(ft, treeList[l]->root, treeList[l]);
		fprintf(ft, "\n\n");
		fprintSubtreeNewick(ft, tree->root, tree);
		fprintf(ft, "\n\n");
		fclose(ft);
	}
for(i=0; i<nstop; i++)
	printf("stop %d = %lf\n", i, stop[i]);
printf("final number of stops %d\n", nstop);
	printf("min Time %lf\n",treeList[l]->minTime);
	error("Problem tree %d %le\n", l, logProbSplittedX(treeList[l], (TypeFBDInfo*)treeList[l]->info, stop, nstop, &coeff));
}
			res += logProbSplitted(treeList[l], (TypeFBDInfo*)treeList[l]->info, stop, nstop, &coeff);
		}
		if(treeList[l]->info!=NULL)
			free((void*)treeList[l]->info);
		freeTree(treeList[l]);
	}
	free((void*)stop);
	free((void*)treeList);
	free((void*)coeff.row);
	return res;
}





double logLikelihoodTreeFossilX(TypeTree *tree, TypeFossilFeature *fos, TypeFBDTimeConstraint *constraint, TypePiecewiseModelParam *model) {
	TypeTree **treeList;
	int i, n, l, size;
	double *stop, res = 0.;
	TypeFBDCoeff coeff = getCoeff(model);
	TypeFBDInfo *info;
	info = (TypeFBDInfo*) malloc(tree->size*sizeof(TypeFBDInfo));
{
int f;
for(f=0; f<fos->size; f++)
	if(isnan(fos->fossilList[f].time))
		error("Problem F %d\n", f);
}
	for(n=0; n<tree->size; n++) {
if(tree->node[n].child == NOSUCH) {
	if(tree->time[n] > model->startTime[model->size]) {
		warning("node %d time %le > %le (set to %le)\n", n, tree->time[n], model->startTime[model->size], model->startTime[model->size]);
		tree->time[n] = model->startTime[model->size];
	}
	if(tree->time[n] < model->startTime[0]) {
		warning("node %d (%d) time %le < %le (set to %le)\n", n, tree->node[n].child, tree->time[n], model->startTime[0], model->startTime[0]);
		tree->time[n] = model->startTime[model->size];
	}
} else
	tree->time[n] = NO_TIME;

		if(tree->node[n].child != NOSUCH && tree->time[n] != NO_TIME) {
			warning("Internal node %d has time %.2lf (time discarded)\n", n, tree->time[n]);
			tree->time[n] = NO_TIME;
		}
		info[n].min = NO_TIME;
		info[n].max = NO_TIME;
	}
	if(constraint != NULL)
		for(i=0; constraint[i].time != NO_TIME; i++) {
			if(constraint[i].type == FBD_negative_constraint_type)
				info[constraint[i].node].max = constraint[i].time;
			if(constraint[i].type == FBD_positive_constraint_type)
				info[constraint[i].node].min = constraint[i].time;
		}
//printf("\nTree complet min time %lf\n", tree->minTime);
//printTreeDebug(stdout, tree->root, tree, NULL);
//{
	//int n;
	//for(n=0; n<tree->size; n++)
//if(tree->node[n].child != NOSUCH && tree->time[n] != NO_TIME && tree->time[n] > 0)
	//exit(0);
//}

	splitTreeFossilFBD(tree, info, fos, &treeList, &size);
	free((void*)info);
	stop = (double*) malloc((3*tree->size+fos->size+model->size+2)*sizeof(double));
	for(l=0; l<size; l++) {
{int n;
for(n=0; n<treeList[l]->size; n++)
	if(treeList[l]->node[n].child == NOSUCH && treeList[l]->time[n]<treeList[l]->minTime) {
		printf("Problem tree %d node %d time %lf %lf\n", l, n, treeList[l]->time[n], treeList[l]->minTime);
	FILE *ft;
	if((ft = fopen("Problem_Basic_TreeXX.newick", "w"))) {
		if(treeList[l]->name != NULL)
			free((void*)treeList[l]->name);
		treeList[l]->name = nameAll("Tip_", "Int_",treeList[l]);
		if(treeList[l]->comment != NULL)
			free((void*)treeList[l]->comment);
		treeList[l]->comment = NULL;
		fprintSubtreeNewick(ft, treeList[l]->root, treeList[l]);
		fprintf(ft, "\n");
		fclose(ft);
	}
	printf("min Time %lf\n",treeList[l]->minTime);
	error("Problem tree %d\n", l);
	exit(0);
}
}	
//printf("\nTree %d\n", l);
//printTreeDebug(stdout, treeList[l]->root, treeList[l], NULL);
		if(treeList[l]->size == 0) {
			res += logProbSplitted(treeList[l], (TypeFBDInfo*)treeList[l]->info, NULL, 0, &coeff);
		} else {
			int j, k, nstop;
			double timeMin, timeMax;
			timeMin = treeList[l]->minTime;
			timeMax = timeMin;
			nstop=0;
			stop[nstop++] = timeMin;
			for(n=0; n<treeList[l]->size; n++) { /*add all the times from tree[l] to the stop*/
				if(treeList[l]->time[n] != NO_TIME) {
					stop[nstop++] = treeList[l]->time[n];
//printf("stop tree time %d = %lf\n", n, treeList[l]->time[n]);
					if(treeList[l]->time[n] > timeMax)
						timeMax = treeList[l]->time[n];
				}
//printf("stop tree time %d min max %lf\n", l, treeList[l]->time[n);
				if(((TypeFBDInfo*)treeList[l]->info)[n].min != NO_TIME)
					stop[nstop++] = ((TypeFBDInfo*)treeList[l]->info)[n].min;
				if(((TypeFBDInfo*)treeList[l]->info)[n].max != NO_TIME) {
					stop[nstop++] = ((TypeFBDInfo*)treeList[l]->info)[n].max;
					if(((TypeFBDInfo*)treeList[l]->info)[n].max > timeMax)
						timeMax = ((TypeFBDInfo*)treeList[l]->info)[n].max;
				}
			}
			for(i=0; i<=model->size; i++) /*add the times of the model between timeMin and timeMax to the stops*/
				if(model->startTime[i]>timeMin && model->startTime[i]<timeMax)
					stop[nstop++] = model->startTime[i];
//for(i=0; i<nstop; i++)
	//printf("stop %d = %lf\n", i, stop[i]);
//printf("final number of stops %d\n", nstop);

			if(nstop>0) { /*sort and suppress doublons among the stops*/
				double tmp;
				qsort(stop, nstop, sizeof(double), compareDouble);
				tmp = stop[0];
				j = 1;
				for(k=1; k<nstop; k++)
					if(stop[k] != tmp) {
						stop[j++] = stop[k];
						tmp = stop[k];
					}
				nstop = j;
			}
//for(i=0; i<nstop; i++)
	//printf("stop %d = %lf\n", i, stop[i]);
//printf("final number of stops %d\n", nstop);
//printPiecewiseModel(stdout, model);
//exit(0);
			fillFBDInfo(treeList[l]->root, treeList[l]->minTime, treeList[l], (TypeFBDInfo*)treeList[l]->info); /*Fill min and max times and number of tips of treeList[l] if not empty*/
//{
		//if(treeList[l]->name != NULL)
			//free((void*)treeList[l]->name);
		//treeList[l]->name = nameAll("Tip_", "Int_",treeList[l]);
		//if(treeList[l]->comment != NULL)
			//free((void*)treeList[l]->comment);
		//treeList[l]->comment = NULL;
		//fprintSubtreeNewick(stdout, treeList[l]->root, treeList[l]);
		//fprintf(stdout, "\n");
//}
if(isinf(logProbSplitted(treeList[l], (TypeFBDInfo*)treeList[l]->info, stop, nstop, &coeff))) {
	FILE *ft;
	if((ft = fopen("Problem_Basic_Tree.newick", "w"))) {
		if(treeList[l]->name != NULL)
			free((void*)treeList[l]->name);
		treeList[l]->name = nameAll("Tip_", "Int_",treeList[l]);
		if(treeList[l]->comment != NULL)
			free((void*)treeList[l]->comment);
		treeList[l]->comment = NULL;
		fprintSubtreeNewick(ft, treeList[l]->root, treeList[l]);
		fprintf(ft, "\n\n");
		fprintSubtreeNewick(ft, tree->root, tree);
		fprintf(ft, "\n\n");
		fclose(ft);
	}
for(i=0; i<nstop; i++)
	printf("stop %d = %lf\n", i, stop[i]);
printf("XX final number of stops %d\n", nstop);
	printf("min Time %lf\n",treeList[l]->minTime);
	error("Problem tree %d\n", l);
}
printf("tree %d %le mintime %lf\n", l, logProbSplitted(treeList[l], (TypeFBDInfo*)treeList[l]->info, stop, nstop, &coeff), treeList[l]->minTime);
//printTreeDebug(stdout, treeList[l]->root, treeList[l], NULL);
fprintTreeNewick(stdout, treeList[l]);
printf("\n");
			res += logProbSplitted(treeList[l], (TypeFBDInfo*)treeList[l]->info, stop, nstop, &coeff);
		}
		if(treeList[l]->info!=NULL)
			free((void*)treeList[l]->info);
		freeTree(treeList[l]);
	}
	free((void*)stop);
	free((void*)treeList);
	free((void*)coeff.row);
	return res;
}
#define PRECISION_LOG_LIKE -5*log(10.)

int fillBasicTree(TypeTree *basic, int n, TypeTree *tree, TypeFBDInfo *info, TypeFossilTab *ftab, int toCompute, int *newIndex) {
	if(n == toCompute)
		*newIndex = basic->size;
	int curInd = basic->size++;
	if(ftab[n].size>0) {
		basic->time[curInd] = ftab[n].time[0];
		basic->node[curInd].child = NOSUCH;
		((TypeFBDInfo*)basic->info)[curInd].type = FBD_fossil_type;
		((TypeFBDInfo*)basic->info)[curInd].min = NO_TIME;
		((TypeFBDInfo*)basic->info)[curInd].max = NO_TIME;
	} else {
		basic->time[curInd] = tree->time[n];
		if(tree->node[n].child != NOSUCH) {
			((TypeFBDInfo*)basic->info)[curInd].type = FBD_none_type;
			((TypeFBDInfo*)basic->info)[curInd].min = info[n].min;
			((TypeFBDInfo*)basic->info)[curInd].max = info[n].max;
			basic->node[curInd].child = fillBasicTree(basic, tree->node[n].child, tree, info, ftab, toCompute, newIndex);
			basic->node[basic->node[curInd].child].sibling = fillBasicTree(basic, tree->node[tree->node[n].child].sibling, tree, info, ftab, toCompute, newIndex);
			basic->node[basic->node[basic->node[curInd].child].sibling].sibling = NOSUCH;
		} else {
			((TypeFBDInfo*)basic->info)[curInd].min = NO_TIME;
			((TypeFBDInfo*)basic->info)[curInd].max = NO_TIME;
			if(tree->time[n] == tree->maxTime) {
				((TypeFBDInfo*)basic->info)[curInd].type = FBD_end_type;
			} else {
				((TypeFBDInfo*)basic->info)[curInd].type = FBD_stop_type;
			}
			basic->node[curInd].child = NOSUCH;
		}
	}
	return curInd;
}

int fillSplitTreeFossilFBD(TypeTree *tcur, int n, TypeTree *tree, TypeFBDInfo *info, TypeFossilTab *ftab, TypeTree **treeList, int *size) {
	int curInd = tcur->size++;
	if(ftab[n].size>0) {
		tcur->time[curInd] = ftab[n].time[0];
if(tcur->time[curInd] <= tcur->minTime) {
	printf("Problem A node %d (%s) %lf %lf\n", n, tree->name[n], tcur->time[curInd], tcur->minTime);
}
		tcur->node[curInd].child = NOSUCH;
		((TypeFBDInfo*)tcur->info)[curInd].type = FBD_fossil_type;
		((TypeFBDInfo*)tcur->info)[curInd].min = NO_TIME;
		((TypeFBDInfo*)tcur->info)[curInd].max = NO_TIME;
		if(ftab[n].size == 1 && tree->node[n].child == NOSUCH && (tree->time[n] == ftab[n].time[0] || tree->time[n]==NO_TIME)) { /*empty tree*/
			treeList[(*size)] = newTree(0);
			treeList[(*size)]->root = 0;
			treeList[(*size)]->minTime = ftab[n].time[0];
			treeList[(*size)]->maxTime = tree->maxTime;
			(*size)++;
			ftab[n].size--;
			ftab[n].time++;
		} else { /*tree not empty*/
			treeList[(*size)] = newTree(tree->size);
			treeList[(*size)]->root = 0;
			treeList[(*size)]->minTime = tcur->time[curInd];
			treeList[(*size)]->maxTime = tree->maxTime;
			treeList[(*size)]->info = (void*) malloc(tree->size*sizeof(TypeFBDInfo));
			((TypeFBDInfo*)treeList[(*size)]->info)[0].type = FBD_none_type;
			((TypeFBDInfo*)treeList[(*size)]->info)[0].min = NO_TIME;
			((TypeFBDInfo*)treeList[(*size)]->info)[0].max = NO_TIME;
			(*size)++;
			ftab[n].size--;
			ftab[n].time++;
			fillSplitTreeFossilFBD(treeList[(*size)-1], n, tree, info, ftab, treeList, size);
		}
	} else {
		tcur->time[curInd] = tree->time[n];
if(tcur->time[curInd]!=NO_TIME && tcur->time[curInd] <= tcur->minTime) {
	printf("Problem B node %d (%s) %lf %lf\n\n\n", n, tree->name[n], tcur->time[curInd], tcur->minTime);
	exit(0);
}
		if(tree->node[n].child != NOSUCH) {
			int c, prec;
			((TypeFBDInfo*)tcur->info)[curInd].type = FBD_none_type;
			((TypeFBDInfo*)tcur->info)[curInd].min = info[n].min;
			((TypeFBDInfo*)tcur->info)[curInd].max = info[n].max;
			tcur->node[curInd].child = fillSplitTreeFossilFBD(tcur, tree->node[n].child, tree, info, ftab, treeList, size);
			prec = tcur->node[curInd].child;
			for(c=tree->node[tree->node[n].child].sibling; c!=NOSUCH; c=tree->node[c].sibling) {
				tcur->node[prec].sibling = fillSplitTreeFossilFBD(tcur, c, tree, info, ftab, treeList, size);
				prec = tcur->node[prec].sibling;
			}
			tcur->node[prec].sibling = NOSUCH;
		} else {
			if(tree->time[n] == NO_TIME)
				warning("Warning in fillSplitTreeFossilFBD function tip %d (%s) has no time\n", n, tree->name[n]);
			((TypeFBDInfo*)tcur->info)[curInd].min = NO_TIME;
			((TypeFBDInfo*)tcur->info)[curInd].max = NO_TIME;
			if(tree->time[n] == tree->maxTime) {
				((TypeFBDInfo*)tcur->info)[curInd].type = FBD_end_type;
			} else {
				((TypeFBDInfo*)tcur->info)[curInd].type = FBD_stop_type;
			}
			tcur->node[curInd].child = NOSUCH;
		}
	}
	return curInd;
}

/*split the tree "tree" each time a fossill occurs*/
void splitTreeFossilFBD(TypeTree *tree, TypeFBDInfo *info, TypeFossilFeature *fos, TypeTree ***treeList, int *size) {
    int i;
    TypeFossilTab *ftab, *ftmp;	
{
	int f;
for(f=0; f<fos->size; f++)
	if(isnan(fos->fossilList[f].time) || fos->fossilList[f].time>0)
		warning("Problem D %d %lf\n", f, fos->fossilList[f].time);
}
{
int n;
	for(n=0; n<tree->size; n++)
if(tree->node[n].child != NOSUCH && tree->time[n] != NO_TIME && tree->time[n] > 0)
	exit(0);
}
//printf("\n");
//fprintFossilFeature(stdout, fos, tree->name, tree->size);
//printf("\n");
    ftab = listToFossilTab(fos, tree->size);
    ftmp = (TypeFossilTab*) malloc(tree->size*sizeof(TypeFossilTab));
    for(i=0; i<tree->size; i++) {
        ftmp[i].size = ftab[i].size;
        ftmp[i].time = ftab[i].time;
    }
    *treeList = (TypeTree **) malloc((fos->size+1)*sizeof(TypeTree*));
    (*treeList)[0] = newTree(tree->size);
    (*treeList)[0]->root = 0;
    (*treeList)[0]->node[0].sibling = NOSUCH;
    (*treeList)[0]->minTime = tree->minTime;
//printf("node %d new tree %lf (total size %d)\n", tree->root, tree->minTime, (fos->size+1));
    //(*treeList)[0]->minTimeInt.inf = tree->minTimeInt.inf;
    //(*treeList)[0]->minTimeInt.sup = tree->minTimeInt.sup;
    (*treeList)[0]->maxTime = tree->maxTime;
    (*treeList)[0]->info = (void*) malloc(tree->size*sizeof(TypeFBDInfo));
	((TypeFBDInfo*)(*treeList)[0]->info)[0].type = info[tree->root].type;
	((TypeFBDInfo*)(*treeList)[0]->info)[0].min = info[tree->root].min;
	((TypeFBDInfo*)(*treeList)[0]->info)[0].max = info[tree->root].max;
	*size = 1;
    fillSplitTreeFossilFBD((*treeList)[0], tree->root, tree, info, ftmp, *treeList, size);
    for(i=0; i<tree->size; i++)
        if(ftab[i].time != NULL)
            free((void*)ftab[i].time);
    free((void*)ftab);
    free((void*)ftmp);
    for(i=0; i<*size; i++)
        reallocTree((*treeList)[i]->size, (*treeList)[i]);
}


int fillBasicTreeFBD(TypeTree *basic, int n, TypeTree *tree, TypeFossilTab *ftab, int toCompute, int *newIndex) {
	if(n == toCompute)
		*newIndex = basic->size;
	int curInd = basic->size++;
	if(ftab[n].size>0) {
		basic->time[curInd] = ftab[n].time[0];
		basic->node[curInd].child = NOSUCH;
		((TypeNodeStatus*)basic->info)[curInd]=FBD_fossil_type;
	} else {
		basic->time[curInd] = tree->time[n];
		if(tree->node[n].child != NOSUCH) {
			((TypeNodeStatus*)basic->info)[curInd] = FBD_none_type;
			basic->node[curInd].child = fillBasicTreeFBD(basic, tree->node[n].child, tree, ftab, toCompute, newIndex);
			basic->node[basic->node[curInd].child].sibling = fillBasicTreeFBD(basic, tree->node[tree->node[n].child].sibling, tree, ftab, toCompute, newIndex);
			basic->node[basic->node[basic->node[curInd].child].sibling].sibling = NOSUCH;
		} else {
			if(tree->time[n] == tree->maxTime)
				((TypeNodeStatus*)basic->info)[curInd]=FBD_end_type;
			else
				((TypeNodeStatus*)basic->info)[curInd]=FBD_stop_type;
			basic->node[curInd].child = NOSUCH;
		}
	}
	return curInd;
}

TypeTree *getBasicTreeDivergenceFBD(int n, TypeTree *tree, TypeFossilTab *ftab, int *newIndex) {
	TypeTree *basic;
	int root;
	if(tree->parent == NULL)
		setParent(tree);
	for(root=n; root!=tree->root && ftab[root].size==0; root=tree->parent[root])
		;
	basic = newTree(tree->size);
	basic->root = 0;
	if(ftab[root].size>0) {
		basic->minTime = ftab[root].time[ftab[root].size-1];
		//basic->minTimeInt.inf = ftab[root].time[ftab[root].size-1];
		//basic->minTimeInt.sup = ftab[root].time[ftab[root].size-1];
	} else {
		basic->minTime = tree->minTime;
		//basic->minTimeInt.inf = tree->minTimeInt.inf;
		//basic->minTimeInt.sup = tree->minTimeInt.sup;
	}
	basic->maxTime = tree->maxTime;
	basic->info = (TypeNodeStatus*) malloc(tree->size*sizeof(TypeNodeStatus));
	ftab[root].size = 0;
	fillBasicTreeFBD(basic, root, tree, ftab, n, newIndex);
	reallocTree(basic->size, basic);
	basic->info = (TypeNodeStatus*) realloc((void*)basic->info, basic->size*sizeof(TypeNodeStatus));
	return basic;
}
