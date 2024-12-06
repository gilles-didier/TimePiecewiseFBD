#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#include "Utils.h"
#include "MyR.h"
#include "MyRandom.h"
#include "Tree.h"
#include "Fossil.h"
#include "FossilInt.h"
#include "PiecewiseModel.h"
#include "FBDDensity.h"
//#include "SimulPiecewise.h"
#include "MCMCImportanceSamplingFossilInt.h"
#include "MinimizeNLOpt.h"



//./xfbd -u 40 -p -20. 0.5 0.1 0.5 1. -p -10 0.05 0.2 0.1 1. -s 5000 -b 10000 -g 50 -f 1. -a 0.33 0.33 -w 0.5 0.5 0.5 -i 5. 5. 5.  -r 0000 
//valgrind ./dfbd -s 200 ~/Dropbox/FossilPiecewise/dataMS/TreeEupely2.phy ~/Dropbox/FossilPiecewise/dataMS/CotylosauriaAgesNew.csv ~/Dropbox/FossilPiecewise/modelMS/schemeM0.txt

#define NAME_CODA "mcmc_sample"

#define MAX_PARAM 10

#define STRING_SIZE 300
#define HELP_MESSAGE "\nusage: sample [options] [<output file>]\n\nsample simulates random trees and fossils finds and saves them in Newick format\n\nOptions:\n\t-h : display help\n\t-b <birth>\t: set birth rate\n\t-d <death>\t: set death rate\n\t-f <fossil>\t: set fossil find rate\n\t-m <min>\t: set minimum number of contemporary species of a simulation to be considered\n\t-M <size>\t: set maximum size of a simulation to be considered\n\t-i <niter>\t: set the number of simulations\n\t-t <time> : the end time of the diversification (start is always 0)\n"


//./dfbd -s 20000 ../../data/TreeEupely40.phy ../../data/CotylosauriaAgesNew.csv model/schemeExtKR.txt 

int main(int argc, char **argv) {	
	char *inputFileNameTree, *inputFileNameFossil, *inputFileNameScheme, outputName[STRING_SIZE], outputFileName[STRING_SIZE], *outputPrefix, option[256];
	FILE *fit, *fif, *fis, *fo;
	int i, niter = 30, nSamp = 5000, nBurn = 10000, nGap = 10;
	//double al = 0.75, probTime = 0.05, probSpe = 0.24, probExt = 0.24, probFos = 0.24, propParam = 0.2, stopTime = DBL_MAX;
	//TypeModelParam windSize = {.birth=0.05, .death = 0.05, .fossil = 0.05, .sampling = 0.5}, init = {.birth=0.5, .death = 0.5, .fossil = 0.1, .sampling = 1.};
	double al = 0.99, probTime = 0.1, probSpe = 0.2, probExt = 0.2, probFos = 0.2, propParam = 0.9, stopTime = DBL_MAX;
	TypeModelParam windSize = {.birth=10., .death = 10., .fossil = 5., .sampling = 2.}, init = {.birth=0.5, .death = 0.5, .fossil = 0.1, .sampling = 1.};
	TypeNLOptOption nloptOption;

	nloptOption.infSpe = 0.;
	nloptOption.supSpe = 1.;
	nloptOption.infExt = 0.;
	nloptOption.supExt = 1.;
	nloptOption.infFos = 0.;
	nloptOption.supFos = 1.;
	nloptOption.trials = 10;
	nloptOption.tolOptim = 0.0001;
	nloptOption.maxIter = 10000;
	
	for(i=0; i<256; i++)
		option[i] = 0;
	for(i=1; i<argc && *(argv[i]) == '-'; i++) {
		int j;
		for(j=1; argv[i][j] != '\0'; j++)
			option[argv[i][j]] = 1;
		if(option['a']) {
			option['a'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%le", &probSpe) == 1)
				i++;
			else
				error("3 values are expected after -a");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &probExt) == 1)
				i++;
			else
				error("3 values are expected after -a");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &probFos) == 1)
				i++;
			else
				error("3 values are expected after -a");
		}
		if(option['w']) {
			option['w'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%le", &windSize.birth) == 1)
				i++;
			else
				error("4 values are expected after -w");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &windSize.death) == 1)
				i++;
			else
				error("4 values are expected after -w");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &windSize.fossil) == 1)
				i++;
			else
				error("4 values are expected after -w");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &windSize.sampling) == 1)
				i++;
			else
				error("4 values are expected after -w");
		}
		if(option['i']) {
			option['i'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%le", &init.birth) == 1)
				i++;
			else
				error("4 values are expected after -w");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &init.death) == 1)
				i++;
			else
				error("4 values are expected after -w");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &init.fossil) == 1)
				i++;
			else
				error("4 values are expected after -w");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &init.sampling) == 1)
				i++;
			else
				error("4 values are expected after -w");
		}
		if(option['s']) {
			option['s'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &nSamp) == 1)
				i++;
			else
				error("a number is expected after -s");
		}
		if(option['b']) {
			option['b'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &nBurn) == 1)
				i++;
			else
				error("a number is expected after -b");
		}
		if(option['g']) {
			option['g'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &nGap) == 1)
				i++;
			else
				error("a number is expected after -b");
		}
		if(option['u']) {
			option['u'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &niter) == 1)
				i++;
		}
		if(option['S']) {
			option['S'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &stopTime) == 1)
				i++;
		}
		if(option['h']) {
			printf("%s\n", HELP_MESSAGE);
			return 0;
		}
	}
	if(i<argc) {
		inputFileNameTree = argv[i++];
	} else {
		fprintf(stderr, "Please provide the name of a file containing a phylogenetic tree in Newick format\n");
		exit(1);
	}
	if(i<argc) {
		inputFileNameFossil = argv[i++];
	} else {
		fprintf(stderr, "Please provide the name of a file containing a fossil list\n");
		exit(1);
	}
	if(i<argc) {
		inputFileNameScheme = argv[i++];
	} else {
		fprintf(stderr, "Please provide the name of a file containing a  model scheme\n");
		exit(1);
	}
	strcpy(outputName, inputFileNameScheme);
	if((outputPrefix = strrchr(outputName, '.')) != NULL)
		outputPrefix[0] = '\0';
	//if((outputPrefix=strrchr(outputName, '/')) == NULL)
		//outputPrefix = outputName;
	//else
		//outputPrefix++;
	
	if((fit = fopen(inputFileNameTree, "r")) && (fif = fopen(inputFileNameFossil, "r")) && (fis = fopen(inputFileNameScheme, "r"))) {
		TypeTree **tree;
		TypeFossilIntFeature *fint;
		TypeFossilFeature *foss, *fsave;
		TypeSamplingScheme scheme;
		TypeSamplingCurrent sampCurrent;
		TypePiecewiseModelParam maxParam;
		int sizeTree, nParam, n;
		double logLike;
		void *rand_data;
		TypeEstimation estim;
		rand_data = my_rand_get_data();
		scheme = readSamplingScheme(fis);
        fclose(fis);
        fsave = NULL;
printSamplingScheme(stdout, scheme);
//exit(0);
		tree = readTrees(fit);
        fclose(fit);
		sizeTree = 0;
		for(i=0; tree[i]!=NULL; i++) {
			int n;
			toBinary(tree[i]);
			if(tree[i]->name!=NULL)
				for(n=0; n<tree[i]->size; n++)
					if(tree[i]->name[n]!=NULL)
						fixSpace(tree[i]->name[n]);
			for(n=0; n<tree[i]->size; n++)
				if(tree[i]->time[n] < scheme.param.startTime[scheme.param.size])
					tree[i]->time[n] = NO_TIME;
				else
					tree[i]->time[n] = scheme.param.startTime[scheme.param.size];
			tree[i]->minTime = scheme.param.startTime[0];
			tree[i]->maxTime = scheme.param.startTime[scheme.param.size];
			sizeTree++;
		}
		if(sizeTree == 0)
			error("No tree!\n");
		reindexTreesFromName(tree, sizeTree);
		fint = getFossilIntFeature(fif, tree[0]->name, tree[0]->size);
        fclose(fif);
		if(getMaxFossilIntTime(fint) > 0.)
			negateFossilInt(fint);
		fixStatus(tree[0], fint);
		for(n=0; n<tree[0]->size; n++)
			if(fint->status[n] == unknownNodeStatus)
				for(i=0; tree[i]!=NULL; i++)
					tree[i]->time[n] = fint->endTimeTable[n].inf;
		nParam = scheme.sizeParamTime+scheme.sizeParamBirth+scheme.sizeParamDeath+scheme.sizeParamFossil+scheme.sizeParamSampling;
		maxParam.size = scheme.param.size;
		maxParam.param = (TypeModelParam*) malloc(scheme.param.size*sizeof(TypeModelParam));
		maxParam.startTime = (double*) malloc((scheme.param.size+1)*sizeof(double));
		estim.param.size = scheme.param.size;
		estim.param.param = (TypeModelParam*) malloc(scheme.param.size*sizeof(TypeModelParam));
		estim.param.startTime = (double*) malloc((scheme.param.size+1)*sizeof(double));
		
//		maxParam.startTime = scheme.param.startTime;
		//logLike = MCMCFillMaxParamSampleSC(tree, sizeTree, fint, scheme, al, nSamp, propParam, &windSize, &init, probTime, probSpe,probExt, probFos, &maxParam, &foss, &sampCurrent, rand_data);
		//fprintFossilFeature(stdout, foss, tree[0]->name, tree[0]->size);
		//printf("log-Likelihood\t%le\nAIC\t%le\n", logLike, 2.*(((double)nParam)-logLike));
		//printPiecewiseModel(stdout, &(maxParam));
		
//exit(0);
int iter;
double logMax = -DBL_MAX;
for(iter = 0; iter <niter; iter++) {
		if(scheme.sizeParamTime>0) {
			int i;
			sampCurrent.currentTime = (double*) malloc(scheme.sizeParamTime*sizeof(double));
			for(i=0; i<scheme.sizeParamTime; i++) {
				//sampCurrent.currentTime[i] = scheme.param.startTime[scheme.indexTime[i]-1]+(scheme.param.startTime[scheme.indexTime[i]+1]-scheme.param.startTime[scheme.indexTime[i]-1])*my_rand_unif_unit(rand_data);
				//printf("%.2lf\t%.2lf\t%.2lf\n", scheme.param.startTime[scheme.indexTime[i]-1], scheme.param.startTime[scheme.indexTime[i]+1], sampCurrent.currentTime[i]);
				sampCurrent.currentTime[i] = scheme.param.startTime[scheme.indexTime[i]-1]+(scheme.param.startTime[scheme.indexTime[i]+1]-scheme.param.startTime[scheme.indexTime[i]-1])*my_rand_unif_unit(rand_data);
				//printf("%.2lf\t%.2lf\t%.2lf\n", scheme.param.startTime[scheme.indexTime[i]-1], scheme.param.startTime[scheme.indexTime[i]+1], sampCurrent.currentTime[i]);
//				sampCurrent.currentTime[i] = -260;
				printf("\n\ntime %.2lf\n", sampCurrent.currentTime[i]);
			}
		} else
			sampCurrent.currentTime = NULL;
		if(scheme.sizeParamBirth>0) {
			int i;
			sampCurrent.currentBirth = (double*) malloc(scheme.sizeParamBirth*sizeof(double));
			for(i=0; i<scheme.sizeParamBirth; i++)
				sampCurrent.currentBirth[i] = my_rand_unif_unit(rand_data)*init.birth+0.0001;
		} else
			sampCurrent.currentBirth = NULL;
		if(scheme.sizeParamDeath>0) {
			int i;
			sampCurrent.currentDeath = (double*) malloc(scheme.sizeParamDeath*sizeof(double));
			for(i=0; i<scheme.sizeParamDeath; i++)
				sampCurrent.currentDeath[i] = my_rand_unif_unit(rand_data)*init.death+0.0001;
		} else
			sampCurrent.currentDeath = NULL;
		if(scheme.sizeParamFossil>0) {
			int i;
			sampCurrent.currentFossil = (double*) malloc(scheme.sizeParamFossil*sizeof(double));
			for(i=0; i<scheme.sizeParamFossil; i++)
				sampCurrent.currentFossil[i] = my_rand_unif_unit(rand_data)*init.fossil+0.0001;
		} else
			sampCurrent.currentFossil = NULL;
		if(scheme.sizeParamSampling>0) {
			int i;
			sampCurrent.currentSampling = (double*) malloc(scheme.sizeParamSampling*sizeof(double));
			for(i=0; i<scheme.sizeParamSampling; i++)
				sampCurrent.currentSampling[i] = my_rand_unif_unit(rand_data)+0.0001;
		} else
			sampCurrent.currentSampling = NULL;
		foss = sampleFossilIntSpecial(fint, tree[0]->size);
//		fprintFossilFeature(stdout, foss, tree[0]->name, tree[0]->size);
		minimizePiecewiseParamFromSetTreeFossilInt(getLogDensitySum, tree, sizeTree, fint, foss, sampCurrent, scheme, &nloptOption, &estim);
		printf("estim %d log-Likelihood\t%le\nAIC\t%le\n", iter, estim.logLikelihood, 2.*(((double)nParam)-estim.logLikelihood));
//		printf("model size %d\n", estim.param.size);
		printPiecewiseModel(stdout, &(estim.param));
//		fprintFossilFeature(stdout, foss, tree[0]->name, tree[0]->size);
				printf("\n\n");

		//printf("log-Likelihood\t%le\nAIC\t%le\n", logLike, 2.*(((double)nParam)-logLike));
		//sprintf(outputFileName, "%s_ML_model.txt", outputName);
		//printf("outputFile %s\n", outputFileName);
		//if((fo = fopen(outputFileName, "w"))) {
			//fprintf(fo, "#%d iterations\n#log-likelihood max\t%le\n#%d parameters\n#AIC\t%le\n", nSamp, logLike, nParam, 2.*(((double)nParam+fint->sizeFossil)-logLike));
			//printPiecewiseModel(fo, &maxParam);
			//fclose(fo);
		//} else
			//error("Can't open %s\n", outputFileName);
			//fprintf(stdout, "#%d iterations\n#log-likelihood max\t%le\n#%d parameters\n#AIC\t%le\n", nSamp, estim.logLikelihood, nParam, 2.*(((double)nParam+fint->sizeFossil)-estim.logLikelihood));
			//printPiecewiseModel(stdout, &(estim.param));
			if(estim.logLikelihood>logMax) {
				int k;
				logMax = estim.logLikelihood;
				if(fsave != NULL)
					freeFossilFeature(fsave);
				fsave = cpyFossilFeature(foss, tree[0]->size);
				for(k=0; k<maxParam.size; k++)
					maxParam.param[k] = estim.param.param[k];
				for(k=0; k<=maxParam.size; k++)
					maxParam.startTime[k] = estim.param.startTime[k];
			}
			printf("log max %le\n", logMax);
			freeFossilFeature(foss);
}
	sprintf(outputFileName, "%s_ML_model.txt", outputName);
	printf("outputFile %s\n", outputFileName);
	if((fo = fopen(outputFileName, "w"))) {
		fprintf(fo, "#%d iterations\n#log-likelihood max\t%le\n#%d parameters\n#AIC\t%le\n", niter, logMax, nParam, 2.*(((double)nParam+fint->sizeFossil)-logMax));
		printPiecewiseModel(fo, &maxParam);
		fclose(fo);
	} else
		error("Can't open %s\n", outputFileName);
	sprintf(outputFileName, "%s_fossil_ages.txt", outputName);
	printf("outputFile %s\n", outputFileName);
	if((fo = fopen(outputFileName, "w")) && fsave != NULL) {
		fprintFossilFeature(fo, fsave, tree[0]->name, tree[0]->size);
		fclose(fo);
	} else
		error("Can't open %s\n", outputFileName);
		freeSamplingScheme(&scheme);
		free((void*)maxParam.param);
	}
	
	return 0;
}
