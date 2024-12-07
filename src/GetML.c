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
#include "MCMCImportanceSamplingFossilInt.h"
#include "MinimizeNLOpt.h"


#define STRING_SIZE 300
#define HELP_MESSAGE "\nusage: sample [options] [<output file>]\n\nsample simulates random trees and fossils finds and saves them in Newick format\n\nOptions:\n\t-h : display help\n\t-b <birth>\t: set birth rate\n\t-d <death>\t: set death rate\n\t-f <fossil>\t: set fossil find rate\n\t-m <min>\t: set minimum number of contemporary species of a simulation to be considered\n\t-M <size>\t: set maximum size of a simulation to be considered\n\t-i <niter>\t: set the number of simulations\n\t-t <time> : the end time of the diversification (start is always 0)\n"


int main(int argc, char **argv) {	
	char *inputFileNameTree, *inputFileNameFossil, *inputFileNameScheme, *outputName, outputFileName[STRING_SIZE+50], option[256];
	FILE *fit, *fif, *fis, *fo;
	int i, niter = 30;
	TypeModelParam init = {.birth=0.5, .death = 0.5, .fossil = 0.1, .sampling = 1.};
	TypeNLOptOption nloptOption;

	nloptOption.infSpe = 0.;
	nloptOption.supSpe = 1.;
	nloptOption.infExt = 0.;
	nloptOption.supExt = 1.;
	nloptOption.infFos = 0.;
	nloptOption.supFos = 1.;
	nloptOption.trials = 10;
	nloptOption.tolOptim = 0.00001;
	nloptOption.maxIter = 10000;
	
	for(i=0; i<256; i++)
		option[i] = 0;
	for(i=1; i<argc && *(argv[i]) == '-'; i++) {
		int j;
		for(j=1; argv[i][j] != '\0'; j++)
			option[argv[i][j]] = 1;
		if(option['o']) {
			FILE *fopt;
			option['o'] = 0;
			if((i+1)<argc) {
				printf("file nlopt %s\n", argv[i+1]);
				if((fopt = fopen(argv[++i], "r"))) {
					fscanNLoptOptionTag(fopt, &nloptOption);
					fclose(fopt);
				} else
					error("Can't open NLopt option file %s\n", argv[i]);
			} else
				error("File name missing after option -o\n");

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
		if(option['r']) {
			option['r'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &niter) == 1)
				i++;
		}
		if(option['h']) {
			printf("%s\n", HELP_MESSAGE);
			return 0;
		}
	}
	if(i<argc)
		inputFileNameTree = argv[i++];
	else
		error("Please provide the name of a file containing a phylogenetic tree in Newick format\n");
	if(i<argc)
		inputFileNameFossil = argv[i++];
	else
		error("Please provide the name of a file containing a fossil list\n");
	if(i<argc)
		inputFileNameScheme = argv[i++];
	else
		error("Please provide the name of a file containing a  model scheme\n");
	if(i<argc)
		outputName = argv[i++];
	else {
		char STMP[STRING_SIZE], *outputPrefix;
		strcpy(STMP, inputFileNameScheme);
		if((outputPrefix = strrchr(STMP, '.')) != NULL)
			outputPrefix[0] = '\0';
		outputName = STMP;
		//if((outputPrefix=strrchr(STMP, '/')) == NULL)
			//outputPrefix = STMP;
		//else
			//outputPrefix++;
	}
	
	if((fit = fopen(inputFileNameTree, "r")) && (fif = fopen(inputFileNameFossil, "r")) && (fis = fopen(inputFileNameScheme, "r"))) {
		TypeTree **tree;
		TypeFossilIntFeature *fint;
		TypeFossilFeature *foss, *fsave;
		TypeSamplingScheme scheme;
		TypeSamplingCurrent sampCurrent;
		TypePiecewiseModelParam maxParam;
		int sizeTree, n, iter;
		double logMax = -DBL_MAX;
		void *rand_data;
		TypeEstimation estim;
		rand_data = my_rand_get_data();
		scheme = readSamplingScheme(fis);
        fclose(fis);
        fsave = NULL;

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
		maxParam.size = scheme.param.size;
		maxParam.param = (TypeModelParam*) malloc(scheme.param.size*sizeof(TypeModelParam));
		maxParam.startTime = (double*) malloc((scheme.param.size+1)*sizeof(double));
		estim.param.size = scheme.param.size;
		estim.param.param = (TypeModelParam*) malloc(scheme.param.size*sizeof(TypeModelParam));
		estim.param.startTime = (double*) malloc((scheme.param.size+1)*sizeof(double));
		for(iter = 0; iter <niter; iter++) {
			if(scheme.sizeParamTime>0) {
				int i;
				sampCurrent.currentTime = (double*) malloc(scheme.sizeParamTime*sizeof(double));
				for(i=0; i<scheme.sizeParamTime; i++) {
					sampCurrent.currentTime[i] = scheme.param.startTime[scheme.indexTime[i]-1]+(scheme.param.startTime[scheme.indexTime[i]+1]-scheme.param.startTime[scheme.indexTime[i]-1])*my_rand_unif_unit(rand_data);
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
			minimizePiecewiseParamFromSetTreeFossilInt(getLogDensitySum, tree, sizeTree, fint, foss, sampCurrent, scheme, &nloptOption, &estim);
			printf("estim %d log-Likelihood\t%le\n", iter, estim.logLikelihood);
			printPiecewiseModel(stdout, &(estim.param));
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
			printf("log max %le\n\n\n", logMax);
			freeFossilFeature(foss);
		}
		snprintf(outputFileName, STRING_SIZE+49, "%s_ML_model.txt", outputName);
		printf("outputFile model: %s\n", outputFileName);
		if((fo = fopen(outputFileName, "w"))) {
			fprintf(fo, "#%d iterations\n#log-likelihood max\t%le\n", niter, logMax);
			printPiecewiseModel(fo, &maxParam);
			fclose(fo);
		} else
			error("Can't open %s\n", outputFileName);
		snprintf(outputFileName, STRING_SIZE+49, "%s_fossil_ages.txt", outputName);
		printf("outputFile fossils: %s\n", outputFileName);
		if((fo = fopen(outputFileName, "w")) && fsave != NULL) {
			fprintFossilFeature(fo, fsave, tree[0]->name, tree[0]->size);
			fclose(fo);
		} else
			error("Can't open %s\n", outputFileName);
		freeSamplingScheme(&scheme);
		free((void*)maxParam.param);
	} else
		error("Can't open %s, %s or %s\n", inputFileNameTree, inputFileNameFossil, inputFileNameScheme);
	return 0;
}
