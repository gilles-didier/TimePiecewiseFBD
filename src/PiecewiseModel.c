#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "MyR.h"
#include "PiecewiseModel.h"

/*return a piecewise model with a single piece and parameters param*/
TypePiecewiseModelParam simple2piecewise(TypeModelParam *param, double startTime, double endTime) {
	TypePiecewiseModelParam res;
	res.size = 1;
	res.startTime = (double*) malloc(2*sizeof(double));
	res.startTime[0] = startTime;
	res.startTime[1] = endTime;
	res.param = (TypeModelParam*) malloc(sizeof(TypeModelParam));
	res.param[0] = *param;
	return res;
}

/*return index i such that v in [param->startTime[i], param->startTime[i+1][*/
int getPieceIndex(double v, TypePiecewiseModelParam *param) {
	int a, b, c;
	if(v<param->startTime[0]) {
		error("Error in 'getPieceIndex': value %.2lf too small (start time = %.2lf)\n", v, param->startTime[0]);
		return -1;
	}
	if(v>param->startTime[param->size]) {
		warning("Error in 'getPieceIndex': value %.2lf too high (end time = %.2lf)\n", v, param->startTime[param->size]);
		return -1;
	}
	a = 0;
	b = param->size;
	while(b > a+1) {
		c = (a+b)/2;
		if(param->startTime[c] < v)
			a = c;
		else
			b = c;
	}
	return a;
}
			
/*just print the model param*/
void printPiecewiseModel(FILE *f, TypePiecewiseModelParam *param) {
	int i;
//	fprintf(f, "size %d\n", param->size);
	for(i=0; i<param->size; i++)
		fprintf(f, "%.17g\n%lf %lf %lf %lf\n", param->startTime[i], param->param[i].birth, param->param[i].death, param->param[i].fossil, param->param[i].sampling);
	fprintf(f, "%.17g\n", param->startTime[param->size]);
}

#define INC_SIZEBUF_SC 10
#define MAX_NAME_SIZE 100

char readSep(FILE *f, char c) {
	while(c == '#' || issepline(c)) {
		for(; c != EOF && issepline(c); c = fgetc(f))
			;
		if(c == '#')
			for(c = fgetc(f); c != EOF && c != '\n' && c != '\r'; c = fgetc(f))
				;
	}
	return c;
}

char readNumber(FILE *f, double *numb) {
	int i;
	char c;
	char tmp[MAX_NAME_SIZE+1];
	//for(c = fgetc(f); c != EOF && issepline(c); c = fgetc(f))
		//;
	c = readSep(f, fgetc(f));
	for(i=0; i<MAX_NAME_SIZE && c !=EOF && !issepline(c); i++) {
		tmp[i] = c;
		c = fgetc(f);
	}
	c = readSep(f, fgetc(f));
	//for(c = fgetc(f); c != EOF && issepline(c); c = fgetc(f))
		//;
	if(c != EOF)
		ungetc(c, f);
	if(i == MAX_NAME_SIZE) {
		tmp[i] = '\0';
		error("number too much long:\n'%s'", tmp);
	}
	tmp[i++] = '\0';
	if(i>1) {
		*numb = atof(tmp);
		//if(isdigit(tmp[0])) {
			//*numb = atof(tmp);
		//} else {
			//error("Reading %s in place of a number\n", tmp);
		//}
//printf("item '%s' scheme %d\n", tmp, *scheme);
	} else
		error("Something goes wrong while reading the model (l 74)\n");
	return c;
}

char readParamLine(FILE *f, double *start, TypeModelParam *param) {
	char c;
	c = readNumber(f, start);
	if(c != EOF) {
		c = readNumber(f, &(param->birth));
		if(c != EOF) {
			c = readNumber(f, &(param->death));
			if(c != EOF) {
				c = readNumber(f, &(param->fossil));
				if(c != EOF) {
					c = readNumber(f, &(param->sampling));
				} else
					error("Something goes wrong while reading the model A\n");
			} else
				error("Something goes wrong while reading the model B\n");
		} else
			error("Something goes wrong while reading the model C\n");
	}
	return c;
}

TypePiecewiseModelParam readPiecewiseModelParam(FILE *f) {
	char c;
	int sizeBuf;
	TypePiecewiseModelParam param;
	
	sizeBuf = INC_SIZEBUF_SC;
	param.size = 0;
	param.startTime = (double*) malloc((sizeBuf+1)*sizeof(double));
	param.param = (TypeModelParam*) malloc(sizeBuf*sizeof(TypeModelParam));
	//for(c = fgetc(f); c != EOF && c != '\n' && c != '\r'; c = fgetc(f))
		//;
    do {
		if(param.size >= sizeBuf) {
			sizeBuf += INC_SIZEBUF_SC;
			param.startTime = (double*) realloc((void*) param.startTime, (sizeBuf+1)*sizeof(double));
			param.param = (TypeModelParam*) realloc((void*) param.param, sizeBuf*sizeof(TypeModelParam));
		}
		c = readParamLine(f, &(param.startTime[param.size]), &(param.param[param.size]));
		if(c != EOF)
			param.size++;
	} while(c != EOF);
	param.startTime = (double*) realloc((void*) param.startTime, (param.size+1)*sizeof(double));
	param.param = (TypeModelParam*) realloc((void*) param.param, param.size*sizeof(TypeModelParam));
	return param;
}

void printSamplingScheme(FILE *f, TypeSamplingScheme sampScheme) {
	int i;
	//for(i=0; i<sampScheme.param.size; i++)
		//printf("sampling %d -> %d\n", i, sampScheme.schemeSampling[i]);

	for(i=0; i<sampScheme.param.size; i++) {
		if(sampScheme.schemeTime[i]>=0)
			fprintf(f, "t_%c\n", 'a'+sampScheme.schemeTime[i]);
		else
			fprintf(f, "%.17g\n", sampScheme.param.startTime[i]);
		if(sampScheme.schemeBirth[i]>=0)
			fprintf(f, "b_%c", 'a'+sampScheme.schemeBirth[i]);
		else
			fprintf(f, "%le", sampScheme.param.param[i].birth);
		if(sampScheme.schemeDeath[i]>=0)
			fprintf(f, "\td_%c", 'a'+sampScheme.schemeDeath[i]);
		else
			fprintf(f, "\t%le", sampScheme.param.param[i].death);
		if(sampScheme.schemeFossil[i]>=0)
			fprintf(f, "\tf_%c", 'a'+sampScheme.schemeFossil[i]);
		else
			fprintf(f, "\t%le", sampScheme.param.param[i].fossil);
		if(sampScheme.schemeSampling[i]>=0)
			fprintf(f, "\ts_%c", 'a'+sampScheme.schemeSampling[i]);
		else
			fprintf(f, "\t%le", sampScheme.param.param[i].sampling);
		fprintf(f, "\n");
	}
	if(sampScheme.schemeTime[sampScheme.param.size]>=0)
		fprintf(f, "t_%c\n", 'a'+sampScheme.schemeTime[sampScheme.param.size]);
	else
		fprintf(f, "%.17g\n", sampScheme.param.startTime[sampScheme.param.size]);
}

char setSamplingSchemeTime(FILE *f, int current, int *index, int *size, int *scheme, double *value) {
	int i;
	char c;
	char tmp[MAX_NAME_SIZE+1];
	//for(c = fgetc(f); c != EOF && issepline(c); c = fgetc(f))
		//;
	c = readSep(f, fgetc(f));
	if(c == '\'' || c == '"') {
		c = fgetc(f);
		for(i=0; i<MAX_NAME_SIZE && c != EOF && c != '\'' && c != '"'; i++) {
			tmp[i] = c;
			c = fgetc(f);
		}
		if(c == '\'' || c == '"')
		c = fgetc(f);
	} else {
		for(i=0; i<MAX_NAME_SIZE && c !=EOF && !issepline(c); i++) {
			tmp[i] = c;
			c = fgetc(f);
		}
	}
	if(c!=EOF)
		ungetc(c, f);
	if(i == MAX_NAME_SIZE) {
		tmp[i] = '\0';
		error("Name too much long:\n'%s'", tmp);
	}
	tmp[i++] = '\0';
	if(i>1) {
		if(isdigit(tmp[0]) || tmp[0] == '-' || tmp[0] == '.') {
			*scheme = -1;
			*value = atof(tmp);
		} else {
			index[(*size)] = current;
			*scheme = (*size)++;
			*value = DBL_MAX;
		}
	} else
		error("Something goes wrong while reading the sampling scheme (l 143)\n");
	c = readSep(f, fgetc(f));
	if(c != EOF)
		ungetc(c, f);
	return c;
}

char setSamplingScheme(FILE *f, TypeLexiTree *dict, int *size, int *scheme, double *value) {
	int i;
	char c;
	char tmp[MAX_NAME_SIZE+1];
	for(c = fgetc(f); c != EOF && issepline(c); c = fgetc(f))
		;
	if(c == '\'' || c == '"') {
		c = fgetc(f);
		for(i=0; i<MAX_NAME_SIZE && c != EOF && c != '\'' && c != '"'; i++) {
			tmp[i] = c;
			c = fgetc(f);
		}
		if(c == '\'' || c == '"')
		c = fgetc(f);
	} else {
		for(i=0; i<MAX_NAME_SIZE && c !=EOF && !issepline(c); i++) {
			tmp[i] = c;
			c = fgetc(f);
		}
	}
	if(c!=EOF)
		ungetc(c, f);
	if(i == MAX_NAME_SIZE) {
		tmp[i] = '\0';
		error("Name too much long:\n'%s'", tmp);
	}
	tmp[i++] = '\0';
	if(i>1) {
		if(isdigit(tmp[0])) {
			*scheme = -1;
			*value = atof(tmp);
		} else {
			int index = findWordLexi(tmp, dict);
			if(index < 0) {
				addWordLexi(tmp, *size, dict);
				index = *size;
				(*size)++;
			}
			*scheme = index;
			*value = DBL_MAX;
		}
	} else
		error("Something goes wrong while reading the sampling scheme (in function setSamplingScheme)\n");
	return c;
}

TypeSamplingScheme readSamplingScheme(FILE *f) {
	char c;
	int sizeBuf;
	TypeLexiTree *dictBirth, *dictDeath, *dictFossil, *dictSampling;
	TypeSamplingScheme sampScheme;

	dictBirth = newLexiTree();
	dictDeath = newLexiTree();
	dictFossil = newLexiTree();
	dictSampling = newLexiTree();
	sizeBuf = INC_SIZEBUF_SC;
	sampScheme.sizeParamTime = 0;
	sampScheme.sizeParamBirth = 0;
	sampScheme.sizeParamDeath = 0;
	sampScheme.sizeParamFossil = 0;
	sampScheme.sizeParamSampling = 0;
	sampScheme.indexTime = (int*) malloc(sizeBuf*sizeof(int));
	sampScheme.schemeTime = (int*) malloc(sizeBuf*sizeof(int));
	sampScheme.schemeBirth = (int*) malloc(sizeBuf*sizeof(int));
	sampScheme.schemeDeath = (int*) malloc(sizeBuf*sizeof(int));
	sampScheme.schemeFossil = (int*) malloc(sizeBuf*sizeof(int));
	sampScheme.schemeSampling = (int*) malloc(sizeBuf*sizeof(int));
	sampScheme.param.size = 0;
	sampScheme.param.param = (TypeModelParam*) malloc(sizeBuf*sizeof(TypeModelParam));
	sampScheme.param.startTime = (double*) malloc((sizeBuf+1)*sizeof(double));
    do {
		c = setSamplingSchemeTime(f, sampScheme.param.size, sampScheme.indexTime, &(sampScheme.sizeParamTime), &(sampScheme.schemeTime[sampScheme.param.size]), &(sampScheme.param.startTime[sampScheme.param.size]));
		if(c!=EOF) {
			c = setSamplingScheme(f, dictBirth, &(sampScheme.sizeParamBirth), &(sampScheme.schemeBirth[sampScheme.param.size]), &(sampScheme.param.param[sampScheme.param.size].birth));
			if(c!=EOF)
				c = setSamplingScheme(f, dictDeath, &(sampScheme.sizeParamDeath), &(sampScheme.schemeDeath[sampScheme.param.size]), &(sampScheme.param.param[sampScheme.param.size].death));
			else
				error("Something goes wrong while reading the sampling scheme (l 173)\n");
			if(c!=EOF)
				c = setSamplingScheme(f, dictFossil, &(sampScheme.sizeParamFossil), &(sampScheme.schemeFossil[sampScheme.param.size]), &(sampScheme.param.param[sampScheme.param.size].fossil));
			else
				error("Something goes wrong while reading the sampling scheme (l 177)\n");
			c = setSamplingScheme(f, dictSampling, &(sampScheme.sizeParamSampling), &(sampScheme.schemeSampling[sampScheme.param.size]), &(sampScheme.param.param[sampScheme.param.size].sampling));
			sampScheme.param.size++;
			if(sampScheme.param.size>=sizeBuf) {
				sizeBuf += INC_SIZEBUF_SC;
				sampScheme.indexTime = (int*) realloc(sampScheme.indexTime, (sizeBuf+1)*sizeof(int));
				sampScheme.schemeTime = (int*) realloc(sampScheme.schemeTime, (sizeBuf+1)*sizeof(int));
				sampScheme.schemeBirth = (int*) realloc(sampScheme.schemeBirth, sizeBuf*sizeof(int));
				sampScheme.schemeDeath = (int*) realloc(sampScheme.schemeDeath, sizeBuf*sizeof(int));
				sampScheme.schemeFossil = (int*) realloc(sampScheme.schemeFossil, sizeBuf*sizeof(int));
				sampScheme.schemeSampling = (int*) realloc(sampScheme.schemeSampling, sizeBuf*sizeof(int));
				sampScheme.param.param = (TypeModelParam*) realloc(sampScheme.param.param, sizeBuf*sizeof(TypeModelParam));
				sampScheme.param.startTime = (double*) realloc(sampScheme.param.startTime, (sizeBuf+1)*sizeof(double));
			}
		}
		for(c = fgetc(f); c != EOF && issepline(c); c = fgetc(f))
			;
		if(c!=EOF)
			ungetc(c, f);	
	} while(c != EOF);
	freeLexiTree(dictBirth);
	freeLexiTree(dictDeath);
	freeLexiTree(dictFossil);
	freeLexiTree(dictSampling);
	sampScheme.indexTime = (int*) realloc(sampScheme.indexTime, (sizeBuf+1)*sizeof(int));
	sampScheme.schemeTime = (int*) realloc(sampScheme.schemeTime, (sizeBuf+1)*sizeof(int));
	sampScheme.schemeBirth = (int*) realloc(sampScheme.schemeBirth, sizeBuf*sizeof(int));
	sampScheme.schemeDeath = (int*) realloc(sampScheme.schemeDeath, sizeBuf*sizeof(int));
	sampScheme.schemeFossil = (int*) realloc(sampScheme.schemeFossil, sizeBuf*sizeof(int));
	sampScheme.schemeSampling = (int*) realloc(sampScheme.schemeSampling, sizeBuf*sizeof(int));
	sampScheme.param.param = (TypeModelParam*) realloc(sampScheme.param.param, sampScheme.param.size*sizeof(TypeModelParam));
	sampScheme.param.startTime = (double*) realloc(sampScheme.param.startTime, (sampScheme.param.size+1)*sizeof(double));
	return sampScheme;
}

void freeSamplingScheme(TypeSamplingScheme *sampScheme) {
	if(sampScheme->indexTime != NULL)
		free((void*) sampScheme->indexTime);
	if(sampScheme->schemeTime != NULL)
		free((void*) sampScheme->schemeTime);
	if(sampScheme->schemeBirth != NULL)
		free((void*) sampScheme->schemeBirth);
	if(sampScheme->schemeDeath != NULL)
		free((void*) sampScheme->schemeDeath);
	if(sampScheme->schemeFossil != NULL)
		free((void*) sampScheme->schemeFossil);
	if(sampScheme->schemeSampling != NULL)
		free((void*) sampScheme->schemeSampling);
	if(sampScheme->param.param != NULL)
		free((void*) sampScheme->param.param);
	if(sampScheme->param.startTime != NULL)
		free((void*) sampScheme->param.startTime);
	sampScheme->indexTime = NULL;
	sampScheme->schemeTime = NULL;
	sampScheme->schemeBirth = NULL;
	sampScheme->schemeDeath = NULL;
	sampScheme->schemeFossil = NULL;
	sampScheme->schemeSampling = NULL;
	sampScheme->param.param = NULL;
	sampScheme->param.startTime = NULL;
}

TypeSamplingCurrent getSamplingCurrent(TypeSamplingScheme sampScheme) {
	TypeSamplingCurrent sampCurrent;
	if(sampScheme.sizeParamTime>0) {
		sampCurrent.currentTime = (double*) malloc(sampScheme.sizeParamTime*sizeof(double));
	} else
		sampCurrent.currentTime = NULL;
	if(sampScheme.sizeParamBirth>0) {
		sampCurrent.currentBirth = (double*) malloc(sampScheme.sizeParamBirth*sizeof(double));
	} else
		sampCurrent.currentBirth = NULL;
	if(sampScheme.sizeParamDeath>0) {
		sampCurrent.currentDeath = (double*) malloc(sampScheme.sizeParamDeath*sizeof(double));
	} else
		sampCurrent.currentDeath = NULL;
	if(sampScheme.sizeParamFossil>0) {
		sampCurrent.currentFossil = (double*) malloc(sampScheme.sizeParamFossil*sizeof(double));
	} else
		sampCurrent.currentFossil = NULL;
	if(sampScheme.sizeParamSampling>0) {
		sampCurrent.currentSampling = (double*) malloc(sampScheme.sizeParamSampling*sizeof(double));
	} else
		sampCurrent.currentSampling = NULL;
	return sampCurrent;
}

void freeSamplingCurrent(TypeSamplingCurrent *sampCurrent) {
	if(sampCurrent->currentTime != NULL)
		free((void*)sampCurrent->currentTime);
	if(sampCurrent->currentBirth != NULL)
		free((void*)sampCurrent->currentBirth);
	if(sampCurrent->currentDeath != NULL)
		free((void*)sampCurrent->currentDeath);
	if(sampCurrent->currentFossil != NULL)
		free((void*)sampCurrent->currentFossil);
	if(sampCurrent->currentSampling != NULL)
		free((void*)sampCurrent->currentSampling);
	sampCurrent->currentTime = NULL;
	sampCurrent->currentBirth = NULL;
	sampCurrent->currentDeath = NULL;
	sampCurrent->currentFossil = NULL;
	sampCurrent->currentSampling = NULL;
}


void setParamFromSamplingCurrent(TypePiecewiseModelParam *param, TypeSamplingCurrent sampCurrent, TypeSamplingScheme sampScheme) {
	int i;
	for(i=0; i<=param->size; i++)
		if(sampScheme.schemeTime[i] >= 0)
			param->startTime[i] = sampCurrent.currentTime[sampScheme.schemeTime[i]];
		else
			param->startTime[i]  = sampScheme.param.startTime[i];
	for(i=0; i<param->size; i++) {
		if(sampScheme.schemeBirth[i] >= 0)
			param->param[i].birth = sampCurrent.currentBirth[sampScheme.schemeBirth[i]];
		else
			param->param[i].birth = sampScheme.param.param[i].birth;
		if(sampScheme.schemeDeath[i] >= 0)
			param->param[i].death = sampCurrent.currentDeath[sampScheme.schemeDeath[i]];
		else
			param->param[i].death = sampScheme.param.param[i].death;
		if(sampScheme.schemeFossil[i] >= 0)
			param->param[i].fossil = sampCurrent.currentFossil[sampScheme.schemeFossil[i]];
		else
			param->param[i].fossil = sampScheme.param.param[i].fossil;
		if(sampScheme.schemeSampling[i] >= 0)
			param->param[i].sampling = sampCurrent.currentSampling[sampScheme.schemeSampling[i]];
		else
			param->param[i].sampling = sampScheme.param.param[i].sampling;
	}
}

TypePiecewiseModelParam getPiecewiseParamFromSamplingScheme(TypeSamplingScheme sampScheme) {
	TypePiecewiseModelParam param;
	int i;
	param.size = sampScheme.param.size;
	param.startTime = (double*) malloc((param.size+1)*sizeof(double));
	for(i=0; i<=param.size; i++)
		param.startTime[i] = sampScheme.param.startTime[i];
	param.param = (TypeModelParam*) malloc(param.size*sizeof(TypeModelParam));
	return param;
}

void freePiecewiseParam(TypePiecewiseModelParam *param) {
	free((void*)param->startTime);
	free((void*)param->param);
	param->startTime = NULL;
	param->param = NULL;
}
