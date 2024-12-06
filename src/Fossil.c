#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "MyR.h"
#include "Utils.h"
#include "Fossil.h"

#define INC_FOSSIL_ITEM 50
#define MAX_NAME_SIZE 3000




int getMinMaxFossil(int n, TypeFossilFeature *fos, double *min, double *max) {

		if(fos && fos->fossil[n] != NOSUCH) {
			int f;
			*min = fos->fossilList[fos->fossil[n]].time;
			*max = fos->fossilList[fos->fossil[n]].time;
			for(f=fos->fossilList[fos->fossil[n]].prec; f!=NOSUCH; f=fos->fossilList[f].prec) {
				if(fos->fossilList[f].time<*min)
					*min = fos->fossilList[f].time;
				if(fos->fossilList[f].time>*max)
					*max = fos->fossilList[f].time;
			}
			return 1;
		} else
			return 0;
}

void checkTimeConsistency(int n, double min, TypeTree *tree, TypeFossilFeature *ff) {
	int f, c;
	for(f=ff->fossil[n]; f!=NOSUCH; f=ff->fossilList[f].prec)
		if(ff->fossilList[f].time < min)
			warning("Time consistency problem in fossil %d from node %d fos %lf < min %lf\n", n, f, ff->fossilList[f].time, min);
	if(ff->fossil[n]!=NOSUCH && ff->fossilList[ff->fossil[n]].time>min)
		min = ff->fossilList[ff->fossil[n]].time;
	if(tree->time[n] != NO_TIME && tree->time[n] < min) {
		warning("Time consistency problem node %d (status %d, child %d) time %lf < min %lf\n", n, ff->status[n], tree->node[n].child, tree->time[n], min);
		for(f=ff->fossil[n]; f!=NOSUCH; f=ff->fossilList[f].prec)
			printf("fossil %lf\n", ff->fossilList[f].time);
	}
	for(c=tree->node[n].child; c!=NOSUCH; c=tree->node[c].sibling)
		checkTimeConsistency(c, min, tree, ff);
}


void checkConsistencyFossil(TypeTree *tree,  TypeFossilFeature *fos) {
	int n;
	double *min, *max;
	
	printf("check consistency Fossil\n");
	min = (double*) malloc(tree->size*sizeof(double));
	max = (double*) malloc(tree->size*sizeof(double));
	for(n=0; n<tree->size; n++) {
		if(!getMinMaxFossil(n, fos, &(min[n]), &(max[n]))) {
			min[n] = sqrt(-1.);
			max[n] = sqrt(-1.);
		}
	}
	for(n=0; n<tree->size; n++) {
		if(!isnan(max[n])) {
			int c;
			for(c=tree->node[n].child; c!=NOSUCH; c=tree->node[c].sibling) {
				if(!isnan(min[c]) && min[c]<max[n]) {
					if(tree->name != NULL)
						printf("issue %d (%s) max %.2lf -> %d (%s) min %.2lf\n", n, tree->name[n], max[n], c, tree->name[c], min[c]);
					else
						printf("issue %d max %.2lf -> %d min %.2lf\n", n, max[n], c, min[c]);
				}
				if(tree->time[n] != NO_TIME && tree->time[n]<max[n]) {
					if(tree->name != NULL)
						printf("issue %d (%s) max %.2lf -> %d (%s) min %.2lf\n", n, tree->name[n], max[n], c, tree->name[c], tree->time[n]);
					else
						printf("issue %d max %.2lf -> %d min %.2lf\n", n, max[n], c, tree->time[n]);
				}				
			}
		}
	}
}

int checkConsistency(TypeFossilFeature *fp, int size) {
	int i, k;
	for(i=0; i<size; i++) {
		if(fp->fossil[i] != NOSUCH) {
			double last = fp->fossilList[fp->fossil[i]].time;
			for(k=fp->fossilList[fp->fossil[i]].prec; k!=NOSUCH; k=fp->fossilList[k].prec) {
				if(last < fp->fossilList[k].time) {
					printf("problem %.2lf %.2lf\n", last, fp->fossilList[k].time);
					return 0;
				}
				last = fp->fossilList[k].time;
			}
		}
	}
	return 1;
}
		
/*fully duplicate "feat"*/
TypeFossilFeature *cpyFossilFeature(TypeFossilFeature *feat, int n) {
    TypeFossilFeature *res;
    int i;
    if(feat == NULL)
        return NULL;
    res = (TypeFossilFeature*) malloc(sizeof(TypeFossilFeature));
    res->size = feat->size;
    res->sizeBuf = feat->size;
    res->fossil = (int*) malloc(n*sizeof(int));
    for(i=0; i<n; i++)
        res->fossil[i] = feat->fossil[i];
	res->status = (TypeNodeStatus*) malloc(res->sizeBuf*sizeof(TypeNodeStatus));
    res->fossilList = (TypeFossilList*) malloc(res->sizeBuf*sizeof(TypeFossilList));
    for(i=0; i<res->size; i++) {
        res->fossilList[i] = feat->fossilList[i];
        res->status[i] = feat->status[i];
    }
    return res;
}
void negateFossil(TypeFossilFeature* feat) {
    int i;
    for(i=0; i<feat->size; i++) {
        double tmp = feat->fossilList[i].time;
        feat->fossilList[i].time = -feat->fossilList[i].time;
    }
    //if(feat->endTimeTable != NULL) {
		//for(i=0; i<feat->sizeNode; i++) {
			//double tmp = feat->endTimeTable[i].inf;
			//feat->endTimeTable[i].inf = -feat->endTimeTable[i].sup;
			//feat->endTimeTable[i].sup = -tmp;
		//}
    //}		
}

double getMinFossilTime(TypeFossilFeature* feat) {
	int i;
	double min;
	if(feat->size == 0)
		return 0.;
	min = feat->fossilList[0].time;
	for(i=1; i<feat->size; i++)
		if(feat->fossilList[i].time<min)
			min = feat->fossilList[i].time;
	return min;
}

double getMaxFossilTime(TypeFossilFeature* feat) {
	int i;
	double max;
	if(feat->size == 0)
		return 0.;
	max = feat->fossilList[0].time;
	for(i=1; i<feat->size; i++)
		if(feat->fossilList[i].time>max)
			max = feat->fossilList[i].time;
	return max;
}

int compareFossilList(const void* a, const void* b) {
	if(((TypeFossilList*)a)->time>((TypeFossilList*)b)->time)
		return 1;
	if(((TypeFossilList*)a)->time<((TypeFossilList*)b)->time)
		return -1;
	return 0;
}

void freeFossilFeature(TypeFossilFeature *fos) {
	if(fos == NULL)
		return;
	if(fos->fossil != NULL)
		free((void*)fos->fossil);
	if(fos->status != NULL)
		free((void*)fos->status);
	if(fos->fossilList != NULL)
		free((void*)fos->fossilList);
	free((void*)fos);
}


/*print fossilInt*/	
void fprintFossilFeature(FILE *f, TypeFossilFeature *feat, char **name, int size) {
	int i;
	for(i=0; i<size; i++) {
		if(feat->fossil[i]!=NOSUCH) {
			if(name != NULL && name[i] != NULL) {
				fprintf(f, "%s", name[i]);
				fprintFossilList(f, feat->fossil[i], feat->fossilList);
				fprintf(f, "\n");
			}
		}
	}
}
///*print fossilInt*/	
//void fprintFossilFeature(FILE *f, TypeFossilFeature *feat, char **name, int size) {
	//int i;
	//for(i=0; i<size; i++) {
		//if(feat->fossil[i]!=NOSUCH) {
			//if(name != NULL && name[i] != NULL)
				//fprintf(f, "%s (%d)", name[i], i);
			//else
				//fprintf(f, "--- (%d)", i);
			//fprintFossilList(f, feat->fossil[i], feat->fossilList);
			//fprintf(f, "\n");
		//}
	//}
//}

void fixTreeFossil(TypeTree *tree, TypeFossilFeature *fos) {
	int n;
	if(fos == NULL)
		return;
	for(n=0; n<tree->size; n++)
        if(tree->node[n].child == NOSUCH && (tree->time[n] != NO_TIME || fos->fossil[n] == NOSUCH))
            tree->time[n] = tree->maxTime;
        else
            tree->time[n] = NO_TIME;
}

void fixTreeBis(TypeTree *tree, TypeFossilFeature *fos) {
	int n;
	if(fos == NULL)
		return;
	for(n=0; n<tree->size; n++)
        if(tree->node[n].child == NOSUCH && tree->time[n] == NO_TIME && fos->fossil[n] != NOSUCH) {
			int f;
			tree->time[n] = fos->fossilList[fos->fossil[n]].time;
			for(f=fos->fossilList[fos->fossil[n]].prec; f!=NOSUCH; f=fos->fossilList[f].prec)
				if(fos->fossilList[f].time>tree->time[n])
					tree->time[n] = fos->fossilList[f].time;
		}
}

void fixTreeTer(TypeTree *tree, TypeFossilFeature *fos) {
	int n;
	if(fos == NULL)
		return;
	for(n=0; n<tree->size; n++)
        if(tree->node[n].child == NOSUCH && tree->time[n] == NO_TIME) {
            if(fos->fossil[n] != NOSUCH) {
				int f;
				tree->time[n] = fos->fossilList[fos->fossil[n]].time;
				for(f=fos->fossilList[fos->fossil[n]].prec; f!=NOSUCH; f=fos->fossilList[f].prec)
					if(fos->fossilList[f].time>tree->time[n])
						tree->time[n] = fos->fossilList[f].time;
			} else {
				tree->time[n] = tree->maxTime;
			}
		}
}


/*print fossilInt table*/	
void fprintFossilList(FILE *fo, int e, TypeFossilList *list) {
	int f;
	for(f=e; f!=NOSUCH; f=list[f].prec)
		fprintf(fo, "\t %lf", list[f].time);
}


/*returns the number of internal fossil finds of "tree", i.e. the ones which are not just before an extinction*/
int countInternalFossils(double *time, TypeFossilFeature *fos, int size) {
	int n, count = 0;
	for(n=0; n<size; n++) {
		int f;
		for(f=fos->fossil[n]; f!=NOSUCH; f = fos->fossilList[f].prec)
			if(fos->fossilList[f].time<time[n])
				count++;
	}
	return count;
}




TypeInterval readInterval(char *s) {
	int i, ind;
	TypeInterval fossil;
	char tmp[MAX_NAME_SIZE];
	for(i=0; s[i] != '\0' && issep(s[i]); i++);
	if(s[i] == '(') {
		i++;
		for(; s[i] != '\0' && issep(s[i]); i++);
		ind = 0;
		for(; s[i] != '\0' && s[i] != ':' && !issep(s[i]); i++)
			tmp[ind++] = s[i];
		tmp[ind] = '\0';
		fossil.inf =	atof(tmp);
		for(; s[i] != '\0' && issep(s[i]); i++);
		if(s[i] != ':')
			exitProg(ErrorReading, "Missing ':' while reading a fossil time interval.");
		i++;
		for(; s[i] != '\0' && issep(s[i]); i++);
		ind = 0;
		for(; s[i] != '\0' && s[i] != ')' && !issep(s[i]); i++)
			tmp[ind++] = s[i];
		tmp[ind] = '\0';
		fossil.sup =	atof(tmp);
		for(; s[i] != '\0' && issep(s[i]); i++);
		if(s[i] != ')')
			exitProg(ErrorReading, "Missing ')' while reading a fossil time interval.");
	} else {
		fossil.inf = atof(s);		
		fossil.sup = fossil.inf;
	}
	return fossil;
}

/*turn a fossil list into a fossilTab*/
TypeFossilTab *listToFossilTab(TypeFossilFeature *fos, int size) {
	int i;
	TypeFossilTab *res;
	res = (TypeFossilTab *) malloc(size*sizeof(TypeFossilTab));
	for(i=0; i<size; i++)
		if(fos->fossil[i] >= 0) {
			int f, n = 0;
			for(f=fos->fossil[i]; f!=NOSUCH; f=fos->fossilList[f].prec)
				n++;
			res[i].size = n;
			res[i].time = (double*) malloc(n*sizeof(double));
			for(f=fos->fossil[i]; f!=NOSUCH; f=fos->fossilList[f].prec) {
				res[i].time[--n] = fos->fossilList[f].time;
				if(isnan(res[i].time[n])) {
					for(f=fos->fossil[i]; f!=NOSUCH; f=fos->fossilList[f].prec)
						printf("%lf\t", fos->fossilList[f].time);
					error("Problem C node %d\n", i);
				}
			}
		} else {
			res[i].size = 0;
			res[i].time = NULL;
		}
	return res;
}

/*free a fossil tab*/
void freeFossilTab(TypeFossilTab *fosTab, int size) {
	int i;
	if(fosTab != NULL) {
		for(i=0; i<size; i++)
			if(fosTab[i].time != NULL)
				free((void*)fosTab[i].time);
		free((void*)fosTab);
	}
}
			
		

int compareInterval(const void* a, const void* b) {
	if(((TypeInterval*)a)->inf>((TypeInterval*)b)->inf)
		return 1;
	if(((TypeInterval*)a)->inf<((TypeInterval*)b)->inf)
		return -1;
	if(((TypeInterval*)a)->sup>((TypeInterval*)b)->sup)
		return 1;
	if(((TypeInterval*)a)->sup<((TypeInterval*)b)->sup)
		return -1;
	return 0;
}


/*print fossil*/	
void fprintInterval(FILE *f, TypeInterval fos) {
	if(fos.inf == fos.sup)
		fprintf(f, "%lf", fos.inf);
	else
		fprintf(f, "(%lf:%lf)", fos.inf, fos.sup);
}



/*print fossil table*/	
void sprintFossilNHX(char *s, char *prefix, int e, TypeFossilList *list) {
	int f;
	if(prefix != NULL)
		s += sprintf(s, "%s", prefix);
	for(f=e; f!=NOSUCH; f = list[f].prec)
		s += sprintf(s, ":FOS=%lf", list[f].time);
}

char **getFossilComment(char **comment, TypeFossilFeature *fos, int size) {
	char **res, buffer[2000];
	int i;

	if((res = (char**) malloc(size*sizeof(char*))) == NULL)
		return NULL;
	buffer[0] = '\0';
	for(i=0; i<size; i++) {
		if(comment != NULL)
			sprintFossilNHX(buffer, comment[i], fos->fossil[i], fos->fossilList);
		else
			sprintFossilNHX(buffer, NULL, fos->fossil[i], fos->fossilList);
		res[i] = (char*) malloc((strlen(buffer)+1)*sizeof(char));
		strcpy(res[i], buffer);
	}
	return res;
}

/*print tree in newick format*/	
void fprintFossilTreeNewick(FILE *f, TypeTree *tree,  TypeFossilFeature *fos) {
	int i;
    char **commentSave;
    commentSave = tree->comment;
    getFossilComment(commentSave, fos, tree->size);
    fprintTreeNewick(f, tree);
	for(i=0; i<tree->size; i++)
        free((void*)tree->comment[i]);
    free((void*)tree->comment);
    tree->comment = commentSave;
}


/*prune "tree" to that can be observed from contemporary lineages and fossil finds*/
TypeTree *pruneFossil(TypeTree *tree,  TypeFossilFeature *fos) {
	TypeTree *resT;
	TypeFossilFeature *resF;
	int n, f, *new_feat, *parent;
	if(tree->time == NULL)
		return NULL;
	resT = newTree(tree->size);
	resT->maxTime = tree->maxTime;
	resT->maxTimeInt = tree->maxTimeInt;
	resT->minTime = tree->minTime;
	resT->minTimeInt = tree->minTimeInt;
	if(fos)
		resF = (TypeFossilFeature*) malloc(sizeof(TypeFossilFeature));
	else
		resF = NULL;
    resT->info = (void*) resF;
	if(tree->size == 0) {
		resT->size = 0;
		resT->sizeBuf = 0;
		free((void*)resT->node);
		resT->node = NULL;
		free((void*)resT->time);
		resT->time = NULL;
		if(resF) {
			resF->size = 0;
			resF->sizeBuf = 0;
			resF->fossilList = NULL;
			resF->fossil = NULL;
		}
		return resT;
	}
	resT->size = 0;
	if(tree->name)
		resT->name = (char**) malloc(tree->size*sizeof(char*));
	else
		resT->name = NULL;
	if(tree->comment)
		resT->comment = (char**) malloc(tree->size*sizeof(char*));
	else
		resT->comment = NULL;
	if(fos) {
		resF->fossil = (int*) malloc(tree->size*sizeof(int));
		resF->status = NULL;
		resF->size = fos->size;
		resF->sizeBuf = fos->size;
        resF->fossilList = (TypeFossilList*) malloc(fos->size*sizeof(TypeFossilList));
		for(f=0; f<fos->size; f++)
			resF->fossilList[f] = fos->fossilList[f];
	}
	new_feat = (int*) malloc(tree->size*sizeof(int));
	for(n=0; n<tree->size; n++)
        new_feat[n] = NOSUCH;
	parent = getParent(tree);
	for(n=0; n<tree->size; n++) {
		if(tree->time[n]>=tree->maxTime || (fos && fos->fossil[n]>=0)) {
			int m;
			new_feat[n] = resT->size;
            resT->node[resT->size].child = NOSUCH;
            resT->node[resT->size].sibling = NOSUCH;
			resT->time[resT->size] = tree->time[n];
            if(tree->name)
                resT->name[resT->size] = strdpl(tree->name[n]);
            if(tree->comment)
                resT->comment[resT->size] = strdpl(tree->comment[n]);
			if(fos)
				resF->fossil[resT->size] = fos->fossil[n];
			if(tree->time[n]>=tree->maxTime)
				resT->time[resT->size] = tree->time[n];
			else
				resT->time[resT->size] = fos->fossilList[fos->fossil[n]].time;
			resT->size++;
			for(m=parent[n]; m>=0 && new_feat[m]<0; m=parent[m]) {
				new_feat[m] = resT->size;
				resT->node[resT->size].child = resT->size-1;
                resT->node[resT->size].sibling = NOSUCH;
				resT->time[resT->size] = tree->time[m];
                if(tree->name)
                    resT->name[resT->size] = strdpl(tree->name[m]);
                if(tree->comment)
                    resT->comment[resT->size] = strdpl(tree->comment[m]);
				if(fos)
					resF->fossil[resT->size] = fos->fossil[m];
				resT->size++;
			}
			if(m>=0) {
				resT->time[new_feat[m]] = tree->time[m];
				if(resT->node[new_feat[m]].child>=0) {
					resT->node[resT->node[new_feat[m]].child].sibling = resT->size-1;
				} else {
					resT->node[new_feat[m]].child = resT->size-1;
				}
			}
		}
	}
	resT->node = (TypeNode*) realloc(resT->node, resT->size*sizeof(TypeNode));
	if(resT->time)
		resT->time = (double*) realloc(resT->time, resT->size*sizeof(double));
	if(tree->name)
		resT->name = (char**) realloc(resT->name, resT->size*sizeof(char*));
	if(tree->comment)
		resT->comment = (char**) realloc(resT->comment, resT->size*sizeof(char*));
	if(resF)
		resF->fossil = (int*) realloc(resF->fossil, resT->size*sizeof(int));
	free((void*)new_feat);
	free((void*)parent);
	parent = getParent(resT);
	for(n=0; n<resT->size && parent[n]>=0; n++);
	free((void*)parent);
	resT->root = n;
	return resT;
}

/*prune "tree" to that can be observed from contemporary lineages and fossil finds*/
TypeTree *pruneFossilBis(TypeTree *tree,  TypeFossilFeature *fos) {
	TypeTree *resT;
	TypeFossilFeature *resF;
	int n, f, *parent, *index;
	if(tree->time == NULL)
		return NULL;
	resT = newTree(tree->size);
	resT->maxTime = tree->maxTime;
	resT->minTime = tree->minTime;
	if(fos)
		resF = (TypeFossilFeature*) malloc(sizeof(TypeFossilFeature));
	else
		resF = NULL;
    resT->info = (void*) resF;
	if(tree->size == 0) {
		resT->size = 0;
		resT->sizeBuf = 0;
		free((void*)resT->node);
		resT->node = NULL;
		free((void*)resT->time);
		resT->time = NULL;
		if(resF) {
			resF->size = 0;
			resF->sizeBuf = 0;
			resF->fossilList = NULL;
			resF->fossil = NULL;
		}
		return resT;
	}
	resT->size = 0;
	if(fos) {
		resF->fossil = (int*) malloc(tree->size*sizeof(int));
		resF->status = (TypeNodeStatus*) malloc(tree->size*sizeof(TypeNodeStatus));
		resF->size = fos->size;
		resF->sizeBuf = fos->size;
        resF->fossilList = (TypeFossilList*) malloc(fos->size*sizeof(TypeFossilList));
		for(f=0; f<fos->size; f++)
			resF->fossilList[f] = fos->fossilList[f];
	}
	index = (int*) malloc(tree->size*sizeof(int));
	for(n=0; n<tree->size; n++)
		index[n] = NOSUCH;
	parent = getParent(tree);
	for(n=0; n<tree->size; n++)
		if(index[n] == NOSUCH && (tree->time[n]>=tree->maxTime || (fos && fos->fossil[n]>=0))) {
			int m;
			index[n] = resT->size++;
			for(m=parent[n]; m != NOSUCH && index[m] == NOSUCH; m = parent[m])
				index[m] = resT->size++;
		}
	for(n=0; n<tree->size; n++) {
		if(index[n] != NOSUCH) {
			int *prec, c;
			resT->time[index[n]] = tree->time[n];
            if(tree->name)
                resT->name[index[n]] = strdpl(tree->name[n]);
            if(tree->comment)
                resT->comment[index[n]] = strdpl(tree->comment[n]);
			if(fos) {
				resF->fossil[index[n]] = fos->fossil[n];
			}
			prec = &(resT->node[index[n]].child);
			for(c=tree->node[n].child; c!=NOSUCH; c=tree->node[c].sibling)
				if(index[c] != NOSUCH) {
					*prec = index[c];
					prec =  &(resT->node[index[c]].sibling);
				}
			*prec = NOSUCH;
			if(resT->node[index[n]].child == NOSUCH) {
				//if(tree->time[n] == tree->maxTime) {
					//resT->time[index[n]] = tree->time[n];
					//resF->status[index[n]] = contempNodeStatus;
				//} else {
					//if(tree->time[n] == NO_TIME) {
						//resT->time[index[n]] = fos->fossilList[fos->fossil[n]].time;
						//resF->status[index[n]] = extinctNodeStatus;
					//} else {
						//warning("bizarre node %d time %lf (status %d)\n", index[n], tree->time[n], fos->status[n]);
						//resT->time[index[n]] = tree->time[n];
						//resF->status[index[n]] = unknownNodeStatus;
					//}	
				//}
				if(tree->node[n].child == NOSUCH) {
					if(fos->status[n] == contempNodeStatus) {
						resT->time[index[n]] = tree->maxTime;
						resF->status[index[n]] = contempNodeStatus;
					} else {
						if(fos->status[n] == extinctNodeStatus) {
							resT->time[index[n]] = fos->fossilList[fos->fossil[n]].time;
							resF->status[index[n]] = extinctNodeStatus;
						} else {
							if(fos->status[n] == unknownNodeStatus) {
								if(tree->time[n] == NO_TIME)
									error("Error while pruning fossil, unknown node without time\n");
								resT->time[index[n]] = tree->time[n];
								resF->status[index[n]] = unknownNodeStatus;
							} else
								error("tip %d time %lf (status %d)\n", index[n], tree->time[n], fos->status[n]);
						}	
					}
				} else {
					if(fos->fossil[n] == NOSUCH)
						error("internal node %d time %lf (status %d)\n", index[n], tree->time[n], fos->status[n]);
					resT->time[index[n]] = fos->fossilList[fos->fossil[n]].time;
					resF->status[index[n]] = extinctNodeStatus;
				}
			} else {
				resT->time[index[n]] = tree->time[n];
				resF->status[index[n]] = noneNodeStatus;
			}
		}
	}
	free((void*)index);
	resT->node = (TypeNode*) realloc(resT->node, resT->size*sizeof(TypeNode));
	if(resT->time)
		resT->time = (double*) realloc(resT->time, resT->size*sizeof(double));
	if(tree->name)
		resT->name = (char**) realloc(resT->name, resT->size*sizeof(char*));
	if(tree->comment)
		resT->comment = (char**) realloc(resT->comment, resT->size*sizeof(char*));
	if(resF) {
		resF->fossil = (int*) realloc(resF->fossil, resT->size*sizeof(int));
		resF->status = (TypeNodeStatus*) realloc(resF->status, resT->size*sizeof(TypeNodeStatus));
	}
	resT->sizeBuf = resT->size;
	free((void*)parent);
	parent = getParent(resT);
	for(n=0; n<resT->size && parent[n]>=0; n++);
	free((void*)parent);
	resT->root = n;
	return resT;
}


/*stop "tree"at stopTime*/
TypeTree *stopTreeFossil(TypeTree *tree,  TypeFossilFeature *fos, double stopTime) {
	TypeTree *resT;
	TypeFossilFeature *resF;
	int n, f, *parent, *indexT;
	if(tree->time == NULL)
		return NULL;
	if(tree->time == NULL)
		return NULL;
	for(n=0; n<tree->size; n++)
		if(tree->time[n] == NO_TIME)
			return NULL;
	resT = newTree(tree->size);
	if(tree->name != NULL)
		resT->name = (char**) malloc(tree->size*sizeof(char*));
	if(tree->comment != NULL)
		resT->comment = (char**) malloc(tree->size*sizeof(char*));
	resT->maxTime = tree->maxTime;
	resT->minTime = tree->minTime;
	if(fos)
		resF = (TypeFossilFeature*) malloc(sizeof(TypeFossilFeature));
	else
		resF = NULL;
    resT->info = (void*) resF;
	if(tree->size == 0) {
		resT->size = 0;
		resT->sizeBuf = 0;
		free((void*)resT->node);
		resT->node = NULL;
		free((void*)resT->time);
		resT->time = NULL;
		if(resF) {
			resF->size = 0;
			resF->sizeBuf = 0;
			resF->fossilList = NULL;
			resF->fossil = NULL;
		}
		return resT;
	}
	resT->size = 0;
	if(fos) {
		resF->fossil = (int*) malloc(tree->size*sizeof(int));
		resF->status = (TypeNodeStatus*) malloc(tree->size*sizeof(TypeNodeStatus));
		resF->size = 0;
		resF->sizeBuf = fos->size;
        resF->fossilList = (TypeFossilList*) malloc(fos->size*sizeof(TypeFossilList));
	}
	parent = getParent(tree);
	indexT = (int*) malloc(tree->size*sizeof(int));
	for(n=0; n<tree->size; n++) {
		if(tree->time[n] < stopTime || parent[n] == NOSUCH || (parent[n] != NOSUCH && tree->time[parent[n]] < stopTime))
			indexT[n] = resT->size++;
		else
			indexT[n] = NOSUCH;
	}
//printf("stoptime %lf size %d\n", stopTime, resT->size);
	for(n=0; n<tree->size; n++) {
		if(indexT[n] != NOSUCH) {
			int *prec, c;
			if(tree->time[n] > stopTime) {
				resT->time[indexT[n]] = stopTime;
				resT->node[indexT[n]].child = NOSUCH;
			} else {
				resT->time[indexT[n]] = tree->time[n];
				prec = &(resT->node[indexT[n]].child);
				for(c=tree->node[n].child; c!=NOSUCH; c=tree->node[c].sibling)
					if(indexT[c] != NOSUCH) {
						*prec = indexT[c];
						prec =  &(resT->node[indexT[c]].sibling);
					}
				*prec = NOSUCH;
			}
            if(tree->name)
                resT->name[indexT[n]] = strdpl(tree->name[n]);
            if(tree->comment)
                resT->comment[indexT[n]] = strdpl(tree->comment[n]);
			if(fos) {
				if(fos->fossil[n] != NOSUCH) {
					for(f=fos->fossil[n]; f!=NOSUCH && fos->fossilList[f].time >= stopTime; f=fos->fossilList[f].prec)
						;
					if(f != NOSUCH) {
						resF->fossil[indexT[n]] = resF->size;
						resF->fossilList[resF->size].time = fos->fossilList[f].time;
						resF->size++;
						for(f=fos->fossilList[f].prec; f!=NOSUCH; f=fos->fossilList[f].prec) {
							resF->fossilList[resF->size].time = fos->fossilList[f].time;
							resF->fossilList[resF->size-1].prec = resF->size;
							resF->size++;
						}
						resF->fossilList[resF->size-1].prec = NOSUCH;
					} else
						resF->fossil[indexT[n]] = NOSUCH;
//{
//printf("node %d -> %d (%s %s)\n", n, indexT[n], tree->name[n], resT->name[indexT[n]]);
//for(f=fos->fossil[n]; f!=NOSUCH; f=fos->fossilList[f].prec)
	//printf("\t%lf", fos->fossilList[f].time);
//printf("\n");
//for(f=resF->fossil[indexT[n]]; f!=NOSUCH; f=resF->fossilList[f].prec)
	//printf("\t%lf", resF->fossilList[f].time);
//printf("\n\n");
//}
				} else
					resF->fossil[indexT[n]] = NOSUCH;
			}
			if(resT->node[indexT[n]].child == NOSUCH) {
				if(tree->node[n].child == NOSUCH) {
					if(tree->time[n] >= stopTime)
						resF->status[indexT[n]] = unknownNodeStatus;
					else {
						resT->time[indexT[n]] = resF->fossilList[resF->fossil[indexT[n]]].time;
						resF->status[indexT[n]] = extinctNodeStatus;
					}
				} else
					resF->status[indexT[n]] = unknownNodeStatus;
			} else {
				resT->time[indexT[n]] = tree->time[n];
				resF->status[indexT[n]] = noneNodeStatus;
			}
		}
	}
	free((void*)indexT);
	resT->node = (TypeNode*) realloc(resT->node, resT->size*sizeof(TypeNode));
	if(resT->time)
		resT->time = (double*) realloc(resT->time, resT->size*sizeof(double));
	if(tree->name)
		resT->name = (char**) realloc(resT->name, resT->size*sizeof(char*));
	if(tree->comment)
		resT->comment = (char**) realloc(resT->comment, resT->size*sizeof(char*));
	if(resF) {
		resF->fossil = (int*) realloc((void*) resF->fossil, resT->size*sizeof(int));
		resF->status = (TypeNodeStatus*) realloc((void*) resF->status, resT->size*sizeof(TypeNodeStatus));
		resF->fossilList = (TypeFossilList*) realloc((void*) resF->fossilList, resF->size*sizeof(TypeFossilList));
		resF->sizeBuf = resF->size;
	}
	resT->sizeBuf = resT->size;
	free((void*)parent);
	parent = getParent(resT);
	for(n=0; n<resT->size && parent[n]>=0; n++);
	free((void*)parent);
	resT->root = n;
	return resT;
}


int iterateBinaryFossil(int n, TypeTree *resT,  TypeFossilFeature *resF, TypeTree *tree,  TypeFossilFeature *fos) {
    int m, ftmp = NOSUCH;
	if(n<0)
		return -1;
	for(m=n; tree->node[m].child!=NOSUCH && tree->node[tree->node[m].child].sibling==NOSUCH; m=tree->node[m].child) {
		if(fos && fos->fossil[m]!=NOSUCH) {
			int f;
			for(f=fos->fossil[m]; resF->fossilList[f].prec>=0; f=resF->fossilList[f].prec);
			resF->fossilList[f].prec = ftmp;
			ftmp = fos->fossil[m];
		}
	}
	if(fos && fos->fossil[m]!=NOSUCH) {
		int f;
		for(f=fos->fossil[m]; resF->fossilList[f].prec!=NOSUCH; f=resF->fossilList[f].prec);
		resF->fossilList[f].prec = ftmp;
		ftmp = fos->fossil[m];
	}
	if(tree->node[m].child!=NOSUCH) {
		int c1, c2, c;
        c1 = iterateBinaryFossil(tree->node[m].child, resT, resF, tree, fos);
		c2 = c1;
		for(c=tree->node[tree->node[m].child].sibling; c >= 0; c = tree->node[c].sibling) {
            resT->node[c2].sibling = iterateBinaryFossil(c, resT, resF, tree, fos);
			c2 = resT->node[c2].sibling;
		}
        resT->node[c2].sibling = NOSUCH;
		resT->node[resT->size].child = c1;
	} else {
        resT->node[resT->size].child = NOSUCH;
	}
	if(fos) {
		resF->status[resT->size] = fos->status[m];
		resF->fossil[resT->size] = ftmp;
	}
	if(tree->time)
		resT->time[resT->size] = tree->time[m];
	if(tree->name)
		resT->name[resT->size] = strdpl(tree->name[m]);
	if(tree->comment)
		resT->comment[resT->size] = strdpl(tree->comment[m]);
	resT->size++;
	return resT->size-1;
}

TypeTree *fixBinaryFossil(TypeTree *tree,  TypeFossilFeature *fos) {
	TypeTree *resT;
	TypeFossilFeature *resF;
    int f;
	
	resT = (TypeTree*) malloc(sizeof(TypeTree));
	if(fos)
		resF = (TypeFossilFeature*) malloc(sizeof(TypeFossilFeature));
	else
		resF = NULL;
    resT->info = resF;
	resT->maxTime = tree->maxTime;
	resT->maxTimeInt = tree->maxTimeInt;
	resT->minTime = tree->minTime;
	resT->minTimeInt = tree->minTimeInt;
	if(tree->size == 0) {
		resT->root = 0;
		resT->size = 0;
		resT->sizeBuf = 0;
		resT->node = NULL;
		resT->time = NULL;
        resT->parent = NULL;
        resT->name = NULL;
        resT->comment = NULL;
		if(resF) {
			resF->size = 0;
			resF->sizeBuf = 0;
			resF->fossilList = NULL;
			resF->fossil = NULL;
		}
		return resT;
	}
	resT->size = 0;
	resT->sizeBuf = tree->size;
	resT->node = (TypeNode*) malloc(tree->size*sizeof(TypeNode));
	resT->time = (double*) malloc(tree->size*sizeof(double));
    resT->parent = NULL;
	if(tree->name)
		resT->name = (char**) malloc(tree->size*sizeof(char*));
	else
		resT->name = NULL;
	if(tree->comment)
		resT->comment = (char**) malloc(tree->size*sizeof(char*));
	else
		resT->comment = NULL;
	if(fos) {
		resF->fossil = (int*) malloc(tree->size*sizeof(int));
		resF->status = (TypeNodeStatus*) malloc(tree->size*sizeof(TypeNodeStatus));
		resF->size = fos->size;
		resF->sizeBuf = fos->size;
        resF->fossilList = (TypeFossilList*) malloc(fos->size*sizeof(TypeFossilList));
		for(f=0; f<fos->size; f++)
			resF->fossilList[f] = fos->fossilList[f];
	}
    iterateBinaryFossil(tree->root, resT, resF, tree, fos);
	resT->node = (TypeNode*) realloc(resT->node, resT->size*sizeof(TypeNode));
	if(resT->time)
		resT->time = (double*) realloc(resT->time, resT->size*sizeof(double));
	if(tree->name)
		resT->name = (char**) realloc(resT->name, resT->size*sizeof(char*));
	if(tree->comment)
		resT->comment = (char**) realloc(resT->comment, resT->size*sizeof(char*));
	if(resF) {
		resF->fossil = (int*) realloc(resF->fossil, resT->size*sizeof(int));		
		resF->status = (TypeNodeStatus*) realloc(resF->status, resT->size*sizeof(TypeNodeStatus));
	}	
	resT->sizeBuf = resT->size;
	resT->root = getRoot(resT);
	return resT;
}



/*print fossilInt table*/
int sprintFossilListNHX(char *s, int e, TypeFossilList *list) {
    int f, ind = 0;
    for(f=e; f!=NOSUCH; f = list[f].prec) {
        ind += sprintf(s+ind, " :FOS=%lf", list[f].time);
    }
    return ind;
}


#define STRING_INTER_LENGTH 60
void fillCommentFossil(TypeTree *tree, TypeFossilFeature *fos) {
    if(tree == NULL || fos == NULL || tree->size == 0)
        return;
    if(tree->comment == NULL)
        tree->comment = (char**) malloc(tree->size*sizeof(char*));
    int n;
    for(n=0; n<tree->root; n++) {
        int f, nb = 0, length;
        for(f=fos->fossil[n]; f!=NOSUCH; f = fos->fossilList[f].prec)
            nb++;
        if(tree->comment[n] != NULL)
            free((void*)tree->comment[n]);
        if(nb>0) {
            tree->comment[n] = (char*) malloc(nb*STRING_INTER_LENGTH*sizeof(char)+1);
            length = sprintFossilListNHX(tree->comment[n], fos->fossil[n], fos->fossilList);
            tree->comment[n][length] = '\0';
            tree->comment[n] = (char*) realloc(tree->comment[n], (length+1)*sizeof(char));
        } else
            tree->comment[n] = NULL;
    }
    int f, nb = 0, ind = 0;
    for(f=fos->fossil[tree->root]; f!=NOSUCH; f = fos->fossilList[f].prec)
        nb++;
    if(tree->comment[tree->root] != NULL)
        free((void*)tree->comment[tree->root]);
    if(tree->minTimeInt.inf!=NO_TIME && tree->minTimeInt.sup!=NO_TIME)
        nb++;
    if(tree->maxTimeInt.inf!=NO_TIME && tree->maxTimeInt.sup!=NO_TIME)
        nb++;
    if(nb>0) {
        tree->comment[tree->root] = (char*) malloc((nb)*STRING_INTER_LENGTH*sizeof(char)+1);
        if(tree->minTimeInt.inf!=NO_TIME && tree->minTimeInt.sup!=NO_TIME) {
            ind += sprintf(tree->comment[tree->root]+ind, " :ORI=%lf", tree->minTimeInt.inf);
        }
        if(tree->maxTimeInt.inf!=NO_TIME && tree->maxTimeInt.sup!=NO_TIME) {
            ind += sprintf(tree->comment[tree->root]+ind, " :END=%lf", tree->maxTimeInt.sup);
        }
        ind += sprintFossilListNHX(tree->comment[tree->root]+ind, fos->fossil[tree->root], fos->fossilList);
        tree->comment[tree->root][ind] = '\0';
        tree->comment[tree->root] = (char*) realloc(tree->comment[tree->root], (ind+1)*sizeof(char));
    } else
        tree->comment[tree->root] = NULL;
    for(n=tree->root+1; n<tree->size; n++) {
        int f, nb = 0, length;
        for(f=fos->fossil[n]; f!=NOSUCH; f = fos->fossilList[f].prec)
            nb++;
        if(tree->comment[n] != NULL)
            free((void*)tree->comment[n]);
        if(nb>0) {
            tree->comment[n] = (char*) malloc(nb*STRING_INTER_LENGTH*sizeof(char)+1);
            length = sprintFossilListNHX(tree->comment[n], fos->fossil[n], fos->fossilList);
            tree->comment[n][length] = '\0';
            tree->comment[n] = (char*) realloc(tree->comment[n], (length+1)*sizeof(char));
        } else
            tree->comment[n] = NULL;
    }
}

/*Save tree with fossil as comments*/
void fprintTreeFossil(FILE *f, TypeTree *tree, TypeFossilFeature *fos) {
    char **comment_saved, **name_saved;
	int n;
    name_saved = tree->name;
    if(tree->name == NULL)
        tree->name = nameLeaves("leaf", tree);
    comment_saved = tree->comment;
    tree->comment = (char**) malloc(tree->sizeBuf*sizeof(char*));
    for(n=0; n<tree->sizeBuf; n++)
        tree->comment[n] = NULL;
    fillCommentFossil(tree, fos);
    tree->name = nameLeaves("leaf", tree);
    fprintSubtreeNewick(f, tree->root, tree);
    for(n=0; n<tree->sizeBuf; n++)
        if(tree->comment[n] != NULL)
            free((void*)tree->comment[n]);
    free((void*)tree->comment);
    tree->comment = comment_saved;
    if(name_saved == NULL) {
        for(n=0; n<tree->sizeBuf; n++)
            if(tree->name[n] != NULL)
                free((void*)tree->name[n]);
        free((void*)tree->name);
    }
    tree->name = name_saved;
}

void fillBoundsFossil(int n, double tmin, double tmax, TypeTree *tree,  TypeFossilFeature *fos, double *min, double *max, int *dmax) {
	int c;
	if(tree->time[n] != NO_TIME) {
		min[n] = tree->time[n];
	} else {
		if(fos && fos->fossil[n] != NOSUCH) {
			min[n] = fos->fossilList[fos->fossil[n]].time;
		} else {
			min[n] = tmin;
		}
	}
	for(c=tree->node[n].child; c!=NOSUCH; c=tree->node[c].sibling)
		fillBoundsFossil(c, min[n], tmax, tree,  fos, min, max, dmax);
	if(tree->time[n] != NO_TIME) {
		max[n] = tree->time[n];
		dmax[n] = 0;
	} else {
		if(tree->node[n].child==NOSUCH) {
			max[n] = tmax;
			dmax[n] = 0;
		} else {	
			max[n] = tmax+1;
			dmax[n] = 0;
			for(c=tree->node[n].child; c!=NOSUCH; c=tree->node[c].sibling) {
				if(fos && fos->fossil[c] != NOSUCH) {
					int f;
					for(f=fos->fossil[c]; f!=NOSUCH; f=fos->fossilList[f].prec)
						if(fos->fossilList[f].time<max[n]) {
							max[n] = fos->fossilList[f].time;
							dmax[n] = 0;
						}
				} else {
					if((max[c])<(max[n])) {
						max[n] = max[c];
						dmax[n] = dmax[c]+1;
					}
				}
			}
		}
	}
}

void fillUnknownTimesFossil(double tmin, double tmax, TypeTree *tree,  TypeFossilFeature *fos) {
	int *dmax;
	double *min, *max;
	
	min = (double*) malloc(tree->size*sizeof(double));
	max = (double*) malloc(tree->size*sizeof(double));
	dmax = (int*) malloc(tree->size*sizeof(int));
	fillBoundsFossil(tree->root, tmin, tmax, tree,  (TypeFossilFeature*) fos, min, max, dmax);
	fillTime(tree->root, tmin, tree, min, max, dmax);
	free((void*)min);
	free((void*)max);
	free((void*)dmax);
}


void fixFossilOrder(TypeFossilFeature* sample, int size) {
    int n, f, *tmp;
    size_t *index;
    index = qsortindex(sample->fossilList, sample->size, sizeof(TypeFossilList), compareFossilList);
    for(n=0; n<size; n++)
        if(sample->fossil[n]>=0)
            sample->fossil[n] = index[sample->fossil[n]];
    for(f=0; f<sample->size; f++)
        if(sample->fossilList[f].prec!=NOSUCH)
            sample->fossilList[f].prec = index[sample->fossilList[f].prec];
    free((void*) index);
    tmp = (int*) malloc((sample->size+1)*sizeof(int));
    for(n=0; n<size; n++) {
        int ind = 0;
        for(f=sample->fossil[n]; f!=NOSUCH; f=sample->fossilList[f].prec)
            tmp[ind++] = f;
        if(ind>0) {
            qsort(tmp, ind, sizeof(int), compareInt);
            sample->fossilList[tmp[0]].prec = NOSUCH;
            for(f=1; f<ind; f++)
                sample->fossilList[tmp[f]].prec = tmp[f-1];
            sample->fossil[n] = tmp[ind-1];
        }
    }
    free((void*) tmp);
}

void printFossilDebug(FILE *fo, TypeFossilFeature* sample, int size) {
    int n, f;
    for(n=0; n<size; n++) {
        if(sample->fossil[n] != NOSUCH) {
			fprintf(fo, "node %d", n);
			for(f=sample->fossil[n]; f!=NOSUCH; f=sample->fossilList[f].prec)
				fprintf(fo, ", %.3lf", sample->fossilList[f].time);
			fprintf(fo, "\n");
		}
     }
 	fprintf(fo, "\n");
}
