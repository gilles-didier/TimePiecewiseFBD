#ifndef FossilF
#define FossilF

#include <stdlib.h>
#include <stdio.h>
#include "Utils.h"
#include "Tree.h"

#define TAG_FOSSIL "FOS"
#define TAG_ORIGIN "ORI"
#define TAG_END "END"
#define TAG_STATUS "STAT"
#define TAG_TIME "TIME"

#define INC_BUF_TIME_TABLE 3

typedef enum {
    noneNodeStatus=0,
    contempNodeStatus,
    extinctNodeStatus,
    unknownNodeStatus,
    divergenceNodeStatus
} TypeNodeStatus;



typedef struct INTERVAL {
	double inf, sup;
} TypeInterval;

typedef struct FOSSIL_LIST {
	double time;
	int prec;
} TypeFossilList;


typedef struct FOSSIL_TAB {
	double *time;
	int size;
} TypeFossilTab;


/*Entries of fossilList are assumed sorted in increasing order with regard to the times while the chain prec have to go in decreasing order again with time*/
typedef struct FOSSIL_FEATURE {
	TypeFossilList *fossilList;
    TypeTimeInterval *endTimeTable;
    TypeNodeStatus *status;
    int *fossil, size, sizeBuf;
} TypeFossilFeature;

typedef struct STANDARD_AND_FOSSIL_FEATURE {
//	TypeStandardFeature *std;
	TypeFossilFeature *fos;
} TypeStandardAndFossilFeature;


typedef struct BASIC_FOSSIL_FEATURE {	char **name, **comment;
	double *time, maxTime, minTime;
	TypeFossilList *fossilList;
	int *fossil, size, sizeBuf;
} TypeBasicFossilFeature;


#ifdef __cplusplus
extern "C" {
#endif

void negateFossil(TypeFossilFeature* feat);
void checkConsistencyFossil(TypeTree *tree,  TypeFossilFeature *fos);
TypeFossilFeature *cpyFossilFeature(TypeFossilFeature *feat, int n);
double getMinFossilTime(TypeFossilFeature* feat);
double getMaxFossilTime(TypeFossilFeature* feat);
void fixTreeBis(TypeTree *tree, TypeFossilFeature *fos);
void fixTreeTer(TypeTree *tree, TypeFossilFeature *fos);
void fillUnknownTimesFossil(double tmin, double tmax, TypeTree *tree,  TypeFossilFeature *fos);
int checkConsistency(TypeFossilFeature *fp, int size);
void checkTimeConsistency(int n, double min, TypeTree *tree, TypeFossilFeature *ff);

/*compare two fossils a and b with regard to their times*/
int compareFossilList(const void* a, const void* b);
/*desallocate memory space of fos*/
void freeFossilFeature(TypeFossilFeature *fos);
/*returns the number of internal fossil finds of "tree", i.e. the ones which are not just before an extinction*/
int countInternalFossils(double *time, TypeFossilFeature *fos, int size);
int compareInterval(const void* a, const void* b);
/*print fossil*/	
void fprintInterval(FILE *f, TypeInterval fos);
/*print tree in newick format*/	
void fprintFossilTreeNewick(FILE *f, TypeTree *tree,  TypeFossilFeature *fos);
/*prune "tree" to that can be observed from contemporary lineages and fossil finds*/
TypeTree *pruneFossil(TypeTree *tree,  TypeFossilFeature *fos);
int iterateBinaryFossil(int n, TypeTree *resT,  TypeFossilFeature *resF, TypeTree *tree,  TypeFossilFeature *fos);
TypeTree *fixBinaryFossil(TypeTree *tree,  TypeFossilFeature *fos);
/*print fossil table*/	
void sprintFossilNHX(char *s, char *prefix, int e, TypeFossilList *list);
char **getFossilComment(char **comment, TypeFossilFeature *fos, int size);
TypeFossilTab *listToFossilTab(TypeFossilFeature *fos, int size);
/*free a fossil tab*/
void freeFossilTab(TypeFossilTab *fosTab, int size);
void fprintFossilFeature(FILE *f, TypeFossilFeature *feat, char **name, int size);
void fprintFossilList(FILE *fo, int e, TypeFossilList *list);
int sprintFossilListNHX(char *s, int e, TypeFossilList *list);
void fillCommentFossil(TypeTree *tree, TypeFossilFeature *fos);
/*Save tree with fossil as comments*/
void fprintTreeFossil(FILE *f, TypeTree *tree, TypeFossilFeature *fos);
void fixTreeFossil(TypeTree *tree, TypeFossilFeature *fos);
TypeTree *pruneFossilBis(TypeTree *tree,  TypeFossilFeature *fos);
void fixFossilOrder(TypeFossilFeature* sample, int size);
void printFossilDebug(FILE *f, TypeFossilFeature* sample, int size);
/*stop "tree"at stopTime*/
TypeTree *stopTreeFossil(TypeTree *tree,  TypeFossilFeature *fos, double stopTime);

#ifdef __cplusplus
}
#endif

#endif
