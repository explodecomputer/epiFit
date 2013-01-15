#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <dlfcn.h>

//#include "CMemLeak.h"

#define PI 3.1415927
#define NPATS 36

// globals
void *sofile;
void (*fprob)(int*,int*,double*,double*,double*,int*);

int USEED, UNID, UINCREMENT, UNQTL, *UPAT, UGENERATIONS, UREPEATS, UCONDITION, NENUM;
float UHSQ, **URSQ, UPROPSAMPLE, *UEFFECT, **UFREQ, **allpatterns;
char UFILENAME[50];

typedef struct genotype {
        char a[2];
        char b[2];
        char X;
        char alink[2];
        char blink[2];
        char Xlink;
} genotype;

typedef struct record {
        float phen;
        char sex;
        genotype *qtl;
} record;

typedef struct stats {
	float A;
	float B;
	float Vcomp[8];
	float maf[2];
	float F[2];
	float F2[2];
	float Fmarg;
	float Fint;
	float Ffull;
} stats;

typedef struct allstats {
	stats report[2];
} allstats;

void rnorm(float *x, int n);
void wam(int n, double *p, int *a, int nans, int *ans);
void sample2vals(float freq, int N, char *destination);
void addnoise(record *pop, int N, float vare);
void tables(record *pop, int N, int nqtl, float **pattern);
float basepop(int N, int nqtl, float **p, float **pattern, float hsq, record *pop);
float basepopmutant(int N, int nqtl, float **p, float **pattern, float hsq, record *pop);
void nextgen(int N, int nqtl, float **pattern, float vare, record *pop);
float calcstats(record *pop, int N, int nqtl, allstats *stat);
float calcstats1(record *pop, int N, int nqtl, allstats *stat);
float calcstatsnoia(record *pop, int N, int nqtl, allstats *stat);
void calcmaf(record *pop, int loc, float maf[2]);
void calchapfreq(float maf, float Rsq, float hapfreq[4]);
void makeld(record *pop, float **Rsq);
void calcld(record *pop);
void ftest8df(record *pop, int NID, int loc, allstats *stat, int linked);
void ftest2df(record *pop, int NID, int loc, allstats *stat, int linked);
void ftestlin(record *pop, int NID, int loc, allstats *stat, int linked);
void printstats(allstats *stat, int nqtl);
void printpop(record *pop, int N, int nqtl);
void printpattern(float **pattern, float **p, int nqtl);
int checkfixation1(allstats *stat, int nqtl, int repeat);
int checkfixation(allstats *stat, int nqtl);
void countgenerations(allstats *stat, int nqtl, int generations, int *ngen);
int runsimul(record *pop, int N, int nqtl, float hsq, int *pat, int rep, float **pattern, float vare, int generations, allstats *stat, FILE *out, int condition);
void inputsummary(char *inputfile);
void ShowError();

int ***make_cube(int x, int y, int z);
void free_cube(int ***array);
float ***make_cube_d(int x, int y, int z);
void free_cube_d(float ***array);
int **get_space(int m, int n);
void release_space(int **X);
float **get_space_f(int m, int n);
void release_space_f(float **X);
void make_ssi(float **Ssi, float *p);
void make_ss(float **Ss, float *p);
void solve(float **Ssi, float *e, float *g, int ngen);
void kronecker(float **A, float **B, float **AB, int n1, int n2, int m1, int m2);
void noiapar(record *pop, int M, int loc, allstats *stat, char test, int linked);

