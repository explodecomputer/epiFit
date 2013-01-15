// Gib Hemani
// Created - 03/08/2010
// Converting popsim<x>.R to C, to allow larger simulations
// Warning - svn is very out of date

#include "pop.h"
//#define _DEBUG
/*
float allpatterns[56][9] = {{0,0,0,0,0,0,0,0,0}, // NULL
{0,0,0,0,0,0,0,0,1},
{0,0,0,0,0,0,0,1,0},
{0,0,0,0,0,0,0,1,1},
{0,0,0,0,0,0,1,0,1},
{0,0,0,0,0,0,1,1,1}, // A+D
{0,0,0,0,0,1,0,1,0},
{0,0,0,0,0,1,0,1,1},
{0,0,0,0,0,1,1,0,0},
{0,0,0,0,0,1,1,0,1},
{0,0,0,0,0,1,1,1,0},
{0,0,0,0,0,1,1,1,1},
{0,0,0,0,1,0,0,0,0},
{0,0,0,0,1,0,0,0,1},
{0,0,0,0,1,0,0,1,0},
{0,0,0,0,1,0,0,1,1},
{0,0,0,0,1,0,1,0,1},
{0,0,0,0,1,0,1,1,1},
{0,0,0,0,1,1,0,1,0},
{0,0,0,0,1,1,0,1,1},
{0,0,0,0,1,1,1,0,0},
{0,0,0,0,1,1,1,0,1},
{0,0,0,0,1,1,1,1,0},
{0,0,0,1,0,1,0,0,0},
{0,0,0,1,0,1,0,0,1},
{0,0,0,1,0,1,0,1,0},
{0,0,0,1,0,1,0,1,1},
{0,0,0,1,0,1,1,0,1},
{0,0,0,1,1,1,0,0,0}, // D
{0,0,0,1,1,1,0,0,1},
{0,0,0,1,1,1,0,1,0},
{0,0,0,1,1,1,0,1,1},
{0,0,0,1,1,1,1,0,1},
{0,0,1,0,0,0,1,0,0},
{0,0,1,0,0,0,1,0,1},
{0,0,1,0,0,0,1,1,0},
{0,0,1,0,0,1,1,1,0},
{0,0,1,0,1,0,1,0,0},
{0,0,1,0,1,0,1,0,1},
{0,0,1,0,1,0,1,1,0},
{0,0,1,0,1,1,1,1,0},
{0,0,1,1,0,0,0,0,1},
{0,0,1,1,0,0,0,1,0},
{0,0,1,1,0,0,0,1,1},
{0,0,1,1,0,0,1,0,1},
{0,0,1,1,0,1,0,1,0},
{0,0,1,1,0,1,1,0,0},
{0,0,1,1,1,0,0,0,1},
{0,0,1,1,1,0,0,1,0},
{0,1,0,1,0,1,0,1,0}, // DxD
{0,1,0,1,1,1,0,1,0},
{1,0.5,0,0.5,0.5,0.5,0,0.5,1}, // AxA
{0,1,0,0.5,0.5,0.5,1,0,1},
{1,1,1,0,0,0,1,1,1}, // D
{1,1,1,1,0,0,1,0,0}, // AxD
{0,0,0,0.5,0.5,0.5,1,1,1}}; // A
*/


/*
float allpatterns[40][9] = {{0,0,0,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0,1},
{0,0,0,0,0,0,0,1,0},
{0,0,0,0,0,0,0,1,1},
{0,0,0,0,0,0,1,0,1},
//{0,0,0,0,0,0,1,1,1},
{0,0,0,0,0,1,0,1,0},
{0,0,0,0,0,1,0,1,1},
{0,0,0,0,0,1,1,0,0},
{0,0,0,0,0,1,1,0,1},
{0,0,0,0,0,1,1,1,0},
{0,0,0,0,0,1,1,1,1},
{0,0,0,0,1,0,0,0,0},
{0,0,0,0,1,0,0,0,1},
{0,0,0,0,1,0,0,1,0},
{0,0,0,0,1,0,0,1,1},
{0,0,0,0,1,0,1,0,1},
{0,0,0,0,1,0,1,1,1},
{0,0,0,0,1,1,0,1,0},
{0,0,0,0,1,1,0,1,1},
{0,0,0,0,1,1,1,0,0},
{0,0,0,0,1,1,1,0,1},
{0,0,0,0,1,1,1,1,0},
{0,0,0,1,0,1,0,0,0},
{0,0,0,1,0,1,0,0,1},
{0,0,0,1,0,1,0,1,0},
{0,0,0,1,0,1,0,1,1},
{0,0,0,1,0,1,1,0,1},
//{0,0,0,1,1,1,0,0,0},
//{0,0,0,1,1,1,0,0,1},
//{0,0,0,1,1,1,0,1,0},
{0,0,0,1,1,1,0,1,1},
{0,0,0,1,1,1,1,0,1},
//{0,0.5,1,0.5,0.5,0.5,1,0.5,0},
//{1,0.5,0,0,0.5,1,1,0.5,0},
{0,0,0,0.5,0.5,0.5,1,1,1},
//{1,0,1,0,0,0,1,0,1},
//{1,0,1,0,1,0,1,0,1},
{1,1,1,0,0,0,1,1,1}};
*/


void calcmaf(record *pop, int loc, float maf[2])
{
	int i;
	maf[0] = 0;
	maf[1] = 0;
	for(i = 0; i < UNID; i++)
	{
		maf[0] += (pop[i].qtl[loc].a[0] + pop[i].qtl[loc].a[1]);
		maf[1] += (pop[i].qtl[loc].b[0] + pop[i].qtl[loc].b[1]);
	}
	maf[0] /= (UNID*2);
	maf[1] /= (UNID*2);
	maf[0] = 1 - maf[0];
	maf[1] = 1 - maf[1];
}


void calchapfreq(float maf, float Rsq, float hapfreq[4])
{
	float pq;
	pq = maf * maf * (1-maf) * (1-maf);
	hapfreq[0] = sqrt(Rsq * pq) + maf * maf;
	hapfreq[1] = maf * (1-maf) - sqrt(Rsq * pq);
	hapfreq[2] = (1-maf) * maf - sqrt(Rsq * pq);
	hapfreq[3] = sqrt(Rsq * pq) + (1-maf) * (1-maf);
}

void makeld(record *pop, float **Rsq)
{
	int i,j;
	float hapfreq[4]; // x11, x12, x21, x22
	float maf[2], temp, temp0, temp1;

	for(i = 0; i < UNQTL; i++)
	{
		if(Rsq[i][0] != 1)
		{
			calcmaf(pop,i,maf);
			// locus 1
			calchapfreq(maf[0],Rsq[i][0],hapfreq);
			temp0 = RAND_MAX * hapfreq[0] / (hapfreq[0]+hapfreq[1]);
			temp1 = RAND_MAX * hapfreq[2] / (hapfreq[2]+hapfreq[3]);
			for(j = 0; j < UNID; j++)
			{
				if(pop[j].qtl[i].a[0])
				{
					pop[j].qtl[i].alink[0] = rand() < temp1 ? 0 : 1;
				} else {
					pop[j].qtl[i].alink[0] = rand() < temp0 ? 0 : 1;
				}
				if(pop[j].qtl[i].a[1])
				{
					pop[j].qtl[i].alink[1] = rand() < temp1 ? 0 : 1;
				} else {
					pop[j].qtl[i].alink[1] = rand() < temp0 ? 0 : 1;
				}
			}
		} else {
			for(j = 0; j < UNID; j++)
			{
				pop[j].qtl[i].alink[0] = pop[j].qtl[i].a[0];
				pop[j].qtl[i].alink[1] = pop[j].qtl[i].a[1];
			}
		}
		
		// locus 2
		if(Rsq[i][0] != 1)
		{
			calchapfreq(maf[1],Rsq[i][1],hapfreq);
			temp0 = (float)RAND_MAX * hapfreq[0] / (hapfreq[0]+hapfreq[1]);
			temp1 = (float)RAND_MAX * hapfreq[2] / (hapfreq[2]+hapfreq[3]);
			for(j = 0; j < UNID; j++)
			{
				if(pop[j].qtl[i].b[0])
				{
					pop[j].qtl[i].blink[0] = rand() < temp1 ? 0 : 1;
				} else {
					pop[j].qtl[i].blink[0] = rand() < temp0 ? 0 : 1;
				}
				temp = rand();
				if(pop[j].qtl[i].b[1])
				{
					pop[j].qtl[i].blink[1] = rand() < temp1 ? 0 : 1;
				} else {
					pop[j].qtl[i].blink[1] = rand() < temp0 ? 0 : 1;
				}
			}
		} else {
			for(j = 0; j < UNID; j++)
			{
				pop[j].qtl[i].blink[0] = pop[j].qtl[i].b[0];
				pop[j].qtl[i].blink[1] = pop[j].qtl[i].b[1];
			}
		}			
		
		// X
		for(j = 0; j < UNID; j++)
		{
			pop[j].qtl[i].Xlink = pop[j].qtl[i].alink[0]+pop[j].qtl[i].alink[1]+3*(pop[j].qtl[i].blink[0]+pop[j].qtl[i].blink[1]);
		}
	}
}

void calcld(record *pop)
{
	int i,j,k;
	float hapfreq[4], maf[2], Rsq;
	
	for(i = 0; i < UNQTL; i++)
	{
		maf[0] = 0; maf[1] = 0;
		for(j = 0; j < 4; j++) hapfreq[j] = 0;
		for(j = 0; j < UNID; j++)
		{
			maf[0] += (pop[j].qtl[i].a[0] + pop[j].qtl[i].a[1]);
			maf[1] += (pop[j].qtl[i].alink[0] + pop[j].qtl[i].alink[1]);
			k = pop[j].qtl[i].a[0] + 2 * pop[j].qtl[i].alink[0];
			hapfreq[k]++;
			k = pop[j].qtl[i].a[1] + 2 * pop[j].qtl[i].alink[1];
			hapfreq[k]++;
		}
		maf[0] /= (UNID*2);
		maf[1] /= (UNID*2);
		for(j = 0; j < 4; j++) hapfreq[j] /= (2*UNID);
		
		Rsq = (hapfreq[0] - (1-maf[0])*(1-maf[1]))*(hapfreq[0] - (1-maf[0])*(1-maf[1])) / (maf[0]*maf[1]*(1-maf[0])*(1-maf[1]));
		printf("a - %f - %f %f\n",Rsq,maf[0],maf[1]);
		
		maf[0] = 0; maf[1] = 0;
		for(j = 0; j < 4; j++) hapfreq[j] = 0;
		for(j = 0; j < UNID; j++)
		{
			maf[0] += (pop[j].qtl[i].b[0] + pop[j].qtl[i].b[1]);
			maf[1] += (pop[j].qtl[i].blink[0] + pop[j].qtl[i].blink[1]);
			k = pop[j].qtl[i].b[0] + 2 * pop[j].qtl[i].blink[0];
			hapfreq[k]++;
			k = pop[j].qtl[i].b[1] + 2 * pop[j].qtl[i].blink[1];
			hapfreq[k]++;
		}
		maf[0] /= (UNID*2);
		maf[1] /= (UNID*2);
		for(j = 0; j < 4; j++) hapfreq[j] /= (2*UNID);
		
		Rsq = (hapfreq[0] - (1-maf[0])*(1-maf[1]))*(hapfreq[0] - (1-maf[0])*(1-maf[1])) / (maf[0]*maf[1]*(1-maf[0])*(1-maf[1]));
		printf("b - %f - %f %f\n",Rsq,maf[0],maf[1]);
	}
}
	

void wam(int n, double *p, int *a, int nans, int *ans)
{
	double *q, rU;
	int i,j,k;
	int *HL,*H,*L;
	HL = calloc(n,sizeof(int));
	q = calloc(n,sizeof(double));
	H = HL - 1; L = HL + n;
	double sum = 0;
	for(i = 0; i < n; i++)
	{
		sum += p[i];
	}
	for(i = 0; i < n; i++)
	{
		p[i] /= sum;
	}
	for(i = 0; i < n; i++)
	{
		q[i] = p[i] * n;
		if(q[i] < 1.) *++H = i; else *--L = i;
	}
	if(H >= HL && L < HL +n)
	{
		for(k = 0; k < n-1; k++)
		{
			i = HL[k];
			j = *L;
			a[i] = j;
			q[j] += q[i] - 1;
			if(q[j] < 1.) L++;
			if(L >= HL + n) break;
		}
	}
	for(i = 0; i < n; i++) q[i] += i;
	for(i = 0; i < nans; i++)
	{
		rU = (double) rand() / RAND_MAX * n;
		k = (int) rU;
		ans[i] = (rU < q[k]) ? k : a[k];
	}
	free(HL);
	free(q);
}

// this requires minimum size of 20 columns to work. I don't know why
/*void **get_space(int rows, int cols, size_t size, size_t psize)
{
	int i;
	void *p, **a;
	p = malloc(rows * cols * size);
	a = malloc(rows * psize);
	for(i = 0; i < rows; i++)
	{
		a[i] = p + (i * cols);
	}
	return a;
}*/

// Box-Muller method approximation. mean = 0; var = 1
// Marsaglia polar method may be faster
void rnorm(float *x, int n)
{
	int
		i;
	float
		tmp1f, tmp2f,
		x1,x2;

	for(i = 0; i < n; i++)
	{
		x1 = (float)rand() / RAND_MAX;
		x2 = (float)rand() / RAND_MAX;
		tmp1f = sqrtf(-2 * logf(x1));
		tmp2f = cosf(2*PI*x2);
		x[i] = tmp1f * tmp2f;
	}
}

// if 0 <= a < b <= 1, then p random uniform sample lies between a and b is b-a;
void sample2vals(float freq, int N, char *destination)
{
	int i, thresh = (int) RAND_MAX * freq;
	for(i = 0; i < N; i++)
	{
		destination[i] = (rand() > thresh) ? 1 : 0;
	}
}

void addnoise(record *pop, int N, float vare)
{
	float *temp = malloc(sizeof(float) * N);
	rnorm(temp, N);
	int i;
	float min = 0;
	for(i = 0; i < N; i++)
	{
//		pop[i].phen += (fabsf(temp[i]) * sqrt(vare));
		pop[i].phen += (temp[i] * sqrt(vare));
		if(min > pop[i].phen) min = pop[i].phen;
	}
	for(i = 0; i < N; i++)
	{
		pop[i].phen -= min;
		//printf("%f ",pop[i].phen);fflush(stdout);
	}
	free(temp);
}

void tables(record *pop, int N, int nqtl, float **pattern)
{
	int i,j,k;
	float gmeans[9];
	int gcount[9];
	for(i = 0; i < nqtl; i++)
	{
		for(j=0;j<9;j++){ gmeans[j] = 0; gcount[j] = 0; }
		for(j = 0; j < N; j++)
		{
			for(k = 0; k < 9; k++)
			{
				if(pop[j].qtl[i].X == k)
				{
					gmeans[k] += pop[j].phen;
					gcount[k]++;
					break;
				}
			}
		}
		for(j = 0; j < 9; j++) printf("%0.3f\t",pattern[i][j]);
		printf("\n");
		for(j = 0; j < 9; j++)
		{
			if(gcount[j] > 0)
				gmeans[j] /= gcount[j];
			printf("%0.3f\t",gmeans[j]);
		}
		printf("\n\n");
	}
}


float basepop(int N, int nqtl, float **p, float **pattern, float hsq, record *pop)
{
	int i, j, k;
	char *values = malloc(sizeof(char)*2);
	values[0] = 0; values[1] = 1;
	char *sample = malloc(sizeof(char)*N);
	float *freq = malloc(sizeof(float)*2);
	for(i = 0; i < N; i++)
	{
		pop[i].sex = ((i % 2) == 0) ? 0 : 1;
		pop[i].phen = 0;
	}
	for(i = 0; i < nqtl; i++)
	{

		if(p[i][0] != 0 && p[i][0] != 1)
		{
			sample2vals(p[i][0],N,sample);
			for(j = 0; j < N; j++) pop[j].qtl[i].a[0] = sample[j];
			sample2vals(p[i][0],N,sample);
			for(j = 0; j < N; j++) pop[j].qtl[i].a[1] = sample[j];
		} else {
			// mutant
			for(j = 0; j < N; j++) pop[j].qtl[i].a[0] = p[i][0];
			for(j = 0; j < N; j++) pop[j].qtl[i].a[1] = p[i][0];
			j = (int)rand() % N;
			pop[j].qtl[i].a[1] = 1-p[i][0];
		}

		if(p[i][1] != 0 && p[i][1] != 1)
		{
			sample2vals(p[i][1],N,sample);
			for(j = 0; j < N; j++) pop[j].qtl[i].b[0] = sample[j];
			sample2vals(p[i][1],N,sample);
			for(j = 0; j < N; j++) pop[j].qtl[i].b[1] = sample[j];
		} else {
			for(j = 0; j < N; j++) pop[j].qtl[i].b[0] = p[i][1];
			for(j = 0; j < N; j++) pop[j].qtl[i].b[1] = p[i][1];
			j = (int)rand() % N;
			pop[j].qtl[i].b[1] = 1-p[i][1];
		}
		
		for(j = 0; j < N; j++)
		{
			pop[j].qtl[i].X = pop[j].qtl[i].a[0] + pop[j].qtl[i].a[1] + 3*(pop[j].qtl[i].b[0] + pop[j].qtl[i].b[1]);
		}
	}

	for(i = 0; i < nqtl; i++)
	{
		for(j = 0; j < N; j++)
		{
			for(k = 0; k < 9; k++)
			{
				if (pop[j].qtl[i].X == k)
					pop[j].phen += pattern[i][k];
			}
		}
	}
	float varg = 0, mean = 0, vare;
	for(i = 0; i < N; i++)
	{
		mean += pop[i].phen;
	}
	mean /= N;
	for(i = 0; i < N; i++)
	{
		varg += (pop[i].phen - mean) * (pop[i].phen - mean);
	}
	varg /= N;
	vare = varg == 0 ? 1 : (1 / hsq) * varg - varg;
	if(vare > 0)
	{
		addnoise(pop,N,vare);
	}
	// normalise and scale
//	float sd = sqrt(varg+vare);
//	for(i = 0; i < N; i++)
//	{
//		pop[i].phen = (pop[i].phen - mean) / sd;
//	}
//	vare /= sd;
	free(values);
	free(sample);
	free(freq);
	
	makeld(pop,URSQ);
	
	return vare;
}

void nextgen(int N, int nqtl, float **pattern, float vare, record *pop)
{
	int *maleso = malloc(sizeof(int)*N/2);
	int *males = malloc(sizeof(int)*N/2);
	int *selectedmales = malloc(sizeof(int)*N/4);
	double *mphen = malloc(sizeof(double)*N/2);
	int *femaleso = malloc(sizeof(int)*N/2);
	int *females = malloc(sizeof(int)*N/2);
	int *selectedfemales = malloc(sizeof(int)*N/4);
	double *fphen = malloc(sizeof(double)*N/2);
	int i,j,k,l,gam, m = 0, f = 0;

	// select for the next generation
	for(i = 0; i < N; i++)
	{
		if(pop[i].sex == 0)
		{
			maleso[m] = i;
			males[m] = i;
			mphen[m++] = (double)pop[i].phen;
		} else {
			femaleso[f] = i;
			females[f] = i;
			fphen[f++] = (double)pop[i].phen;
		}
	}
	wam(f, fphen, females, f/2, selectedfemales);
	wam(m, mphen, males, m/2, selectedmales);
	for(i = 0; i < N/4; i++)
	{
		selectedmales[i] = maleso[(selectedmales[i])];
		selectedfemales[i] = femaleso[(selectedfemales[i])];
	}

	// Create new temporary structure
	record *pop2 = malloc(N * sizeof(record));
	genotype *block = malloc(N*nqtl*sizeof(genotype));
	for(i = 0; i < N; i++)
	{
		pop2[i].qtl = block + i*nqtl;
	}

	// sample gametes for each new child at each locus
	for(i = 0; i < nqtl; i++)
	{
		l = 0;
		for(j = 0; j < N/4; j++) // each couple has 4 children, 2 male 2 female
		{
			m = selectedmales[j];
			f = selectedfemales[j];
			for(k = 0; k < 4; k++)
			{
				gam = rand() % 2;
				pop2[l].qtl[i].a[0] = pop[m].qtl[i].a[gam];
				gam = rand() % 2;
				pop2[l].qtl[i].a[1] = pop[f].qtl[i].a[gam];
				gam = rand() % 2;
				pop2[l].qtl[i].b[0] = pop[m].qtl[i].b[gam];
				gam = rand() % 2;
				pop2[l].qtl[i].b[1] = pop[f].qtl[i].b[gam];

				pop2[l].qtl[i].X = 
					pop2[l].qtl[i].a[0] + 
					pop2[l].qtl[i].a[1] + 3*(
					pop2[l].qtl[i].b[0] + 
					pop2[l].qtl[i].b[1]);
				l++;
			}
		}
	}

	// phenotype
	for(i = 0; i < N; i++)
	{
		pop[i].phen = 0;
	}
	for(i = 0; i < nqtl; i++)
	{
		for(j = 0; j < N; j++)
		{
			for(k = 0; k < 9; k++)
			{
				if (pop2[j].qtl[i].X == k)
					pop[j].phen += pattern[i][k];
			}
		}
	}
	if(vare > 0)
	{
		addnoise(pop,N,vare);
	}

	// free original genotypes and make pop2 point to new block
	genotype *pointer = pop[0].qtl;
	free(pointer);	
	for(i = 0; i < N; i++)
	{
		pop[i].qtl = block + i*nqtl;
	}
	free(pop2);
	free(maleso);
	free(males);
	free(selectedmales);
	free(mphen);
	free(femaleso);
	free(females);
	free(selectedfemales);
	free(fphen);
	
	// now generate loci in LD with causal SNPs
	makeld(pop,URSQ);
	// calcld(pop);

}

float calcstatsnoia(record *pop, int N, int nqtl, allstats *stat)
{
	int i,n;
	for(i = 0; i < nqtl; i++)
	{
		noiapar(pop,N,i,stat,0,0);
		n = N * UPROPSAMPLE;
		ftestlin(pop,n,i,stat,0);
		ftest8df(pop,n,i,stat,0);
		ftest2df(pop,n,i,stat,0);
		if((URSQ[i][0]+URSQ[i][1]) < 2.0)
		{
			noiapar(pop,N,i,stat,0,1);
			ftestlin(pop,n,i,stat,1);
			ftest8df(pop,n,i,stat,1);
			ftest2df(pop,n,i,stat,1);
		} else {
			stat[i].report[1] = stat[i].report[0];
		}
			
		if(UPROPSAMPLE != 0)
		{
			noiapar(pop,n,i,stat,1,0);
			if((URSQ[i][0]+URSQ[i][1]) < 2.0) noiapar(pop,n,i,stat,1,1);
		}
	}

	float var = 0, mean = 0;
	for(i = 0; i < N; i++)
	{
		mean += pop[i].phen;
	}
	mean /= N;
	for(i = 0; i < N; i++)
	{
		var += (pop[i].phen - mean) * (pop[i].phen - mean);
	}
	var /= N;
	return(var);
}

/*
void printstats(stats *stat, int nqtl)
{
	int i;
	printf("QTL\tA\tB\tRsq\tVt\tVg\tVe\taf1\taf2\n");
	for(i = 0; i < nqtl; i++)
	{
//		printf("%d:\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\n",i,stat[i].A,stat[i].B,stat[i].Rsq,stat[i].Vt,stat[i].Vg,stat[i].Ve,stat[i].maf[0],stat[i].maf[1]);
	}
}
*/

void printpop(record *pop, int N, int nqtl)
{
	int i,j;
	for(i = 0; i < N; i++)
	{
		printf("%d: %d %f\t",i,pop[i].sex,pop[i].phen);
		for(j = 0; j < nqtl; j++)
		{
			printf("%d ",pop[i].qtl[j].X);
		}
		printf("\n");
	}
}

void printpattern(float **pattern, float **p, int nqtl)
{
	int i,j;
	for(i = 0; i < nqtl; i++)
	{
		printf("%0.3f\t%0.3f\t",p[i][0],p[i][1]);
		for(j = 0; j < 9; j++)
		{
			printf("%0.3f ",pattern[i][j]);
		}
		printf("\n");
	}
}

// if any allele goes to fixation return 1, else return 0
int checkfixation1(allstats *stat, int nqtl, int repeat)
{
	int i, check = 0;
	for(i = 0; i < nqtl; i++)
	{
	//	printf(" - %f %f\n",stat[i].maf[0],stat[i].maf[1]);
		check = (stat[i].report[0].maf[0] == 0 || stat[i].report[0].maf[0] == 1) ? 1 : 0;
		check += (stat[i].report[0].maf[1] == 0 || stat[i].report[0].maf[1] == 1) ? 1 : 0;
		if(check) return 1;
	}
	return 0;
}
// if any allele is still unfixed return 1, else return 0
int checkfixation(allstats *stat, int nqtl)
{
	int i, check = 0;
	for(i = 0; i < nqtl; i++)
	{
		check = (stat[i].report[0].maf[0] == 0 || stat[i].report[0].maf[0] == 1) ? 1 : 0;
		check += (stat[i].report[0].maf[1] == 0 || stat[i].report[0].maf[1] == 1) ? 1 : 0;
		if(check) return 1;
	}
	return 0;
}

void countgenerations(allstats *stat, int nqtl, int generations, int *ngen)
{
	int i,j,k;
	for(i = 0; i < nqtl; i++)
	{
		for(j = 0; j < generations; j++)
		{
			k = nqtl*j + i;
			if(stat[k].report[0].maf[0] == 0 || stat[k].report[0].maf[0] == 1) break;
			if(stat[k].report[0].maf[1] == 0 || stat[k].report[0].maf[1] == 1) break;
			ngen[i]++;
		}
	}
}

int runsimul(record *pop, int N, int nqtl, float hsq, int *pat, int rep, float **pattern, float vare, int generations, allstats *stat, FILE *out, int condition)
{
	float *var = malloc(sizeof(float)*generations);
//	float chipsize = 300000;
//	float thresh1 = -log10f(0.05/chipsize);
//	float thresh2 = -log10f(0.05/(chipsize*(chipsize-1)*0.5));
	var[0] = calcstatsnoia(pop, N, nqtl, stat);

	int ad, i, j, k;
	for(i = 1; i < generations; i++)
	{
	//	if((i % 10) == 0) printf("Generation %d\n",i);
		ad = i*nqtl;
		nextgen(N, nqtl, pattern, vare, pop);
		var[i] = calcstatsnoia(pop, N, nqtl, stat+ad);
		if(checkfixation1(stat+ad,nqtl,rep)) break;
	}
	int *ngen = malloc(sizeof(int)*nqtl);
	for(j = 0; j < nqtl; j++) ngen[j] = 0;
	countgenerations(stat,nqtl,generations,ngen);


// Condition Repeat Pattern Generation RSQ1 RSQ2 MAF1 MAF2 Va1 Va2 Vd1 Vd2 Vaa Vad Vda Vdd Vt pa1 pa2 pad1 pad2 pmarg pepi pfull
	for(j = 0; j < nqtl; j++)
	{
		for(k = 0; k < ngen[j]; k+=UINCREMENT)
		{
			fprintf(out,"%d %d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f ",
				condition,rep,pat[j],k,URSQ[j][0],URSQ[j][1],
				stat[(j+k*nqtl)].report[0].maf[0],
				stat[(j+k*nqtl)].report[0].maf[1],
				stat[(j+k*nqtl)].report[0].Vcomp[0],
				stat[(j+k*nqtl)].report[0].Vcomp[2],
				stat[(j+k*nqtl)].report[0].Vcomp[1],
				stat[(j+k*nqtl)].report[0].Vcomp[5],
				stat[(j+k*nqtl)].report[0].Vcomp[3],
				stat[(j+k*nqtl)].report[0].Vcomp[6],
				stat[(j+k*nqtl)].report[0].Vcomp[4],
				stat[(j+k*nqtl)].report[0].Vcomp[7],
				var[k],
				stat[(j+k*nqtl)].report[0].F[0],
				stat[(j+k*nqtl)].report[0].F[1],
				stat[(j+k*nqtl)].report[0].F2[0],
				stat[(j+k*nqtl)].report[0].F2[1],
				stat[(j+k*nqtl)].report[0].Fmarg,
				stat[(j+k*nqtl)].report[0].Fint,
				stat[(j+k*nqtl)].report[0].Ffull);
			fprintf(out,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
				stat[(j+k*nqtl)].report[1].Vcomp[0],
				stat[(j+k*nqtl)].report[1].Vcomp[2],
				stat[(j+k*nqtl)].report[1].Vcomp[1],
				stat[(j+k*nqtl)].report[1].Vcomp[5],
				stat[(j+k*nqtl)].report[1].Vcomp[3],
				stat[(j+k*nqtl)].report[1].Vcomp[6],
				stat[(j+k*nqtl)].report[1].Vcomp[4],
				stat[(j+k*nqtl)].report[1].Vcomp[7],
				var[k],
				stat[(j+k*nqtl)].report[1].F[0],
				stat[(j+k*nqtl)].report[1].F[1],
				stat[(j+k*nqtl)].report[1].F2[0],
				stat[(j+k*nqtl)].report[1].F2[1],
				stat[(j+k*nqtl)].report[1].Fmarg,
				stat[(j+k*nqtl)].report[1].Fint,
				stat[(j+k*nqtl)].report[1].Ffull);
		}
	}

	free(ngen);
	free(var);
	return i;
}

void ShowError()
{
	char *dlError = dlerror();
	if(dlError) printf("Error with dl function: %s\n", dlError);
}


/*
Parameters:
N
Nsampled for tests
nqtl
for each qtl:
	pattern
	effect
	p1 - if p1 = 0 or 1 then simulate mutant
	p2 - if p2 = 0 or 1 then simulate mutant
	rsq1
	rsq2
hsq
generations
repeats
filename
record increment
seed
condition
*/

void inputsummary(char *inputfile)
{
	int i, j;

	FILE *in = fopen(inputfile, "r");
	
	fscanf(in,"%d\n",&NENUM);
	allpatterns = get_space_f(NENUM,9);
	for(i = 0; i < NENUM; i++)
	{
		for(j = 0; j < 9; j++)
		{
			fscanf(in,"%f ",&allpatterns[i][j]);
		}
		fscanf(in,"\n");
	}
	fscanf(in,"%d\n",&UNID);
	fscanf(in,"%f\n",&UPROPSAMPLE);
	fscanf(in,"%d\n",&UNQTL);
	UPAT = malloc(sizeof(int) * UNQTL);
	UEFFECT = malloc(sizeof(float) * UNQTL);
	UFREQ = get_space_f(UNQTL, 2);
	URSQ = get_space_f(UNQTL, 2);
	for(i = 0; i < UNQTL; i++)
	{
		fscanf(in,"%d %f %f %f %f %f\n",&UPAT[i],&UEFFECT[i],&UFREQ[i][0],&UFREQ[i][1],&URSQ[i][0],&URSQ[i][1]);
	}
	fscanf(in,"%f\n",&UHSQ);
	fscanf(in,"%d\n",&UGENERATIONS);
	fscanf(in,"%d\n",&UREPEATS);
	fscanf(in,"%s\n",UFILENAME);
	fscanf(in,"%d\n",&UINCREMENT);
	fscanf(in,"%d\n",&USEED);
	USEED ? srand(USEED) : srand(time(NULL));
	fscanf(in,"%d\n",&UCONDITION);


	for(i = 0; i < NENUM; i++)
	{
		for(j = 0; j < 9; j++)
		{
			printf("%f ",allpatterns[i][j]);
		}
		printf("\n");
	}

	printf("NID = %d\n",UNID);
	printf("PROPSAMPLE = %f\n",UPROPSAMPLE);
	printf("nQTL = %d\n\n",UNQTL);
	for(i = 0; i < UNQTL; i++)
	{
		printf("QTL 1:\n  Pattern: %d\n  Effect: %f\n  Freq: %f %f\n  Rsq: %f %f\n",UPAT[i],UEFFECT[i],UFREQ[i][0],UFREQ[i][1],URSQ[i][0],URSQ[i][1]);
	}
	
	printf("\nHsq = %f\nGenerations = %d\nRepeats = %d\nFilename = %s\nRec incr = %d\nSeed = %d\nCondition = %d\n", UHSQ, UGENERATIONS, UREPEATS, UFILENAME, UINCREMENT, USEED, UCONDITION);
}
	
	

int main(int argc, char **argv)
{
	if(argc != 2)
	{
		printf("\n%s input.dat\n\n",argv[0]);
		exit(1);
	}

	int i, j;
	allstats *stat;
	float *var;
	FILE *out;

	inputsummary(argv[1]);

	sofile = dlopen("fprob.so",RTLD_LAZY);
	ShowError();
	fprob = dlsym(sofile,"fprob_");
	ShowError();

	stat = (allstats *)malloc(sizeof(allstats) * UNQTL * UGENERATIONS);
	var = (float *)malloc(sizeof(float) * UGENERATIONS);

	// create population structures
	record *pop = malloc(UNID * sizeof(record));
	genotype *block = malloc(UNID*UNQTL*sizeof(genotype));
	for(i = 0; i < UNID; i++)
	{
		pop[i].qtl = block + i*UNQTL;
	}

	// create QTL frequencies and patterns
	float **pattern = get_space_f(UNQTL,9);
	for(i = 0; i < UNQTL; i++)
	{
		for(j = 0; j < 9; j++)
		{
			//pattern[i][j] = 1 - allpatterns[UPAT[i]][j] * (1 - UEFFECT[i]);
			pattern[i][j] = allpatterns[UPAT[i]][j] * UEFFECT[i];
		}
	}
	// generate initial population
	char line1;
	out = fopen(UFILENAME,"r");
	if(out)
	{
		fscanf(out,"%c\t",&line1);
		fclose(out);
		if(line1 == 'C')
		{	
			out = fopen(UFILENAME,"a");
		}
	} else {
		out = fopen(UFILENAME,"w");
		fprintf(out,"Condition Repeat Pattern Generation Rsq1 Rsq2 maf1 maf2 Va1 Va2 Vd1 Vd2 Vaa Vad Vda Vdd Vt pa1 pa2 pad1 pad2 pmarg pint pfull Va1l Va2l Vd1l Vd2l Vaal Vadl Vdal Vddl Vtl pa1l pa2l pad1l pad2l pmargl pintl pfulll\n");
	}


	float vare;
	int *totgen = malloc(sizeof(int) * UREPEATS);
	for(i = 0; i < UREPEATS; i++)
	{
//		vare = basepopmutant(UNID,UNQTL,UFREQ,pattern,UHSQ,pop);
		vare = basepop(UNID,UNQTL,UFREQ,pattern,UHSQ,pop);
//		tables(pop,N,nqtl,pattern);
		totgen[i] = runsimul(pop, UNID, UNQTL, UHSQ, UPAT, i, pattern, vare, UGENERATIONS, stat, out, UCONDITION);
		printf("Repeat\t%d:\t%d\tgenerations\n",i,totgen[i]);
	}


	fclose(out);
	block = pop[0].qtl;
	free(block);
	free(pop);
	free(totgen);
	free(UPAT);
	free(UEFFECT);
	release_space_f(pattern);
	release_space_f(UFREQ);
	release_space_f(URSQ);
	free(stat);
	free(var);
	dlclose(sofile);
	ShowError();

	return 0;
}

