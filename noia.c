#include "pop.h"

void noiapar(record *pop, int m, int loc, allstats *stat, char test, int linked)
{
    int
		n = 2,
		i,j,k,ii,
		ngen,size,
		*X;

    float
		***SsM,
		***SsiM,
		**Ss,**Ssi,
		**pM, **TEMP, **TEMPi,
		*p, *g, mY,temp1,temp2,
		*yhat;

	ngen = (int)pow(3,(float)n); // always 9

	pM = get_space_f(n,3);
	Ss = get_space_f(ngen,ngen);
	Ssi = get_space_f(ngen,ngen);
	TEMP = get_space_f(ngen,ngen);
	TEMPi = get_space_f(ngen,ngen);
	
	SsM = make_cube_d(n,3,3);
	SsiM = make_cube_d(n,3,3);
	
	X = malloc(m * sizeof(int));
	p = malloc(ngen * sizeof(float));	
	g = malloc(ngen * sizeof(float));	
	
	// initialise
	for(i = 0; i < ngen; i++)
	{
		p[i] = 0;
		g[i] = 0;
	}
	for(i = 0; i < n; i++)
		for(j = 0; j < 3; j++)
			pM[i][j] = 0;
	
	// find frequencies and means
	mY = 0;
	if(linked)
	{
		for(i = 0; i < m; i++)
		{
			X[i] = pop[i].qtl[loc].Xlink;
			k = pop[i].qtl[loc].alink[0] + pop[i].qtl[loc].alink[1];
			(pM[0][k])++;
			k = pop[i].qtl[loc].blink[0] + pop[i].qtl[loc].blink[1];
			(pM[1][k])++;
			p[X[i]]++;
			g[X[i]] += pop[i].phen;
			mY += pop[i].phen;
		}
	} else {
		for(i = 0; i < m; i++)
		{
			X[i] = pop[i].qtl[loc].X;
			k = pop[i].qtl[loc].a[0] + pop[i].qtl[loc].a[1];
			(pM[0][k])++;
			k = pop[i].qtl[loc].b[0] + pop[i].qtl[loc].b[1];
			(pM[1][k])++;
			p[X[i]]++;
			g[X[i]] += pop[i].phen;
			mY += pop[i].phen;
		}
	}
	
	mY /= m;
	for(i = 0; i < ngen; i++)
	{
		if(p[i] > 0) g[i] /= p[i];
		p[i] /= m;
	}
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < 3; j++)
		{
			pM[i][j] /= m;
		}
	}

	temp1 = pM[0][2] + pM[0][1] / 2;
	temp2 = pM[1][2] + pM[1][1] / 2;

	if(!test)
	{
		stat[loc].report[linked].maf[0] = temp1;
		stat[loc].report[linked].maf[1] = temp2;

		for(i = 0; i < n; i++)
		{
			if(temp1 == 1 || temp1 == 0) return;
			if(temp2 == 1 || temp2 == 0) return;
		}
	}

	// Parameterise according to NOIA
	// This should make each variance component orthogonal
	for(i = 0; i < n; i++)
	{
		make_ss(SsM[i],pM[i]);
		make_ssi(SsiM[i],pM[i]);
	}
	// general kronecker
	// need scratch space to store intermediates
	// need general space to store intermediate/final answers and use for next steps
	// initialise general with final locus, work backwards filling scratch and then replacing general with scratch

	ii = n-1;
	for(i = 0; i < (n-1); i++)
	{
		size = (int)pow(3,(float)(i+1));
		if(i == 0)
		{
			for(j=0;j<3;j++)
			{
				for(k=0;k<3;k++)
				{
					TEMP[j][k]=SsM[ii][j][k];
					TEMPi[j][k]=SsiM[ii][j][k];
				}
			}
		} else {
			for(j=0;j<size;j++)
			{
				for(k=0;k<size;k++)
				{
					TEMP[j][k]=Ss[j][k];
					TEMPi[j][k]=Ssi[j][k];
				}
			}
		}
//		printf("%d : in... ",i);
		kronecker(TEMPi,SsiM[ii-1],Ssi,size,3,size,3);
		kronecker(TEMP,SsM[ii-1],Ss,size,3,size,3);
//		printf("out, %d\n",ii);
		ii--;
	}

/*	for(i = 0; i < ngen; i++)
	{
		for(j = 0; j < ngen; j++)
		{
			printf("%f ",Ss[i][j]);
		}
		printf("\n");
	}
	printf("\n");
*/
	float *e = malloc(sizeof(float)*ngen);
	solve(Ssi,e,g,ngen);
	
	// calculate variances for each variance component
	float *vars = malloc(sizeof(float)*ngen);
	for(i = 0; i < 9; i++) vars[i] = 0;
	float **Sstemp = get_space_f(ngen,ngen);
	for(i = 0; i < ngen; i++) // rows
	{
		for(j = 0; j < ngen; j++) // cols
		{
			Sstemp[i][j] = Ss[i][j] * Ss[i][j] * p[i];
			vars[j] += Sstemp[i][j];
		}
	}
	for(i = 0; i < ngen; i++)
	{
		vars[i] *= (e[i]*e[i]);
	}
	vars[0] = 0;

	// .. a. d. .a aa da .d ad dd 
	if(!test)
	{
		stat[loc].report[linked].A = e[1];
		stat[loc].report[linked].B = e[3];

		// genetic variances
		for(i = 1; i < ngen; i++)
		{
			stat[loc].report[linked].Vcomp[i-1] = vars[i];
		}

	}
	// .. a. d. .a aa da .d ad dd
	// F test for epistatic components

	yhat = (float *)malloc(sizeof(float) * m);
	int nparameter = 0;
	int intpar[5] = {0,4,5,7,8};
//	int intpar[5] = {0,1,2,3,6};
//	int intpar[9] = {0,1,2,3,4,5,6,7,8};
	int ngen4 = 5;
	for(i = 0; i < ngen4; i++)
	{
		if(e[intpar[i]] != 0) nparameter++;
	}
	if(nparameter < 2)
	{
		stat[loc].report[linked].Fint = 0;
		return;
	}

	for(i = 0; i < m; i++)
	{
		yhat[i] = 0;
		for(j = 0; j < ngen4; j++)
		{
			yhat[i] += Ss[X[i]][intpar[j]] * e[intpar[j]];
		}
	}

	float SST = 0, SSE = 0, SSM = 0;
	for(i = 0; i < m; i++)
	{
		SST += (pop[i].phen - mY) * (pop[i].phen - mY);
		SSE += (pop[i].phen - yhat[i]) * (pop[i].phen - yhat[i]);
		SSM += (yhat[i] - mY) * (yhat[i] - mY);
	}

	int dfm = nparameter - 1, dfe = m - nparameter, err;
	float MSM = SSM / dfm;
	float MSE = SSE / dfe;
	double F = (double)MSM/MSE,pval,qval;
	if(isnan(F)) F = 1;
	(*fprob)(&dfm,&dfe,&F,&pval,&qval,&err);
	stat[loc].report[linked].Fint = (float)-log10(qval);
	

	// F test for non epistatic components
	nparameter = 0;
	int intpar2[5] = {0,1,2,3,6};
	for(i = 0; i < ngen4; i++)
	{
		if(e[intpar2[i]] != 0) nparameter++;
	}
	if(nparameter < 2)
	{
		stat[loc].report[linked].Fmarg = 0;
		return;
	}

	for(i = 0; i < m; i++)
	{
		yhat[i] = 0;
		for(j = 0; j < ngen4; j++)
		{
			yhat[i] += Ss[X[i]][intpar2[j]] * e[intpar2[j]];
		}
	}

	SST = 0; SSE = 0; SSM = 0;
	for(i = 0; i < m; i++)
	{
		SST += (pop[i].phen - mY) * (pop[i].phen - mY);
		SSE += (pop[i].phen - yhat[i]) * (pop[i].phen - yhat[i]);
		SSM += (yhat[i] - mY) * (yhat[i] - mY);
	}

	dfm = nparameter - 1; dfe = m - nparameter;
	MSM = SSM / dfm;
	MSE = SSE / dfe;
	F = (double)MSM/MSE,pval,qval;
	if(isnan(F)) F = 1;
	(*fprob)(&dfm,&dfe,&F,&pval,&qval,&err);
	stat[loc].report[linked].Fmarg = (float)-log10(qval);
	
	release_space_f(Ss);
	release_space_f(Ssi);
	release_space_f(TEMP);
	release_space_f(TEMPi);
	release_space_f(pM);
	release_space_f(Sstemp);
	
	free_cube_d(SsM);
	free_cube_d(SsiM);

	free(vars);
	free(yhat);
	free(p);
	free(g);
	free(e);
	free(X);
}
// 1 df F values are nonsense - they act as independent tests. However, for full (e.g. 2df 1 locus or 8df 2 loci) tests, the F value results are identical to full 8df tests in anova(lm) in R
// The sum of variances of each individual parameter doesn't sum the total variance explained because everything tested has been correlated, so departure from the assumption of orthogonality.
// Speed - ~200 times faster that R/noia, ~20 times faster than R/anova
	
int **get_space(int m, int n)
{
	int i, *p, **a;
	p = malloc(m * n * sizeof(int));
	a = malloc(m * sizeof(int *));
	for(i = 0; i < m; i++) a[i] = p + (i * n);
	return a;
}

void release_space(int** X)
{
	int * p;
	p = (int *) X[0];
	free(p);
	free(X);
}
	
float **get_space_f(int m, int n)
{
	int i;
	float *p, **a;
	p = malloc(m * n * sizeof(float));
	a = malloc(m * sizeof(float *));
	for(i = 0; i < m; i++) a[i] = p + (i * n);
	return a;
}

void release_space_f(float **X)
{
	float * p;
	p = (float *) X[0];
	free(p);
	free(X);
}
	
int ***make_cube(int x, int y, int z)
{
	int i,j;
	int *p, **p1, ***a;
	p = malloc(x * y * z * sizeof(int));
	p1 = malloc(x * y * sizeof(int *));
	a = malloc(x * sizeof(int **));
	
	for(i = 0; i < x; i++)
	{
		a[i] = p1 + i*y;
		for(j = 0; j < y; j++)
		{
			a[i][j] = p + i*y*z + j*y;
		}
	}
	return(a);
}
	
void free_cube(int ***array)
{
	int **p1, *p2;
	
	p1 = (int **) array[0];
	p2 = (int *) array[0][0];
	free(p2);
	free(p1);
	free(array);
}
	
float ***make_cube_d(int x, int y, int z)
{
	int i,j;
	float *p, **p1, ***a;
	p = malloc(x * y * z * sizeof(float));
	p1 = malloc(x * y * sizeof(float *));
	a = malloc(x * sizeof(float **));
	
	for(i = 0; i < x; i++)
	{
		a[i] = p1 + i*y;
		for(j = 0; j < y; j++)
		{
			a[i][j] = p + i*y*z + j*y;
		}
	}
	return(a);
}
	
void free_cube_d(float ***array)
{
	float **p1, *p2;
	
	p1 = (float **) array[0];
	p2 = (float *) array[0][0];
	free(p2);
	free(p1);
	free(array);
}
	
	
void make_ssi(float** Ssi, float *p)
{
	float temp;
	Ssi[0][0] = p[0]; Ssi[0][1] = p[1]; Ssi[0][2] = p[2];
	temp = p[0]+p[2]-(p[0]-p[2])*(p[0]-p[2]);
	Ssi[1][0] = (-p[0]*(p[1]+2*p[2]))/temp;
	Ssi[1][1] = (p[1]*(p[0]-p[2]))/temp;
	Ssi[1][2] = (p[2]*(p[1]+2*p[0]))/temp;
	Ssi[2][0] = -0.5; Ssi[2][1] = 1; Ssi[2][2] = -0.5;
}
	
void make_ss(float** Ss, float *p)
{
	float temp;
	temp = p[0]+p[2]-(p[0]-p[2])*(p[0]-p[2]);
	Ss[0][0] = Ss[1][0] = Ss[2][0] = 1;
	// additive
	Ss[0][1] = -p[1]-2*p[2];
	Ss[1][1] = 1 + Ss[0][1];
	Ss[2][1] = 2 + Ss[0][1];
	
	// dominance
	Ss[0][2] = -2*p[1]*p[2] / temp;
	Ss[1][2] = 4*p[0]*p[2] / temp;
	Ss[2][2] = -2*p[0]*p[1] / temp;
}
	
void solve(float **Ssi, float *e, float *g, int ngen)
{
	int i,j;
	for(i = 0; i < ngen; i++) e[i] = 0;
	
	for(i = 0; i < ngen; i++)
		for(j = 0; j < ngen; j++)
			e[i] += Ssi[i][j] * g[j];
}

// kronecker of SsiA SsiB = Ssi and SsA SsB = Ss
void kronecker(float **A, float **B, float **AB, int n1, int n2, int m1, int m2)
{
	int i,j,k,l,jj,kk;
	
	jj = 0;
	for(i = 0; i < m1; i++) // row mat 1
	{
		for(j = 0; j < m2; j++) // row mat 2
		{
			kk = 0;
			for(k = 0; k < n1; k++) // col mat 1
			{
				for(l = 0; l < n2; l++) // col mat 2
				{
					AB[jj][kk++] = A[i][k] * B[j][l];
				}
			}
			jj++;
		}
	}
}

