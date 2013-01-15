#include "pop.h"

void ftestlin(record *pop, int n, int loc, allstats* stat, int linked)
{
	int i, df2;
	char *Xa = malloc(sizeof(char) * n);
	char *Xb = malloc(sizeof(char) * n);

	float 
		Sxxa = 0, Sxxb = 0,
		Sxya = 0, Sxyb = 0,
		Syy = 0,
		Sa = 0, Sb = 0,
		bhata, bhatb,
		ahata, ahatb,
		mXa = 0, mXb = 0,
		mY = 0;

	float *Yhata = malloc(n * sizeof(float));
	float *Yhatb = malloc(n * sizeof(float));

// CHECK
	if(linked)
	{
		for(i = 0; i < n; i++)
		{
			Xa[i] = pop[i].qtl[loc].alink[0] + pop[i].qtl[loc].alink[1];
			Xb[i] = pop[i].qtl[loc].blink[0] + pop[i].qtl[loc].blink[1];
		}	
	} else {
		for(i = 0; i < n; i++)
		{
			Xa[i] = pop[i].qtl[loc].a[0] + pop[i].qtl[loc].a[1];
			Xb[i] = pop[i].qtl[loc].b[0] + pop[i].qtl[loc].b[1];
		}
	}

// linear regression and anova

	for(i = 0; i < n; i++)
	{
		mXa += Xa[i];
		mXb += Xb[i];
		mY += pop[i].phen;
	}
	mXa /= n;
	mXb /= n;
	mY /= n;

	for(i = 0; i < n; i++)
	{
		Sxxa += (Xa[i] - mXa) * (Xa[i] - mXa);
		Sxxb += (Xb[i] - mXb) * (Xb[i] - mXb);
		Sxya += (Xa[i] - mXa) * (pop[i].phen - mY);
		Sxyb += (Xb[i] - mXb) * (pop[i].phen - mY);
		Syy += (pop[i].phen - mY) * (pop[i].phen - mY);
	}

	bhata = Sxya / Sxxa;
	bhatb = Sxyb / Sxxb;
	ahata = mY - bhata * mXa;
	ahatb = mY - bhatb * mXb;

	for(i = 0; i < n; i++)
	{
		Yhata[i] = ahata + bhata * Xa[i];
		Yhatb[i] = ahatb + bhatb * Xb[i];
	}

	for(i = 0; i < n; i++)
	{
		Sa += (pop[i].phen - Yhata[i]) * (pop[i].phen - Yhata[i]);
		Sb += (pop[i].phen - Yhatb[i]) * (pop[i].phen - Yhatb[i]);
	}	

	df2 = n - 2;
	Sa /= df2;
	Sb /= df2;

	double Fa = (Sxya*Sxya/Sxxa) / Sa;
	double Fb = (Sxyb*Sxyb/Sxxb) / Sb;
	double p,q;
	int df1 = 1, err;
	if(isnan(Fa)) Fa = 1;
	if(isnan(Fb)) Fb = 1;

	(*fprob)(&df1,&df2,&Fa,&p,&q,&err);
	stat[loc].report[linked].F[0] = (float)-log10(q);
	(*fprob)(&df1,&df2,&Fb,&p,&q,&err);
	stat[loc].report[linked].F[1] = (float)-log10(q);


	free(Yhata);
	free(Yhatb);
	free(Xa);
	free(Xb);
}


void ftest8df(record *pop, int n, int loc, allstats *stat, int linked)
{
	int i, tempx, nfac, factor_count[9], dfm, dfe;
	float MSM, MSE, SSM = 0, SST = 0, SSE = 0, mY = 0;
	float mY_fac[9];

// find genotype groups (nfac), counts (factor_count) and names (factors)
// only works with 3 genotypes - for epistasis change array length declaration
	for(i = 0; i < 9; i++)
	{
		factor_count[i] = 0;
		mY_fac[i] = 0;
	}

	for (i = 0; i < n; i++)
	{
		if(linked)
		{
			tempx = pop[i].qtl[loc].Xlink;
		} else {
			tempx = pop[i].qtl[loc].X;
		}
		mY += pop[i].phen;
		switch (tempx){
			case 0:
				factor_count[0]++;
				mY_fac[0] += pop[i].phen;
				break;
			case 1:
				factor_count[1]++;
				mY_fac[1] += pop[i].phen;
				break;
			case 2:
				factor_count[2]++;
				mY_fac[2] += pop[i].phen;
				break;
			case 3:
				factor_count[3]++;
				mY_fac[3] += pop[i].phen;
				break;
			case 4:
				factor_count[4]++;
				mY_fac[4] += pop[i].phen;
				break;
			case 5:
				factor_count[5]++;
				mY_fac[5] += pop[i].phen;
				break;
			case 6:
				factor_count[6]++;
				mY_fac[6] += pop[i].phen;
				break;
			case 7:
				factor_count[7]++;
				mY_fac[7] += pop[i].phen;
				break;
			case 8:
				factor_count[8]++;
				mY_fac[8] += pop[i].phen;
				break;
		}
	}
	mY = mY / n;

	nfac = 0;
	for(i = 0; i < 9; i++)
	{
		if (factor_count[i]>0)
		{
			nfac++;
			mY_fac[i] /= factor_count[i];
			SSM += (float)factor_count[i] * (mY_fac[i] - mY) * (mY_fac[i] - mY);
		}
	}
	for (i = 0; i < n; i++)
	{
		SST += (pop[i].phen - mY) * (pop[i].phen - mY);
	}
	dfm = nfac - 1;
	dfe = n - nfac;
	SSE = SST - SSM; // within group sum squares

	MSM = SSM / dfm;
	MSE = SSE / dfe;
	
	double F = (double)MSM / MSE;
	if(isnan(F)) F = 1;

	double p,q;
	int err;
	(*fprob)(&dfm,&dfe,&F,&p,&q,&err);
	stat[loc].report[linked].Ffull = (float)-log10(q);

}


void ftest2df(record *pop, int n, int loc, allstats *stat, int linked)
{
	int i,j, nfacA, nfacB, factor_countA[3], factor_countB[3], dfm, dfe;
	float MSM, MSE, SSMA = 0, SSMB = 0, SST = 0, SSE = 0, mY = 0;
	float mY_facA[3], mY_facB[3];

	for(i = 0; i < 3; i++)
	{
		factor_countA[i] = 0;
		mY_facA[i] = 0;
		factor_countB[i] = 0;
		mY_facB[i] = 0;
	}

	for (i = 0; i < n; i++)
	{
		mY += pop[i].phen;
		if(linked)
		{
			j = pop[i].qtl[loc].alink[0] + pop[i].qtl[loc].alink[1];
		} else {
			j = pop[i].qtl[loc].a[0] + pop[i].qtl[loc].a[1];
		}
		switch (j){
			case 0:
				factor_countA[0]++;
				mY_facA[0] += pop[i].phen;
				break;
			case 1:
				factor_countA[1]++;
				mY_facA[1] += pop[i].phen;
				break;
			case 2:
				factor_countA[2]++;
				mY_facA[2] += pop[i].phen;
				break;
		}
		if(linked)
		{
			j = pop[i].qtl[loc].blink[0] + pop[i].qtl[loc].blink[1];
		} else {
			j = pop[i].qtl[loc].b[0] + pop[i].qtl[loc].b[1];
		}
		switch (j){
			case 0:
				factor_countB[0]++;
				mY_facB[0] += pop[i].phen;
				break;
			case 1:
				factor_countB[1]++;
				mY_facB[1] += pop[i].phen;
				break;
			case 2:
				factor_countB[2]++;
				mY_facB[2] += pop[i].phen;
				break;
		}
	}
	mY = mY / n;


	nfacA = 0;
	nfacB = 0;
	SSMA = 0; SSMB = 0; SST = 0;
	for(i = 0; i < 3; i++)
	{
		if (factor_countA[i]>0)
		{
			nfacA++;
			mY_facA[i] /= factor_countA[i];
			SSMA += (float)factor_countA[i] * (mY_facA[i] - mY) * (mY_facA[i] - mY);
		}
		if (factor_countB[i]>0)
		{
			nfacB++;
			mY_facB[i] /= factor_countB[i];
			SSMB += (float)factor_countB[i] * (mY_facB[i] - mY) * (mY_facB[i] - mY);
		}
	}

	for (i = 0; i < n; i++)
	{
		SST += (pop[i].phen - mY) * (pop[i].phen - mY);
	}
	dfm = nfacA - 1;
	dfe = n - nfacA;
	SSE = SST - SSMA; // within group sum squares

	MSM = SSMA / dfm;
	MSE = SSE / dfe;
	
	double F = (double)MSM / MSE;
	if(isnan(F)) F = 1;

	double p,q;
	int err;
	(*fprob)(&dfm,&dfe,&F,&p,&q,&err);
	stat[loc].report[linked].F2[0] = (float)-log10(q);

	dfm = nfacB - 1;
	dfe = n - nfacB;
	SSE = SST - SSMB; // within group sum squares

	MSM = SSMB / dfm;
	MSE = SSE / dfe;
	
	F = (double)MSM / MSE;
	if(isnan(F)) F = 1;
	p = 0; q = 0;
	(*fprob)(&dfm,&dfe,&F,&p,&q,&err);
	stat[loc].report[linked].F2[1] = (float)-log10(q);

}



