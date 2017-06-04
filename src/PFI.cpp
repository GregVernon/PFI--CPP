#include <iostream>
#include <cstdio>
#include <chrono>
#include <stdlib.h>
#include <string>
#include "mkl.h"
#include "mkl_blas.h"
typedef std::chrono::high_resolution_clock Clock;

// function declarations
double * linspace(double vMin, double vMax, MKL_INT nInterval);
struct CartesianCoordinates Para2Cart(struct ParabolicCoordinates *pCOORD);

struct ParabolicCoordinates
{
	int NX, MY;
	double xmin, xmax, ymin, ymax;
	double *x, *y;
	double dx, dy;
	double dx2, dy2;
	double dxx, dyy;

};

struct CartesianCoordinates
{
	
	double *bX, *bY; 
	double *pX, *pY;
	double *mu, *dmux, *dmuy;
	double *eta, *detax, *detay;

};

void main()
{
	auto allStart = Clock::now();

	int	NX = 201, MY = 401;
	int i, j;
	double xmin = -20, xmax = 20, ymin = 1, ymax = 11;
	double dx, dy;
	double dx2, dy2;
	double dxx, dyy;
	double	*x, *y;
	struct ParabolicCoordinates *pCOORD = new ParabolicCoordinates;
	struct CartesianCoordinates *cCOORD = new CartesianCoordinates;

	// Define Parabolic Coordinates
	x = linspace(xmin, xmax, NX);
	y = linspace(ymin, ymax, MY);

	dx = x[1] - x[0];
	dx2 = 2 * dx;
	dxx = dx * dx;

	dy = y[1] - y[0];
	dy2 = 2 * dy;
	dyy = dy * dy;

	pCOORD->NX = NX;
	pCOORD->MY = MY;
	pCOORD->xmin = xmin;
	pCOORD->xmax = xmax;
	pCOORD->ymin = ymin;
	pCOORD->ymax = ymax;
	pCOORD->x = x;
	pCOORD->y = y;
	pCOORD->dx = dx;
	pCOORD->dx2 = dx2;
	pCOORD->dxx = dxx;
	pCOORD->dy = dy;
	pCOORD->dy2 = dy2;
	pCOORD->dyy = dyy;

	// Compute Mapping: Parabolic -> Cartesian
	*cCOORD = Para2Cart(pCOORD);
	
	auto allEnd = Clock::now();
	std::cout << "Total Runtime: " << std::chrono::duration_cast<std::chrono::nanoseconds>(allEnd - allStart).count() << " nanoseconds" << std::endl;
}


double * linspace(double vMin, double vMax, MKL_INT nInterval)
{
	int i;
	double dv;
	double *v;
	
	v = (double *)mkl_malloc(nInterval * sizeof(double), 64);

	dv = (vMax - vMin) / (nInterval - 1);
	v[0] = vMin;
	for (i = 1; i < nInterval; ++i)
	{
		v[i] = v[i - 1] + dv;
	}

	return v;
}



CartesianCoordinates Para2Cart(struct ParabolicCoordinates *pCOORD)
{
	int i, j;
	pCOORD->x;
	int NX = pCOORD->NX;
	int MY = pCOORD->MY;
	MKL_INT nNODE = NX * MY;
	double xmin = pCOORD->xmin;
	double xmax = pCOORD->xmax;
	double ymin = pCOORD->ymin;
	double ymax = pCOORD->ymax;
	double *x = pCOORD->x;
	double xskew, yskew;

	double *pX, *pY;
	double *bX, *bY;
	double *mu, *eta;
	double *dmux, *dmuy;
	double *detax, *detay;
	double *r;
	
	struct CartesianCoordinates * cCOORD = new CartesianCoordinates;
	xskew = 1;
	yskew = 1;
	pX    =	(double *)mkl_malloc(NX*MY * sizeof(double), 64);
	pY    =	(double *)mkl_malloc(NX*MY * sizeof(double), 64);
	dmux  =	(double *)mkl_malloc(NX*MY * sizeof(double), 64);
	dmuy  =	(double *)mkl_malloc(NX*MY * sizeof(double), 64);
	detax = (double *)mkl_malloc(NX*MY * sizeof(double), 64);
	detay = (double *)mkl_malloc(NX*MY * sizeof(double), 64);
	r	  = (double *)mkl_malloc(NX*MY * sizeof(double), 64);
	mu    =	(double *)mkl_malloc(NX * sizeof(double), 64);
	eta   =	(double *)mkl_malloc(MY * sizeof(double), 64);
	bX    =	(double *)mkl_malloc(NX * sizeof(double), 64);
	bY    =	(double *)mkl_malloc(NX * sizeof(double), 64);

	for (i = 0; i < NX; ++i)
	{
		if (x[i] <= 0.)
		{
			mu[i] = pow((-1 / 2) * abs(xmax - xmin)*abs((i - ((NX - 1) / 2) - 1) / ((NX - 1) / 2)), xskew);
		}
		else if (x[i] > 0.)
		{
			mu[i] = pow(( 1 / 2) * abs(xmax - xmin)*abs((i - ((NX - 1) / 2) - 1) / ((NX - 1) / 2)), xskew);
		}
	}

	for (j = 0; j < MY; ++j)
	{
		eta[j] = ymin + abs(ymax - ymin) * pow(((j - 1) / (MY - 1)), yskew);
	}

	for (j = 0; j < MY; ++j)
		for (i = 0; i < NX; ++i)
		{
			pX[j*NX + i] = 1 / 2 * (pow(mu[i], 2) - pow(eta[j], 2));
			pY[j*NX + i] = mu[i] * eta[j];
			r[j*NX + i]  = pow(pX[j*NX + i], 2) + pow(pY[j*NX + i], 2);
		}
	
	for (i = 0; i < NX; ++i)
	{
		bX[i] = 1 / 2 * (pow(mu[i], 2) - 1);
		bY[i] = mu[i];
	}
	
	for (j = 1; j < MY; ++j)
		for (i = 0; i < NX; ++i)
		{
			if (x[i] == 0)
			{
				continue; 
				// SKIP to avoid a "divide by zero" case
				// Note that this should be an edge case as we would usually expect 
				// x[i] to be *close to* but not exactly equal to zero
				// However, x[i] should only be zero when NX is ODD and we check for this case later
			}
			else
			{
				dmux[j*NX + i]  = (1 / 2) * ((pX[j*NX + i] / r[j*NX + i]) + 1) / (mu[i]);
				dmuy[j*NX + i]  = (1 / 2) * ((pY[j*NX + i] / r[j*NX + i]) + 0) / (mu[i]);
				detax[j*NX + i] = (1 / 2) * ((pX[j*NX + i] / r[j*NX + i]) - 1) / (eta[j]);
				detay[j*NX + i] = (1 / 2) * ((pY[j*NX + i] / r[j*NX + i]) - 0) / (eta[j]);
			}
		}
	
	if (fmod(NX, 2) != 0)
	{
		// NX is ODD
		i = (NX - 1) / 2;
		for (j = 1; j < MY; ++j)
		{
			dmux[j*NX + i] = (dmux[j*NX + (i - 1)] + dmux[j*NX + (i + 1)]) / 2;
			dmuy[j*NX + i] = (dmuy[j*NX + (i - 1)] + dmuy[j*NX + (i + 1)]) / 2;
			detax[j*NX + i] = (detax[j*NX + (i - 1)] + detax[j*NX + (i + 1)]) / 2;
			detay[j*NX + i] = (detay[j*NX + (i - 1)] + detay[j*NX + (i + 1)]) / 2;
		}
	}
	
	
	cCOORD->bX = bX; 
	cCOORD->bY = bY; 
	cCOORD->pX = pX;
	cCOORD->pY = pY;
	cCOORD->mu = mu;
	cCOORD->dmux = dmux;
	cCOORD->dmuy = dmuy;
	cCOORD->eta = eta;
	cCOORD->detax = detax;
	cCOORD->detay = detay;

	return *cCOORD;
}