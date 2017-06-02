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

void main()
{
	auto allStart = Clock::now();

	MKL_INT	NX = 201, MY = 401;
	int i, j;
	double xmin = -20, xmax = 20, ymin = 1, ymax = 11;
	double dx, dy;
	double dx2, dy2;
	double dxx, dyy;
	double	*x, *y;
	
	// Define Parabolic Coordinates
	x = (double *)mkl_malloc(NX * sizeof(double), 64);
	y = (double *)mkl_malloc(MY * sizeof(double), 64);
	
	x = linspace(xmin, xmax, NX);
	y = linspace(ymin, ymax, MY);

	dx = x[1] - x[0];
	dx2 = 2 * dx;
	dxx = dx * dx;

	dy = y[1] - y[0];
	dy2 = 2 * dy;
	dyy = dy * dy;

	// Compute Mapping: Parabolic -> Cartesian

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

std::tuple<double *, double *, double *, double *, double *, double *, double *, double *, double *, double *> Para2Cart()
{
	int i, j;
	double *pX, *pY;

	pX = (double *)mkl_malloc(NX * sizeof(double), 64);
	pY = (double *)mkl_malloc(MY * sizeof(double), 64);

	for (i = 0; i < NX; ++i)
	{
		if (i <= NX)
		{
			pX[i] = pow((-1 / 2) * abs(xmax - xmin)*abs((ii - ((NX - 1) / 2) - 1) / ((NX - 1) / 2)), xskew);
		}
		else if (i > NX)
		{
			pX[i] = pow(( 1 / 2) * abs(xmax - xmin)*abs((ii - ((NX - 1) / 2) - 1) / ((NX - 1) / 2)), xskew);
		}
	}
}