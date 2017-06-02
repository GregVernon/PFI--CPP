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
	double	*x, *y;
	double	*XX, *YY;

	x = (double *)mkl_malloc(NX * sizeof(double), 64);
	y = (double *)mkl_malloc(MY * sizeof(double), 64);

	XX = (double *)mkl_malloc(NX*MY * sizeof(double), 64);
	YY = (double *)mkl_malloc(NX*MY * sizeof(double), 64);
	
	// LINSPACE ROUTINE -- Construct x-vector
	x = linspace(xmin, xmax, NX);
	y = linspace(ymin, ymax, MY);

	// Construct XX array
	for (j = 0; j < MY; ++j)
		for (i = 0; i < NX; ++i)
		{
			XX[j * NX + i] = x[i];
		}

	// Construct YY array
	for (j = 0; j < MY; ++j)
		for (i = 0; i < NX; ++i)
		{
			YY[j * NX + i] = y[j];
		}

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
