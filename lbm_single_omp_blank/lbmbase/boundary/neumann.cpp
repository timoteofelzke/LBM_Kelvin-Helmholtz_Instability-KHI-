/**
 *  @file   neumann.cpp
 *  @author Diogo Nardelli Siebert
 *  @brief  Implementation of functions to impose Neumann type boundary conditions.
 */

#include "neumann.h"
#include "../geometry.h"
#include "../lattice/lattice.h"

stBoundaryNeumann defineNeumannYZ( stGeometry& geo, int x , int deltaX )
{
	stBoundaryNeumann boundary;
	for ( int z = 0; z < geo.nz; z++ )
	{
		for ( int y = 0; y < geo.ny; y++ )
		{
			unsigned int idLocal = getIndex(geo,x,y,z);
			if ( idLocal )
			{
				unsigned int idNext = getIndex(geo,x + deltaX, y , z);
	        	        boundary.locBoundary.push_back(idLocal-1);
        	       		boundary.locNext.push_back(idNext-1);
			}
		}
	}
	return boundary;
}

stBoundaryNeumann defineNeumannXY( stGeometry& geo, int z , int deltaZ )
{
	stBoundaryNeumann boundary;
	for ( int y = 0; y < geo.ny; y++ )
	{
		for ( int x = 0; x < geo.nx; x++ )
		{
			unsigned int idLocal = getIndex(geo,x,y,z);
			if ( idLocal )
			{
				unsigned int idNext = getIndex(geo,x, y , z  + deltaZ);
	        	        boundary.locBoundary.push_back(idLocal-1);
        	       		boundary.locNext.push_back(idNext-1);
			}
		}
	}
	return boundary;
}

stBoundaryNeumann defineNeumannXZ( stGeometry& geo, int y , int deltaY )
{
	stBoundaryNeumann boundary;
	for ( int z = 0; z < geo.nz; z++ )
	{
		for ( int x = 0; x < geo.nx; x++ )
		{
			unsigned int idLocal = getIndex(geo,x,y,z);
			if ( idLocal )
			{
				unsigned int idNext = getIndex(geo,x, y   + deltaY , z);
       			        boundary.locBoundary.push_back(idLocal-1);
                		boundary.locNext.push_back(idNext-1);
            		}
		}
	}
	return boundary;
}

void neumann(stBoundaryNeumann& boundary, double * iniN, double rho)
{
    int numberOfPoints = boundary.locBoundary.size();

    #pragma omp parallel for
    for (int k=0; k < numberOfPoints; k++)
    {
        int64_t  locBoundary = boundary.locBoundary[k];
        int64_t  locNext = boundary.locNext[k];

        double* fBoundary = iniN + locBoundary*NUM_OF_VEL;

        if (locNext != locBoundary)
        {
            double* fNext =  iniN + locNext * NUM_OF_VEL;
            for (int i  = 0 ; i< NUM_OF_VEL; i++) fBoundary[i] = fNext[i];
        }
        else
        {
            equilibriumDistribution( 0.0 , 0.0, 0.0 , rho, fBoundary);
        }
    }
}

void neumannEquilibrium(stBoundaryNeumann& boundary, double * iniN, double rho)
{
    int numberOfPoints = boundary.locBoundary.size();

    #pragma omp parallel for
    for (int k=0; k < numberOfPoints; k++)
    {
        int64_t  locBoundary = boundary.locBoundary[k];
        int64_t  locNext = boundary.locNext[k];

	double vxA, vyA, vzA, rhoA;
        double *f;

        if ( locNext != locBoundary )
        {
            f = iniN + locNext*NUM_OF_VEL;
            calculateMacro( f, vxA, vyA, vzA, rhoA );
        }
        else
        {
            vxA = vyA = vzA = 0.0;
        }

        f = iniN + locBoundary *NUM_OF_VEL ;
        equilibriumDistribution( vxA, vyA, vzA, rho, f);
    }
}

void neumannVelocity( stBoundaryNeumann& boundary, double * iniN, double rho )
{
	int numberOfPoints =  boundary.locBoundary.size();

	#pragma omp parallel for
	for (int k=0; k < numberOfPoints; k++)
	{
	        int64_t  locBoundary = boundary.locBoundary[k];
	   	int64_t  locNext = boundary.locNext[k];

        	double* fLocal = iniN + locBoundary*NUM_OF_VEL;

		if ( locNext >= 0)
		{
            		double* fNext =  iniN + locNext*NUM_OF_VEL;
		        double rhoNext = calculateMass( fNext );
            		double factor = rho / rhoNext;

			for (int i  = 0 ; i< NUM_OF_VEL; i++) fLocal[i] = fNext[i] * factor;
		}
		else
		{
			equilibriumDistribution( 0.0 , 0.0, 0.0 , rho, fLocal);
		}
	}
}
