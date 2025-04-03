/**
 *  @file   boundary.cpp
 *  @author Diogo Nardelli Siebert
 *
 *  Declaration of functions to impose boundary conditions.
 */

#include <iostream>
#include <cstdlib>
#include "convective.h"
#include "../geometry.h"
#include "../lattice/lattice.h"

using namespace std;

stBoundaryConvective defineConvectiveYZ( stGeometry& geo, int x , int deltaX )
{
	stBoundaryConvective boundary;
    boundary.axis = 0;

	for ( int z = 0; z < geo.nz; z++ )
	{
		for ( int y = 0; y < geo.ny; y++ )
		{
            unsigned int idZero = getIndex(geo,x,y,z);
            if ( idZero )
			{
                unsigned int idFirst  = getIndex(geo,x +   deltaX, y , z);

                if ( (idFirst==0) )
                {
                    cout << "Error, geometry do not allow for Convective Boundary Condition" << endl;
                    exit(3);
                }

                boundary.zeroLayer.push_back(idZero-1);
                boundary.firstLayer.push_back(idFirst-1);
			}
		}
	}

    boundary.direction = 0;  
    boundary.backup = NULL;

	return boundary;
}

void allocate(stBoundaryConvective& boundary,double* iniN)
{
    int numberOfPoints = boundary.zeroLayer.size();

    // if it is the first time the function is executed the buffer will no be allocated
    // so is necessary to allocate it and to copy the information from the adjecnt sites

    boundary.backup  = new double[  numberOfPoints * NUM_OF_VEL ];

    for (int k=0; k < numberOfPoints; k++)
    {
        int nFirst = boundary.firstLayer[k];

        double* fBackup = boundary.backup + k*NUM_OF_VEL;
        double* f     = iniN + nFirst * NUM_OF_VEL;

        for (int i  = 0 ; i< NUM_OF_VEL; i++)
        {
            fBackup[i]   = f[i];
        }
    }
}

void convective(stBoundaryConvective& boundary, double * iniN)
{
    if (boundary.backup == NULL) allocate(boundary,iniN);

    int numberOfPoints = boundary.zeroLayer.size();

    double U = 0.;

    for (int k=0; k < numberOfPoints; k++)
    {
        int nFirst = boundary.firstLayer[k];

        double* f  = iniN +  nFirst * NUM_OF_VEL;
        double vx,vy,vz,rho;
        calculateMacro(f,vx,vy,vz,rho);
        U+= vx;
    }

    U = U / numberOfPoints;

    #pragma omp parallel for
    for (int k=0; k < numberOfPoints; k++)
    {
        int nZero   = boundary.zeroLayer[k];
        int nFirst  = boundary.firstLayer[k];

        double* fBoundary = iniN + nZero  * NUM_OF_VEL;
        double* fNeighbor = iniN + nFirst * NUM_OF_VEL;
        double* fBackup   = boundary.backup + k*NUM_OF_VEL;

        for (int i  = 0 ; i< NUM_OF_VEL; i++)
        {
            fBoundary[i] = (fBackup[i] + U * fNeighbor[i])/(1+U);
            fBackup[i] = fBoundary[i];
        }

    }
}

