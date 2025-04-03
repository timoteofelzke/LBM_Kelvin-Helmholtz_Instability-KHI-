/**
 *  @file   convective_multi.cpp
 *  @author Diogo Nardelli Siebert
 *
 *  Declaration of functions to impose boundary conditions.
 */

#include <iostream>
#include <cstdlib>
#include "convective_multi.h"
#include "../geometry.h"
#include "../lattice/lattice.h"

using namespace std;

namespace multi
{

stBoundaryConvective defineConvective(stGeometry& geo, double* iniR, double* iniB, double* iniW, boundaryPlane plane, velocityType uType)
{
    stBoundaryConvective boundary;

    // Compute the axis which the condition will be applied;
    boundary.axis = plane / 2;

    // Compute the direction (+1 to BEGIN) e (-1 to END) to which the boundary condition will be applied
    // and coordinate value of the boundary

    int end[3] = {geo.nx,geo.ny,geo.nz};
    if (plane % 2 == 0)
    {
        boundary.direction = +1;
        boundary.coordinate = 0;
    }
    else
    {
        boundary.direction = -1;
        boundary.coordinate = end[ boundary.axis] - 1;
    }

    boundary.uType = uType;

    for (unsigned int id = 1; id <= (unsigned int) geo.numberOfPoints; id++)
    {
        int r[3];
        getPosition(geo,id,r[0],r[1],r[2]);
        if ( boundary.coordinate == r[ boundary.axis] )
        {
            boundary.boundaryLayer.push_back(id);
            // Computes the position of the neighbouring site
            r[boundary.axis] += boundary.direction;
            unsigned int idNeighbour  = getIndex(geo, r[0], r[1] , r[2]);
            if ( (idNeighbour==0) )
            {
                cout << "Error, geometry do not allow for Convective Boundary Condition" << endl;
                exit(3);
            }
            boundary.neighborLayer.push_back(idNeighbour);
        }
    }

    allocate(boundary,iniR,iniB,iniW);

	return boundary;
}

void updateConvectiveVelocity(stBoundaryConvective& boundary,double* iniR,double* iniB)
{
    double meanU = 0;

    #pragma omp parallel for reduction(+:meanU)
    for (size_t k=0; k < boundary.neighborLayer.size(); k++)
    {
        unsigned int idN = boundary.neighborLayer[k];
        double* R  = iniR + ( idN - 1 ) * NUM_OF_VEL;
        double* B  = iniB + ( idN - 1 ) * NUM_OF_VEL;
        double f[NUM_OF_VEL];

        for (int i  = 0 ; i< NUM_OF_VEL; i++) f[i] = R[i] + B[i];

        double rho;
        double v[3];
        calculateMacro(f,v[0],v[1],v[2],rho);

        boundary.U[k] = v[ boundary.axis ] ;
        meanU += v[ boundary.axis ];
    }

    boundary.meanU = meanU / boundary.neighborLayer.size() ;
}

void allocate(stBoundaryConvective& boundary,double* iniR,double* iniB ,double* iniW)
{
    int boundaryPoints = boundary.boundaryLayer.size();

    boundary.wLast = new double[  boundaryPoints * NUM_OF_VEL ];
    boundary.blueLast = new double[  boundaryPoints * NUM_OF_VEL ];
    boundary.redLast = new double[  boundaryPoints * NUM_OF_VEL ];
    boundary.U = new double[  boundaryPoints ];

    #pragma omp parallel for
    for (int k=0; k < boundaryPoints; k++)
    {
        unsigned int id  = boundary.boundaryLayer[k];

        double* R  = iniR + ( id - 1 )  * NUM_OF_VEL;
        double* Rl = boundary.redLast + k * NUM_OF_VEL;

        double* B  = iniB + ( id - 1 )  * NUM_OF_VEL;
        double* Bl = boundary.blueLast + k * NUM_OF_VEL;

        double* W  = iniW + ( id - 1 )  * NUM_OF_VEL;
        double* Wl = boundary.wLast + k * NUM_OF_VEL;

        for (int i  = 0 ; i< NUM_OF_VEL; i++)
        {
            Rl[i] = R[i];
            Bl[i] = B[i];
            Wl[i] = W[i];
        }
    }
}

void applyConvectiveRB(stBoundaryConvective& boundary, double * iniR, double * iniB)
{
    updateConvectiveVelocity(boundary,iniR,iniB);

    int boundaryPoints = boundary.boundaryLayer.size();
    #pragma omp parallel for
    for (int k=0; k < boundaryPoints; k++)
    {
        unsigned int id  = boundary.boundaryLayer[k];
        unsigned int idN = boundary.neighborLayer[k];

        double* R  = iniR + ( id - 1 )  * NUM_OF_VEL;
        double* Rn = iniR + ( idN - 1 ) * NUM_OF_VEL;
        double* Rl = boundary.redLast + k * NUM_OF_VEL;

        double* B  = iniB + ( id - 1 )  * NUM_OF_VEL;
        double* Bn = iniB + ( idN - 1 ) * NUM_OF_VEL;
        double* Bl = boundary.blueLast + k * NUM_OF_VEL;

        double U;
        if ( boundary.uType == MEAN_VELOCITY)
        {
            U = - boundary.direction * boundary.meanU;
        }
        else
        {
            U = - boundary.direction * boundary.U[k];
        }

        for (int i  = 0 ; i< NUM_OF_VEL; i++)
        {
            R[i] = (Rl[i] + U * Rn[i] ) / (1+U);
            B[i] = (Bl[i] + U * Bn[i] ) / (1+U);
            Rl[i] = R[i];
            Bl[i] = B[i];
        }
    }
}

void applyConvectiveW(stBoundaryConvective& boundary, double * iniW)
{
    int boundaryPoints = boundary.boundaryLayer.size();
    #pragma omp parallel for
    for (int k=0; k < boundaryPoints; k++)
    {
        unsigned int id  = boundary.boundaryLayer[k];
        unsigned int idN = boundary.neighborLayer[k];

        double* W  = iniW + ( id - 1 )  * NUM_OF_VEL;
        double* Wn = iniW + ( idN - 1 ) * NUM_OF_VEL;
        double* Wl = boundary.wLast + k * NUM_OF_VEL;

        double U;
        if ( boundary.uType == MEAN_VELOCITY)
        {
            U = - boundary.direction * boundary.meanU;
        }
        else
        {
            U = - boundary.direction * boundary.U[k];
        }

        for (int i  = 0 ; i< NUM_OF_VEL; i++)
        {
            W[i] = (Wl[i] + U * Wn[i] ) / (1+U);
            Wl[i] = W[i];
        }
    }
}



}
