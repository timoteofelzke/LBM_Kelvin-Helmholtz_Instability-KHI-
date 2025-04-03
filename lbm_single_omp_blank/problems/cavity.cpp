/**
 *  \file default.cpp
 *
 *  Description: A empty problem (default) to be used as template to implement new problems
 *  into the code
 */

#include <iostream>
#include <fstream>
#include <cmath>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include "../lbmbase/geometry.h"
#include "../parameters.h"
#include "../lbmbase/boundary/zouhe.h"
#include "../lbmbase/lattice/lattice.h"

#include "cavity.h"

using namespace std;

namespace cavityProblem
{

stProblem problemParameters( string filename )
{
    stProblem prm;
    using boost::property_tree::ptree;
    ptree pt;
    read_ini(filename.c_str(),pt);

    prm.Re  =  pt.get<double>("Cavity.re");

    return prm;
}

void problemInitialize(double* iniN, stCollision& col , stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force)
{
	sp.u_x = sp.Re * col.kinematicViscosity / geo.nx;

	for (int x = 0 ; x < geo.nx; x++)
	{
		unsigned int id = getIndex(geo,x,geo.ny-1,0);
		if (id) sp.northBoundary.push_back(id);
	}

        for (int x = 0 ; x < geo.nx; x++)
        {
                unsigned int id = getIndex(geo,x,0,0);
		if (id) sp.southBoundary.push_back(id);
        }

        for (int y = 1 ; y < geo.ny-1; y++)
        {
                unsigned int id = getIndex(geo,0,y,0);
		if (id) sp.westBoundary.push_back(id);
        }

        for (int y = 1; y < geo.ny-1; y++)
        {
                unsigned int id = getIndex(geo,geo.nx-1,y,0);
		if (id) sp.eastBoundary.push_back(id);
        }
}


void problemFinalize(double* iniN, stCollision& col , stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force, int step )
{
}

void problemBoundary(double* iniN, stCollision& col , stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force, int step )
{
	#pragma omp parallel for
	for (auto & id : sp.northBoundary) 
	{
		double *f = iniN + (id-1)*NUM_OF_VEL;
                zouHeVelocityBack(f,sp.u_x,sp.u_y,sp.u_z);
	}

	double u_x = 0.0;
	#pragma omp parallel for
	for (auto & id : sp.southBoundary) 
	{
		double *f = iniN + (id-1)*NUM_OF_VEL;
                zouHeVelocityFront(f,u_x,sp.u_y,sp.u_z);
	}

	#pragma omp parallel for
	for (auto & id : sp.westBoundary) 
        {
                double *f = iniN + (id-1)*NUM_OF_VEL;
		zouHeVelocityLeft(f,u_x,sp.u_y,sp.u_z);
        }

	#pragma omp parallel for
	for (auto & id : sp.eastBoundary) 
        {
                double *f = iniN + (id-1)*NUM_OF_VEL;
		zouHeVelocityRight(f,u_x,sp.u_y,sp.u_z);
        }

}

void problemPosCollision(double* iniN, stCollision& col , stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force, int step )
{

}

bool problemOutput(double* iniN, stCollision& col , stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force, int step )
{
    return true;
}

}

