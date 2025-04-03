/**
 *  @file   permeability.cpp
 *  @author Diogo Nardelli Siebert
 *  @brief  Implementation of functions to compute the intrinsic permeability
 *
 *  Implementation of functions to compute the intrinsic permeability of a porous media image
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include "../lbmbase/geometry.h"
#include "../parameters.h"
#include "../lbmbase/lattice/lattice.h"

#include "permeability.h"

using namespace std;

namespace permeabilityProblem
{

stProblem problemParameters( string filename )
{
	stProblem prm;
	using boost::property_tree::ptree;
	ptree pt;
	read_ini(filename.c_str(),pt);

	prm.filename  =  pt.get<string>("Permeability.filename","permeability.csv");
	prm.axis      =  pt.get<int>("Permeability.axis");
	prm.err 	  =  pt.get<double>("Permeability.error",0);
	
    cout << endl;
    cout << left << setw(40) << "Using module:";
    cout << right << setw(30) << "Permeability" << endl;

    cout << left << setw(40) << "Permeability computed on axis: ";
    cout << right<< setw(30) << prm.axis << endl;

    cout << left << setw(40) << "Saving permeability data on the file: ";
    cout << right<< setw(30) << prm.filename << endl << endl;

    return prm;
}

void problemInitialize(double* iniN, stCollision& col , stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force)
{

}

void problemFinalize(double* iniN, stCollision& col , stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force, int step )
{

}

void problemBoundary(double* iniN, stCollision& col , stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force, int step )
{

}

void problemPosCollision(double* iniN, stCollision& col , stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force, int step )
{

}

bool problemOutput(double* iniN, stCollision& col , stGeometry &geo, stParameters& prm, stProblem& per, stForce& force, int step )
{
    string filename = per.filename;
    static double permOld = -1;
    bool keepRunning = 1;

    double permNew;

    if (step == 1)
    {
        ofstream output(filename.c_str());
        output << "Step\tPermeability\n";
        output.close();
    }

    if ( (step%100 == 0 ) )
    {
        int axis = per.axis;
        double g[3] = {prm.accX , prm.accY , prm.accZ};
        double u = phaseVelocity(iniN, geo, axis, g[axis] );

        permNew = u * prm.viscosity / g[axis];
        ofstream output(filename.c_str(), ofstream::app );
        output << step << "\t" << permNew << endl;
        output.close();

        if ( abs(permOld - permNew) < per.err )
        {
            keepRunning = false;
        }

        permOld = permNew;
    }

    return keepRunning;
}

double phaseVelocity(double* iniN, stGeometry &geo, int axis, double g)
{
	double flowRate = 0;
	
	#pragma omp parallel for reduction(+:flowRate)
	for (int n=0; n < geo.numberOfPoints ; n++)
	{
		double  rho;
		double  v[3];
                double* f = iniN + mapMemory(n);
		calculateMacro(f,v[0],v[1],v[2],rho);
		flowRate += *(v+axis) +0.5*g;
	}
	
	return flowRate / static_cast<double> (geo.nx * geo.ny * geo.nz);
}



}
//
