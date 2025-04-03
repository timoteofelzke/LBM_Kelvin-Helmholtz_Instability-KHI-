/**
 *  @file   mrt.cpp
 *  @author Diogo Nardelli Siebert
 *  @brief  Implementation of  functions to perform the collision step using the MRT operator.
 */

#include <iostream>
#include <string>
#include <iomanip>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include "mrt.h"
#include "../lattice/lattice.h"

using namespace  std;

namespace mrt
{

stCollision collisionParameters(std::string filename)
{
    using boost::property_tree::ptree;
    stCollision prm;

    ptree pt;
    read_ini(filename.c_str(),pt);
    prm.tau  =  pt.get<double>("MRT.tau");
    prm.s    =  relaxationMRT(prm.tau);
    prm.kinematicViscosity = c_s2 * (prm.tau - 0.5 );

    return prm;
}

void reportCollision(stCollision& prm)
{
    cout << endl;
    cout << "Collision Model: MRT";
    cout << endl;
    cout << left << setw(40) << "Relaxation Time:";
    cout << right << setw(30) << prm.tau << endl;

    for (int n = 0 ; n < NUM_OF_VEL ; n++ )
    {
        cout << "s_" << n << left  << setw(40) << ":";
        cout << right << setw(30) << prm.s[n] << endl;
    }

    cout << left << setw(40) << "Viscosity:";
    cout << right << setw(30) << prm.kinematicViscosity << endl;

    cout << endl;
}

void collision( double* iniN, stForce& force,  stCollision& prm, int numberOfPoints)
{
	#pragma omp parallel for
	for (int n = 0; n < numberOfPoints; n++)
	{
        double* f = iniN + n*NUM_OF_VEL;
		double m[NUM_OF_VEL];
		double meq[NUM_OF_VEL];

		calculateMomentsMRT(f,m);
		equilibriumMomentsMRT(meq,m);

		for (int i = 4; i< NUM_OF_VEL; i++)
		{
            m[i] = m[i] - prm.s[i] * (m[i] - meq[i]);
		}

		reconstructedDistributionMRT(f,m);
	}
}

void collisionSwap( double* iniN, stForce& force,  stCollision& prm, int numberOfPoints)
{


	#pragma omp parallel for
	for (int n = 0; n < numberOfPoints; n++)
	{

        double* f = iniN + n*NUM_OF_VEL;
		double m[NUM_OF_VEL];
		double meq[NUM_OF_VEL];

		calculateMomentsMRT(f,m);
		equilibriumMomentsMRT(meq,m);

		for (int i = 4; i< NUM_OF_VEL; i++)
		{
            m[i] = m[i] - prm.s[i] * (m[i] - meq[i]);
		}

		reconstructedDistributionMRT(f,m);

		double temp;
		for (int i = 1; i< NUM_OF_VEL; i = i + 2)
		{
			temp = f[i];
			f[i] = f[i+1];
			f[i+1] = temp;
		}

	}
}

}
