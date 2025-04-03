/**
 *  @file   regularized.cpp
 *  @author Diogo Nardelli Siebert
 *  @brief  Implementation of  functions to perform the collision step using the BGK operator with a regularization procedure.
 */
#include <iostream>
#include <string>
#include <iomanip>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include "regularized.h"
#include "../lattice/lattice.h"

using namespace  std;

namespace regularized
{

stCollision collisionParameters(std::string filename)
{
    using boost::property_tree::ptree;
    stCollision prm;

    ptree pt;
    read_ini(filename.c_str(),pt);
    prm.tau      =  pt.get<double>("Regularized.tau");
    prm.alphaEq  =  1.0/prm.tau;
    prm.alphaNon =  (1.0 - 1.0/prm.tau);
    prm.kinematicViscosity = c_s2 * (prm.tau - 0.5);
    return prm;
}

void reportCollision(stCollision& prm)
{
    cout << endl;
    cout << "Collision Model: Regularized (BGK) ";
    cout << endl;
    cout << left << setw(40) << "Relaxation Time:";
    cout << right << setw(30) << prm.tau << endl;
    cout << left << setw(40) << "Viscosity:";
    cout << right << setw(30) << prm.kinematicViscosity << endl;
    cout << endl;
}


void collision( double* iniN, stForce& force,  stCollision& prm, int numberOfPoints)
{
	#pragma omp parallel for
	for (int n = 0; n < numberOfPoints; n++)
	{
        double* f = iniN + n * NUM_OF_VEL;
		double m[10];

		calculateMoment2nd(f,m);

		double ux = m[1] / m[0];
		double uy = m[2] / m[0];
		double uz = m[3] / m[0];

        /* Computes the the moments (until second order) after collision  */
        m[4] = prm.alphaNon * m[4] +  prm.alphaEq * m[0] * (ux * ux + c_s2);
        m[5] = prm.alphaNon  * m[5] +  prm.alphaEq * m[0] * (uy * uy + c_s2);
        m[6] = prm.alphaNon  * m[6] +  prm.alphaEq * m[0] * (uz * uz + c_s2);
        m[7] = prm.alphaNon  * m[7] +  prm.alphaEq * m[0] * (ux * uy);
        m[8] = prm.alphaNon  * m[8] +  prm.alphaEq * m[0] * (ux * uz);
        m[9] = prm.alphaNon  * m[9] +  prm.alphaEq * m[0] * (uy * uz);

		reconstructedDistribution2nd(f, m);
	}
}

void collisionSwap( double* iniN, stForce& force,  stCollision& prm, int numberOfPoints)
{
	#pragma omp parallel for
	for (int n = 0; n < numberOfPoints; n++)
	{
        double* f = iniN + n * NUM_OF_VEL;
		double m[10];

		calculateMoment2nd(f,m);

		double ux = m[1] / m[0];
		double uy = m[2] / m[0];
		double uz = m[3] / m[0];

        /* Computes the the moments (until second order) after collision  */
        m[4] = prm.alphaNon * m[4] +  prm.alphaEq * m[0] * (ux * ux + c_s2);
        m[5] = prm.alphaNon  * m[5] +  prm.alphaEq * m[0] * (uy * uy + c_s2);
        m[6] = prm.alphaNon  * m[6] +  prm.alphaEq * m[0] * (uz * uz + c_s2);
        m[7] = prm.alphaNon  * m[7] +  prm.alphaEq * m[0] * (ux * uy);
        m[8] = prm.alphaNon  * m[8] +  prm.alphaEq * m[0] * (ux * uz);
        m[9] = prm.alphaNon  * m[9] +  prm.alphaEq * m[0] * (uy * uz);

		reconstructedDistribution2nd(f, m);

		double temp;
		for (int i =1; i < NUM_OF_VEL; i=i+2)
		{
			temp = f[i];
			f[i] = f[i+1];
			f[i+1] = temp;
		}

	}
}

}
