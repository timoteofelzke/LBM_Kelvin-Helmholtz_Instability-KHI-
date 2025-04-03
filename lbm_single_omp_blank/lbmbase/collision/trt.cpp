/**
 *  @file   trt.cpp
 *  @author Diogo Nardelli Siebert
 *  @brief  Implementation of  functions to perform the collision step using the TRT operator.
 */

#include <iostream>
#include <string>
#include <iomanip>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include "trt.h"
#include "../force/force.h"
#include "../lattice/lattice.h"

using namespace std;

namespace trt
{

stCollision collisionParameters(std::string filename)
{
    using boost::property_tree::ptree;
    stCollision prm;

    ptree pt;
    read_ini(filename.c_str(),pt);
    prm.tauSym      =  pt.get<double>("TRT.tau_symmetric");
    prm.tauAnt      =  pt.get<double>("TRT.tau_antisymmetric",  ( 8.0 -1.0 /prm.tauSym ) / ( 8.0 * (  2.0 - 1.0 / prm.tauSym ) ) );

    //  Stores the division of 0.5 by tau for optimization. This value is used in all collision.

    prm.alphaSym = 0.5 / prm.tauSym;
    prm.alphaAnt = 0.5 / prm.tauAnt;
    prm.kinematicViscosity = c_s2 * (prm.tauSym - 0.5);
    return prm;
}

void reportCollision(stCollision& prm)
{
    cout << endl;
    cout << "Collision Model: TRT";
    cout << endl;
    cout << left << setw(40) << "Relaxation Time (Sym):";
    cout << right << setw(30) << prm.tauSym << endl;
    cout << left << setw(40) << "Relaxation Time (Ant):";
    cout << right << setw(30) << prm.tauAnt << endl;
    cout << left << setw(40) << "Viscosity:";
    cout << right << setw(30) << prm.kinematicViscosity << endl;
    cout << endl;
}

void collision(double* iniN, stForce& force,  stCollision& prm, int numberOfPoints)
{
	//  Stores (for optimization) the ammount of velocity to be added in each direction because of the external force field

	#pragma omp parallel for
	for (int n = 0; n < numberOfPoints; n++)
	{
		// Point f to the distribution function at the nth site.
        double* f = iniN + n * NUM_OF_VEL;

		// Computes the symmetric and antisymmetric distribution functions.
		double vx, vy, vz, rho;

        double* F = force.constant ? force.constantForce : force.forceField + n * NUM_OF_DIM;

        calculateMacroEquilibrium( f, vx, vy, vz, rho, F[0] , F[1]  , F[2]  );

		double f_eq_sim[NUM_OF_VEL];
		double f_eq_ant[NUM_OF_VEL];
        double F_sim[NUM_OF_VEL];
        double F_ant[NUM_OF_VEL];

        equilibriumDistributionSym( vx  , vy  , vz , rho, f_eq_sim);
        equilibriumDistributionAnt( vx  , vy  , vz , rho, f_eq_ant);

        double trt_sim = prm.alphaSym *  ( f_eq_sim[0]- 2.0 * f[0]);
		double trt_ant = 0;

        forceTermSym(F_sim, F[0] ,  F[1] ,  F[2]  ,vx,vy,vz, prm.tauSym);
        forceTermAnt(F_ant, F[0] ,  F[1] ,  F[2]  ,vx,vy,vz, prm.tauAnt);

        f[0] = f[0] + trt_sim + F_sim[0] ;

		//  Performe the collision for all the directions. The collision is performed explicitly only for odd directions (1,3,5 ,...)
		//  For even directions the symmtrie property is used.

        for (int i = 1; i< NUM_OF_VEL; i = i + 2)
		{
            trt_sim = prm.alphaSym *  ( f_eq_sim[i]   - f[i] - f[i+1] );
            trt_ant = prm.alphaAnt *  ( f_eq_ant[i]   - f[i] + f[i+1] );
            f[i] = f[i] + trt_sim + trt_ant + F_sim[i] + F_ant[i];
            f[i+1] = f[i+1] + trt_sim - trt_ant + F_sim[i] - F_ant[i];
		}

	}

	// Add the force and momentum to the external variables.

}

void collisionSwap( double* iniN, stForce& force,  stCollision& prm, int numberOfPoints)
{
	#pragma omp parallel for
	for (int n = 0; n < numberOfPoints; n++)
	{
		// Point f to the distribution function at the nth site.
        double* f = iniN + n * NUM_OF_VEL;

		// Computes the symmetric and antisymmetric distribution functions.
		double vx, vy, vz, rho;

        double* F = force.constant ? force.constantForce : force.forceField + n * NUM_OF_DIM;

        calculateMacroEquilibrium( f, vx, vy, vz, rho, F[0] , F[1]  , F[2]  );

		double f_eq_sim[NUM_OF_VEL];
		double f_eq_ant[NUM_OF_VEL];

        equilibriumDistributionSym( vx  , vy  , vz , rho, f_eq_sim);
        equilibriumDistributionAnt( vx  , vy  , vz , rho, f_eq_ant);

        double trt_sim = prm.alphaSym *  ( f_eq_sim[0]- 2.0 * f[0]);
		double trt_ant = 0;

        double F_sim[NUM_OF_VEL];
        double F_ant[NUM_OF_VEL];

        forceTermSym(F_sim, F[0] ,  F[1] ,  F[2]  ,vx,vy,vz, prm.tauSym);
        forceTermAnt(F_ant, F[0] ,  F[1] ,  F[2]  ,vx,vy,vz, prm.tauAnt);

        f[0] = f[0] + trt_sim + F_sim[0] ;

		//  Performe the collision for all the directions. The collision is performed explicitly only for odd directions (1,3,5 ,...)
		//  For even directions the symmtrie property is used. The final value is store in the opposite direction (pre-streaming process).

		double temp;
		for (int i = 1; i< NUM_OF_VEL; i = i + 2)
		{
            trt_sim = prm.alphaSym *  ( f_eq_sim[i]   - f[i] - f[i+1] );
            trt_ant = prm.alphaAnt *  ( f_eq_ant[i]   - f[i] + f[i+1] );
            temp = f[i] + trt_sim + trt_ant + F_sim[i] + F_ant[i];
            f[i] = f[i+1] + trt_sim - trt_ant + F_sim[i] - F_ant[i];
            f[i+1] = temp;
		}
	}

}

}
