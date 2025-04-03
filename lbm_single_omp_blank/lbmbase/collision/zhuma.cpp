

#include <iostream>
#include <string>
#include <iomanip>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include "zhuma.h"
#include "../force/force.h"
#include "../lattice/lattice.h"

using namespace std;

namespace zhuma
{

stCollision collisionParameters(std::string filename)
{
    using boost::property_tree::ptree;
    stCollision prm;

    ptree pt;
    read_ini(filename.c_str(),pt);
    prm.tau          =  pt.get<double>("ZHUMA.tau");
    prm.alphaEq  =  1.0/prm.tau;
    prm.alphaNon =  (1.0 - 1.0/prm.tau);
    prm.kinematicViscosity = c_s2 * (prm.tau - 0.5);
    return prm;
}

void reportCollision(stCollision& prm)
{
    cout << endl;
    cout << "Collision Model: Zhu and Ma (gLBM)";
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
        double* f =  iniN + n*NUM_OF_VEL;
        double vx, vy, vz, rho;
        calculateMacro( f, vx, vy, vz, rho );

        double f_eq[NUM_OF_VEL];

        equilibriumDistribution( vx, vy , vz, rho, f_eq);

        for (int i = 0; i< NUM_OF_VEL; i = i + 1)
        {
            f[i] = prm.alphaNon * f[i] + prm.alphaEq * f_eq[i];
        }

    }
}

void collisionSwap( double* iniN, stForce& force,  stCollision& prm, int numberOfPoints)
{
    double* ns = (double *) prm.data;

    #pragma omp parallel for
    for (int n = 0; n < numberOfPoints; n++)
    {
        double* f = iniN + n*NUM_OF_VEL;
        double vx, vy, vz, rho;

        double* F = force.constant ? force.constantForce : force.forceField + n * NUM_OF_DIM;

        calculateMacroEquilibrium( f, vx, vy, vz, rho, F[0] , F[1]  , F[2]  );

        double f_eq[NUM_OF_VEL];
        double F_sim[NUM_OF_VEL];
        double F_ant[NUM_OF_VEL];

        forceTermSym(F_sim, F[0] ,  F[1] ,  F[2]  ,vx,vy,vz, prm.tau);
        forceTermAnt(F_ant, F[0] ,  F[1] ,  F[2]  ,vx,vy,vz, prm.tau);

        equilibriumDistribution( vx, vy , vz, rho, f_eq);

        f[0] = prm.alphaNon  * f[0] + prm.alphaEq * f_eq[0] + F_sim[0] ;

        double temp;
        for (int i = 1; i< NUM_OF_VEL; i = i + 2)
        {
            double fa = prm.alphaNon * f[i] + prm.alphaEq * f_eq[i] + F_sim[i] + F_ant[i];
            double fb = prm.alphaNon * f[i+1] +  prm.alphaEq  * f_eq[i+1] + F_sim[i] - F_ant[i];
            double df = (*ns) *(fb - fa);
            temp = fa + df;
            f[i] = fb - df;
            f[i+1] = temp;
        }
        ns++;
    }
}

}

