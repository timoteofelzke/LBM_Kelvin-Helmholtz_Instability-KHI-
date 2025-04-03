/**
 *  \file default.cpp
 *
 *  Description: A empty problem (default) to be used as template to implement new problems
 *  into the code
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <stdexcept>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include "../lbmbase/geometry.h"
#include "../parameters.h"
#include "../lbmbase/boundary/neumann.h"
#include "../lbmbase/lattice/lattice.h"

#include "blank.h"

using namespace std;

namespace blankProblem
{

stProblem problemParameters(string filename) {
    
    stProblem prm;
    prm.filename = "dados.dat";

    using boost::property_tree::ptree;
    ptree pt;
    read_ini(filename.c_str(),pt);
    prm.U0 = pt.get<double>("KelvinHelmholtz.U0");
    prm.lbd = pt.get<double>("KelvinHelmholtz.lbd");
    prm.eps = pt.get<double>("KelvinHelmholtz.eps");
    return prm;
}

void problemInitialize(double* iniN, stCollision& col, stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force)
{
       int L = geo.nx;     
       double Re = sp.U0 * L / col.kinematicViscosity;
       double rho = gp.initialRho;
       cout << "Reynolds Number: "  << Re << endl;

       double y_quarter = geo.ny / 4.0;  
       double y_3quarter = 3.0 * geo.ny / 4.0; 
   
       for (int y = 0; y < geo.ny; ++y) {
           for (int x = 0; x < geo.nx; ++x) {
               int idx = (y * geo.nx + x) * NUM_OF_VEL;
               
               double u;
               if (y <= geo.ny/2) {
                   u = sp.U0 * tanh(sp.lbd * (y/(double)geo.ny - 0.25));  
               } else {
                   u = sp.U0 * tanh(sp.lbd * (0.75 - y/(double)geo.ny));  
               }
               
               double phase = (x/(double)geo.nx) + 0.25;          
               double v_perturb = sp.eps * sp.U0 * sin(2.0 * M_PI * phase); 
   
               equilibriumDistribution(u, v_perturb, 0.0, rho, iniN + idx);
           }
       }
   }

void problemFinalize(double* iniN, stCollision& col , stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force, int step )
{
    // Finalização do problema (pode ser usado para cálculos finais ou logs)
}

void problemBoundary(double* iniN, stCollision& col, stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force, int step)
{
   
}

void problemPosCollision(double* iniN, stCollision& col , stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force, int step )
{
    // Atualizações pós-colisão, se necessário
}

bool problemOutput(double* iniN, stCollision& col , stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force, int step )
{
    return true;
}


}

