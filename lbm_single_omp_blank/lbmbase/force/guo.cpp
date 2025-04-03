/**
 *  @file   d3q19.cpp
 *  @author Diogo Nardelli Siebert
 *  @brief  Implementation of Functions for the D3Q19 lattice
 *
 *  Implementation of functions to compute the discrite equilibrium  distribuition and the macroscopic variables for the D3Q19 lattice
 */

#include "../lattice/lattice.h"
#include "guo.h"

namespace d3q19
{

namespace guoForce
{

void forceTermSym(double *F, double& Fx, double&Fy, double&Fz, double& vx, double& vy, double& vz, double& tau)
{
    double lambda = ( 1. - 1/(2*tau) ) *one_over_c_s4;

    double zeroOrder = -  c_s2 * ( Fx*vx + Fy * vy + Fz*vz);
    F[0] = w0 * lambda  * zeroOrder;

    double lambda1 = w1 * lambda;

    F[1] = lambda1 * ( zeroOrder + Fx*vx );
    F[3] = lambda1 * ( zeroOrder + Fy*vy );
    F[5] = lambda1 * ( zeroOrder + Fz*vz );

    double lambda2 = w2 * lambda;
    F[7] = lambda2 * ( zeroOrder + (Fx + Fy)*(vx + vy) );
    F[9] = lambda2 * ( zeroOrder + (Fx - Fy)*(vx - vy) );
    F[11] = lambda2 * ( zeroOrder + (Fx + Fz)*(vx + vz) );
    F[13] = lambda2 * ( zeroOrder + (Fx - Fz)*(vx - vz));
    F[15] = lambda2 * ( zeroOrder + (-Fy - Fz)*(-vy - vz));
    F[17] = lambda2 * ( zeroOrder + (-Fy + Fz)*(-vy + vz) );
}

void forceTermAnt(double *F, double& Fx, double&Fy, double&Fz, double& vx, double& vy, double& vz, double& tau)
{
    double lambda = ( 1. - 1/(2*tau) ) *one_over_c_s2;

    F[0] = 0;

    double lambda1 = w1 * lambda;
    F[1] = lambda1 * ( Fx );
    F[3] = lambda1 * ( Fy );
    F[5] = lambda1 * ( Fz );

    double lambda2 = w2 * lambda;
    F[7] =  lambda2 * (Fx + Fy);
    F[9] =  lambda2 * (Fx - Fy);
    F[11] = lambda2 * (Fx + Fz);
    F[13] = lambda2 * (Fx - Fz);
    F[15] = lambda2 * (-Fy - Fz);
    F[17] = lambda2 * (-Fy + Fz);
}

void forceTerm(double *F, double& Fx, double&Fy, double&Fz, double& vx, double& vy, double& vz, double& tau)
{
    double lambda = ( 1. - 1/(2*tau) ) *one_over_c_s4;

    double zeroOrder = -  c_s2 * ( Fx*vx + Fy * vy + Fz*vz);
    F[0] = w0 * lambda  * zeroOrder;

    double lambda1 = w1 * lambda;

    double sym = zeroOrder + Fx*vx;
    double ant = Fx;
    F[1] = lambda1 * ( sym + ant );
    F[2] = lambda1 * ( sym - ant );

    sym = zeroOrder + Fy*vy;
    ant = Fy;
    F[3] = lambda1 * (sym + ant);
    F[4] = lambda1 * (sym - ant);

    sym = zeroOrder + Fz*vz;
    ant = Fz;
    F[5] = lambda1 * (sym + ant );
    F[6] = lambda1 * (sym - ant );

    double lambda2 = w2 * lambda;
    sym = zeroOrder + (Fx + Fy)*(vx + vy);
    ant = Fx + Fy;
    F[7] = lambda2 * (sym + ant);
    F[8] = lambda2 * (sym - ant);

    sym = zeroOrder + (Fx - Fy)*(vx - vy);
    ant =(Fx - Fy);
    F[9] = lambda2 * (sym + ant);
    F[10] = lambda2 * (sym - ant);

    sym = zeroOrder + (Fx + Fz)*(vx + vz);
    ant =(Fx + Fz);
    F[11] = lambda2 * (sym + ant);
    F[12] = lambda2 * (sym - ant);

    sym = zeroOrder + (Fx - Fz)*(vx - vz);
    ant = (Fx - Fz);
    F[13] = lambda2 * (sym + ant);
    F[14] = lambda2 * (sym - ant);

    sym = zeroOrder + (-Fy - Fz)*(-vy - vz);
    ant = (-Fy - Fz);
    F[15] = lambda2 * (sym + ant);
    F[16] = lambda2 * (sym - ant);

    sym = zeroOrder + (-Fy + Fz)*(-vy + vz);
    ant = (-Fy + Fz);
    F[17] = lambda2 * (sym + ant);
    F[18] = lambda2 * (sym - ant);
}


void calculateMacroEquilibrium( double *f, double& vx, double& vy, double& vz, double& rho , double& Fx, double &Fy, double& Fz)
{
    calculateMacro(f,vx,vy,vz,rho);
    vx += 0.5*Fx/rho;
    vy += 0.5*Fy/rho;
    vz += 0.5*Fz/rho;
}

void calculateMacroOutput( double *f, double& vx, double& vy, double& vz, double& rho, double& Fx, double &Fy, double& Fz)
{
    calculateMacro(f,vx,vy,vz,rho);
    vx += 0.5*Fx/rho;
    vy += 0.5*Fy/rho;
    vz += 0.5*Fz/rho;
}


}

}
