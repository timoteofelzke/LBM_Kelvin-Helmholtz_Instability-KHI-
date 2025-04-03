/**
 *  @file   zhuma.h
 *  @author Diogo Nardelli Siebert
 *  @brief  Header of functions needed in the collision step for the gray Lattice boltzmann model of Zhu and Ma
 *
 *  [1] Zhu, J., & Ma, J. (2013). An improved gray lattice Boltzmann model for simulating fluid flow in multi-scale porous media. Advances in Water Resources, 56, 61–76. https://doi.org/10.1016/j.advwatres.2013.03.001
 */


#ifndef __ZHUMA_H__
#define __ZHUMA_H__

#include <string>

struct stForce; // Foward Declare

namespace zhuma
{

struct stCollision
{
    double  tau;      /* Relaxation Time */
    double  alphaEq;
    double  alphaNon;
    double  kinematicViscosity;
    void* data;       /* Fraction of reflected particles */
};

stCollision collisionParameters(std::string filename);
void reportCollision(stCollision& prm);

void collision( double* iniN, stForce& force,  stCollision& prm, int numberOfPoints);
void collisionSwap( double* iniN, stForce& force,  stCollision& prm, int numberOfPoints);

}

#endif // ZHUMA_H
