/**
 *  @file   force.h
 *  @author Diogo Nardelli Siebert
 *  @brief  Declaration of constant and functions for the lattice D3Q19
 *
 *  This file contains explicite formulas for the force term proposed by Guo, Zheng and Shi.
 *
 */

#ifndef __FORCE_H_INCLUDED__
#define __FORCE_H_INCLUDED__

#include "../lattice/lattice.h"
#include "guo.h"

struct stGeometry; // Foward Declared

using namespace guoForce;

struct stForce
{
    bool constant;
    double* forceField;
    double  constantForce[3];
};

stForce setConstantForce(double& forceX, double& forceY, double& forceZ );

#endif
