/**
 *  @file   bgk.h
 *  @author Diogo Nardelli Siebert
 *  @brief  Header of functions to calculate the collision using BGK operator
 */

#ifndef __BGK_H_INCLUDED__   // if collision.h hasn't been included yet...
#define __BGK_H_INCLUDED__

#include <string>

struct stForce;

namespace bgk
{

struct stCollision
{
    double tau;
    double alphaEq;
    double alphaNon;
    double kinematicViscosity;
    void* data;
};

stCollision collisionParameters(std::string filename);

void reportCollision(stCollision& prm);

/**
 *  Compute the collision in all domain  using the Bhatnagar–Gross–Krook (BGK) model
 *
 *  \f$ f_i = \left(1-\frac{1}{\tau} \right) f_i^{'} + \frac{1}{\tau} f_i^{eq} \f$
 *
 *  [1] Qian, Y. H., Dominique d'Humières, and Pierre Lallemand. "Lattice BGK models for Navier-Stokes equation." EPL (Europhysics Letters) 17.6 (1992): 479.
 *
 *  @param  pDistribution   pointer to the the distribution function.
 *  @param  tau             relaxation time of the the collision of the symmetric distribution function
 *  @param  gx              accelaration in the x direction as result of an external force field.
 *  @param  gy              accelaration in the x direction as result of an external force field.
 *  @param  gz              accelaration in the x direction as result of an external force field.
 *  @param  numberOfPoints  total number of sites containing fluids in the domain.
 */

void collision( double* iniN, stForce& force,  stCollision& prm, int numberOfPoints);

/**
 *  Compute the collision in all domain  using the Bhatnagar–Gross–Krook (BGK) model and swap the distribution
 *  in oposite direction after the collision. This last step is necessary for the swap streaming method.
 *
 *  \f$ f_i = \left(1-\frac{1}{\tau} \right) f_i^{'} + \frac{1}{\tau} f_i^{eq} \f$
 *
 *  [1] Qian, Y. H., Dominique d'Humières, and Pierre Lallemand. "Lattice BGK models for Navier-Stokes equation." EPL (Europhysics Letters) 17.6 (1992): 479.
 *
 *  [2] Mattila, Keijo, et al. "An efficient swap algorithm for the lattice Boltzmann method." Computer Physics Communications 176.3 (2007): 200-210.
 *
 *  @param  pDistribution   pointer to the the distribution function.
 *  @param  tau             relaxation time of the the collision of the symmetric distribution function
 *  @param  gx              accelaration in the x direction as result of an external force field.
 *  @param  gy              accelaration in the x direction as result of an external force field.
 *  @param  gz              accelaration in the x direction as result of an external force field.
 *  @param  numberOfPoints  total number of sites containing fluids in the domain.
 */


void collisionSwap( double* iniN, stForce& force,  stCollision& prm, int numberOfPoints);

}
#endif
