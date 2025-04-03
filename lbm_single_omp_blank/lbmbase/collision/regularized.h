/**
 *  @file   regularized.h
 *  @author Diogo Nardelli Siebert
 *  @brief  Header of functions needed in the collision step usign the BGK collision operator with a kinection projection (regularization) procedure.
 */

#ifndef __REGULARIZED_H_INCLUDED__
#define __REGULARIZED_H_INCLUDED__

struct stForce;

namespace regularized
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
 * @brief Single relaxation time collision (BGK) with regularization[1].
 * @param iniN pointer to the the distribution function.
 * @param tau relaxation time parameter of the BGK collision
 * @param gx external accelaration in the x-axis due to an external force
 * @param gy external accelaration in the y-axis due to an external force
 * @param gz external accelaration in the z-axis due to an external force
 * @param numberOfPoints number of (domain/fluid) sites
 * 
 * [1] Mattila, Keijo K., Luiz A. Hegele Jr, and Paulo C. Philippi. "Investigation of an entropic stabilizer for the lattice-Boltzmann method." Physical Review E 91.6 (2015): 063010.
 *
 * Computes the collision for the given number of sites stored in the ini_n pointer using the BGK collision 
 * operator with regularization [1].
 */

void collision( double* iniN, stForce& force,  stCollision& prm, int numberOfPoints);

/**
 * @brief Single relaxation time collision (BGK) with regularization[1]. Also performs the swap between oposite directions for swap method of streaming.
 * @param iniN pointer to the the distribution function.
 * @param tau relaxation time parameter of the BGK collision
 * @param gx external accelaration in the x-axis due to an external force
 * @param gy external accelaration in the y-axis due to an external force
 * @param gz external accelaration in the z-axis due to an external force
 * @param numberOfPoints number of (domain/fluid) sites
 *
 * [1] Mattila, Keijo K., Luiz A. Hegele Jr, and Paulo C. Philippi. "Investigation of an entropic stabilizer for the lattice-Boltzmann method." Physical Review E 91.6 (2015): 063010.
 * [2] Mattila, Keijo, et al. "An efficient swap algorithm for the lattice Boltzmann method." Computer Physics Communications 176.3 (2007): 200-210.
 *
 * Computes the collision for the given number of sites stored in the ini_n pointer using the BGK collision
 * operator and performing a kinetic profection (or normalization) procedure. After the collision it stores
 * the value of the distribution function in the opposite direction (pre-streaming of the swap algorithm:
 * http://www.lbmethod.org/palabos/downloads/plb-tr1.pdf)
 */
 
void collisionSwap( double* iniN, stForce& force,  stCollision& prm, int numberOfPoints);

}
#endif
