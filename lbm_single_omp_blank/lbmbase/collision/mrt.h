/**
 *  @file   mrt.h
 *  @author Diogo Nardelli Siebert
*   @brief  Header of functions needed in the collision step usign the MRT collision operator .
 */

#ifndef __MRT_H_INCLUDED__   // if collision.h hasn't been included yet...
#define __MRT_H_INCLUDED__

struct stForce;

namespace mrt
{

struct stCollision
{
    double tau;
    double* s;
    double kinematicViscosity;
    void* data;
};


stCollision collisionParameters(std::string filename);
void reportCollision(stCollision& prm);

/**
 *
 * @brief Multiple Relaxation Time collision (MRT).
 * @param iniN pointer to the the distribution function.
 * @param tau relaxation time parameter of the BGK collision
 * @param gx external accelaration in the x-axis due to an external force
 * @param gy external accelaration in the y-axis due to an external force
 * @param gz external accelaration in the z-axis due to an external force
 * @param numberOfPoints number of (domain/fluid) sites
 *
 * Computes the collision for the given number of sites stored in the iniN pointer using the Multiple Relaxation Time (MRT).
 *
 * [1] D’Humières, D., Ginzburg, I., Krafczyk, M., Lallemand, P., & Luo, L.-S. (2002). Multiple-relaxation-time lattice Boltzmann models in three dimensions. Philosophical Transactions. Series A, Mathematical, Physical, and Engineering Sciences, 360(1792), 437–51.
 * http://doi.org/10.1098/rsta.2001.0955
 *
 */

void collision( double* iniN, stForce& force,  stCollision& prm, int numberOfPoints);

/**
 *
 * @brief Multiple Relaxation Time collision (MRT) with with Pre-Swap Streaming.
 * @param iniN pointer to the the distribution function.
 * @param tau relaxation time parameter of the BGK collision
 * @param gx external accelaration in the x-axis due to an external force
 * @param gy external accelaration in the y-axis due to an external force
 * @param gz external accelaration in the z-axis due to an external force
 * @param numberOfPoints number of (domain/fluid) sites
 *
 * Computes the collision for the given number of sites stored in the iniN pointer using the Multiple Relaxation Time (MRT)
 * collision and swap the distribution in oposite direction after the collision. This last step is necessary for the swap streaming method.
 *
 * [1] D’Humières, D., Ginzburg, I., Krafczyk, M., Lallemand, P., & Luo, L.-S. (2002). Multiple-relaxation-time lattice Boltzmann models in three dimensions. Philosophical Transactions. Series A, Mathematical, Physical, and Engineering Sciences, 360(1792), 437–51.
 * http://doi.org/10.1098/rsta.2001.0955
 *
 * [2] Mattila, Keijo, et al. "An efficient swap algorithm for the lattice Boltzmann method." Computer Physics Communications 176.3 (2007): 200-210.
 *
 */

void collisionSwap( double* iniN, stForce& force,  stCollision& prm, int numberOfPoints);

}

#endif
