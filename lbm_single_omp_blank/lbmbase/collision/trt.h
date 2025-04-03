/**
 *  @file   trt.h
 *  @author Diogo Nardelli Siebert
 *  @brief  Header of functions needed in the collision step.
 */

#ifndef __TRT_H_INCLUDED__   // if collision.h hasn't been included yet...
#define __TRT_H_INCLUDED__

struct stForce; // Foward Declare

namespace trt
{

struct stCollision
{
    double tauSym;
    double tauAnt;
    double alphaSym;
    double alphaAnt;
    double kinematicViscosity;
    void* data;
};

stCollision collisionParameters(std::string filename);
void reportCollision(stCollision& prm);

/**
 *  Compute the collision in all domain  using the two relaxation time (TRT) model [1]. This model relay on the concept of decomposition
 *  of the distribuntion function in a symetric and and antisymetric part. Then the collision is performed using a BGK operator for each part.
 *
 *  \f{equation}{
 *  f_{s,i} = \frac{f_i + f_{-i}}{2}
 *  \f}
 *  \f{equation}{
 *  f_{a,i} = \frac{f_i - f_{-i}}{2}
 *  \f}
 *  \f{equation}{
 *  f^{'}_i = \frac{f_{s,i} - f^{eq}_{s,i} }{\tau_s} + \frac{f_{a,i} - f^{eq}_{a,i} }{\tau_a}
 *  \f}
 *
 *  In this model viscosity is calculated by \f$\nu = \frac{1}{c_s^2} (\tau_s - 1/2) \f$. The \f$ \tau_a \f$ is a function of \f$ \tau_s \f$ in order to vanish high order spurious terms. For more details [1].
 *
 *  [1] Ginzburg, Irina. "Equilibrium-type and link-type lattice Boltzmann models for generic advection and anisotropic-dispersion equation." Advances in Water resources 28.11 (2005): 1171-1195.
 *
 *  @param  iniN            pointer to the the distribution function.
 *  @param  tau             relaxation time of the the collision of the symmetric distribution function
 *  @param  gx              accelaration in the x direction as result of an external force field.
 *  @param  gy              accelaration in the x direction as result of an external force field.
 *  @param  gz              accelaration in the x direction as result of an external force field.
 *  @param  numberOfPoints  total number of sites containing fluids in the domain.
 */

void collision( double* iniN, stForce& force,  stCollision& prm, int numberOfPoints);

/**
 *  Compute the collision in all domain  using the two relaxation time (TRT) model [1]. This model relay on the concept of decomposition
 *  of the distribuntion function in a symetric and and antisymetric part. Then the collision is performed using a BGK operator for each part.
 *
 *  \f{equation}{
 *  f_{s,i} = \frac{f_i + f_{-i}}{2}
 *  \f}
 *  \f{equation}{
 *  f_{a,i} = \frac{f_i - f_{-i}}{2}
 *  \f}
 *  \f{equation}{
 *  f^{'}_i = \frac{f_{s,i} - f^{eq}_{s,i} }{\tau_s} + \frac{f_{a,i} - f^{eq}_{a,i} }{\tau_a}
 *  \f}
 *
 *  In this model viscosity is calculated by \f$\nu = \frac{1}{c_s^2} (\tau_s - 1/2) \f$. The \f$ \tau_a \f$ is a function of \f$ \tau_s \f$ in order to vanish high order spurious terms. For more details [1].
 *
 *  [1] Ginzburg, Irina. "Equilibrium-type and link-type lattice Boltzmann models for generic advection and anisotropic-dispersion equation." Advances in Water resources 28.11 (2005): 1171-1195.
 *
 *  @param  iniN            pointer to the the distribution function.
 *  @param  tau             relaxation time of the the collision of the symmetric distribution function
 *  @param  gx              accelaration in the x direction as result of an external force field.
 *  @param  gy              accelaration in the x direction as result of an external force field.
 *  @param  gz              accelaration in the x direction as result of an external force field.
 *  @param  numberOfPoints  total number of sites containing fluids in the domain.
 */

void collisionSwap( double* iniN, stForce& force,  stCollision& prm, int numberOfPoints);

}

#endif
