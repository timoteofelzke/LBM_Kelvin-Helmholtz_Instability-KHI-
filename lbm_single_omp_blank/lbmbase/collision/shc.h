/**
 *  @file   shc.h
 *  @author Diogo Nardelli Siebert
 *
 *  Header of functions needed to apply the collision operator
 *  proposed by Spencer, Halliday & Care. The models is base
 *  on the Rothmann and Keler colocar gradient model.
 *
 *  For more details.
 *
 *  [1] Spencer, T. J., Halliday, I., & Care, C. M. (2010). Lattice Boltzmann equation method for multiple immiscible continuum fluids. Physical Review. E, Statistical, Nonlinear, and Soft Matter Physics, 82(6 Pt 2), 66701. http://doi.org/10.1103/PhysRevE.82.066701
 *
 *  [2] Spencer, T. J., Halliday, I., & Care, C. M. (2011). A local lattice Boltzmann method for multiple immiscible fluids and dense suspensions of drops. Philosophical Transactions. Series A, Mathematical, Physical, and Engineering Sciences, 369(1944), 2255–63. http://doi.org/10.1098/rsta.2011.0029
 */

#ifndef __SHC_H_INCLUDED__   // if collision.h hasn't been included yet...
#define __SHC_H_INCLUDED__

#include <string>

struct stGeometry;  // Foward declare
struct stParameters; // Forward declare
struct stForce; // Foward declare

namespace d3q19
{

/**
 *  @brief Computes the local gradient using the lattice neighbors
 *
 *
 *
 *  @param  grad      pointer where the values of the gradient will be stored.
 *  @param  f         pointer to the array containing the value of the quantity in the lattice neighbors.
 */

void getGradient(double* grad, double* f);

/**
 *  @brief Computes the force term due to the interfacial tension
 *
 *  \f{equation}{
 *  F_i = \rho w_i \left( \frac{ \alpha | \nabla w_r | }{c_s^2 \tau} \right) ( n_\alpha n_\beta - \delta_{\alpha \beta} ) (c_{i,\alpha} c_{i,\beta} - c_s^2 \delta_{\alpha \beta} )
 *  \f}
 *
 *  Using this force term the surface tension will be given by
 *
 *  \f{equation}{
 *  \gamma = 2 \rho c_s^2 \alpha
 *  \f}
 *
 *  @param  F         pointer where the force distribution will be stored
 *  @param  normal    normalized red concetrantion gradient
 *  @param  rho       density
 *  @param  coeff     the term between parentheses in the above equation
 */

void shcForce(double* F, double* normal,double& rho,  double coeff );
void shcRecoloring(double *f, double *R, double *B, double* normal, double rho, double w_r, double beta);
}

/**
 *  Compute the collision in all domain  using the two relaxation time (TRT) model
 *  (I. Ginzburg, Adv Water Resour, 28, 2005)
 *
 *  @param  pDistribution   pointer to the the distribution function.
 *  @param  symTau          relaxation time of the the collision of the symmetric distribution function
 *  @param  antTau          relaxation time of the the collision of the antisymmetric distribution function
 *  @param  gx              accelaration in the x direction as result of an external force field.
 *  @param  gy              accelaration in the x direction as result of an external force field.
 *  @param  gz              accelaration in the x direction as result of an external force field.
 *  @param  numberOfPoints  total number of sites containing fluids in the domain.
 *  @param  sumMx           reference to the variable storing the total momentum in the x direction
 *  @param  sumForce        reference to the variable storing the total force in the x direction
 */

namespace shcTrt
{

    void generateWr(double* ini_W,double *ini_R, double* ini_B, int numberOfPoints);

struct stCollision
{
    double tauSymRed;
    double tauAntRed;
    double alphaSymRed;
    double alphaAntRed;
    double tauSymBlue;
    double tauAntBlue;
    double alphaSymBlue;
    double alphaAntBlue;
    double wrWall;

    double coeffRed;
    double coeffBlue;

    double alpha;
    double beta;

    double surfaceTension;
    double interfaceWidth;
    double kinematicViscosityRed;
    double kinematicViscosityBlue;
    void* data;
};


stCollision collisionParameters(std::string filename);

void reportCollision(stCollision& prm);

void collisionSwap(double* ini_R, double* ini_B, double* ini_W,  stForce& Fred, stForce& Fblue,  stCollision& prm, int numberOfPoints);

}




#endif
