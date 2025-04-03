/**
 *  @file   d3q19.h
 *  @author Diogo Nardelli Siebert
 *  @brief  Declaration of constant and functions for the lattice D3Q19
 *
 *  This file contains the definition of D3Q19 velocity vectors and
 *  weights and the header of functions to compute the discrite equilibrium
 *  distribuition and the macroscopic variables for the D3Q19 lattice
 *
 *  @image html d3q19.png
 *
 */

#ifndef __LATTICE_D3Q19_H_INCLUDED__
#define __LATTICE_D3Q19_H_INCLUDED__

namespace d3q19
{
	
const int NUM_OF_VEL = 19;    /*!< Number of velocites of the lattice */
const int NUM_OF_DIM = 3;     /*!< Dimension of the lattice  */

const double w0 = 0.3333333333333333333333333333333333333333333333333333;       /*!< Weight of first  set of velocities i = {0} with \f$ |\vec{c}_i| = 0 \f$ */
const double w1 = 0.0555555555555555555555555555555555555555555555555555;       /*!< Weight of second set of velocities i = {1, .. ,6}   with \f$ |\vec{c}_i| = 1 \f$ */
const double w2 = 0.0277777777777777777777777777777777777777777777777777;       /*!< Weight of third  set of velocities i = {7, .. ,18}  with \f$ |\vec{c}_i| = \sqrt{2} \f$ */

const double one_over_c_s2 = 3.0;           /*!< Inverse of the lattice sound speed squared*/
const double one_over_c_s4 = 9.0;           /*!< Inverse of 4th power the lattice sound speed*/
const double c_s2 = 1.0/one_over_c_s2;      /*!< Sound speed lattice squared */
const double c_s4 = 1.0/one_over_c_s4;      /*!< 4th power of the Sound speed lattice */

//                  			 0   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15  16  17  18
const int i_op[NUM_OF_VEL] =  {  0,  2,  1,  4,  3,  6,  5,  8,  7,  10,  9, 12, 11, 14, 13, 16, 15, 18, 17};     /*!< Number of the opposite lattice velocity */

//                             0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18
const int cx[NUM_OF_VEL] =  {  0,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0};     /*!< X component of the lattice velocities */
const int cy[NUM_OF_VEL] =  {  0,  0,  0,  1, -1,  0,  0,  1, -1, -1,  1,  0,  0,  0,  0, -1,  1, -1,  1};     /*!< Y component of the lattice velocities */
const int cz[NUM_OF_VEL] =  {  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1, -1, -1,  1, -1,  1,  1, -1};     /*!< Z component of the lattice velocities */

const double w[NUM_OF_VEL]  =  {w0, w1, w1, w1, w1, w1, w1, w2, w2, w2, w2, w2, w2, w2, w2, w2, w2, w2, w2};

/**
 * @brief Compute the equilibrium Distribution for D3Q19 lattice
 *
 * @param vx macroscopic velocity in the x direction
 * @param vy macroscopic velocity in the y direction
 * @param vz macroscopic velocity in the z direction
 * @param rho macroscopic density
 * @param feq pointer where the equilibrium distribution function computed will written.
 *
 * Computes the second order equilibrium distribution for the D3Q19 lattice using the given macroscopic
 * variables.
 *
 *  \f{equation}{
 *  f^{eq}_i = \rho w_i \left[ 1 + 3 (\vec{c_i} \cdot \vec{u}) + \frac{9}{2} (\vec{c_i} \cdot \vec{u})^2 - \frac{3}{2} u^2 \right]
 *  \f}
 *
 * In order to get a faster algorithm, this function do not use loops in
 * in the different velocities using explicit formulas instead in order to avoid unnecessary computation.
 */

void equilibriumDistribution   ( double vx, double vy, double vz, double rho, double* feq);

/**
 * @brief Compute the Symetric Equilibrium Distribution (TRT) for D3Q19 lattice
 * 
 * @param vx macroscopic velocity in the x direction
 * @param vy macroscopic velocity in the y direction
 * @param vz macroscopic velocity in the z direction
 * @param rho macroscopic density
 * @param feq pointer where the equilibrium distribution function computed will written.
 * 
 * Computes the symmetric equilibrium distribution for the D3Q19 lattice using the given macroscopic
 * variables. This function is used in the Two Relaxation Collision Model (TRT) [1].
 *
 *  \f{equation}{
 *  f^{eq}_{s}(\vec{c}_i) = \frac{f^{eq}(\vec{c}_i)  + f^{eq}(-\vec{c}_i)}{2}
 *  \f}
 *
 * In order to get a faster algorithm, this function do not use loops in the discrite velocities to avoid unnecessary computation using instead
 * explicit formulas to each direction. Be carefull, since this distribution is symmetric \f$ f^{eq}(\vec{c}_i) = f^{eq}(-\vec{c}_i) \f$ this function is only computed
 * for directions with even indexes for optimization.
 *
 * [1] Ginzburg, Irina. "Equilibrium-type and link-type lattice Boltzmann models for generic advection and anisotropic-dispersion equation." Advances in Water resources 28.11 (2005): 1171-1195.
 *
 **/

void equilibriumDistributionSym( double vx, double vy, double vz, double rho, double* feq);

/**
 * @brief Compute the Antisymmetric Equilibrium Distribution (TRT) for D3Q19 lattice
 *
 * @param vx macroscopic velocity in the x direction
 * @param vy macroscopic velocity in the y direction
 * @param vz macroscopic velocity in the z direction
 * @param rho macroscopic density
 * @param feq pointer where the equilibrium distribution function computed will written.
 *
 * Computes the antisymmetric equilibrium distribution for the D3Q19 lattice using the given macroscopic
 * variables. This function is used in the Two Relaxation Collision Model (TRT) [1].
 *
 *  \f{equation}{
 *  f^{eq}_{a}(\vec{c}_i) = \frac{f^{eq}(\vec{c}_i) -  f^{eq}(-\vec{c}_i)}{2}
 *  \f}
 *
 * In order to get a faster algorithm, this function do not use loops in the discrite velocities to avoid unnecessary computation using instead
 * explicit formulas to each direction. Be carefull, since this distribution is antisymmetric \f$ f^{eq}(\vec{c}_i) = - f^{eq}(-\vec{c}_i) \f$ this function is only computed
 * for directions with even indexes for optimization.
 *
 * [1] Ginzburg, Irina. "Equilibrium-type and link-type lattice Boltzmann models for generic advection and anisotropic-dispersion equation." Advances in Water resources 28.11 (2005): 1171-1195.
 **/

void equilibriumDistributionAnt( double vx, double vy, double vz, double rho, double* feq);

/**
 *  @brief  Compute macroscopic variables from the distribution function for the D3Q19 lattice.
 *
 *  @param  f         ṕointer to the place in memory where the distribution function
 *                      is stored.
 *  @param  vx        reference to double where the computed x component of the velocity will be stored.
 *  @param  vy        reference to double where the computed y component of the velocity will be stored.
 *  @param  vz        reference to double where the computed z component of the velocity will be stored.
 *  @param  rho       reference to double where the computed density will be stored.
 *
 *  Compute macroscopic variables from the distribution function for the D3Q19 lattice, i. e.,
 *
 *  \f{equation}{
 *  \rho = \sum_i f_i
 *  \f}
 *
 *  \f{equation}{
 *  \rho \vec{u} = \sum_i f_i \vec{c}_i
 *  \f}
 **/

void calculateMacro( double *f, double& vx, double& vy, double& vz, double& rho );

/**
 *  @brief  Computes the mass density from a given distribution function for the D3Q19 lattice.
 *
 *  @param  f           pointer to the (D3Q19) distribution function 
 *
 *  @return             the fluid mass density
 *
 *  Computes the mass density from a given distribution function for the D3Q19 lattice.
 *
 *  \f{equation}{
 *  \rho = \sum_i f_i
 *  \f}
 *
 **/

double calculateMass ( double* f );

/**
 *  @brief  Computes all the moments of the distribution function up to 2nd Order for D3Q19 lattice.
 *
 *  @param  f           pointer to the (D3Q19) distribution function
 *  @param  m           pointer to the memory region where the values of the moment will be stored.
 *
 *  Computes all the momentos of the distribution function up to 2nd Order for D3Q19 lattice, i. e, \f$ (m^{(0)}, m^{(1)}_x,m^{(1)}_y,m^{(1)}_z,m^{(2)}_{xx},m^{(2)}_{yy}, m^{(2)}_{zz},m^{(2)}_{xy},m^{(2)}_{xz}, m^{(2)}_{yz}) \f$.
 *  This Function is needed for the regularized collision.
 *
 *  \f{equation}{
 *  m^{(0)} =  \sum_i f_i
 *  \f}
 *
 *  \f{equation}{
 *  m^{(1)}_{\alpha} =  \sum_i f_i c_{i,\alpha}
 *  \f}
 *
 *  \f{equation}{
 *  m^{(2)}_{\alpha,\beta} =  \sum_i f_i c_{i,\alpha} c_{i,\beta}
 *  \f}
 **/

void calculateMoment2nd( double *f, double* m);

/**
 *  @brief  Computes the distribution function up to 2nd Order for the D3Q19 lattice with the moments up this order
 *
 *  @param  f           pointer to the position of memory where the (D3Q19) distribution function will be stored
 *  @param  m           pointer to array with all moments up to 2nd order.
 *
 *  Computes all the moments of the distribution function up to 2nd Order for D3Q19 lattice, i. e, \f$ (m^{(0)}, m^{(1)}_x,m^{(1)}_y,m^{(1)}_z,m^{(2)}_{xx},m^{(2)}_{yy}, m^{(2)}_{zz},m^{(2)}_{xy},m^{(2)}_{xz}, m^{(2)}_{yz}) \f$.
 *  This Function is needed for the regularized collision.
 *
 *  \f{equation}{
 *  f_i = w_i \left( -\frac{1}{6} m^{(2)}_{\alpha \alpha} + \frac{1}{18} c_{i,\alpha} m^{(2)}_{\alpha \beta} {c}_{i,\beta} -\frac{1}{6} m^{(0)} {c}_i^2+ \frac{1}{2}\text{D} m^{(0)} + \frac{1}{3} m^{(1)}_\alpha c_{i,\alpha}+m^{(0)} \right)
 *  \f}
 **/

void reconstructedDistribution2nd( double *f, double* m);

/**
 * @brief   Computes the moments used with the MRT collision operator for the D3Q19
 *
 * @param  f           pointer to the (D3Q19) distribution function
 * @param  m           pointer to the memory region where the values of the moments will be stored.
 *
 * Computes the moments used with the MRT collision operator for the D3Q19. The explicit formulas where obtained by Diogo Siebert using mathematica.
 *
 * https://www.dropbox.com/s/dm9pe1mkqot32xe/MRT_D3Q19.nb?dl=0
 *
 * For details see [1].
 *
 * [1]  D’Humières, D., Ginzburg, I., Krafczyk, M., Lallemand, P., & Luo, L.-S. (2002). Multiple-relaxation-time lattice Boltzmann models in three dimensions. Philosophical Transactions. Series A, Mathematical, Physical, and Engineering Sciences, 360(1792), 437–51.
 * http://doi.org/10.1098/rsta.2001.0955
 *
 **/

void calculateMomentsMRT( double *f, double* m);

/**
 * @brief  Convert the values of the moments in to a distribution function for the MRT method (D3Q19 lattice)
 *
 * @param  f           pointer to the position of memory where the (D3Q19) distribution function will be stored.
 * @param  m           pointer to array with all 19 (MRT) moments of the D3Q19.
 *
 * Computes the moments used with the MRT collision operator for the D3Q19. For details see [1]. The explicit formulas where obtained by Diogo Siebert using mathematica.
 *
 * https://www.dropbox.com/s/dm9pe1mkqot32xe/MRT_D3Q19.nb?dl=0
 *
 * For details see [1].
 *
 * [1]  D’Humières, D., Ginzburg, I., Krafczyk, M., Lallemand, P., & Luo, L.-S. (2002). Multiple-relaxation-time lattice Boltzmann models in three dimensions. Philosophical Transactions. Series A, Mathematical, Physical, and Engineering Sciences, 360(1792), 437–51.
 * http://doi.org/10.1098/rsta.2001.0955
 *
 **/

void reconstructedDistributionMRT( double *f, double* m);

/**
 * @brief  Computes the values of the all equilibrium moments used in MRT for a given set of non equilibrium moments for the D3Q19 lattice.
 *
 * @param  meq           pointer to the position of memory where the 19 equilibrium moments of the D3Q19 lattice will be stored. (output)
 * @param  m             pointer to the position of memory with 19 non-equilibrium moments. (input)
 *
 * Computes the moments used with the MRT collision operator for the D3Q19. In MRT, the number of necessary moments is equal to the number of velocities
 * of the lattice. For details see [1].
 *
 * The explicit formulas where obtained by Diogo Siebert using mathematica.
 *
 * https://www.dropbox.com/s/dm9pe1mkqot32xe/MRT_D3Q19.nb?dl=0
 *
 * For details see [1].
 *
 * [1]  D’Humières, D., Ginzburg, I., Krafczyk, M., Lallemand, P., & Luo, L.-S. (2002). Multiple-relaxation-time lattice Boltzmann models in three dimensions. Philosophical Transactions. Series A, Mathematical, Physical, and Engineering Sciences, 360(1792), 437–51.
 * http://doi.org/10.1098/rsta.2001.0955
 *
 **/

void equilibriumMomentsMRT(double* meq,double* m);

/**
 * @brief  Computes the inverse of the relaxation time for each of the 19 moments used in the MRT for D3Q19 lattice.
 *
 * @param  tau           The standart relaxation time which is related with the viscosity
 *
 * Computes the inverse of the relaxation time for each of the 19 moments used in the MRT for D3Q19 lattice. For details see [1].
 *
 * [1]  D’Humières, D., Ginzburg, I., Krafczyk, M., Lallemand, P., & Luo, L.-S. (2002). Multiple-relaxation-time lattice Boltzmann models in three dimensions. Philosophical Transactions. Series A, Mathematical, Physical, and Engineering Sciences, 360(1792), 437–51.
 * http://doi.org/10.1098/rsta.2001.0955
 *
 **/

double* relaxationMRT(double tau);

/**
 *  @fn     calculateMomentum
 *
 *  Computes the the momentum density from a given distribution function for the D3Q19 lattice.
 *
 *  @param  f           pointer to a (D3Q19) distribution function 
 *  @param  m   pointer to place where the momentum values will be stored
 */

void calculateMomentum ( double* f , double* m );
/**
 *  @fn     calculateMomentum
 *
 *  Computes the the momentum density from a given distribution function for the D3Q19 lattice.
 *
 *  @param  f           pointer to a (D3Q19) distribution function 
 *  @param  m   pointer to place where the momentum values will be stored
 */

void calculateMomentum ( double* f , double* m );

}
#endif
