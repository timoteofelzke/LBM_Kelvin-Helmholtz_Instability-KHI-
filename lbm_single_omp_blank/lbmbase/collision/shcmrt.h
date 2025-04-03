/**
 *  @file   shcmrt.h
 *  @author Lauro José Biesek Junior
 *  @brief  Header of functions needed in the collision step.
 */

#ifndef __SHCMRT_H_INCLUDED__   // if collision.h hasn't been included yet...
#define __SHCMRT_H_INCLUDED__

#include <string>

struct stGeometry;   // Foward Declare
struct stParameters; // Foward Declare
struct stForce;      // Foward Declare

namespace shcMrt
{

struct stCollision
{
    double delta;
    double Beta;
    double A;
    double surfaceTension;
    double kinematicViscosityRed;
    double tauRed;
    double tauBlue;
    double kinematicViscosityBlue;
    double alphaBlue;
    double alphaRed;
    double drWall;
    double s1;
    double s2;
    double s3;
    double t1;
    double t2;
    double t3;
};

stCollision collisionParameters(std::string filename);

void reportCollision(stCollision& prm);

/**
 * @brief
 * @param   ini_R               Pointer to the Red fluid distribution.
 * @param   ini_B               Pointer to the Blue fluid distribution.
 * @param   ini_W               Pointer to the Ghost distribution, which it will be used to calculate the Gradient.
 * @param   numberOfPoints      Number of (domain/fluid) sites.
 *
 * This collision model consists in solving the following equation [1]:
 *      \f{equation}{
 *          f_{i}^{k*} = f_{i}^{k} + (\Omega_{i}^{k})^{1} + (\Omega_{i}^{k})^{2}
 *      \f}
 *
 * Where the first collision operator corresponds to the Multi Relaxation Time Model[2], applied in the blind distribution, given by:
 *      \f{equation}{
 *          f_i = f_{i}^{r} + f_{i}^{b}
 *      \f}
 *
 * And the second collision operator proposed by Reis and Phillips[3], is responsible for inducing the appropriate surface tension term in the macroscopic equations.
 *
 *      \f{equation}{
 *          ( \Omega_{i}^{k} )^{2} = \frac{A}{2}|\textbf{f}|\left[  w_{i} \frac{ (e_{i} \cdot \textbf{f} )^{2} }{| \textbf{f} |^{2}} - B_i \right]
 *      \f}
 *
 * [1]HUANG, Haibo; HUANG, Jun-jie; LU, Xi-yun. Study of immiscible displacements in porous media using a color-gradient-based multiphase lattice Boltzmann method. Computers & Fluids, [s.l.], v. 93, p.164-172, abr. 2014. Elsevier BV. http://dx.doi.org/10.1016/j.compfluid.2014.01.025.
 * http://dx.doi.org/10.1016/j.compfluid.2014.01.025
 *
 * [2] D’Humières, D., Ginzburg, I., Krafczyk, M., Lallemand, P., & Luo, L.-S. (2002). Multiple-relaxation-time lattice Boltzmann models in three dimensions. Philosophical Transactions. Series A, Mathematical, Physical, and Engineering Sciences, 360(1792), 437–51.
 * http://doi.org/10.1098/rsta.2001.0955
 *
 * [3] Reis T, Phillips TN. Lattice Boltzmann model for simulating immiscible two-phase flows. J Phys A: Math Theor 2007;40:4033–53.
 * https://doi.org/10.1088/1751-8113/40/14/018
 */

void collisionSwap(double* ini_R, double* ini_B, double* ini_W,  stForce& forceRed, stForce& forceBlue,  stCollision& sMRT, int numberOfPoints);

/**
*@brief Function used to interpolate the relaxation parameter \tau and change it smoothly at the interfaces between two fluids[1].
*@param    tauRed      relaxation time parameter of the Red fluid
*@param    tauBlue     relaxation time parameter of the Blue fluid
*@param    rhoR        Red density in the site where is being calculated the interpolation
*@param    rhoB        Blue density in the site where is being calculated the interpolation
*
*
*[1]Grunau D, Chen S, Eggert K. A lattice Boltzmann model for multiphase fluid flows. Phys Fluids A: Fluid Dyn 1993;5(10):2557–62.
https://doi.org/10.1063/1.858769
*
*/

double computeTau(stCollision& sMRT,double& rhoR,double& rhoB);


void generateWr(double *ini_R, double *ini_B, double *ini_W,int numberOfPoints);


}

namespace d3q19
{

/**
 * @brief   Function used to the recolloring step[1]. It consists in the maximization of the work done against the color gradient[2]
 * @param   alphaRed            parameter related to the density of the red fluid, used in the equilibrium distribution.
 * @param   alphaBlue           parameter related to the density of the blue fluid, used in the equilibrium distribution.
 * @param   Beta                parameter related to the interface width.
 * @param   grad                pointer to the color gradient.
 * @param   R                   pointer to the array of Red distribution function of the site.
 * @param   B                   pointer to the array of Blue distribution function of the site.
 * @param   blindDistribution   pointer to the array of the blind distribution function of the site.
 *
 *
 * [1]Latva-Kokko M, Rothman DH. Static contact angle in lattice Boltzmann models of immiscible fluids. Phys Rev E 2005;72(4):046701.
 * https://doi.org/10.1103/PhysRevE.72.046701
 *
 * [2]Reis T, Phillips TN. Lattice Boltzmann model for simulating immiscible two-phase flows. J Phys A: Math Theor 2007;40:4033–53.
 * https://doi.org/10.1088/1751-8113/40/14/018
 */
void shcMrtRecolloring(double *f, double *R, double *B, double* normal,double module, double beta, double alphaRed, double alphaBlue);

/**
 * @brief   Function Used to calculate the Second collision term.
 * @param   A           parameter of the fluid related to the superficial tension.
 * @param   grad        pointer to the color gradient.
 * @param   Omega2      pointer that stores the second collision term.
 * 
 * The second collision term is calculated by :
 *      \f{equation}{
 *          ( \Omega_{i}^{k} )^{2} = \frac{A}{2}|\textbf{f}|\left[  w_{i} \frac{ (e_{i} \cdot \textbf{f} )^{2} }{| \textbf{f} |^{2}} - B_i \right]
 *      \f}
 * Where[1]:
 * B_0 = - \frac{1}{3} for i = 0
 * B_i = \frac{1}{18} for i = 1,2...,6
 * B_i = \frac{1}{36} for i = 7,8,...18
 * 
 * [1]SAITO, Shimpei; ABE, Yutaka; KANEKO, Akiko; KANAGAWA, Tetsuya; IWASAWA, Yuzuru; KOYAMA, Kazuya. Simulation of a Liquid Jet using the Lattice Boltzmann Model for Immiscible Two-Phase Flow. Japanese Journal Of Multiphase Flow, [s.l.], v. 29, n. 5, p.433-441, 2016. The Japanese Society for Multiphase Flow. 
 * http://dx.doi.org/10.3811/jjmf.29.433.
 *
 */

void shcMrtForce(double* omg,  double* normal,double module, double A);
    
}

#endif
