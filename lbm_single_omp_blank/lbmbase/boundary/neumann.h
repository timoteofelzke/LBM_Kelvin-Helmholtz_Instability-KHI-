/**
 *  @file   neumann.h
 *  @author Diogo Nardelli Siebert
 *  @brief  Header of functions to impose the Neumann type boundary condition.
 */

#ifndef __NEUMANN_H_INCLUDED__   // if boundary.h hasn't been included yet...
#define __NEUMANN_H_INCLUDED__

#include <vector>
#include <stdint.h>

struct stGeometry;

struct stBoundaryNeumann
{
    std::vector<int64_t> locBoundary;  /*!< List of the location in the distribution function of the  neumann type boundary sites   */
    std::vector<int64_t> locNext;      /*!< List of the location in the distribution function of neighboring sites where the Neumann boundary Condition must be applied */
};


/**
 * @brief Gather information about the point where to apply in Neumann Boundary Condition.
 * @param geo struct containing all the geometry infomation about the domain.
 * @param z position where the boundary condition will be applied.
 * @param deltaZ increment in the position to locate the plain where the information will be copied.
 * @return It returns a StBoundaryNeummann structure containg all information about the points where to apply the neumann boundary condition.
 *
 * It performs a setup or configuration at the beginning of the simulation to apply the neumann boundary condition \f$ \partial_z f_i = 0 \f$, or in the discrete form \f$  f_i(z) = f_i(z+\Delta z) \f$
 * The resulting stBoundaryNeumann struct will containg all points where the boundary condition will be applied  and all the points where the information will come from
 * (z + deltaZ position). This Functions is to be used if the boundary condition needs to be apllied in a  z = constat plane.
 *
 */

stBoundaryNeumann defineNeumannXY( stGeometry& geo, int z , int deltaZ );

/**
 * @brief Gather information about the point where to apply in Neumann Boundary Condition.
 * @param geo struct containing all the geometry infomation about the domain.
 * @param x position  where the boundary condition will be applied.
 * @param deltaX increment in the position to locate the plain where the information will be copied.
 * @return It returns a StBoundaryNeummann structure containg all information about the points where to apply the neumann boundary condition.
 *
 * It performs a setup or configuration at the beginning of the simulation to apply the neumann boundary condition \f$ \partial_x f_i = 0 \f$, or in the discrete form \f$  f_i(x) = f_i(x+\Delta x) \f$
 * The resulting stBoundaryNeumann struct will containg all points where the boundary condition will be applied  and all the points where the information will come from
 * (x + deltaX position). This Functions is to be used if the boundary condition needs to be apllied in a  x = constat plane.
 *
 */

stBoundaryNeumann defineNeumannYZ( stGeometry& geo, int x , int deltaX );

/**
 * @brief Gather information about the point where to apply in Neumann Boundary Condition.
 * @param geo struct containing all the geometry infomation about the domain.
 * @param y position  where the boundary condition will be applied.
 * @param deltaY increment in the position to locate the plain where the information will be copied.
 * @return It returns a StBoundaryNeummann structure containg all information about the points where to apply the neumann boundary condition.
 *
 * It performs a setup or configuration at the beginning of the simulation to apply the neumann boundary condition \f$ \partial_y f_i = 0 \f$, or in the discrete form \f$  f_i(y) = f_i(y+\Delta y) \f$
 * The resulting stBoundaryNeumann struct will containg all points where the boundary condition will be applied  and all the points where the information will come from
 * (y + deltaY position). This Functions is to be used if the boundary condition needs to be apllied in a  y = constat plane.
 *
 */

stBoundaryNeumann defineNeumannXZ( stGeometry& geo, int y , int deltaY );

/**
 * @brief Apply the Neumann Boundary Condition to a set of sites
 * @param boundary stBoundaryNeumann including information about the points where to apply the boundary condition.
 * @param iniN pointer to the distribution function of all domain
 * @param rho density used in the boundary contition 
 * 
 * Apply the Neumann Boundary Condition (zero gradient) to a set of sites specified in the stBoundaryNeumann structure. It actually copies
 * the distribuntion from to the previous plain.
 */
 
void neumann(stBoundaryNeumann& boundary, double * iniN, double rho);

/**
 * @brief Apply the Neumann Boundary Condition to a set of sites. This function only uses the macroscopic information of the neighboard plain to apply the equilibrium distribution.
 * @param boundary stBoundaryNeumann including information about the points where to apply the boundary condition.
 * @param iniN pointer to the distribution function of all domain
 * @param rho density used in the boundary contition
 *
 * Apply the Neumann Boundary Condition (zero gradient) to a set of sites specified in the stBoundaryNeumann structure. It copies
 * macroscopic variables of the previous plain and set the equilibrium distribution of local plain using this variables.
 *
 *  \f{equation}{
 *  f_i( x) = f_i^{eq} (\rho(x+\Delta x), \vec u (x+\Delta x) )
 *  \f}
 */

void neumannEq(stBoundaryNeumann& boundary, double * iniN, double rho);


/**
 * @brief Apply the Neumann Boundary Condition to the velocity field of a set of sites. The density is imposed.
 * @param boundary stBoundaryNeumann including information about the points where to apply the boundary condition.
 * @param iniN pointer to the distribution function of all domain
 * @param rho density imposed in the boundary contition 
 * 
 * Apply the Neumann Boundary Condition (zero gradient) to the velocity field of a set of sites specified in the stBoundaryNeumann structure. 
 * It actually copies the distribuntion from to the previous plain, multiplying by a factor to imposed the desired rho density.
 */
 
void neumannVelocity(stBoundaryNeumann& boundary, double * iniN, double rho);

#endif
