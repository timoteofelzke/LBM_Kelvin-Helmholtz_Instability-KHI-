/**
 *  @file   convective.h
 *  @author Diogo Nardelli Siebert
 *
 *  Header of functions to impose the convective boundary conditions for details
 *  please refer to the paper doi: 10.1103/PhysRevE.87.063301
 */

#ifndef __CONVECTIVE_H_INCLUDED__
#define __CONVECTIVE_H_INCLUDED__

#include <vector>

struct stGeometry;

struct stBoundaryConvective
{
    std::vector<unsigned int>  zeroLayer;     /*!< List of IDs of the sites in the boundary plane */
    std::vector<unsigned int> firstLayer;     /*!< List of IDs of the sites in the neighboring plane */
    int direction;                            /*!< Direction of the the neighboring plane (-1 or + 1) */
    int axis;                                 /*!< The axis where the boundary condition needs to be applied */
    double* backup;                      /*!< Pointer to position in memory where the distribution of boundary plane will be recorded to be used in the next step */
};

/**
 * @brief Gather information about the point where to apply in Convective Boundary Condition in XY plane
 * @param geo struct containing all the geometry infomation about the domain.
 * @param z position (in the z-axis) where the boundary condition will be applied.
 * @param deltaZ can either +1 or -1, to indicate the direction where the informations will be convected from.
 * @return It returns a StBoundaryConvective structure containg all information about the points where to apply the convective boundary condition.
 *
 * It performs a setup/configuration at the beginning of the simulation to apply the convective boundary condition[1] in constant z plane (or XY plane):
 *
 * \f$ \partial_t f_i + U_z \partial_z f_i = 0 \f$
 *
 * or in the discrete form
 *
 * \f$  f_i(z,t+ \delta_t) = \frac{ f_i(z,t) + U_z f_i(z+\Delta z,t+\Delta)}{1+U_z} \f$.
 *
 * The resulting stBoundaryConvective struct will containg all points where the boundary condition will be applied (x position points) and all the points where the information will come from
 * (x + deltaX position). It also will contain a buffer to store the information of x point, since  this will be lost in propagation step.
 *
 * [1] Lou, Q., Guo, Z., & Shi, B. (2013). Evaluation of outflow boundary conditions for two-phase lattice Boltzmann equation. Physical Review E - Statistical, Nonlinear, and Soft Matter Physics, 87(6), 1–16. http://doi.org/10.1103/PhysRevE.87.063301
 *
 **/

stBoundaryConvective defineConvectiveXY( stGeometry& geo, int z , int deltaZ );

/**
 * @brief Gather information about the point where to apply in Convective Boundary Condition in XY plane
 * @param geo struct containing all the geometry infomation about the domain.
 * @param x position (in the x-axis) where the boundary condition will be applied.
 * @param deltaX can either +1 or -1, to indicate the direction where the informations will be convected from.
 * @return It returns a StBoundaryConvective structure containg all information about the points where to apply the convective boundary condition.
 *
 * It performs a setup/configuration at the beginning of the simulation to apply the convective boundary condition[1] in constant x plane (or YZ plane):
 *
 * \f$ \partial_t f_i + U_x \partial_x f_i = 0 \f$
 *
 * or in the discrete form
 *
 * \f$  f_i(x,t+ \delta_t) = \frac{ f_i(x,t) + U_x f_i(x+\Delta x,t+\Delta)}{1+U_x} \f$.
 *
 * The resulting stBoundaryConvective struct will containg all points where the boundary condition will be applied (x position points) and all the points where the information will come from
 * (x + deltaX position). It also will contain a buffer to store the information of x point, since  this will be lost in propagation step.
 *
 * [1] Lou, Q., Guo, Z., & Shi, B. (2013). Evaluation of outflow boundary conditions for two-phase lattice Boltzmann equation. Physical Review E - Statistical, Nonlinear, and Soft Matter Physics, 87(6), 1–16. http://doi.org/10.1103/PhysRevE.87.063301
 *
 **/

stBoundaryConvective defineConvectiveYZ( stGeometry& geo, int x , int deltaX );

/**
 * @brief Gather information about the point where to apply in Convective Boundary Condition in XY plane
 * @param geo struct containing all the geometry infomation about the domain.
 * @param y position (in the y-axis) where the boundary condition will be applied.
 * @param deltaY can either +1 or -1, to indicate the direction where the informations will be convected from.
 * @return It returns a StBoundaryConvective structure containg all information about the points where to apply the convective boundary condition.
 *
 * It performs a setup/configuration at the beginning of the simulation to apply the convective boundary condition[1] in constant y plane (or XY plane):
 *
 * \f$ \partial_t f_i + U_y \partial_y f_i = 0 \f$
 *
 * or in the discrete form
 *
 * \f$  f_i(y,t+ \delta_t) = \frac{ f_i(y,t) + U_y f_i(y+\Delta y,t+\Delta)}{1+U_y} \f$.
 *
 * The resulting stBoundaryConvective struct will containg all points where the boundary condition will be applied (x position points) and all the points where the information will come from
 * (x + deltaX position). It also will contain a buffer to store the information of x point, since  this will be lost in propagation step.
 *
 * [1] Lou, Q., Guo, Z., & Shi, B. (2013). Evaluation of outflow boundary conditions for two-phase lattice Boltzmann equation. Physical Review E - Statistical, Nonlinear, and Soft Matter Physics, 87(6), 1–16. http://doi.org/10.1103/PhysRevE.87.063301
 *
 **/

stBoundaryConvective defineConvectiveXZ( stGeometry& geo, int y , int deltaY );

/**
 * @brief Allocate memory for backing up the distribution at the boundary site to be used in the next step.
 * @param boundary stBoundaryConvective structure including information about the points where the boundary condition will be apllyied
 * @param iniN pointer to the distribution function of all domain
 **/
 
void allocate(stBoundaryConvective& boundary,double* iniN);

/**
 * @brief Apply the Convective Boundary Condition to a set of sites define in a stBoundaryConvective struct.
 * @param boundary stBoundaryConvective structure including information about the points where the boundary condition will be apllyied
 * @param iniN pointer to the distribution function of all domain
 * @param rho density used in the boundary condition
 *
 * Apply the Convective Boundary Condition (zero gradient) to a set of sites specified in the stBoundaryConvectivestructure. In the X axis,
 *
 * \f$  f_i(x,t+ \delta_t) = \frac{ f_i(x,t) + U_x f_i(x+\Delta x,t+\Delta)}{1+U_x} \f$, which is the discrete form of the relation \f$ \partial_t f_i + U_x \partial_x f_i = 0 \f$.
 *
 */

void convective(stBoundaryConvective& boundary, double * iniN);

#endif
