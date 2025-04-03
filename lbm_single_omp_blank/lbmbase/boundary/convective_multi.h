/**
 *  @file   convective_multi.h
 *  @author Diogo Nardelli Siebert
 *
 *  Header of functions to impose the convective boundary conditions for details
 *  please refer to the paper doi: 10.1103/PhysRevE.87.063301
 */

#ifndef __CONVECTIVE_MULTI_H_INCLUDED__
#define __CONVECTIVE_MULTI_H_INCLUDED__

#include <vector>

struct stGeometry;

namespace multi
{

enum boundaryPlane { X_BEGIN = 0, X_END = 1, Y_BEGIN = 2, Y_END = 3, Z_BEGIN = 4 , Z_END = 5};
enum velocityType  { MEAN_VELOCITY , LOCAL_VELOCITY };

struct stBoundaryConvective
{
    std::vector<unsigned int> neighborLayer;  /*!< List of IDs of the sites in the first neighboring plane */
    std::vector<unsigned int> boundaryLayer;  /*!< List of IDs of the sites in the boundary plane */
    int direction;                            /*!< Direction of the the neighboring plane (-1 or + 1) */
    int axis;                                 /*!< The axis where the boundary condition needs to be applied */
    int coordinate;
    double* wLast;                         /*!< Pointer to position in memory where the red distribution of the plane next to the boundary will be recorded to be used in the next step */
    double* blueLast;                      /*!< Pointer to position in memory where the red distribution of the boundary plane will be recorded to be used in the next step */
    double* redLast;                        /*!< Pointer to position in memory where the blue distribution of the plane next to the boundary will be recorded to be used in the next step */
    double* U;
    double meanU;
    velocityType uType;
};

/**
 * @brief Gathers information about the points where to apply the Convective Boundary Condition in XY plane
 * @param geo struct containing all the geometry infomation about the domain.
 * @param plane one of the six planes of the domain (X_BEGIN, X_END , Y_BEGIN, Y_END , Z_BEGIN , Z_END )
 * @param uType wheter to use local or mean velocity in the material derivative.
 * @return It returns a StBoundaryConvective structure containg all information about the points where to apply the convective boundary condition.
 *
 * It performs a setup/configuration at the beginning of the simulation to apply the convective boundary condition[1]
 *
 * \f$ \partial_t f_i + u_n \partial_n f_i = 0 \f$
 *
 * where $u_n$ and $\partial_n$ are the normal (pointing outwards of the domain) component of the velocity and of the gradient. Using
 * first order finit difference:
 *
 * \f$  f_i(\vec{x}_b,t+1) = \frac{ f_i(\vec{x}_b,t) + u_n f_i(\vec{x_b} - \hat{n} ,t+1)}{1+u_n} \f$.
 *
 * The resulting stBoundaryConvective struct will contain all ids of the \f$ \vec{x}_b \f$ points where the boundary condition will be applied and all the \f$ \vec{x_b} - \hat{n} \f$ neighboring points.
 * It also contains a buffer to store \f$ f_i(\vec{x}_b,t+1) \f$ to apply the boundary condition in the next step, since this value will be lost in the streaming step.
 *
 * [1] Lou, Q., Guo, Z., & Shi, B. (2013). Evaluation of outflow boundary conditions for two-phase lattice Boltzmann equation. Physical Review E - Statistical, Nonlinear, and Soft Matter Physics, 87(6), 1–16. http://doi.org/10.1103/PhysRevE.87.063301
 *
 **/

stBoundaryConvective defineConvective(stGeometry& geo, double* iniR, double* iniB, double* iniW, boundaryPlane plane, velocityType uType);

/**
 * @brief Computes the velocity in each neighboring point to be used by the apply applyConvectiveRB and applyConvectiveW functions
 * @param boundary the stBounndaryConvective storing the information about the boundary region
 * @param iniR pointer to the distribuntion functions for the red fluid
 * @param iniB pointer to the distribuntion functions for the blue fluid
 **/

void updateConvectiveVelocity(stBoundaryConvective& boundary,double* iniR,double* iniB);

/**
 * @brief Allocates the memory for the buffer to store  \f$ R_i(\vec{x}_b,t+1) \f$, \f$ B_i(\vec{x}_b,t+1) \f$ and \f$ W_i(\vec{x}_b,t+1) \f$ and sets the buffer to the value
 *        of the zero step.
 * @param boundary the stBounndaryConvective storing the information about the boundary region
 * @param iniR pointer to the distribuntion functions for the red fluid
 * @param iniB pointer to the distribuntion functions for the blue fluid
 * @param iniW pointer to the distribuntion functions for the mediator (the auxiliar distribution to compute the concentration gradient)
 **/

void allocate(stBoundaryConvective& boundary,double* iniR,double* iniB ,double* iniW);

/**
 * @brief Applies the convective boundary condition for the distributions (R,B) for all points in the selected plane
 * @param iniR pointer to the distribuntion functions for the red fluid
 * @param iniB pointer to the distribuntion functions for the blue fluid
 *
 * \f$  R_i(\vec{x}_b,t+1) = \frac{ R_i(\vec{x}_b,t) + u_n R_i(\vec{x_b} - \hat{n} ,t+1)}{1+u_n} \f$.
 * \f$  B_i(\vec{x}_b,t+1) = \frac{ B_i(\vec{x}_b,t) + u_n B_i(\vec{x_b} - \hat{n} ,t+1)}{1+u_n} \f$.
 *
 * [1] Lou, Q., Guo, Z., & Shi, B. (2013). Evaluation of outflow boundary conditions for two-phase lattice Boltzmann equation. Physical Review E - Statistical, Nonlinear, and Soft Matter Physics, 87(6), 1–16. http://doi.org/10.1103/PhysRevE.87.063301
 *
 **/

void applyConvectiveRB(stBoundaryConvective& boundary, double * iniR, double * iniB);

/**
 * @brief Applies the convective boundary condition for the distribution W for all points in the selected plane
 * @param iniR pointer to the distribuntion functions for the red fluid
 * @param iniB pointer to the distribuntion functions for the blue fluid
 *
 * \f$  W_i(\vec{x}_b,t+1) = \frac{ W_i(\vec{x}_b,t) + u_n W_i(\vec{x_b} - \hat{n} ,t+1)}{1+u_n} \f$.
 *
 **/

void applyConvectiveW(stBoundaryConvective& boundary, double * iniW);
}

#endif
