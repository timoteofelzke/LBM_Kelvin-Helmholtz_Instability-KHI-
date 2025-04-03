/**
 *  @file   bfs.h
 *  @author Diogo Nardelli Siebert
 *  @brief  Declaration of functions used in the Backward Facing Step problem
 *
 *  Header of functions and data structers to run and analise
 *  the Backward Facing Step (BFS) problem
 *
 *  @image html bfs.png
 *
 */

#ifndef __BFS_H_INCLUDED__   // if datain.h hasn't been included yet...
#define __BFS_H_INCLUDED__

#include <string>
#include "../lbmbase/collision/collision.h"
struct stGeometry;
struct stParameters;

namespace bfsProblem
{

struct stProblem
{
    int h1;       /*!< Height of the step */
    int h2;       /*!< Height of the channel before the step  */
    int l1;       /*!< Length of the channel before the step */
    int l2;       /*!< Length of the channel after the step */
    double U;     /*!< Mean velocity at the inlet */
    int tavg;     /*!< Time used to comput average */
    int tstb;     /*!< Minimal time for starting computing results  */
    std::string filename;   /*!< Name of the file that will store  the output in csv format  */

};

/**
 * @brief Create a parabolic velocity profile in the y direction at given x position with a given mean velocity.
 * @param _x adimensional position where the velocity must be evaluated.
 *
 * This functions computes the adimensional velocity in a given position for a parabolic profile.
 *
 *  \f{equation}{
 *  \tilde{u}(\tilde x) = \frac{u}{\bar u} = -6*\left(\tilde{x} - \frac{1}{2}\right)*\left(\tilde{x} + \frac{1}{2}\right)
 *  \f}
 *
 *  Since this profile is normalized with the mean velocity to use this profile the user needs to multiple the result of this function by the desired mean velocity . The dimensionaless position variable haS values in the interval [-1/2,+1/2] and can be computed with
 *
 *  \f{equation}{
 *  \tilde{x} = \frac{x - x_0}{2L}
 *  \f}
 *
 *  where \f$x_0\f$ is the central position of the profile and \f$ L \f$ is the length of the profile.
 *
 * @return the normalized velocity in the given point.
 *
 */
 
double parabolicProfile(double _x);

/**
 * @brief Reads the Backward Facing Step parameters from the INI file.
 * @param Name of the INI file;
 * @return The stProblem Structure with the bfs parameters
 */

stProblem problemParameters( std::string filename );

/**
 * @brief Set up any additional initial configuration to the problem being studied.
 * @param iniN pointer to the memory region where the distribution function of all simulation is stored.
 * @param geo The stGeometry structure with information about the simulation geometry
 * @param generalParameters The stParameters structure with the simulations parameters
 * @param specificParameters The stProblem structure with problem parameters.
 *
 * This function applies some initial conditions or setups to the problem. For the bfs problem the implementation is blank, i.e.,
 * no additional initial condition is needed.
 */

void problemInitialize(double* iniN, stCollision& col , stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force);

/**
 * @brief Performe different tasks before the program exits.
 * @param iniN pointer to the memory region where the distribution function of all simulation is stored.
 * @param geo The stGeometry structure with information about the simulation geometry
 * @param generalParameters The stParameters structure with the simulations parameters
 * @param specificParameters The stProblem structure with problem parameters.
 * @param step The current simulation time step.
 *
 * This function performes some operatiosn before the programs finally exists. For the BFS problem the implementation
 * of this function is blank.
 */

void problemFinalize(double* iniN, stCollision& col , stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force, int step );

/**
 * @brief Applies the boundary condition in the inlet and outlet in accordance with BFS problem.
 * @param iniN pointer to the memory region where the distribution function of all simulation is stored.
 * @param geo The stGeometry structure with information about the simulation geometry
 * @param generalParameters The stParameters structure with the simulations parameters
 * @param specificParameters The stProblem structure with problem parameters.
 * @param step The current simulation time step.
 *
 * It applies at each step the boundary conditions for the BFS problem at the inlet (parabolicProfile with mean velocity bfs.U) and
 * outlet (neumann boundary condition)
 */

void problemBoundary(double* iniN, stCollision& col , stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force, int step );

/**
 * @brief Function used to perform any operations between collision and streaming.
 * @param iniN pointer to the memory region where the distribution function of all simulation is stored.
 * @param geo The stGeometry structure with information about the simulation geometry
 * @param generalParameters The stParameters structure with the simulations parameters
 * @param specificParameters The stProblem structure with problem parameters.
 * @param step The current simulation time step.
 *
 * This function performes some operatiosn before the collision and streaming at each step. The implementation
 * is blank for the BFS problem.
 */

void problemPosCollision(double* iniN, stCollision& col , stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force, int step );

/**
 * @brief Computes and writes output data for the BFS problem.
 * @param iniN pointer to the memory region where the distribution function of all simulation is stored.
 * @param geo The stGeometry structure with information about the simulation geometry
 * @param prm The stParameters structure with the simulations parameters
 * @param info The stProblem structure with drag around object problem parameters .
 * @param step The current simulation time step.
 *
 * Computes the averaged velocity in time for all points at bottom of the channel after the step and writes the
 * data to a CSV file.
 */

bool problemOutput(double* iniN, stCollision& col , stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force, int step );


}
#endif
