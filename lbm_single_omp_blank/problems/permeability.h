/**
 *  @file   permeability.h
 *  @author Diogo Nardelli Siebert
 *  @brief  Declaration of function to compute the intrinsic permeability
 *
 *  Header of functions and data structers to compute the permeability of a porous media image
 */

#ifndef __PERMEABILITY_H_INCLUDED__   // if the header hasn't been included yet...
#define __PERMEABILITY_H_INCLUDED__

#include "../lbmbase/collision/collision.h"

struct stParameters; // Foward Declartion of stParameters (from parameters.h)
struct stGeometry;   // Foward Declartion of stGeometry   (from geometry.h)

namespace permeabilityProblem
{
	
/**
 * @struct stProblem
 * @author Diogo Nardellli Siebert
 * @date 20/07/15
 * @brief A data structure to encapsulate data for the permeability evaluation problem
 **/
		
struct stProblem
{
    int axis;                   /*!< Axis where the permeability will be computed (0:x, 1:y, 2:z)  */
    std::string filename;       /*!< Name of the CSV file that will store the permeability in time */
    double err;                 /*!< The tolerance in the permeability value */
};

/**
 * @brief Reads the parameters from the INI file.for the drag around an object problem
 * @param Name of the INI file;
 * @return The stProblem structure with the drag problem parameters read.
 **/

stProblem problemParameters( std::string filename );

/**
 * @brief Set up any additional initial configuration to the problem being studied.
 * @param iniN pointer to the memory region where the distribution function of all simulation is stored.
 * @param geo The stGeometry structure with information about the simulation geometry
 * @param generalParameters The stParameters structure with the simulations parameters
 * @param specificParameters The stProblem structure with problem parameters.
 *
 * This function applies some initial conditions or setups to the problem. For the this problem the implementation is blank and
 * no special initial condition is applied, i. e., the implementation of this function is left blank
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
 *
 * This function performes some operatiosn before the programs finally exists. For this problem the implementation
 * is blank.
 */

void problemFinalize(double* iniN, stCollision& col , stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force, int step );

/**
* @brief           Empty function since there is no special boundary conditions in this problem
* @param iniN      pointer to the distribution function
* @param geo       A reference to the stGeometry containing all the geometry information
* @param prm       A struct stParameters containing all the initial parameters of the simulation
* @param info      A struct stProblem containing all the information needed in the permeability determination problem
* @param step      iteration value
**/

void problemBoundary(double* iniN, stCollision& col , stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force, int step );

/**
 * @brief Function used to perform any operations between collision and streaming.
 * @param iniN pointer to the memory region where the distribution function of all simulation is stored.
 * @param geo The stGeometry structure with information about the simulation geometry
 * @param generalParameters The stParameters structure with the simulations parameters
 * @param specificParameters The stProblem structure with problem parameters.
 * @param step The current simulation time step.
 *
 * This function performes some operations before the collision and streaming at each step. The implementation
 * is blank for this problem.
 */

void problemPosCollision(double* iniN, stCollision& col , stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force, int step );

/**
 * @brief Computes and writes the permeability of the porous media image over time into an CSV file.
 * @param iniN pointer to the memory region where the distribution function of all simulation is stored.
 * @param geo The stGeometry structure with information about the simulation geometry
 * @param prm The stParameters structure with the simulations parameters
 * @param info The stProblem structure with drag around object problem parameters .
 * @param step The corrent step of the simulation.
 *
 *  Computes and writes the permeability of the porous media image over time into an CSV file.
 *
 * \f{equation}{
 * k = \frac{\mu \bar u}{\Delta P}
 * \f}
 *
 */

bool problemOutput(double* iniN, stCollision& col , stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force, int step );

/**
 * @brief Computes the mean fluid velocity in all porous media.
 * @param iniN  pointer to the memory region where the distribution function of all simulation is stored.
 * @param geo   The stGeometry structure with information about the simulation geometry
 * @param axis  Axis where the permeability will be computed (0:x, 1:y, 2:z)
 * @return the value of the surface velocity
 *
 * \f{equation}{
 * \bar u  = \frac{1}{V} \int_{V_p} u_{\alpha} dV
 * \f}
 *
 * where V_p is the porous volume and V is the sample volume.
 **/

double phaseVelocity(double* iniN, stGeometry &geo, int axis, double g);


}



#endif
