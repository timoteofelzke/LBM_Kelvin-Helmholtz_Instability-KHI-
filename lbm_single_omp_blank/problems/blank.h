/**
 *  @file   blank.h
 *  @author Diogo Nardelli Siebert
 *  @brief  Declaration of functions used in a problem
 *
 *  This is a template problem file. The headers written
 *  here have a blank implementation on the blank.cpp file.
 *  The ideia is to copy this file with a different name
 *  when a new problem is to be added.
 *
 */

#ifndef __BLANK_H_INCLUDED__   // if the header hasn't been included yet...
#define __BLANK_H_INCLUDED__

#include <string.h>
#include "../lbmbase/collision/collision.h"

struct stParameters; // Foward Declartion of stParameters (from parameters.h)
struct stGeometry;   // Foward Declartion of stGeometry   (from geometry.h)

namespace blankProblem
{
	
/**
 * @struct stProblem
 * @author Diogo Nardellli Siebert
 * @date 20/07/15
 * @brief A data structure to encapsulate
 */
		
struct stProblem
{
	// std::string filename;
	double tau;
    double initial_density;
    int steps;
    int size_x, size_y;
};

/**
 * @brief Reads the problem parameters from the INI file.
 * @param Name of the INI file;
 * @return The stProblem Structure with the problem parameters
 */

stProblem problemParameters(std::string filename );

/**
 * @brief Set up any additional initial configuration to the problem being studied.
 * @param iniN pointer to the memory region where the distribution function of all simulation is stored.
 * @param geo The stGeometry structure with information about the simulation geometry
 * @param generalParameters The stParameters structure with the simulations parameters
 * @param specificParameters The stProblem structure with problem parameters.
 *
 * This function applies some initial conditions or setups to the problem. For the blank problem the implementation is blank and
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
 * This function performes some operatiosn before the programs finally exists. Since this is a blank problem, the implementation
 * is blank, but this function can be used to measure or save specific data and end of the simulation.
 */

void problemFinalize(double* iniN, stCollision& col , stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force, int step );

/**
 * @brief Applies the boundary condition to the problem.
 * @param iniN pointer to the memory region where the distribution function of all simulation is stored.
 * @param geo The stGeometry structure with information about the simulation geometry
 * @param generalParameters The stParameters structure with the simulations parameters
 * @param specificParameters The stProblem structure with problem parameters.
 * @param step The current simulation time step.
 *
 * This function performes some operatiosn before the programs finally exists. Since this is a blank problem, the implementation
 * is blank, but this function can be used to measure or save specific data and end of the simulation.
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
 * is blank.
 */

void problemPosCollision(double* iniN, stCollision& col , stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force, int step );

/**
 * @brief This function performe action and end of each step.
 * @param iniN pointer to the memory region where the distribution function of all simulation is stored.
 * @param geo The stGeometry structure with information about the simulation geometry
 * @param prm The stParameters structure with the simulations parameters
 * @param info The stProblem structure with drag around object problem parameters .
 * @param step The current simulation time step.
 * @return return a bool (true if the stop criteria is not yeat fullfilled and false otherwise).
 *
 * The ideia is that this function is to use this function to compute some specific data from the function distribution and write it to a file.
 * This file can also be used to impose specific stop criteria. The simulation will stop if this function return false. For the black problem, the
 * implemenation only contains a return true statement.
 *
 */

bool problemOutput(double* iniN, stCollision& col , stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force, int step );

}


#endif
